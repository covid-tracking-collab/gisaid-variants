"""
Processes GISAID metadata to generate a digested set of data files, and merges with population
and case counts from Our World Of Data.

These are the rules we're using to filter GISAID sequences:
    - at least 20000 in length
    - collection date must not be in the future and must be at the granularity of year/month/day,
      earliest Dec 2019
    - excluding sequences from the Nextstrain exclude list
    - only human samples

Filters for future consideration (not doing these yet):
    - excluding sequences with greater than 5% ambiguous base calls (N)
"""

import datetime
from datetime import date, timedelta
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import numpy as np
import pandas as pd
import requests


parser = ArgumentParser(
    description=__doc__,
    formatter_class=RawDescriptionHelpFormatter)


parser.add_argument('--gisaid-metadata-file', default='', help='GISAID metadata file export path')
parser.add_argument('--merged-gisaid-owid-out', default='',
    help='Path where to write merged GISAID/OWID CSV')
parser.add_argument('--make-weekly-file', action='store_true', default=False,
                    dest='make_weekly',
                    help='Create weekly file with _weekly appended to filename. Daily file is created regardless.')


##############################################################################################
######################################   GISAID data load    #################################
##############################################################################################


def filter_nextstrain_exclude_sequences(gisaid_df):
    # filter out sequences from the exclusion list

    # TODO: do this less awkwardly
    url = 'https://raw.githubusercontent.com/nextstrain/ncov/master/defaults/exclude.txt'
    exclude_sequences_response = requests.get(url)
    exclude_sequences = set(
        [x for x in exclude_sequences_response.text.split('\n') if x != '' \
            and not x.startswith('#')])

    # prepend "hCoV-19/" to match the GISAID metadata column "Virus name"
    exclude_sequences = set(["hCoV-19/%s" % x for x in exclude_sequences])
    gisaid_df = gisaid_df.loc[(~gisaid_df['Virus name'].isin(exclude_sequences))]
    return gisaid_df


def filter_gisaid_sequences(gisaid_df):
    """
    Filter GISAID sequences. These are the rules we're using:
    - at least 20000 in length
    - collection date must not be in the future and must be at the granularity of year/month/day,
      earliest Dec 2019
    - excluding sequences from the Nextstrain exclude list
    - only human samples
    """
    today_str = date.today().strftime('%Y-%m-%d')

    # full date string, between Dec 2019 and today (no future samples allowed)
    def is_legit_date(collection_date):
        return len(collection_date) == 10 and \
            collection_date > '2019-12-01' and collection_date <= today_str
    gisaid_df = gisaid_df[gisaid_df['Collection date'].apply(is_legit_date)]

    gisaid_df = gisaid_df.loc[
        (gisaid_df['Sequence length'] > 20000) &
        (gisaid_df['Collection date'] < today_str) &
        # (gisaid_df['N-Content'] < 0.05) &
        (gisaid_df['Host'] == 'Human')
    ]

    gisaid_df = filter_nextstrain_exclude_sequences(gisaid_df)
    return gisaid_df

def get_weekstartdate(dt_value):
    start = dt_value - timedelta(days=dt_value.weekday())
    return start

def annotate_sequences(gisaid_df):
    gisaid_df['region'] = gisaid_df.Location.apply(lambda x: x.split('/')[0].strip())
    gisaid_df['country'] = gisaid_df.Location.apply(lambda x: x.split('/')[1].strip())
    gisaid_df['division'] = gisaid_df.Location.apply(
        lambda x: x.split('/')[2].strip() if len(x.split('/'))>2 else '')

    # replace 'USA' string with 'United States' etc in location, to match OWID location name
    gisaid_df.loc[gisaid_df[gisaid_df['country'] == 'USA'].index, 'country'] = 'United States'
    gisaid_df.loc[gisaid_df[gisaid_df['country'] == 'Czech Republic'].index, 'country'] = 'Czechia'

    gisaid_df['collect_date'] = pd.to_datetime(gisaid_df['Collection date'])
    gisaid_df['submit_date'] = pd.to_datetime(gisaid_df['Submission date'])

    gisaid_df['lag_days'] = gisaid_df['submit_date'] - gisaid_df['collect_date']
    gisaid_df['lag_days'] = gisaid_df['lag_days'].dt.days.astype('int')

    # using ISO 8601 year and week (Monday as the first day of the week. Week 01 is the week containing Jan 4)
    gisaid_df['collect_yearweek'] = gisaid_df['collect_date'].apply(lambda x: datetime.datetime.strftime(x, "%G-W%V"))
    gisaid_df['submit_yearweek'] = gisaid_df['submit_date'].apply(lambda x: datetime.datetime.strftime(x, "%G-W%V"))

    gisaid_df['collect_weekstartdate'] = gisaid_df['collect_date'].apply(get_weekstartdate)
    gisaid_df['submit_weekstartdate'] = gisaid_df['submit_date'].apply(get_weekstartdate)

    return gisaid_df


def subset_gisaid_df(gisaid_df):
    cols = ['Collection date','Accession ID','Pango lineage',
            'Location','region','country','division',
            'collect_date', 'submit_date', 'lag_days',
            'collect_yearweek','collect_weekstartdate',
            'submit_yearweek','submit_weekstartdate']
    return gisaid_df[cols]


def load_and_filter_gisaid_df(args):
    gisaid_df = pd.read_csv(args.gisaid_metadata_file, sep='\t')
    gisaid_df = filter_gisaid_sequences(gisaid_df)
    gisaid_df = annotate_sequences(gisaid_df)
    gisaid_df = subset_gisaid_df(gisaid_df)
    return gisaid_df


def aggregate_with_lineage(gisaid_df):
    country_variants_df = gisaid_df.groupby(
        ['collect_date','collect_yearweek','collect_weekstartdate','country','Pango lineage']).count()[['Accession ID']].reset_index()
    
    # TODO: do lineage assignment differently
    vocs = ['B.1.1.7', 'B.1.429', 'B.1.427', 'P.1', 'B.1.351']
    vois = ['B.1.526', 'B.1.526.1', 'B.1.526.2', 'B.1.525', 'P.2', 'B.1.617']
    other_important = ['B.1.617.1', 'B.1.617.2', 'B.1.617.3']

    country_variants_df['key_lineages'] = country_variants_df['Pango lineage'].apply(
        lambda x: x if x in vocs + vois + other_important else 'Other')
    
    all_sequences = gisaid_df.groupby(
        ['collect_date','collect_yearweek','collect_weekstartdate','country']).count()[['Accession ID']].reset_index()
    all_sequences['key_lineages'] = 'All lineages'
    country_variants_df = pd.concat([country_variants_df, all_sequences])
 
    # rename columns a bit
    country_variants_df.columns = ['_'.join(c.lower().split()) for c in country_variants_df.columns]
    return country_variants_df

def calc_lagstats(gisaid_df, group_cols=['collect_date','country']):
    # precalculate summary stats about lag time from date_collect to date_submit for all filtered sequences per day and country
    sumstats_df = gisaid_df.groupby(group_cols).describe()['lag_days'].reset_index()
    sumstats_df.rename(columns={'count':'seq_count',
                       '50%':'lagdays_median',
                       '25%':'lagdays_q1',
                       '75%':'lagdays_q3',
                       'min':'lagdays_min',
                       'max':'lagdays_max',
                       }, inplace=True)
    sumstats_df.drop(['seq_count','mean','std'], axis=1, inplace=True)
    return sumstats_df

##############################################################################################
######################################   OWID data load    ###################################
##############################################################################################

def load_owid_df():
    url = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv'
    owid_df = pd.read_csv(url, parse_dates=['date'])
    # only keep data after Dec 2019
    owid_df = owid_df[owid_df['date']>='2019-12-01']

    # drop owid region rows which start with OWID_
    owid_df = owid_df[~owid_df['iso_code'].str.contains('OWID_')]

    # subset of columns, prepend 'owid_' to each column name
    owid_cols = ['date','location','iso_code','new_cases','new_cases_smoothed','population']
    owid_df = owid_df[owid_cols]
    owid_df.columns = ['owid_%s' % x for x in owid_df.columns]
    return owid_df


##############################################################################################
##################################      Merge and pivot data   ###############################
##############################################################################################

def merge_gisaid_owid(country_variants_df, owid_df):    
    merged_df = pd.merge(country_variants_df, owid_df,
        how='outer',
        left_on=['collect_date', 'country'],
        right_on=['owid_date', 'owid_location'],
    )

    # copy missing values from OWID data to GISAID columns and vice versa as necessary (GISAID will lag)
    merged_df['collect_date']= np.where(
        merged_df['collect_date'].isnull(), merged_df['owid_date'], merged_df['collect_date'])
    merged_df['owid_date']= np.where(
        merged_df['owid_date'].isnull(), merged_df['collect_date'], merged_df['owid_date'])
    merged_df['country']= np.where(
        merged_df['country'].isnull(), merged_df['owid_location'], merged_df['country'])
    
    return merged_df


def pivot_merged_df(merged_df):
  # # add a placeholder in order to pivot on key_lineages and not drop empty collect_date rows
    merged_df['key_lineages'] = merged_df['key_lineages'].fillna('placeholder_dropmeplease')
    
    country_variants_pivot = merged_df.pivot_table(
        values=['accession_id'], aggfunc='sum', dropna=False, index=['collect_date','country'],
        columns=['key_lineages']).droplevel(axis=1, level=0).reset_index()

    # merge in owid cases columns which are date-dependent
    cols = ['owid_location', 'owid_date', 'owid_new_cases', 'owid_new_cases_smoothed',]
    country_variants_all_lineages = merged_df[
        merged_df['key_lineages'].isin(['All lineages','placeholder_dropmeplease'])][cols]
    country_variants_pivot = pd.merge(
        country_variants_pivot, country_variants_all_lineages, 
        how='left',
        left_on=['country','collect_date'],
        right_on=['owid_location','owid_date'],
    )

    # merge in owid population regardless of date
    country_variants_pivot = pd.merge(
        country_variants_pivot, 
        merged_df[merged_df['key_lineages'].isin(['All lineages','placeholder_dropmeplease'])][['owid_location', 'owid_population']].drop_duplicates(),
        how='left',
        left_on=['country'],
        right_on=['owid_location'],
        suffixes=('_drop',''),
    )
    # drop the original owid_location col which is sparser across the timeseries
    country_variants_pivot.drop('owid_location_drop', axis=1, inplace=True)

    # copy over collect_date where owid_date is missing, otherwise owid_date is NaT which gets filtered out
    country_variants_pivot['owid_date']= np.where(
        country_variants_pivot['owid_date'].isnull(), country_variants_pivot['collect_date'], country_variants_pivot['owid_date'])
    country_variants_pivot.sort_values(
        ['owid_date','owid_location'], ascending=[False, True], inplace=True)

    # fill out the yearweek and weekstartdate missing value
    country_variants_pivot['collect_yearweek'] = country_variants_pivot['collect_date'].apply(lambda x: datetime.datetime.strftime(x, "%G-W%V"))
    country_variants_pivot['collect_weekstartdate'] = country_variants_pivot['collect_date'].apply(get_weekstartdate)

    country_variants_pivot.drop('placeholder_dropmeplease', axis=1, inplace=True)

    return country_variants_pivot

def add_regions(merged_df, region_path='data/who-regions.csv'):
    who_regions = pd.read_csv(region_path)
    merged_df = pd.merge(merged_df, who_regions[['Entity','WHO region']], how='left', left_on=['owid_location'], right_on=['Entity'])
    merged_df.rename(columns={'WHO region': 'who_region'}, inplace=True)
    merged_df.drop('Entity', axis=1, inplace=True)
    return merged_df

def add_continents(merged_df, region_path='data/continents-according-to-our-world-in-data.csv'):
    continents = pd.read_csv(region_path)
    merged_df = pd.merge(merged_df, continents[['Entity','Continent']], how='left', left_on=['owid_location'], right_on=['Entity'])
    merged_df.rename(columns={'Continent': 'owid_continent'}, inplace=True)
    merged_df.drop('Entity', axis=1, inplace=True)
    return merged_df    

def concat_agglocations(merged_pivoted_df, group_cols=['collect_date','collect_yearweek','collect_weekstartdate','owid_date']):
    continents_df = merged_pivoted_df.groupby(group_cols+['owid_continent']).sum().reset_index()
    whoregions_df = merged_pivoted_df.groupby(group_cols+['who_region']).sum().reset_index()
    global_df = merged_pivoted_df.groupby(group_cols).sum().reset_index()
    
    # create new col to distinguish these from country-level rows
    continents_df.loc[:,'aggregate_location'] = continents_df['owid_continent']
    whoregions_df.loc[:,'aggregate_location'] = 'WHO Region: ' + whoregions_df['who_region']
    global_df.loc[:,'aggregate_location'] = 'Global'
    
    agglocation_df = pd.concat([continents_df, whoregions_df, global_df], sort=False)

    return pd.concat([merged_pivoted_df, agglocation_df], sort=False)

def cleanup_columns(merged_df):
    # prepend gisaid_ to respective columns except for the lineage ones
    renamed_cols = {c:'gisaid_'+c for c in merged_df.columns if ('collect' in c) or ('country' in c) or ('lagdays' in c)}
    merged_df.rename(columns=renamed_cols, inplace=True)
    return merged_df

def calculate_cols(df):
    df['sequences_over_new_cases'] = df['All lineages'] / df['owid_new_cases']
    df.replace(np.inf, np.nan, inplace=True)
    df['new_cases_per_mil'] = df['owid_new_cases'] / (df['owid_population']/1e6)
    return df

def aggregate_weekly(df):
    df = df[[c for c in df.columns if 'lagdays' not in c]] # drop the lagday summary stat cols, need to recompute these
    weekly_agg_df = df.groupby(['gisaid_collect_weekstartdate','gisaid_collect_yearweek','gisaid_country']).sum().reset_index()
    weekly_agg_df.drop(['owid_population','owid_new_cases_smoothed'], axis=1, inplace=True)
    weekly_agg_df = pd.merge(weekly_agg_df, df[['owid_location','owid_population','who_region','owid_continent']].drop_duplicates(), 
                            how='left', left_on=['gisaid_country'], right_on=['owid_location'])
    weekly_agg_df = concat_agglocations(weekly_agg_df, group_cols=['gisaid_collect_weekstartdate','gisaid_collect_yearweek'])
    return weekly_agg_df


def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)

    print('Loading and filtering GISAID data...')
    gisaid_df = load_and_filter_gisaid_df(args)
    print('Done, %d sequences' % gisaid_df.shape[0])

    print('Aggregating GISAID data...')
    gisaid_country_variants_df = aggregate_with_lineage(gisaid_df)
    print('Done.')

    print('Loading OWID data...')
    owid_df = load_owid_df()
    print('Done, %d rows' % owid_df.shape[0])

    print('Merging GISAID and OWID data...')
    merged_df = merge_gisaid_owid(gisaid_country_variants_df, owid_df)
    print('Pivoting merged data...')
    merged_pivoted_df = pivot_merged_df(merged_df)
    print('Add region assignments to countries...')
    merged_pivoted_df = add_regions(merged_pivoted_df)
    merged_pivoted_df = add_continents(merged_pivoted_df)
    print('Aggregate locations and concatenate...')
    merged_pivoted_df = concat_agglocations(merged_pivoted_df)
    print('Add submission lag stats...')
    sumstats_df = calc_lagstats(gisaid_df)  
    merged_pivoted_df = pd.merge(merged_pivoted_df, sumstats_df, how='left')
    print('Final data file cleanup...')
    merged_pivoted_df = cleanup_columns(merged_pivoted_df)
    print('Done.')

    max_gisaid_date = gisaid_df.submit_date.max()
    merged_pivoted_df_latest = merged_pivoted_df.loc[(merged_pivoted_df.owid_date <= max_gisaid_date)]
    merged_pivoted_df_latest.to_csv(args.merged_gisaid_owid_out, index=False)
    print('Wrote output to %s' % args.merged_gisaid_owid_out)

    if args.make_weekly:
        print('Also creating weekly aggregate file...')
        weekly_df = aggregate_weekly(merged_pivoted_df_latest)
        weekly_df = calculate_cols(weekly_df)
        print('Calculating weekly lagtime stats...')
        weekly_sumstats_df = calc_lagstats(gisaid_df, group_cols=['collect_weekstartdate','country'])
        weekly_sumstats_df = cleanup_columns(weekly_sumstats_df)
        weekly_df = pd.merge(weekly_df, weekly_sumstats_df, how='left')
        weekly_df.to_csv(args.merged_gisaid_owid_out.split('.')[0]+'_weekly.csv', index=False)
        print(f"Wrote output to {args.merged_gisaid_owid_out.split('.')[0]+'_weekly.csv'}")

if __name__ == "__main__":
    main()
