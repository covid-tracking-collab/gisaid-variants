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

from datetime import date
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
            collection_date > '2019-12-01' and collection_date < today_str
    gisaid_df = gisaid_df[gisaid_df['Collection date'].apply(is_legit_date)]

    gisaid_df = gisaid_df.loc[
        (gisaid_df['Sequence length'] > 20000) &
        (gisaid_df['Collection date'] < today_str) &
        # (gisaid_df['N-Content'] < 0.05) &
        (gisaid_df['Host'] == 'Human')
    ]

    gisaid_df = filter_nextstrain_exclude_sequences(gisaid_df)
    return gisaid_df


def annotate_sequences(gisaid_df):
    gisaid_df['region'] = gisaid_df.Location.apply(lambda x: x.split('/')[0].strip())
    gisaid_df['country'] = gisaid_df.Location.apply(lambda x: x.split('/')[1].strip())
    gisaid_df['division'] = gisaid_df.Location.apply(
        lambda x: x.split('/')[2].strip() if len(x.split('/'))>2 else '')

    gisaid_df['collect_date'] = pd.to_datetime(gisaid_df['Collection date'])
    gisaid_df['submit_date'] = pd.to_datetime(gisaid_df['Submission date'])

    gisaid_df['lag_days'] = gisaid_df['submit_date'] - gisaid_df['collect_date']
    gisaid_df['lag_days'] = gisaid_df['lag_days'].dt.days.astype('int')

    gisaid_df['collect_week'] = gisaid_df['collect_date'].dt.isocalendar().week
    gisaid_df['submit_week'] = gisaid_df['submit_date'].dt.isocalendar().week
    gisaid_df['collect_week'] = gisaid_df['collect_week'].astype('int')
    gisaid_df['submit_week'] = gisaid_df['submit_week'].astype('int')

    gisaid_df['collect_year'] = gisaid_df['collect_date'].dt.isocalendar().year
    gisaid_df['submit_year'] = gisaid_df['submit_date'].dt.isocalendar().year

    return gisaid_df


def subset_gisaid_df(gisaid_df):
    cols = ['Collection date','Accession ID','Pango lineage',
            'Location','region','country','division',
            'collect_date', 'submit_date','lag_days','submit_week','collect_week']
    return gisaid_df[cols]


def load_and_filter_gisaid_df(args):
    gisaid_df = pd.read_csv(args.gisaid_metadata_file, sep='\t')
    gisaid_df = filter_gisaid_sequences(gisaid_df)
    gisaid_df = annotate_sequences(gisaid_df)
    gisaid_df = subset_gisaid_df(gisaid_df)
    return gisaid_df


def aggregate_with_lineage(gisaid_df):
    country_variants_df = gisaid_df.groupby(
        ['collect_date','country','region','Pango lineage']).count()[['Accession ID']].reset_index()
    all_sequences = gisaid_df.groupby(
        ['collect_date','country','region']).count()[['Accession ID']].reset_index()
    all_sequences['Pango lineage'] = 'All lineages'
    country_variants_df = pd.concat([country_variants_df, all_sequences])

    # TODO: do lineage assignment differently
    vocs = ['B.1.1.7', 'B.1.429', 'B.1.427', 'P.1', 'B.1.351']
    vois = ['B.1.526', 'B.1.526.1', 'B.1.526.2', 'B.1.525', 'P.2', 'B.1.617']
    other_important = ['B.1.617.1', 'B.1.617.2', 'B.1.617.3', 'All lineages']

    country_variants_df['key_lineages'] = country_variants_df['Pango lineage'].apply(
        lambda x: x if x in vocs + vois + other_important else 'Other')

    # rename columns a bit
    country_variants_df.columns = ['_'.join(c.lower().split()) for c in country_variants_df.columns]
    return country_variants_df

##############################################################################################
######################################   OWID data load    ###################################
##############################################################################################

def load_owid_df():
    url = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv'
    owid_df = pd.read_csv(url, parse_dates=['date'])
    # only keep data after Dec 2019
    owid_df = owid_df[owid_df['date']>='2019-12-01']
    
    # replace 'United States' string with 'USA' in location, to match GISAID country name
    owid_df.loc[owid_df[owid_df['location'] == 'United States'].index, 'location'] = 'USA'
    
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

    # TODO: I'm not sure about this... Why say "All lineages" instead of "Other" here?
    merged_df['key_lineages'] = merged_df['key_lineages'].fillna('All lineages')
    
    # copy missing values from OWID data to GISAID columns as necessary (GISAID will lag)
    merged_df['collect_date']= np.where(
        merged_df['collect_date'].isnull(), merged_df['owid_date'], merged_df['collect_date'])
    merged_df['country']= np.where(
        merged_df['country'].isnull(), merged_df['owid_location'], merged_df['country'])
    
    return merged_df


def pivot_merged_df(merged_df):
    country_variants_pivot = merged_df.pivot_table(
        values=['accession_id'], aggfunc='sum', dropna=False, index=['owid_date','owid_location'],
        columns=['key_lineages']).droplevel(axis=1, level=0).reset_index()

    cols = ['country', 'region', 'collect_date', 'owid_new_cases', 'owid_new_cases_smoothed',
            'owid_population']
    country_variants_all_lineages = merged_df[
        merged_df['key_lineages'] == 'All lineages'][cols]
    country_variants_pivot = pd.merge(
        country_variants_pivot, country_variants_all_lineages, 
        how='outer',
        left_on=['owid_location','owid_date'],
        right_on=['country','collect_date']

    )
    country_variants_pivot.sort_values(
        ['owid_date','owid_location'], ascending=[False, False], inplace=True)

    return country_variants_pivot


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
    print('Done.')

    # only output dates up until the latest GISAID available (collection dates will lag a few days)
    max_gisaid_date = gisaid_country_variants_df.collect_date.max()
    merged_pivoted_df_latest = merged_pivoted_df.loc[merged_pivoted_df.owid_date <= max_gisaid_date]
    merged_pivoted_df_latest.to_csv(args.merged_gisaid_owid_out)
    print('Wrote output to %s' % args.merged_gisaid_owid_out)


if __name__ == "__main__":
    main()
