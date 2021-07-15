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
parser.add_argument('--save-filtered-metadata', action='store_true', default=False,
                    dest='save_filtered_local',
                    help='Save filtered GISAID metadata file in local/metadata_filtered.tsv')

##############################################################################################
####################   Designate variants for breakout columns    ############################
##############################################################################################
# every pango lineage in the following lists are broken out as new columns verbatim plus an "All lineages" column as the sum of all sequences and "Other lineages" being "All lineages" minus these designated variant counts
vocs = ['B.1.1.7', # alpha
        'B.1.351', 'B.1.351.2', 'B.1.351.3', # beta
        'P.1', 'P.1.1', 'P.1.2', # gamma
        'B.1.617.2', 'AY.1','AY.2' # delta
        ]
vois = ['B.1.525', # eta
        'B.1.526', # iota
        'B.1.617.1', # kappa
        'C.37', # lambda
        ]
other_important = ['B.1.427', 'B.1.429', 'B.1.427/429', # epsilon
                    'P.2', # zeta
                    'P.3', # theta
                    'B.1.617.3' # CDC VOI
                    ]

# if needed, combine separate lineages under one column. Note that this supersedes the lineage breakout column designations above, i.e. B.1.427 and B.1.429 are collapsed into one column named `B.1.427/429`
lineage_replace_dict = {
    'B.1.427':'B.1.427/429',
    'B.1.429':'B.1.427/429',   
}

# greek naming from WHO at https://www.who.int/en/activities/tracking-SARS-CoV-2-variants/
# names are created as new cols verbatim in the weeky csv with value being the sum of sequence counts for the pango lineages in each corresponding value plus a 'who_other' column that covers the remaining sequences not designated below

greek_dict = {
    'who_alpha' : ['B.1.1.7'],
    'who_beta' : ['B.1.351', 'B.1.351.2', 'B.1.351.3'],
    'who_gamma' : ['P.1', 'P.1.1', 'P.1.2'],
    'who_delta' : ['B.1.617.2', 'AY.1','AY.2'],
    'who_allvois': ['B.1.525', # eta
                    'B.1.526', # iota
                    'B.1.617.1', # kappa
                    'C.37', # lambda
                    ],
}

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

def titlecase_location(location_name, exceptions=['and', 'or', 'the', 'a', 'of', 'in', "d'Ivoire"]):
    word_list = [word if word in exceptions else word.capitalize() for word in location_name.split(' ')]
    return ' '.join(word_list)

def correct_location_names(gisaid_df):
    gisaid_df.loc[:,'country'] = gisaid_df['country'].apply(titlecase_location)
    gisaid_df.loc[gisaid_df['country'].fillna('').str.contains('USA', case=False), 'country'] = 'United States'
    gisaid_df.loc[gisaid_df['country'] == 'Puerto Rico', 'country'] = 'United States'
    gisaid_df.loc[gisaid_df['country'] == 'Guam', 'country'] = 'United States'
    gisaid_df.loc[gisaid_df['country'] == 'Northern Mariana Islands', 'country'] = 'United States'
    gisaid_df.loc[gisaid_df['country'] == 'Czech Republic', 'country'] = 'Czechia'
    gisaid_df.loc[gisaid_df['country'] == 'Antigua', 'country'] = 'Antigua and Barbuda'
    gisaid_df.loc[gisaid_df['country'] == 'Democratic Republic of the Congo', 'country'] = 'Democratic Republic of Congo'
    gisaid_df.loc[gisaid_df['country'] == 'Republic of the Congo', 'country'] = 'Congo'
    gisaid_df.loc[gisaid_df['country'] == 'Faroe Islands', 'country'] = 'Faeroe Islands'
    gisaid_df.loc[gisaid_df['country'] == 'Guinea Bissau', 'country'] = 'Guinea-Bissau'
    gisaid_df.loc[gisaid_df['country'] == 'Niogeria', 'country'] = 'Nigeria'
    # gisaid_df.loc[gisaid_df['country'] == 'mongolia', 'country'] = 'Mongolia'
    # gisaid_df.loc[gisaid_df['country'] == 'morocco', 'country'] = 'Morocco'
    # gisaid_df.loc[gisaid_df['country'] == 'belgium', 'country'] = 'Belgium'
    gisaid_df.loc[gisaid_df['country'] == 'Bosni and Herzegovina', 'country'] = 'Bosnia and Herzegovina'
    return gisaid_df

def annotate_sequences(gisaid_df):
    gisaid_df['region'] = gisaid_df.Location.apply(lambda x: x.split('/')[0].strip())
    gisaid_df['country'] = gisaid_df.Location.apply(lambda x: x.split('/')[1].strip())
    gisaid_df['division'] = gisaid_df.Location.apply(
        lambda x: x.split('/')[2].strip() if len(x.split('/'))>2 else '')

    # replace 'USA' string with 'United States' etc in location, to match OWID location name
    gisaid_df = correct_location_names(gisaid_df)

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
    print(f'Top 20 lineages:\n{gisaid_df["Pango lineage"].value_counts()[:20]}')
    print(f'Pangolin versions present:\n{gisaid_df["Pangolin version"].value_counts(dropna=False)}')
    gisaid_df = filter_gisaid_sequences(gisaid_df)
    gisaid_df = annotate_sequences(gisaid_df)
    gisaid_df = subset_gisaid_df(gisaid_df)
    return gisaid_df


def aggregate_with_lineage(gisaid_df):
    country_variants_df = gisaid_df.groupby(
        ['collect_date','collect_yearweek','collect_weekstartdate','country','Pango lineage']).count()[['Accession ID']].reset_index()

    country_variants_df['key_lineages'] = country_variants_df['Pango lineage'].apply(
        lambda x: x if x in vocs + vois + other_important else 'Other lineages')
    # manually combine lineages, i.e. B.1.427 and B.1.429 under B.1.427/429
    country_variants_df['key_lineages'].replace(lineage_replace_dict, inplace=True)
    
    all_sequences = gisaid_df.groupby(
        ['collect_date','collect_yearweek','collect_weekstartdate','country']).count()[['Accession ID']].reset_index()
    all_sequences['key_lineages'] = 'All lineages'
    country_variants_df = pd.concat([country_variants_df, all_sequences], sort=True)
 
    # rename columns a bit
    country_variants_df.columns = ['_'.join(c.lower().split()) for c in country_variants_df.columns]
    return country_variants_df

def calc_lagstats(gisaid_df, group_cols=['collect_date','country']):
    # precalculate summary stats about lag time from date_collect to date_submit for all filtered sequences per day and country
    sumstats_df = gisaid_df.groupby(group_cols).describe()['lag_days'].reset_index()
    sumstats_df.rename(columns={'count':'seq_count',
                       '50%':'gisaid_lagdays_median',
                       '25%':'gisaid_lagdays_q1',
                       '75%':'gisaid_lagdays_q3',
                       'min':'gisaid_lagdays_min',
                       'max':'gisaid_lagdays_max',
                       }, inplace=True)
    sumstats_df.drop(['seq_count','mean','std'], axis=1, inplace=True)
    return sumstats_df

##############################################################################################
######################################   OWID data load    ###################################
##############################################################################################

def load_owid_df():
    url = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv'
    owid_df = pd.read_csv(url, parse_dates=['date'])
    owid_df.sort_values(['location','date'], ascending=True, inplace=True)
    # only keep data after Dec 2019
    owid_df = owid_df[owid_df['date']>='2019-12-01']

    # drop owid region rows which start with OWID_ except specific locations
    owid_df = owid_df[~owid_df['iso_code'].isin([
                                                'OWID_AFR', # Africa
                                                'OWID_ASI', # Asia
                                                'OWID_EUR', # Europe
                                                'OWID_EUN', # European Union
                                                'OWID_INT', # International
                                                # 'OWID_KOS', # Kosovo
                                                'OWID_NAM', # North America
                                                # 'OWID_CYN', # North Cyprus
                                                'OWID_OCE', # Oceania
                                                'OWID_SAM', # South America
                                                'OWID_WRL', # World
                                                ])]

    # subset of columns, prepend 'owid_' to each column name
    owid_cols = ['date','location','iso_code','continent','new_cases','new_cases_smoothed','population','people_vaccinated','people_fully_vaccinated']
    owid_df = owid_df[owid_cols]
    
    # add vax columns calculating the daily change in people vaccinated to roll up into weekly sums
    for loc in owid_df.location.unique():
        owid_df.loc[owid_df['location']==loc,'new_people_vaccinated'] = owid_df[owid_df['location']==loc]['people_vaccinated'].ffill().fillna(0).diff()
        owid_df.loc[owid_df['location']==loc,'new_people_fully_vaccinated'] = owid_df[owid_df['location']==loc]['people_fully_vaccinated'].ffill().fillna(0).diff()        

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

    # reorganize col order
    country_variants_pivot = country_variants_pivot[['collect_date','country','placeholder_dropmeplease','All lineages']+sorted([c for c in country_variants_pivot.columns if '.' in c])+['Other lineages']]

    # merge in owid cases columns which are date-dependent
    cols = ['owid_location', 'owid_date', 'owid_new_cases', 'owid_new_cases_smoothed','owid_new_people_vaccinated','owid_new_people_fully_vaccinated']
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
        merged_df[merged_df['key_lineages'].isin(['All lineages','placeholder_dropmeplease'])][['owid_location', 'owid_continent', 'owid_population']].drop_duplicates(),
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
        ['owid_date','owid_location'], ascending=[True, True], inplace=True)

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

# OWID dataset already has a continent column so using that instead
# def add_continents(merged_df, region_path='data/continents-according-to-our-world-in-data.csv'):
#     continents = pd.read_csv(region_path)
#     merged_df = pd.merge(merged_df, continents[['Entity','Continent']], how='left', left_on=['owid_location'], right_on=['Entity'])
#     merged_df.rename(columns={'Continent': 'owid_continent'}, inplace=True)
#     merged_df.drop('Entity', axis=1, inplace=True)
#     return merged_df    

def concat_agglocations(merged_pivoted_df, group_cols=['collect_date','collect_yearweek','collect_weekstartdate','owid_date']):
    continents_df = merged_pivoted_df.groupby(group_cols+['owid_continent']).sum().reset_index()
    whoregions_df = merged_pivoted_df.groupby(group_cols+['who_region']).sum().reset_index()
    global_df = merged_pivoted_df.groupby(group_cols).sum().reset_index()
    
    # create new col to distinguish these from country-level rows
    continents_df.loc[:,'aggregate_location'] = continents_df['owid_continent']
    whoregions_df.loc[:,'aggregate_location'] = 'WHO Region: ' + whoregions_df['who_region']
    global_df.loc[:,'aggregate_location'] = 'Global'
    
    agglocation_df = pd.concat([continents_df, whoregions_df, global_df], sort=False)

    return agglocation_df

def calc_regional_lagstats(gisaid_owid_df, group_cols=['collect_weekstartdate']):
    continents_sumstats_df = calc_lagstats(gisaid_owid_df, group_cols=group_cols+['owid_continent'])
    whoregions_sumstats_df = calc_lagstats(gisaid_owid_df, group_cols=group_cols+['who_region'])
    global_sumstats_df = calc_lagstats(gisaid_owid_df, group_cols=group_cols)

    continents_sumstats_df.loc[:,'aggregate_location'] = continents_sumstats_df['owid_continent']
    whoregions_sumstats_df.loc[:,'aggregate_location'] = 'WHO Region: '+ whoregions_sumstats_df['who_region']
    global_sumstats_df.loc[:,'aggregate_location'] = 'Global'

    return pd.concat([whoregions_sumstats_df, continents_sumstats_df, global_sumstats_df], sort=False)

def cleanup_columns(merged_df, gisaid_cols):
    # prepend gisaid_ to respective columns except for the lineage ones
    renamed_cols = {c:'gisaid_'+c for c in merged_df.columns if (c in gisaid_cols)}
    merged_df.rename(columns=renamed_cols, inplace=True)
    return merged_df

def calc_vax_bottomup(vax_df, loc_col = 'owid_location'):
    # adding up the daily or weekly new people vaccinated per region and then calculating the percent of pop
    for loc in vax_df[loc_col].unique():
      vax_df.loc[vax_df[loc_col]==loc,'owid_people_vaccinated'] = vax_df[vax_df[loc_col]==loc]['owid_new_people_vaccinated'].cumsum()
      vax_df.loc[vax_df[loc_col]==loc,'owid_people_fully_vaccinated'] = vax_df[vax_df[loc_col]==loc]['owid_new_people_fully_vaccinated'].cumsum()
    
    vax_df['owid_people_vaccinated_per_hundred'] = np.round(vax_df['owid_people_vaccinated'] / (vax_df['owid_population']/100),2)
    vax_df['owid_people_fully_vaccinated_per_hundred'] = np.round(vax_df['owid_people_fully_vaccinated'] / (vax_df['owid_population']/100),2)
    return vax_df

def get_owid_vax_regional(get_weekly=True):
    iso2loc_dict = {
        'OWID_AFR':'Africa',
        'OWID_ASI':'Asia',
        'OWID_EUR':'Europe',
        'OWID_NAM':'North America',
        'OWID_OCE':'Oceania',
        'OWID_SAM':'South America',
        'OWID_WRL':'Global',
    }
    owid_url = 'https://raw.githubusercontent.com/owid/covid-19-data/master/public/data/owid-covid-data.csv'
    owid_vax_df = pd.read_csv(owid_url, parse_dates=['date'])[['iso_code','date','people_vaccinated_per_hundred','people_fully_vaccinated_per_hundred']]
    owid_vax_df.columns = ['owid_%s' % x for x in owid_vax_df.columns]
    owid_vax_df['aggregate_location'] = owid_vax_df['owid_iso_code'].map(iso2loc_dict)
    owid_vax_df['gisaid_collect_weekstartdate'] = owid_vax_df['owid_date'].apply(get_weekstartdate)
    owid_vax_df.drop('owid_iso_code', axis=1, inplace=True)

    region_vax_df = owid_vax_df[(~owid_vax_df['aggregate_location'].isna())]
    
    # only use data from the last day of each week to represent the week's progress
    if get_weekly: region_vax_df = region_vax_df[(region_vax_df['owid_date'].dt.weekday==6)]
    
    return region_vax_df

def overwrite_vax_regional(df):
    # overwrite these calculated fields at the continent and global level with OWID reported values because of discrepancy vs the bottom-up calcs
    overwrite_cols = ['owid_people_vaccinated_per_hundred','owid_people_fully_vaccinated_per_hundred']
    regional_vax_weekly_df = get_owid_vax_regional()
    overwrite_locations = regional_vax_weekly_df['aggregate_location'].unique()
    
    df.loc[df['aggregate_location'].isin(overwrite_locations), overwrite_cols] = np.nan
    df.set_index(['aggregate_location','gisaid_collect_weekstartdate'], inplace=True)
    regional_vax_weekly_df.set_index(['aggregate_location','gisaid_collect_weekstartdate'], inplace=True)
    for col in overwrite_cols: 
        df[col].update(regional_vax_weekly_df[col])
    df.reset_index(inplace=True)
    return df

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
    return weekly_agg_df

def add_greek_cols(df):
    for k, v in greek_dict.items(): df[k] = df[v].sum(axis=1)
    df['who_other'] = df['All lineages'] - df[[v for v in greek_dict.values() for v in v]].sum(axis=1)
    return df

def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)

    print('Loading and filtering GISAID data...')
    gisaid_df = load_and_filter_gisaid_df(args)
    gisaid_cols = list(gisaid_df.columns)
    print('Done, %d sequences' % gisaid_df.shape[0])

    print('Aggregating GISAID data...')
    print('Break out these pango lineages:', vocs+vois+other_important)
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
    # merged_pivoted_df = add_continents(merged_pivoted_df)
    print('Aggregate locations and concatenate...')
    agglocation_df = concat_agglocations(merged_pivoted_df)
    merged_pivoted_df = pd.concat([merged_pivoted_df, agglocation_df], sort=False)
    print('Add submission lag stats...')
    sumstats_df = calc_lagstats(gisaid_df)  
    merged_pivoted_df = pd.merge(merged_pivoted_df, sumstats_df, how='left')
    print('Final data file cleanup...')
    merged_pivoted_df = cleanup_columns(merged_pivoted_df, gisaid_cols)
    print(f'Locations without OWID join and how many sequences:\n{merged_pivoted_df[(merged_pivoted_df["owid_location"].isna())&(merged_pivoted_df["aggregate_location"].isna())].groupby("gisaid_country").sum()["All lineages"]}')
    print('Done.')

    max_gisaid_date = gisaid_df.submit_date.max()
    merged_pivoted_df_latest = merged_pivoted_df.loc[(merged_pivoted_df.owid_date <= max_gisaid_date)]
    merged_pivoted_df_latest.to_csv(args.merged_gisaid_owid_out, index=False)
    print('Wrote output to %s' % args.merged_gisaid_owid_out)

    if args.make_weekly:
        print('Also creating weekly aggregate file...')
        weekly_df = aggregate_weekly(merged_pivoted_df_latest)
        weekly_df = calc_vax_bottomup(weekly_df, loc_col='owid_location')
        print('Calculating weekly lagtime stats...')
        weekly_sumstats_df = calc_lagstats(gisaid_df, group_cols=['collect_weekstartdate','country'])
        weekly_sumstats_df = cleanup_columns(weekly_sumstats_df, gisaid_cols)
        weekly_df = pd.merge(weekly_df, weekly_sumstats_df, how='left')

        print('Add regional aggregates and calculate weekly regional lagtime stats...')        
        # TODO to refactor
        agglocation_weekly_df = concat_agglocations(weekly_df, group_cols=['gisaid_collect_weekstartdate','gisaid_collect_yearweek'])
        
        print('Add regional vax data...')
        agglocation_weekly_df = calc_vax_bottomup(agglocation_weekly_df, loc_col='aggregate_location')
        agglocation_weekly_df = overwrite_vax_regional(agglocation_weekly_df)

        # drop the aggregated lag stats, need to recalc these
        agglocation_weekly_df.drop([c for c in agglocation_weekly_df.columns if 'lagdays' in c], axis=1, inplace=True)
        
        gisaid_owid_df = add_regions(merge_gisaid_owid(gisaid_df, owid_df))
        regional_sumstats = cleanup_columns(calc_regional_lagstats(gisaid_owid_df), gisaid_cols)
        agglocation_weekly_df = pd.merge(agglocation_weekly_df, regional_sumstats, how='left', 
                                left_on=['aggregate_location','gisaid_collect_weekstartdate'], 
                                right_on=['aggregate_location','gisaid_collect_weekstartdate'],
                                suffixes=('','_drop'))
        agglocation_weekly_df.drop([c for c in agglocation_weekly_df.columns if '_drop' in c], axis=1, inplace=True)

        weekly_df = pd.concat([weekly_df, agglocation_weekly_df], sort=False)
        
        # precalculate cases per mil and percent sequenced for all rows
        weekly_df = calculate_cols(weekly_df)
        
        # add cols for each WHO greek-named VOC + "All VOIs" + "Other"
        weekly_df = add_greek_cols(weekly_df)

        # cut off weekly timeseries at most recent completed week of reporting
        latest_weekstartdate = weekly_df['gisaid_collect_weekstartdate'].max()
        latest_week_numdays = (max_gisaid_date - latest_weekstartdate).days 
        print(f'Latest GISAID submission date: {max_gisaid_date}, latest week starting Monday: {latest_weekstartdate}')
        if latest_week_numdays < 6:
            print('Dropping latest week because of incomplete data')
            weekly_df = weekly_df.loc[(weekly_df.gisaid_collect_weekstartdate < latest_weekstartdate)]

        weekly_df.to_csv(args.merged_gisaid_owid_out.split('.')[0]+'_weekly.csv', index=False)
        print(f"Wrote output to {args.merged_gisaid_owid_out.split('.')[0]+'_weekly.csv'}")

    if args.save_filtered_local:
        print('Saving filtered GISAID metadata to local/metadata_filtered.tsv...')
        gisaid_df.to_csv('./local/metadata_filtered.tsv', index=False, sep='\t')
        print('Done.')

if __name__ == "__main__":
    main()
