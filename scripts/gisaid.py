"""
Processes GISAID metadata to generate a digested set of data files.
"""

from datetime import date
import sys

from argparse import ArgumentParser, RawDescriptionHelpFormatter
import pandas as pd
import requests


parser = ArgumentParser(
    description=__doc__,
    formatter_class=RawDescriptionHelpFormatter)


parser.add_argument('--gisaid-metadata-file', default='', help='GISAID metadata file export path')


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
    - collection date must not be in the future and must be at least year/month, not year alone
    - excluding sequences from the Nextstrain exclude list
    """
    today_str = date.today().strftime('%Y-%m-%d')
    gisaid_df = gisaid_df.loc[
        (gisaid_df['Sequence length'] > 20000) &
        (~gisaid_df['Collection date'].isin(['2020', '2021'])) &
        (gisaid_df['Collection date'] < today_str)
    ]

    gisaid_df = filter_nextstrain_exclude_sequences(gisaid_df)
    return gisaid_df


def annotate_sequences(gisaid_df):
    def get_country(location):
        try:
            return location.split(' / ')[1]
        except IndexError:
            try:
                return location.split('/')[1]
            except IndexError:
                print(location)
                return ''

    gisaid_df['country'] = gisaid_df.Location.apply(get_country)
    return gisaid_df


def load_and_filter_gisaid_df(args):
    gisaid_df = pd.read_csv(args.gisaid_metadata_file, sep='\t')
    gisaid_df = filter_gisaid_sequences(gisaid_df)
    return gisaid_df


def main(args_list=None):
    if args_list is None:
        args_list = sys.argv[1:]
    args = parser.parse_args(args_list)

    gisaid_df = load_and_filter_gisaid_df(args)
    print(gisaid_df.shape)



if __name__ == "__main__":
    main()
