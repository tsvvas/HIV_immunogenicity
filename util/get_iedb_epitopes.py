#!/usr/bin/env python

import os
import shutil
import tempfile
import argparse
import urllib3
import certifi
import pandas as pd

tcell_url = "http://www.iedb.org/downloader.php?file_name=doc/tcell_full_v3.zip"

def retrieve_from_url(url, file_name):
    http = urllib3.PoolManager(cert_reqs='CERT_REQUIRED',
                               ca_certs=certifi.where())

    with http.request('GET', url, preload_content=False) as res, open(file_name, 'wb') as out_file:
        shutil.copyfileobj(res, out_file)
        status = res.status

    return status

def get_epitopes(file_name):
    status = retrieve_from_url(tcell_url, file_name)

    if not status == 200:
        raise Exception("The server returned: {status}")

    if os.stat(file_name).st_size < 1000000:
        raise Exception("The download was incomplete")

def filter_epitopes(file_name):
    df = pd.read_csv(file_name, compression="zip", header=[0, 1])

    # filtering by epitope type, host and disease
    df = df[(df.Epitope['Object Type'] == 'Linear peptide') &
            (df.Epitope['Organism Name'] != 'Homo sapiens') &
            (df.Host.Name.isin(['Homo sapiens', 'Homo sapiens Black', 'Homo sapiens Caucasian'])) &
            (df['1st in vivo Process']['Process Type'] == 'Occurrence of infectious disease')]

    # deleting duplicates in columns (Epitope.Description, Assay['Qualitative Measure'])
    df = df.drop_duplicates(subset=[('Epitope', 'Description'),
                                    ('Assay', 'Qualitative Measure')])

    # filtering duplicates by epitope. Here we have an entire dataset of duplicates with columns
    # Epitope.Description, Assay['Qualitative Measure'] and Assay['Response Frequency']
    dups = df[df.duplicated(subset=[('Epitope', 'Description')], keep=False)]

    dups = dups[[('Epitope', 'Description'),
                 ('Assay', 'Qualitative Measure'),
                 ('Assay', 'Response Frequency')]]

    dups = dups.sort_values(by=[('Epitope', 'Description'),
                                ('Assay', 'Response Frequency')])

    # selecting only rows with epitopes that would be kept with Positive Qualitative Measure:
    # Response Frequency should be greater than 50
    positive = dups[(dups.Assay['Qualitative Measure'] != 'Negative') &
                    (dups.Assay['Response Frequency'] > 50)]

    # dropping the rest of duplicates. Here we keep the first occurrence
    positive = positive.drop_duplicates(subset=[('Epitope', 'Description')])

    # filtering the final dataset:
    # 1. epitope is not duplicated
    # 2. epitope has Positive Qualitative Measure
    # 3. epitope is not in Positive dataset and will be kept with Negative Qualitative Measure
    df = df[(~df.index.isin(dups.index)) |
            (df.index.isin(positive.index)) |
            ((~df.Epitope.Description.isin(positive.Epitope.Description)) &
             (df.Assay['Qualitative Measure'] == 'Negative'))]

    indices = [('Epitope', 'Description'),         # peptide
              ('Assay', 'Qualitative Measure'),    # immunogenicity
              ('MHC', 'Allele Name'),              # HLA allele
              ('Epitope', 'Parent Protein IRI'),   # source protein
              ('Epitope', 'Organism Name')]        # source name

    col_names = ['Peptide','Immunogenicity', 'HLA_Allele', 'Source_Protein', 'Source_Name']

    df = df[indices]
    df.columns = df.columns.droplevel()
    df.columns = col_names
    df.Source_Protein = df.Source_Protein.apply(lambda x: str(x).split('/')[-1])
    df = df.reset_index(drop=True)

    return df

def make_epitope_dataset(folder):
    tmpdir = tempfile.mkdtemp()
    file_name = os.path.join(tmpdir, tcell_url.split("/")[-1])
    path_to_save = os.path.join(folder, 'iedb_immepitopes.csv')

    get_epitopes(file_name)
    df = filter_epitopes(file_name)
    df.to_csv(path_to_save, index=False)

def main():
    parser = argparse.ArgumentParser(description='Utility script to download \
                    epitopes from IEDB (iedb.org) and make training dataset')
    parser.add_argument('--destination', '-d', default='data/',
                        help='Path where to store output')
    args = parser.parse_args()
    dest = args.destination

    make_epitope_dataset(dest)

if __name__ == "__main__":
    main()
