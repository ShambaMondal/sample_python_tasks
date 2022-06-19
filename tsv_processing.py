'''
TSV processing

Load a TSV file. For a column containing an alignment length plot a histogram showing a
frequency of values in this column. If there is more than one alignment for a particular query,
choose a better one. Justify your choice in a short comment within the code. Save the plot to
a PDF file. Save the CSV file with histogram data (the output file’s header is given in the
example).

Unit tests are necessary.
-----------------------------------------------------------
Input:
alignment.b6 -
https://drive.google.com/file/d/1V0f77mgYdyoO7bcl_Gj3OxJseNK818v8/view?usp=sharing
-----------------------------------------------------------

Input
A TSV file containing BLAST results in BLAST6 format (i.e. generated with -outfmt 6
option).
Output
A plot and a CSV file according to the requirements.

Example

Sample input:
read.1
38174173
read.2
30567714
read.3
100091036
read.3
100090908
NC_000008.11
100.000
38174264
7.01e-41
171
NC_000008.11
100.000
30567805
7.01e-41
171
NC_000007.14
100.000
100091087 1.21e-18
97.1
NC_000007.14
100.000
100090950 1.22e-13
80.5

Sample output:
alignment_length,abundance
52,1
92,2

### full names of the columns used below:
   1.  qseqid      query or source (e.g., gene) sequence id

   2.  sseqid      subject  or target (e.g., reference genome) sequence id

   3.  pident      percentage of identical matches

   4.  length      alignment length (sequence overlap)

   5.  mismatch    number of mismatches

   6.  gapopen     number of gap openings

   7.  qstart      start of alignment in query

   8.  qend        end of alignment in query

   9.  sstart      start of alignment in subject

 10.  send        end of alignment in subject

 11.  evalue      expect value

 12.  bitscore    bit score

'''

import pandas as pd
import numpy as np
from matplotlib import pyplot as plt
import seaborn as sns
import sys
import os
import pandas.api.types as ptypes

def preprocess_aln_file(aln_file):
    ''' Takes the alignment tsv file as input.
        Performs simple checks and preprocessing steps on the data.
        Returns a dataframe as output.
    '''
    try:
        data_df = pd.read_table(aln_file, header= None)
    except pd.errors.EmptyDataError:
        raise pd.errors.EmptyDataError('Empty file provided. No columns to parse from file...')

    rows, columns = data_df.shape
    # check number of columns
    try:
        assert columns == 12
    except AssertionError:
        raise AssertionError('The alignment file does not have the standard format of 12 columns...')

    # in case header (string type) already existed, drop first row
    try:
        assert ptypes.is_numeric_dtype(data_df.iloc[0, -1]) # asserting numeric value of bitscore of first row
    except AssertionError:
        #print('Header already exists, dropping existing header...')
        data_df = data_df.iloc[1:, :]
        data_df = data_df.reset_index(drop=True)

    # check data types in columns
    #assert all(ptypes.is_numeric_dtype(data_df[col]) for col in data_df.columns[2:])
    #assert all(ptypes.is_string_dtype(data_df[col]) for col in data_df.columns[:2])
    for col in data_df.columns[2:]:
        try:
            assert ptypes.is_numeric_dtype(data_df[col])
        except AssertionError:
            #print('Non-numeric data type encountered at column: %s , changing to numeric...'%(col))
            data_df[col] = pd.to_numeric(data_df[col])

    # rename the columns to standard format
    data_df = data_df.rename(columns= dict(zip(data_df.columns, ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])))
    
    # drop NaNs
    data_df = data_df.dropna()
    return data_df

def return_best_alignment(preprocessed_dataframe):
    ''' For all alignments against each read, the ones with the highest bitscore indicate the best alignment.
        Even though obtaining the best alignment based on the bitscore only can be achieved by the following single line:
        
        df = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]

        we decided to be conservative, and hierarchically filtered the alignments in case of tie: 
        Highest value of (sequential preference): bitscore > -log10(evalue) > pident > length
        Lowest value of (sequential preference): mismatch > gapopen

        We observed that generally most (even all) ties with the same bitscore had the exact same values in other fields as well.
        For such cases, we just keep the first instance of the identical alignments (or ties. Final step).
    '''
    df = preprocessed_dataframe.copy()

    df.loc[:, '-log10(evalue)'] = -np.log10(df.loc[:, 'evalue'])

    for max_col in ['bitscore', '-log10(evalue)', 'pident', 'length']:
        df = df.groupby('qseqid', as_index=False)[max_col].max().merge(df)
    
    for min_col in ['mismatch', 'gapopen']:
        df = df.groupby('qseqid', as_index=False)[min_col].min().merge(df)
    
    df = df.loc[df.groupby('qseqid')['bitscore'].idxmax()]
    
    # sort reads
    df.loc[:, 'qseqid_num'] = [int(v.split('.')[1]) for v in df.loc[:, 'qseqid'].values]
    df = df.sort_values('qseqid_num')
    df = df.drop(columns=['qseqid_num', '-log10(evalue)'])
    df = df.loc[:, ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']]
    df = df.reset_index(drop=True)
    #df.to_csv('test_ref_dataframe_best_alignments_for_alignment.b6.tsv', sep='\t', index=False) # alignment input: alignment.b6; the tsv file is used in unit testing
    return df

def save_csv_file(best_df):
    values_df = best_df['length'].value_counts().reset_index()
    values_df.columns = ['alignment_length', 'abundance']
    values_df = values_df.sort_values('alignment_length')
    values_df = values_df.reset_index(drop=True)
    values_df.to_csv('tsv_processing_output_histogram_data.csv', index=False)
    return values_df

def plot_histogram(best_df):
    sns.histplot(best_df, x='length', binwidth=1)
    ax = plt.gca()
    patches = ax.patches
    #xy_coords = [p.get_xy() for p in patches]
    all_x_vals = [p.get_xy()[0] for p in patches]
    all_widths = [p.get_width() for p in patches]
    all_heights = [p.get_height() for p in patches]
    ref_hist_df = pd.DataFrame({'x': all_x_vals, 'width': all_widths, 'height': all_heights})
    #ref_hist_df.to_csv('test_ref_dataframe_sns_histogram_values', sep='\t', index=False) # alignment input: alignment.b6; the tsv file is used in unit testing

    plt.savefig('tsv_processing_output_figure_histogram_alignment_length.pdf')
    #plt.clf()
    return ref_hist_df

def check_output_file(out_file):
    ''' Takes the name of a file (the saved pdf or csv file, for example).
        Returns True (boolean) if the file exists and has non-zero size.
        Returns False if the file exists but has size zero bytes.
        Raises an error if file is not found.
    '''
    try:
        statinfo = os.stat('./'+ out_file)
        #print(statinfo)
        if statinfo.st_size > 0:
            return True
        else:
            return False
    except FileNotFoundError:
        raise


def main():
    ''' Usage example: python tsv_processing.py alignment.b6 
    '''
    in_file = sys.argv[-1]
    if in_file == 'tsv_processing.py':
        raise ValueError('No input file provided...')
    preprocessed_df = preprocess_aln_file(in_file)
    best_aln_df = return_best_alignment(preprocessed_df)
    save_csv_file(best_aln_df)
    plot_histogram(best_aln_df)
    assert check_output_file('tsv_processing_output_figure_histogram_alignment_length.pdf') == True
    assert check_output_file('tsv_processing_output_histogram_data.csv') == True

if __name__ == '__main__':
    main()

