'''
Fastq processing
Load a FASTQ file. Generate a plot showing a mean quality and standard deviation on each
position of a read and save it to a PDF file. Save the data to a TSV file (the output file’s
header is given in the example).
Unit tests are necessary.
-----------------------------------------------------------
Input:
reads.fastq -
https://drive.google.com/file/d/1d-LuNTIXbdF4Dgmz5qet2WInYG3gzvbg/view?usp=sharing
-----------------------------------------------------------
Input
A FASTQ file.
Output
A plot and a TSV file according to the requirements.
Example:
Sample input:
@read.1
GTTTGGGGAT
+
BBBFFFFFGH
@read.2
GCCTAGTGGC
+
CCCFFFDFHH
@read.3
GTCAGCGTTT
+
@CCFFDFBHH

Sample output (see the original pdf file of the tasks):

'''


import sys
from Bio import SeqIO
import pandas as pd
from matplotlib import pyplot as plt
import pathlib
import subprocess
import os

def check_zip_status(fq_file):
    ''' Takes the fastq input file as input.
        Decompresses it if the file is gzip compressed.
        Returns the fastq file as output.
    '''
    path = pathlib.Path(fq_file)
    parent = path.parents[0]
    stem = path.stem
    if path.suffix == '.gz':
        subprocess.run('gunzip -d %s'%(path), shell=True)
        return (parent / stem, True)
    else:
        return (fq_file, False)

def parse_fastq(fq_file):
    ''' Takes a decompressed fastq file as input.
        Parses the fastq file.
        Returns Phred quality scores of all reads (list of list)
    '''
    records = SeqIO.parse(fq_file, 'fastq')
    all_phred_scores = [rec.letter_annotations['phred_quality'] for rec in records]
    phred_scores_df = pd.DataFrame(all_phred_scores)
    #phred_scores_df.to_csv('test_ref_phred_scores_df', index=False, sep='\t') # the file is used later for unit testing (for reads.fastq)
    return phred_scores_df

def prepare_stats(phred_df):
    ''' Takes Phred scores of all reads (pandas dataframe) as input.
        Computes descriptive statistics.
        Returns a pandas dataframe for plotting and saving as a tsv file.
    '''
    phred_df = pd.DataFrame(phred_df)
    desc = phred_df.describe()
    desc = desc.T
    #desc.to_csv('test_ref_desc_df_for_reads.fastq.tsv', index=False, sep='\t') # the file is used later for unit testing (for reads.fastq)
    return desc

def prepare_tsv(desc_df):
    ''' Takes the dataframe from prepare_stats() as input.
        Saves the data in a .tsv file formatted as requested.
    '''
    tsv_df = desc_df.drop(columns=[col for col in desc_df.columns if not col in ['mean', 'std']]).copy()
    tsv_df.loc[:, 'read_position'] = tsv_df.index + 1
    tsv_df = tsv_df.rename({'mean':'mean_Phred_qual', 'std':'standard_deviation_Phred_qual'}, axis=1)
    tsv_df = tsv_df.loc[:, ['read_position', 'mean_Phred_qual', 'standard_deviation_Phred_qual']]
    tsv_df.to_csv('fastq_processing_output_dataframe_Phred_mean_std_fastq_reads.tsv', index= False, sep='\t')
    return tsv_df

def plot_figure(desc_df):
    ''' Takes the dataframe from prepare_stats() as input.
        Plots and saves the figure in a pdf file.
    '''
    fig_df = desc_df.copy()
    fig_df.index += 1
    read_length = len(fig_df)
    plt.figure(figsize=(14, 8))
    fig_df.plot(y='mean', yerr='std', kind='bar', legend=None)
    plt.xticks(range(0, read_length+1, 2), fontsize=5, rotation=0)
    plt.xlabel('read position')
    plt.ylabel('mean Phred quality')
    plt.tight_layout()
    plt.savefig('fastq_processing_output_figure_mean_std_fastq_reads.pdf')
    plt.clf()
    del fig_df
    return

def check_output_file(out_file):
    ''' Takes the name of a file (the saved pdf or tsv file, for example).
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
    fastq_file = sys.argv[-1]
    if fastq_file == 'fastq_processing.py':
        raise ValueError('No input file provided...')
    fastq_file, zipped = check_zip_status(fastq_file)
    phred_scores = parse_fastq(fastq_file)
    data_df = prepare_stats(phred_scores)
    plot_figure(data_df)
    prepare_tsv(data_df)
    assert check_output_file('fastq_processing_output_figure_mean_std_fastq_reads.pdf') == True
    assert check_output_file('fastq_processing_output_dataframe_Phred_mean_std_fastq_reads.tsv') == True

if __name__ == '__main__':
    main()

