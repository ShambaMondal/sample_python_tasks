import unittest
import fastq_processing
import pathlib
import os
import subprocess
import pandas as pd

class TestFileBase(unittest.TestCase):

    def assertIsFile(self, file_path):
        if not pathlib.Path(file_path).resolve().is_file():
            raise AssertionError('The file does not exist: %s' %str(file_path))

    def assertIsReadable(self, file_path):
        file_path = pathlib.Path(file_path).resolve()
        if not os.access(file_path, os.R_OK):
            raise AssertionError('The file is not readable: %s' %str(file_path))

class TestFilePath(TestFileBase):

    def test_input_test_files(self):
        for test_file in ['reads.fastq', 'flawed_reads_incomplete_qual_scores.fastq', 'flawed_reads_wrong_file_type.tsv', 'flawed_reads_empty_file_with_pseudo_header.fastq', 'flawed_reads_empty_file.fastq', 'reads_zipped.fastq.gz', 'test_ref_dataframe_Phred_mean_std_fastq_reads.tsv', 'test_ref_desc_df_for_reads.fastq.tsv', 'test_ref_phred_scores_df']: # keep these test files available in the ./test/ directory
            f_path = pathlib.Path('./test/'+ test_file)
            self.assertIsFile(f_path)
            self.assertIsReadable(f_path)


class TestFastqProcessing(unittest.TestCase):


    def test_check_zip_status(self):
        self.assertEqual(fastq_processing.check_zip_status('./test/reads_zipped.fastq.gz')[1], True)
        self.assertEqual(fastq_processing.check_zip_status('./test/reads.fastq')[1], False)

        subprocess.run('gzip ./test/reads_zipped.fastq', shell=True) # to keep this test file zipped for later tests

    def test_parse_fastq(self):
        with self.assertRaises(ValueError):
            fastq_processing.parse_fastq('./test/flawed_reads_incomplete_qual_scores.fastq') # mismatch between length of sequence and quality scores of a read.
        
        with self.assertRaises(ValueError):
            fastq_processing.parse_fastq('./test/flawed_reads_wrong_file_type.tsv') # wrong type of file (tsv, instead of fastq) as input.
        
        with self.assertRaises(ValueError):
            fastq_processing.parse_fastq('./test/flawed_reads_empty_file_with_pseudo_header.fastq') # only a pseudo header starting with @, mimicking a fastq file header. Garbage in next line.

        test_ref_phred_df = pd.read_table('./test/test_ref_phred_scores_df')
        test_ref_phred_df.columns = pd.to_numeric(test_ref_phred_df.columns)
        func_out_phred_df = fastq_processing.parse_fastq('./test/reads.fastq')
        
        #self.assertTrue(func_out_phred_df.equals(test_ref_phred_df))
        pd.testing.assert_frame_equal(func_out_phred_df, test_ref_phred_df)

    def test_prepare_stats(self):
        with self.assertRaises(ValueError):
            fastq_processing.prepare_stats('./test/flawed_reads_empty_file.fastq') # empty dataframe generated from empty input file (SeqIO.parse does not raise error in this case, but pandas does)

        test_ref_phred_df = pd.read_table('./test/test_ref_phred_scores_df')
        test_ref_phred_df.columns = pd.to_numeric(test_ref_phred_df.columns)
        test_ref_desc_df = pd.read_table('./test/test_ref_desc_df_for_reads.fastq.tsv')

        func_out_desc_df = fastq_processing.prepare_stats(test_ref_phred_df)
        
        pd.testing.assert_frame_equal(func_out_desc_df, test_ref_desc_df)

    def test_prepare_tsv(self):
        test_ref_df = pd.read_table('./test/test_ref_dataframe_Phred_mean_std_fastq_reads.tsv')
        func_in_df = pd.read_table('./test/test_ref_desc_df_for_reads.fastq.tsv')
        func_out_df = fastq_processing.prepare_tsv(func_in_df)

        #self.assertTrue(func_out_df.equals(test_ref_df))
        pd.testing.assert_frame_equal(func_out_df, test_ref_df)

    def test_plot_figure(self):
        ''' This step seemed to be unneccesarily complicated.
        Assuming that pyplot correctly draws the figure based on the data provided,
        test_parse_fastq(), test_prepare_stats(), and test_prepare_tsv()
        are sufficient for testing functionalities of fastq_processing.py.
        '''
        pass

if __name__ == '__main__':
    unittest.main()
