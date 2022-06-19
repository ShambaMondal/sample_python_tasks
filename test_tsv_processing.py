import unittest
import tsv_processing
import pathlib
import os
import subprocess
import pandas as pd
import numpy as np

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
        for test_file in ['alignment.b6', 'test_ref_dataframe_best_alignments_for_alignment.b6.tsv', 'test_ref_tsv_processing_output_histogram_data.csv', 'tmp_dataframe_5_columns_no_header.tsv', 'tmp_dataframe_5_columns.tsv', 'tmp_df_max_bitscore_per_read.tsv', 'tmp_empty_file.tsv', 'tmp_df_alignment.b6_with_header', 'test_ref_dataframe_sns_histogram_values']: # keep these 8 test files available in the ./test/ directory
            f_path = pathlib.Path('./test/'+ test_file)
            self.assertIsFile(f_path)
            self.assertIsReadable(f_path)


class TestTsvProcessing(unittest.TestCase):

    def test_preprocess_aln_file(self):

        #print('Testing files with insufficient columns...')
        with self.assertRaises(AssertionError):
            tsv_processing.preprocess_aln_file('./test/tmp_dataframe_5_columns_no_header.tsv')

        with self.assertRaises(AssertionError):
            tsv_processing.preprocess_aln_file('./test/tmp_dataframe_5_columns.tsv')

        with self.assertRaises(pd.errors.EmptyDataError):
            tsv_processing.preprocess_aln_file('./test/tmp_empty_file.tsv')

        # test case: with 'alignment.b6'
        test_out_df_1 = pd.read_table('./test/alignment.b6', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        ##test_out_df_1.to_csv('tmp_df_alignment.b6_with_header', index=False, sep='\t')
        func_out_df_1 = tsv_processing.preprocess_aln_file('./test/alignment.b6')
        #self.assertTrue(func_out_df_1.equals(test_out_df_1))
        pd.testing.assert_frame_equal(func_out_df_1, test_out_df_1)

        #print('Testing alignment file with pre-existing header...')
        # test case: with  'alignment.b6' containing header ('tmp_df_alignment.b6_with_header')
        test_out_df_2 = pd.read_table('./test/tmp_df_alignment.b6_with_header')
        func_out_df_2 = tsv_processing.preprocess_aln_file('./test/tmp_df_alignment.b6_with_header')
        #self.assertTrue(func_out_df_2.equals(test_out_df_2))
        pd.testing.assert_frame_equal(func_out_df_2, test_out_df_2)

    def test_return_best_alignment(self):

        test_in_df = pd.read_table('./test/alignment.b6', names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'])
        test_out_df = pd.read_table('./test/test_ref_dataframe_best_alignments_for_alignment.b6.tsv')
        test_out_df.loc[:, 'evalue'] = -np.log10(test_out_df.loc[:, 'evalue'])
        func_out_df = tsv_processing.return_best_alignment(test_in_df)
        func_out_df.loc[:, 'evalue'] = -np.log10(func_out_df.loc[:, 'evalue'])

        for col in test_out_df.columns: #[2:]:
            self.assertTrue(func_out_df[col].equals(test_out_df[col]))
        
        #self.assertTrue(test_out_df.equals(func_out_df))
        pd.testing.assert_frame_equal(test_out_df, func_out_df)


    def test_save_csv_file(self):

        test_in_df = pd.read_table('./test/test_ref_dataframe_best_alignments_for_alignment.b6.tsv')
        test_ref_hist_df = pd.read_csv('./test/test_ref_tsv_processing_output_histogram_data.csv')

        func_out_df = tsv_processing.save_csv_file(test_in_df)

        #self.assertTrue(func_out_df.equals(test_ref_hist_df))
        pd.testing.assert_frame_equal(func_out_df, test_ref_hist_df)

    def test_plot_histogram(self):
        
        test_ref_hist_df = pd.read_table('./test/test_ref_dataframe_sns_histogram_values')
        func_in_df = pd.read_table('./test/test_ref_dataframe_best_alignments_for_alignment.b6.tsv')
        func_out_hist_df = tsv_processing.plot_histogram(func_in_df)

        #self.assertTrue(func_out_hist_df.equals(test_ref_hist_df))
        pd.testing.assert_frame_equal(func_out_hist_df, test_ref_hist_df)
        


if __name__ == '__main__':
    unittest.main()
