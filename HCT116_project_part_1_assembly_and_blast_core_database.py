#%%
# This project is dedicated specifically for analysis of single TRINITY CCLE FILE with GENCODE v45, mTEC, RiboSeq and immunopeptidomics files.
# Objects created in this project are generating with TRINITY, Kallisto or other additional pipelines.

#%%
# This chunk is to import all libraries we need to analyse data
import os
import pandas as pd
import Bio
import numpy as np
import matplotlib.pyplot as plt

#%%
# importing SeqIO package from Biopython library
# creating an object to store TRINITY CCLE FASTA FILE
from Bio import SeqIO

trinity_ccle_fasta_advanced_blast = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\trinity_CCLE_transcripts.fasta')
 
#%%
# Creating a loop to parse through the whole TRINITY CCLE dataframe using SeqIO to store unique transcript ID and length in dictionary
# converting dictionary into indexed DataFrame object for further work

def get_transcripts_lengths(trinity_ccle_fasta_advanced_blast):
    
    trascripts_lengths = {}

    for record in SeqIO.parse(trinity_ccle_fasta_advanced_blast, 'fasta'):
        transcript_id = record.id
        sequence_length = len(record.seq)

        trascripts_lengths[transcript_id] = sequence_length
    
    return trascripts_lengths

transcripts_lengths = get_transcripts_lengths(trinity_ccle_fasta_advanced_blast)
print(transcripts_lengths)

complex_dataframe_CCLE = pd.DataFrame({'ID' : transcripts_lengths.keys(), 'Transcript Length' : transcripts_lengths.values()})
# %%
# Changing directory to the one with BLAST files outputs

os.chdir("F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\BLAST_OUTPUT")

# %%
# Checking if directory is the proper one

os.getcwd()

# %% 
# FIRST PART OF WHOLE ANALYSIS: TRINITY CCLE GENCODE BLAST ANALYSIS

# Creating an object with column names from OUTPUT (output format type 6) for our BLAST files
# Creating an object in csv format from BLAST.txt file for TRINITY CCLE GENCODE BLAST output

blast_colnames = ['source_gene_id', 'target_gene_id', 'perc_identical', 'query_len','subj_len','alignment_length', 'mismatches', 'gap_openings', 'query_start', 'query_end', 'subj_start', 'subj_end', 'e_value', 'bit_score']
trinity_CCLE_gencode_blast = pd.read_csv('trinity_CCLE_gencode_blast_advanced.txt', sep = '\t', header = None, names = blast_colnames)


#%%
# substracting query_len from subj_len and storing it in a new columnt name "perc_overlap" in trinity_CCLE_gencode_blast dataframe

perc_overlap_results = (trinity_CCLE_gencode_blast['query_len'] / trinity_CCLE_gencode_blast['subj_len']) * 100
perc_overlap_results = perc_overlap_results.round(3)

trinity_CCLE_gencode_blast.insert(loc = 5, column = 'perc_overlap', value = perc_overlap_results)

#%%
# grouping together values in column "source_gene_id" by duplicates in descending order from higher to lower value in 'perc_overlap' column

trinity_CCLE_gencode_blast_sorted = trinity_CCLE_gencode_blast.sort_values(by = ['source_gene_id', 'perc_overlap'], ascending = [True, False])

max_perc_identical = trinity_CCLE_gencode_blast_sorted.groupby('source_gene_id')['perc_identical'].transform('max')

mask = trinity_CCLE_gencode_blast_sorted['perc_identical'] == max_perc_identical

trinity_CCLE_gencode_blast_sorted_mask = trinity_CCLE_gencode_blast_sorted[mask]

trinity_CCLE_gencode_blast_sorted_final = trinity_CCLE_gencode_blast_sorted_mask.drop_duplicates(subset = 'source_gene_id', keep ='first')


# %%
# creating Biotype column for our full-data (without creating subsets by % of identity) trinity CCLE gencode sorted file
# a loop with lambda function to get a biotype stored in unstructured strings in a column 'target_gene_id'


biotype_df_advanced_blast = trinity_CCLE_gencode_blast_sorted_final['target_gene_id'].apply(lambda x: x.split('|')[-2])

#%%
# a loop to add an biotype output to each dataframe   


trinity_CCLE_gencode_blast_sorted_final.insert(loc = 2, column = 'Biotype', value = biotype_df_advanced_blast)


# %%
# collapsing different biotype subtypes into more general categories: protein_coding, pseudogene, lncRNA, decay and other
# part 1: replacing protein_coding
replacing_protein_coding = ['protein_coding_CDS_not_defined']
new_value_protein_coding = ['protein_coding']
trinity_CCLE_gencode_blast_sorted_final['Biotype'] = trinity_CCLE_gencode_blast_sorted_final['Biotype'].replace(replacing_protein_coding, new_value_protein_coding)

# %%
# part 2: replacing pseudogenes
    
replacing_pseudogenes = ['processed_pseudogene', 'transcribed_processed_pseudogene', 'transcribed_unprocessed_pseudogene',
                         'transcribed_unitary_pseudogene', 'unprocessed_pseudogene', 'unitary_pseudogene']
new_value_psedugenes = ['pseudogene', 'pseudogene', 'pseudogene', 'pseudogene', 'pseudogene', 'pseudogene']

trinity_CCLE_gencode_blast_sorted_final['Biotype'] = trinity_CCLE_gencode_blast_sorted_final['Biotype'].replace(replacing_pseudogenes, new_value_psedugenes)


# %%
# part 3: replacing decay
    
replacing_decay = ['nonsense_mediated_decay', 'non_stop_decay']
new_value_decay = ['decay', 'decay']
trinity_CCLE_gencode_blast_sorted_final['Biotype'] = trinity_CCLE_gencode_blast_sorted_final['Biotype'].replace(replacing_decay, new_value_decay)

# %%
# part 3: replacing other
    
replacing_other = ['retained_intron', 'TEC', 'misc_RNA', 'artifact']
new_value_other = ['other', 'other', 'other', 'other']
trinity_CCLE_gencode_blast_sorted_final['Biotype'] = trinity_CCLE_gencode_blast_sorted_final['Biotype'].replace(replacing_other, new_value_other)

# now we have final list of dataframes 'tr_ccle_unique_with_biotypes' which has 5 types of our biotypes: protein_coding, pseudogene, lncRNA, decay and other!


#%%
# merging complex_dataframe_CCLE with trinity_CCLE_gencode_blast_sorted_final

complex_dataframe_CCLE = pd.merge(complex_dataframe_CCLE,trinity_CCLE_gencode_blast_sorted_final, left_on = 'ID', right_on = 'source_gene_id', how = 'left')

# %%
# SECOND PART OF WHOLE ANALYSIS: TRINITY CCLE vs mTEC BLAST ANALYSIS
# import complex_dataframe_mTEC_CCLE

complex_dataframe_mTEC_CCLE = pd.read_excel('complex_dataframe_mTEC_CCLE.xlsx')

# %%
# find IDs in complex_dataframe_mTEC_CCLE that are not matching with IDs in complex_dataframe_CCLE

transcripts_not_in_trinity_CCLE = ~complex_dataframe_mTEC_CCLE['ID'].isin(complex_dataframe_CCLE['ID'])
transcripts_not_in_trinity_CCLE_unique = complex_dataframe_mTEC_CCLE[transcripts_not_in_trinity_CCLE]
# %%
# inner join by ID between complex_dataframe_CCLE and complex_dataframe_mTEC_CCLE

complex_dataframe_CCLE_with_mTEC_CCLE_inner = pd.merge(complex_dataframe_CCLE, complex_dataframe_mTEC_CCLE, on= 'ID', how ='inner')
#%%

complex_dataframe_mTEC_CCLE_ID_renamed = complex_dataframe_mTEC_CCLE.rename(columns = {'ID': 'mTEC_ID'})
#%%
complex_dataframe_CCLE_with_mTEC_CCLE_left = pd.merge(complex_dataframe_CCLE, complex_dataframe_mTEC_CCLE_ID_renamed, left_on= 'ID', right_on = 'mTEC_ID', how ='left')


# %%
# find IDs in complex_dataframe_CCLE that are not matching with IDs in complex_dataframe_mTEC_CCLE

transcripts_not_in_mTEC_CCLE = ~complex_dataframe_CCLE['ID'].isin(complex_dataframe_mTEC_CCLE['ID'])
transcripts_not_in_mTEC_CCLE_unique = complex_dataframe_CCLE[transcripts_not_in_mTEC_CCLE]

# %%
# find IDs in transcripts_not_in_mTEC_CCLE_unique that are not matchin with IDs from GENCODE

trinity_CCLE_cancer_specific_non_gencode_mask = complex_dataframe_CCLE_with_mTEC_CCLE_left['target_gene_id_x'].isna() & complex_dataframe_CCLE_with_mTEC_CCLE_left['mTEC_ID'].isna()
trinity_CCLE_cancer_specific_non_gencode = complex_dataframe_CCLE_with_mTEC_CCLE_left[trinity_CCLE_cancer_specific_non_gencode_mask]
# %%
# create an excel file with complex_dataframe_CCLE_with_mTEC_CCLE_left for further analysis/merging (RNAseq, RiboSeq, translation into AA, immunopeptidomics, plotting)

complex_dataframe_CCLE_with_mTEC_CCLE_left.to_excel(r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\BLAST_OUTPUT\complex_dataframe_CCLE_with_mTEC_CCLE_left.xlsx', index = False)

# %%
