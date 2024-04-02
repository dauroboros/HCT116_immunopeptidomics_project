
# %%
# filtering HCT116_de_novo_assembly_FASTA by merged_df_non_canon IDs to get only non-canonical cancer-specific sequences

min_len = 8


import pandas as pd
import matplotlib as mpl
import numpy as np
import os
import seaborn as sns
import Bio

#change working directory

os.chdir(r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\BLAST_OUTPUT')

#%%
# create a dataframe from RNAseq database

df_full = pd.read_excel(r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\BLAST_OUTPUT\complex_dataframe_CCLE_with_mTEC_CCLE_left.xlsx')

# %%
#create a dataframe from kallisto quant output file

kallisto_df = pd.read_csv(r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\BLAST_OUTPUT\abundance.tsv', sep = '\t')

# %%
# left merge df_full and kallisto_df to see overlaps between RNAseq and RiboSeq data

merged_df = pd.merge(df_full, kallisto_df, left_on = 'ID', right_on = 'target_id', how = 'left')

# %%
# create a dataframe with only non-canonical cancer-specific IDs that have >0 counts on RiboSeq level
merged_df_non_canon =  (
    merged_df
    .query("target_gene_id_x != target_gene_id_x")
    .query("target_gene_id_y != target_gene_id_y")
    .query("est_counts > 0")

)

# %%
# filtering HCT116_de_novo_assembly_FASTA by merged_df_non_canon IDs to get only non-canonical cancer-specific sequences

from Bio import SeqIO

transcripts_ids = set(merged_df_non_canon['ID'])


def filter_fasta(fasta_file, output_file, ids):
    with open(output_file, 'w') as filtered_fasta:
        for record in SeqIO.parse(fasta_file,'fasta'):
            record_id = record.id.split('|')[0]
            if record_id in ids:
                SeqIO.write(record,filtered_fasta,'fasta')

original_fasta = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\trinity_CCLE_transcripts.fasta')
filtered_fasta = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\HCT116_trinity_ccle_non_canonical_kallisto_filtered.fasta')

filter_fasta(original_fasta, filtered_fasta, transcripts_ids)


#%%
# translate NA into 3 AA sequence frames from filtered FASTA file

from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
 
def translate_in_three_frames(seq):
    """
    Translates a nucleotide sequence in three reading frames.
    Returns a list of SeqRecord objects for the translations.
    """
    translations = []
    for frame in range(3):
        translated_seq = seq[frame:].translate()
        translations.append(translated_seq)
    return translations
 
na_seq_fasta = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\HCT116_trinity_ccle_non_canonical_kallisto_filtered.fasta')
amino_acids_fasta = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\HCT116_trinity_ccle_amino_acids.fasta')
 
translated_records = []
 
for record in SeqIO.parse(na_seq_fasta, "fasta"):
    translations = translate_in_three_frames(record.seq)
    for index, translation in enumerate(translations, start=1):
        new_id = f"{record.id}_frame{index}"
        translated_record = SeqRecord(translation, id=new_id, description=f"f{index}")
        translated_records.append(translated_record)
 
SeqIO.write(translated_records, amino_acids_fasta, "fasta")


#%%
# split FASTA file sequences by * symbol


def split_sequences_in_fasta(fasta_file_path, output_file_path):
    
    with open(fasta_file_path, 'r') as fasta_file, open(output_file_path, 'w') as output_file:
        sequence = ''
        header = ''
        for line in fasta_file:
            line = line.strip()  
            if line.startswith('>'):  
                if sequence:  
                    split_and_process_sequence(sequence, header, output_file)
                    sequence = ''  
                
                header = line
                
            else:
                sequence += line 
 
        
        if sequence:
            split_and_process_sequence(sequence, header, output_file)
 
def split_and_process_sequence(sequence, header, output_file):
    
    parts = sequence.split('*')
    
    part_number = 1
    for part in parts:
        if len(part) >= min_len:
            new_header = f"{header}_orf_{part_number}"
            output_file.write(new_header + "\n")
            output_file.write(part + "\n")  
            part_number += 1
 

amino_acids_fasta_split = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\HCT116_trinity_ccle_amino_acids.fasta')
output_fasta_file_path = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\HCT116_RNAseq_split_by_star_kallisto.fasta')
split_sequences_in_fasta(amino_acids_fasta_split, output_fasta_file_path)


# %%
# create a dataframe to analyse our final FASTA

fasta_final = (r'F:\Daur_Meretukov\Dark_Proteome_Data\TRINITY_AND_SPADES_STATISTICS\Raw_files\HCT116_RNAseq_split_by_star_kallisto.fasta')

 
data = []
for record in SeqIO.parse(fasta_final, "fasta"):
    id = record.id
    sequence = str(record.seq)
    data.append([id, sequence])
 
fasta_final_df = pd.DataFrame(data, columns=['ID', 'Sequence'])
 
print(fasta_final_df.head())



