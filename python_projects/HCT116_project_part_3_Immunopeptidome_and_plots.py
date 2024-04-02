#%%
# Part 3 of Immunopeptidomics: converging all dataframes into final one and filtering by final major filters to find best non-canonical peptide candidate
# importing libraries

import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import os
import Bio
import seaborn as sns
import upsetplot

#%%

os.chdir(r'/Users/d/Desktop/Life/UK/Imperial/MRes/Assignments/Project 1/Thesis/data_for_plots')

#%%
# import data to work with: RNAseq-RiboSeq data, Immunopeptidomics (with MHC-I bonding) data, RNAseq Kallisto Output

rna_ribo_seq_df = pd.read_excel(r'/Users/d/Desktop/Life/UK/Imperial/MRes/Assignments/Project 1/Thesis/data_for_plots/complex_dataframe_CCLE_with_mTEC_CCLE_left.xlsx')

rna_seq_kallisto_q_df = pd.read_table(r'/Users/d/Desktop/Life/UK/Imperial/MRes/Assignments/Project 1/Thesis/data_for_plots/HCT116_TRINITY_RNAseq_QUANT.tsv', sep = '\t')

immunopeptides_df = pd.read_table(r'/Users/d/Desktop/Life/UK/Imperial/MRes/Assignments/Project 1/Thesis/data_for_plots/IAT_2_TRINITY_HCT116.txt', sep = '\t')

#%%
# import additional data: RiboSeq Kallisto Output

ribo_seq_kallisto_q_df = pd.read_csv(r'/Users/d/Desktop/Life/UK/Imperial/MRes/Assignments/Project 1/Thesis/data_for_plots/PROTEOFORMER_RIBOSEQ_KALLISTO_QUANT_TRINITY_INDEX.tsv', sep = '\t')


# %%
# left merge of RNAseq-RiboSeq dataframe with RNAseq counts

rna_ribo_seq_counts_merged_df = pd.merge(rna_ribo_seq_df,rna_seq_kallisto_q_df, left_on = 'ID', right_on = 'target_id', how = 'left' )
# %%
# # left merge of RNAseq-RiboSeq-RNAseq counts dataframe with RiboSeq counts

rna_ribo_seq_counts_merged_full_df = pd.merge(rna_ribo_seq_counts_merged_df, ribo_seq_kallisto_q_df, left_on = 'ID', right_on = 'target_id', how = 'left' )

# %%
# transform immunopeptidomics file: take each rows with multiple TRINITY ID and create separate for for each ID and each frame to prerapte for merging

immunopeptides_df['Protein_Accessions'] = immunopeptides_df['Protein_Accessions'].apply(lambda x: x.split(','))

immunopeptides_df_explode = immunopeptides_df.explode('Protein_Accessions')
# %%
# remove '_frame' from each Protein Accession value 

immunopeptides_df_explode['Protein_Accessions'] = immunopeptides_df_explode['Protein_Accessions'].str.replace('_frame\d+', '', regex=True)

#%%
# left merge by TRINIY ID with the previous dataframe to create the final dataframe

final_dataframe = pd.merge(rna_ribo_seq_counts_merged_full_df, immunopeptides_df_explode, left_on = 'ID', right_on = 'Protein_Accessions', how = 'left' )

#%%

final_dataframe = final_dataframe.rename(columns = {'target_id_x' : 'target_id_q_RNAseq', 'length_x' : 'length_q_RNAseq', 'eff_length_x' : 'eff_length_q_RNAseq', 'est_counts_x' : 'est_counts_q_RNAseq', 'tpm_x' : 'tpm_q_RNAseq'})

final_dataframe = final_dataframe.rename(columns = {'target_id_y' : 'target_id_q_RiboSeq', 'length_y' : 'length_q_RiboSeq', 'eff_length_y' : 'eff_length_q_RiboSeq', 'est_counts_y' : 'est_counts_q_RiboSeq', 'tpm_y' : 'tpm_q_q_RiboSeq'})



# %%
# plot the distribution of RiboSeq value to select the threshold

sns.kdeplot(final_dataframe['est_counts_q_RiboSeq'], shade=True)
plt.title('Density Plot of Column a')
plt.xlabel('Value')
plt.ylabel('Density')
plt.show()
# %%
# filtering with the following columns and values: PSMs_High >= 2, RNA-Seq reads >= 5, Ribo-Seq reads >0, BestPSM_PEP <= 0.05, NetMHC Binder

condition_1 = final_dataframe['PSMs_high'] >= 2
condition_2 = final_dataframe['BestPSM_PEP'] <= 0.01
condition_3 = final_dataframe['est_counts_q_RNAseq'] >= 5
condition_4 = final_dataframe['est_counts_q_RiboSeq'] >= 1
condition_5 = final_dataframe['BestBinder'].isin(['WB','SB'])

combined_conditions = condition_1 & condition_2 & condition_3 & condition_4 & condition_5

filtered_final_dataframe = final_dataframe[combined_conditions]

#filtered_final_no_dups = filtered_final_dataframe.drop_duplicates(subset = 'ID', keep = 'first')
# %%
# plot to see distribution of reference and non-canonical

peptides_count_full = filtered_final_dataframe['Peptide_Type'].value_counts()

ax_peptides = peptides_count_full.plot(kind = 'bar')

plt.title('Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)


for p in ax_peptides.patches:
    ax_peptides.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

# %%
# plot distribution of reference and non-canonical without 1 of filters
# filter 1: NetMHC Binder = WB+SB filter removed


combined_conditions_binder = condition_1 & condition_2 & condition_3 & condition_4

filtered_final_dataframe_binder = final_dataframe[combined_conditions_binder]

filtered_final_no_dups = filtered_final_dataframe.drop_duplicates(subset = 'ID', keep = 'first')

peptides_count = filtered_final_no_dups['Peptide_Type'].value_counts()

ax_peptides = peptides_count.plot(kind = 'bar')

plt.title('Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)


for p in ax_peptides.patches:
    ax_peptides.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

# %%
# plot distribution of reference and non-canonical without 1 of filters
# filter 2: PSMs_High >=2 removed


combined_conditions_PSMs_high = condition_2 & condition_3 & condition_4 & condition_5

filtered_final_dataframe_psm = final_dataframe[combined_conditions_PSMs_high]

filtered_final_no_dups = filtered_final_dataframe.drop_duplicates(subset = 'ID', keep = 'first')

peptides_count = filtered_final_no_dups['Peptide_Type'].value_counts()

ax_peptides = peptides_count.plot(kind = 'bar')

plt.title('Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)


for p in ax_peptides.patches:
    ax_peptides.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

# %%
# plot distribution of reference and non-canonical without 1 of filters
# filter 3: BestPSM_PEP <= 0.05 removed


combined_conditions_bestpsm = condition_1 & condition_3 & condition_4 & condition_5

filtered_final_dataframe_bestpsm = final_dataframe[combined_conditions_bestpsm]

filtered_final_no_dups = filtered_final_dataframe.drop_duplicates(subset = 'ID', keep = 'first')

peptides_count = filtered_final_no_dups['Peptide_Type'].value_counts()

ax_peptides = peptides_count.plot(kind = 'bar')

plt.title('Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)


for p in ax_peptides.patches:
    ax_peptides.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

# %%
# plot distribution of reference and non-canonical without 1 of filters
# filter 4: RNA-Seq reads >= 5 removed


combined_conditions_rnaseq = condition_1 & condition_2 & condition_4 & condition_5

filtered_final_dataframe_rnaseq = final_dataframe[combined_conditions_rnaseq]

filtered_final_no_dups = filtered_final_dataframe.drop_duplicates(subset = 'ID', keep = 'first')

peptides_count = filtered_final_no_dups['Peptide_Type'].value_counts()

ax_peptides = peptides_count.plot(kind = 'bar')

plt.title('Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)


for p in ax_peptides.patches:
    ax_peptides.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()
# %%
# plot distribution of reference and non-canonical without 1 of filters
# filter 5: Ribo-Seq reads > 0


combined_conditions_riboseq = condition_1 & condition_2 & condition_3 & condition_5

filtered_final_dataframe_riboseq = final_dataframe[combined_conditions_riboseq]

filtered_final_no_dups = filtered_final_dataframe.drop_duplicates(subset = 'ID', keep = 'first')

peptides_count = filtered_final_no_dups['Peptide_Type'].value_counts()

ax_peptides = peptides_count.plot(kind = 'bar')

plt.title('Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)


for p in ax_peptides.patches:
    ax_peptides.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

# %%
# plot absolute numbers of non-canonical distribution with all filters on or 4/5 on to compare the effect of each filter on protein output



filtered_final_dataframe_nc = filtered_final_dataframe[filtered_final_dataframe['Peptide_Type'].isin(['UNKNOWN'])]

filtered_final_dataframe_binder_nc = filtered_final_dataframe_binder[filtered_final_dataframe_binder['Peptide_Type'].isin(['UNKNOWN'])]

filtered_final_dataframe_psm_nc = filtered_final_dataframe_psm[filtered_final_dataframe_psm['Peptide_Type'].isin(['UNKNOWN'])]

filtered_final_dataframe_bestpsm_nc = filtered_final_dataframe_bestpsm[filtered_final_dataframe_bestpsm['Peptide_Type'].isin(['UNKNOWN'])]

filtered_final_dataframe_rnaseq_nc = filtered_final_dataframe_rnaseq[filtered_final_dataframe_rnaseq['Peptide_Type'].isin(['UNKNOWN'])]

filtered_final_dataframe_riboseq_nc = filtered_final_dataframe_riboseq[filtered_final_dataframe_riboseq['Peptide_Type'].isin(['UNKNOWN'])]

counts_nc = [len(filtered_final_dataframe_nc), len(filtered_final_dataframe_binder_nc), len(filtered_final_dataframe_psm_nc), len(filtered_final_dataframe_bestpsm_nc), len(filtered_final_dataframe_rnaseq_nc), len(filtered_final_dataframe_riboseq_nc)]

df_labels_nc = ['All filters on', 'No BestBinder', 'No PSMs High', 'No BestPSM', 'No RNAseq counts', 'No RiboSeq counts']

# remove bbox_inches: 'tight' to prevent size error in plots with custom plt.text function

#%config InlineBackend.print_figure_kwargs = {'bbox_inches':None}


plt.figure(figsize=(12,8))
y_pos_nc = range(len(df_labels_nc))

bars_nc = plt.bar(y_pos_nc,counts_nc,align = 'center', alpha = 0.5, width = 0.5)

plt.xticks(y_pos_nc, df_labels_nc)
plt.ylabel('TRANSCRIPTS COUNT')
plt.title('Non-canonical proteins counts in relation with filtering on RNAseq, RiboSeq and immunopeptidomic levels')

for bar in bars_nc:

    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, yval, ha = 'center', va = 'bottom')


#%%
# plot % of non-canonical distribution with all filters on or 4/5 on to compare the effect of each filter on protein output

total_count_nc = sum(counts_nc)

percentages_nc = [(count / total_count_nc) * 100 for count in counts_nc]


plt.figure(figsize=(12,8))
y_pos_nc = range(len(df_labels_nc))

bars_nc = plt.bar(y_pos_nc,percentages_nc,align = 'center', alpha = 0.5, width = 0.5)

plt.xticks(y_pos_nc, df_labels_nc)
plt.ylabel('TRANSCRIPTS COUNT PERCENTAGE')
plt.title(r'Non-canonical proteins % in relation with filtering on RNAseq, RiboSeq and immunopeptidomic levels')

for bar, pct in zip (bars_nc, percentages_nc):
    plt.text(bar.get_x() + bar.get_width()/2, pct + 0.3, f'{pct:.2f}%', ha = 'center', va = 'bottom')

#%%
# plot reference proteins distribution with all filters on or 4/5 on to compare the effect of each filter on protein output

filtered_final_dataframe_ref = filtered_final_dataframe[filtered_final_dataframe['Peptide_Type'].isin(['REF'])]

filtered_final_dataframe_binder_ref = filtered_final_dataframe_binder[filtered_final_dataframe_binder['Peptide_Type'].isin(['REF'])]

filtered_final_dataframe_psm_ref = filtered_final_dataframe_psm[filtered_final_dataframe_psm['Peptide_Type'].isin(['REF'])]

filtered_final_dataframe_bestpsm_ref = filtered_final_dataframe_bestpsm[filtered_final_dataframe_bestpsm['Peptide_Type'].isin(['REF'])]

filtered_final_dataframe_rnaseq_ref = filtered_final_dataframe_rnaseq[filtered_final_dataframe_rnaseq['Peptide_Type'].isin(['REF'])]

filtered_final_dataframe_riboseq_ref = filtered_final_dataframe_riboseq[filtered_final_dataframe_riboseq['Peptide_Type'].isin(['REF'])]

counts_ref = [len(filtered_final_dataframe_ref), len(filtered_final_dataframe_binder_ref), len(filtered_final_dataframe_psm_ref), len(filtered_final_dataframe_bestpsm_ref), len(filtered_final_dataframe_rnaseq_ref), len(filtered_final_dataframe_riboseq_ref)]

df_labels_ref = ['All filters on', 'No BestBinder', 'No PSMs High', 'No BestPSM', 'No RNAseq counts', 'No RiboSeq counts']


plt.figure(figsize=(12,8))
y_pos_ref = range(len(df_labels_ref))

bars_nc = plt.bar(y_pos_ref,counts_ref,align = 'center', alpha = 0.5, width = 0.5)

plt.xticks(y_pos_ref, df_labels_ref)
plt.ylabel('TRANSCRIPTS COUNT')
plt.title('Reference proteins counts in relation with filtering on RNAseq, RiboSeq and immunopeptidomic levels')

for bar in bars_nc:

    yval = bar.get_height()
    plt.text(bar.get_x() + bar.get_width()/2, yval + 0.5, yval, ha = 'center', va = 'bottom')


#%%
# plot % of reference distribution with all filters on or 4/5 on to compare the effect of each filter on protein output

total_count_ref = sum(counts_ref)

percentages_ref = [(count / total_count_ref) * 100 for count in counts_ref]


plt.figure(figsize=(12,8))
y_pos_ref = range(len(df_labels_ref))

bars_ref = plt.bar(y_pos_ref,percentages_ref,align = 'center', alpha = 0.5, width = 0.5)

plt.xticks(y_pos_ref, df_labels_ref)
plt.ylabel('TRANSCRIPTS COUNT PERCENTAGE')
plt.title(r'Reference proteins % in relation with filtering on RNAseq, RiboSeq and immunopeptidomic levels')

for bar, pct in zip (bars_ref, percentages_ref):
    plt.text(bar.get_x() + bar.get_width()/2, pct + 0.3, f'{pct:.2f}%', ha = 'center', va = 'bottom')


#%%
#remove duplicated by peptides and specific rules from final_dataframe
    
final_dataframe_upset = final_dataframe.copy()

#final_dataframe_upset['BestBinder'] = final_dataframe_upset['BestBinder'].fillna("-")

#def binder_rule(final_dataframe_upset):

    #return final_dataframe_upset[final_dataframe_upset['BestBinder'] != '-']

#final_dataframe_upset = final_dataframe_upset[final_dataframe_upset['Peptide_Type'].isin(['UNKNOWN'])]
colfil = {
    'BestPSM_PEP': min,
    'PSMs_high': max,
    'est_counts_q_RNAseq': max,
    'est_counts_q_RiboSeq': max
}

grouped = final_dataframe_upset.groupby('Peptide').agg(colfil).reset_index()

final_dataframe_upset_final = pd.merge(grouped, final_dataframe_upset, on = ['BestPSM_PEP', 'PSMs_high', 'est_counts_q_RNAseq', 'est_counts_q_RiboSeq'])

final_dataframe_upset_final.rename(columns = {'Peptide_x': 'Peptide'}, inplace = True)

final_dataframe_upset_final_drop = final_dataframe_upset_final.drop_duplicates(subset = 'Peptide')
# %%
# create an upset plot for data without filtering by reproducibility (5 or more databases sources for peptide)

#Function:

def apply_tests(data, report_all=False):
    
    data_subset = (data[["Peptide", "BestPSM_PEP", "PSMs_high",
                        "BestBinder", "est_counts_q_RNAseq", "est_counts_q_RiboSeq"]]
                    .dropna()           
                    .replace(np.nan, "-")

                    )

    data_tests = (pd
                .DataFrame({'Peptide': data_subset.Peptide,
                            'PEP <= 0.01' : (data_subset['BestPSM_PEP'] <= 0.01 ).astype('int'),
                            'PSM_high >= 2': (data_subset['PSMs_high'] >= 2 ).astype('int'),
                            'NetMHC Binder': (data_subset['BestBinder'].fillna("-") != "-" ).astype('int'),
                            'RNA Reads >= 5': ( data_subset['est_counts_q_RNAseq'].replace("-", 0).astype('float32') >= 5.0 ).astype('int'),
                            'RiboSeq Reads >= 1': ( data_subset['est_counts_q_RiboSeq'].replace("-", 0).astype('float32') >= 1 ).astype('int'),
                            })
                            .set_index("Peptide")
                            .assign(colsum = lambda df: df.apply(sum, axis = 1))
                    )

    if report_all:
        return data_tests
   
    ntests = data_tests.shape[1] - 2
    data_tests = data_tests.query("colsum >= @ntests").drop("colsum", axis=1)
    return data_tests

import upsetplot as upset
upset_df = apply_tests(final_dataframe_upset_final_drop.query("Peptide_Type == 'UNKNOWN'"))

fig = plt.figure(figsize=(5, 5))
upset.plot(upset_df.groupby(upset_df.columns.tolist()).size(), show_counts = "%d")

#%%

reproducible_peptides_count = final_dataframe['File'].str.split(',').str.len()

reproducible_peptides_count_positive = final_dataframe[reproducible_peptides_count >= 5]

reproducible_peptides_count_positive_drop = reproducible_peptides_count_positive.drop_duplicates(subset = 'Peptide')

# %%
# create an upset plot for data

#Function:

def apply_tests(data, report_all=False):
    
    data_subset = (data[["Peptide", "BestPSM_PEP", "PSMs_high",
                        "BestBinder", "est_counts_q_RNAseq", "est_counts_q_RiboSeq"]]
                    .dropna()           
                    .replace(np.nan, "-")

                    )

    data_tests = (pd
                .DataFrame({'Peptide': data_subset.Peptide,
                            'PEP <= 0.01' : (data_subset['BestPSM_PEP'] <= 0.01 ).astype('int'),
                            'PSM_high >= 2': (data_subset['PSMs_high'] >= 2 ).astype('int'),
                            'NetMHC Binder': (data_subset['BestBinder'].fillna("-") != "-" ).astype('int'),
                            'RNA Reads >= 5': ( data_subset['est_counts_q_RNAseq'].replace("-", 0).astype('float32') >= 5.0 ).astype('int'),
                            'RiboSeq Reads >= 1': ( data_subset['est_counts_q_RiboSeq'].replace("-", 0).astype('float32') >= 1 ).astype('int'),
                            })
                            .set_index("Peptide")
                            .assign(colsum = lambda df: df.apply(sum, axis = 1))
                    )

    if report_all:
        return data_tests
   
    ntests = data_tests.shape[1] - 2
    data_tests = data_tests.query("colsum >= @ntests").drop("colsum", axis=1)
    return data_tests

import upsetplot as upset
upset_df = apply_tests(reproducible_peptides_count_positive_drop.query("Peptide_Type == 'UNKNOWN'"))

fig = plt.figure(figsize=(5, 5))
upset.plot(upset_df.groupby(upset_df.columns.tolist()).size(), show_counts = "%d")

# %%
# !!!!!! Part 1 of Plots Creation !!!!!!
# %%
#  Figure 3. Plot 1 (Table): characteristis of RNAseq assembly and blast table (HCT116, GENCODE, mTEC)
 
total_rows_tb1 = final_dataframe.shape[0]
 
a_b_not_nan = final_dataframe.dropna(subset=['ID', 'target_gene_id_x']).shape[0]
 
 
a_c_not_nan = final_dataframe.dropna(subset=['ID', 'target_gene_id_y']).shape[0]
 
abc_not_nan = final_dataframe.dropna(subset=['ID', 'target_gene_id_x', 'target_gene_id_y']).shape[0]
 
a_not_nan_bc_nan = final_dataframe[(final_dataframe['ID'].notna()) & (final_dataframe['target_gene_id_x'].isna()) & (final_dataframe['target_gene_id_y'].isna())].shape[0]
 
summary_table = pd.DataFrame({
    "Description": [
        "Total # of transcripts",
        "Matched with GENCODE",
        "Matched with mTEC",
        "Matched both with GENCODE and mTEC",
        "Not matched"
    ],
    "Count": [
        total_rows_tb1,
        a_b_not_nan,
        a_c_not_nan,
        abc_not_nan,
        a_not_nan_bc_nan
    ]
})
 
styled_table = summary_table.style.set_table_styles([
    {'selector': 'th', 'props': [('font-size', '12pt'), ('text-align', 'center')]},
    {'selector': 'td', 'props': [('text-align', 'center')]},
    {'selector': 'tr:nth-of-type(odd)', 'props': [('background', '#f2f2f2')]},
    {'selector': 'tr:nth-of-type(even)', 'props': [('background', 'white')]},
]).set_properties(**{'border': '1.5px solid black'})
 
styled_table.to_excel("summary_table.xlsx", engine='openpyxl', index=False)
 
 
print(summary_table)
 
# %%
#  Figure 3. Plot 2-4: Transcripts lengths distribution (GENCODE MATCHED, mTEC MATCHED, GENCODE & mTEC MATCHED) density plot
 
#%%
# create different dataframes for further plotting
# filter: to get cancer-specific non-canonical dataframes that do not have matches by ID in GENCODE and mTEC transcripts
 
trinity_CCLE_cancer_specific_non_gencode_mask = final_dataframe['target_gene_id_x'].isna() & final_dataframe['target_gene_id_y'].isna()
trinity_CCLE_cancer_specific_non_gencode = final_dataframe[trinity_CCLE_cancer_specific_non_gencode_mask]
 
#%%
#filter: to get only NOT NaN by GENCODE target gene ID
 
filtered_by_GENCODE_not_nan = final_dataframe[final_dataframe['target_gene_id_x'].notna()]
 
# %%
##filter: to get only NOT NaN by mTEC target gene ID
 
filtered_by_mTEC_not_nan = final_dataframe[final_dataframe['target_gene_id_y'].notna()]
#%%
# filter: to get only NOT NaN both by by GENCODE mTEC target gene ID and mTEC target gene ID
 
density_plot_cond_1 = final_dataframe['target_gene_id_x'].notna()
 
density_plot_cond_2 = final_dataframe['target_gene_id_y'].notna()
 
density_plot_cond_combo = density_plot_cond_1 & density_plot_cond_2
 
filtered_by_GENCODE_and_mTEC_not_nan = final_dataframe[density_plot_cond_combo]
 
# %%
# filter: to get 100%, >95% and <95% matches in GENCODE-positive mTEC positive subset
 
filtered_by_GENCODE_not_nan_100_perc_mTEC_true =  filtered_by_GENCODE_not_nan.loc[(filtered_by_GENCODE_not_nan['perc_identical_x'] == 100) & (filtered_by_GENCODE_not_nan['perc_identical_y'] >=95)] 
filtered_by_GENCODE_not_nan_more_95_perc_mTEC_true=  filtered_by_GENCODE_not_nan.loc[(filtered_by_GENCODE_not_nan['perc_identical_x'] >= 95) & (filtered_by_GENCODE_not_nan['perc_identical_x'] < 100) & (filtered_by_GENCODE_not_nan['perc_identical_y'] >=95)]
filtered_by_GENCODE_not_nan_less_95_perc_mTEC_true =  filtered_by_GENCODE_not_nan.loc[(filtered_by_GENCODE_not_nan['perc_identical_x'] < 95) & (filtered_by_GENCODE_not_nan['perc_identical_y'] >=95)]

filtered_by_GENCODE_not_nan_100_perc_mTEC_false =  filtered_by_GENCODE_not_nan.loc[(filtered_by_GENCODE_not_nan['perc_identical_x'] == 100) & ((filtered_by_GENCODE_not_nan['perc_identical_y'] < 95) | pd.isna(filtered_by_GENCODE_not_nan['perc_identical_y']))]
filtered_by_GENCODE_not_nan_more_95_perc_mTEC_false =  filtered_by_GENCODE_not_nan.loc[(filtered_by_GENCODE_not_nan['perc_identical_x'] >= 95) & (filtered_by_GENCODE_not_nan['perc_identical_x'] < 100) & ((filtered_by_GENCODE_not_nan['perc_identical_y'] < 95) | pd.isna(filtered_by_GENCODE_not_nan['perc_identical_y']))]
filtered_by_GENCODE_not_nan_less_95_perc_mTEC_false =  filtered_by_GENCODE_not_nan.loc[(filtered_by_GENCODE_not_nan['perc_identical_x'] < 95) & ((filtered_by_GENCODE_not_nan['perc_identical_y'] < 95) | pd.isna(filtered_by_GENCODE_not_nan['perc_identical_y']))]
# %%
# create lists with transcripts length: GENCODE positive mTEC positive
 
filtered_by_GENCODE_not_nan_100_perc_mTEC_true_length = filtered_by_GENCODE_not_nan_100_perc_mTEC_true['Transcript Length_x'].tolist()
filtered_by_GENCODE_not_nan_more_95_perc_mTEC_true_length = filtered_by_GENCODE_not_nan_more_95_perc_mTEC_true['Transcript Length_x'].tolist()
filtered_by_GENCODE_not_nan_less_95_perc_mTEC_true_length = filtered_by_GENCODE_not_nan_less_95_perc_mTEC_true['Transcript Length_x'].tolist()
trinity_CCLE_cancer_specific_non_gencode_length = trinity_CCLE_cancer_specific_non_gencode['Transcript Length_x'].tolist()

dp_gc_1 = len(filtered_by_GENCODE_not_nan_100_perc_mTEC_true_length)
dp_gc_2 = len(filtered_by_GENCODE_not_nan_more_95_perc_mTEC_true_length)
dp_gc_3 = len(filtered_by_GENCODE_not_nan_less_95_perc_mTEC_true_length)
dp_commmon = len(trinity_CCLE_cancer_specific_non_gencode_length)

#%%
# create lists with transcripts length: GENCODE positive mTEC negative
 
filtered_by_GENCODE_not_nan_100_perc_mTEC_false_length = filtered_by_GENCODE_not_nan_100_perc_mTEC_false['Transcript Length_x'].tolist()
filtered_by_GENCODE_not_nan_more_95_perc_mTEC_false_length = filtered_by_GENCODE_not_nan_more_95_perc_mTEC_false['Transcript Length_x'].tolist()
filtered_by_GENCODE_not_nan_less_95_perc_mTEC_false_length = filtered_by_GENCODE_not_nan_less_95_perc_mTEC_false['Transcript Length_x'].tolist()
 
dp_mtc_1 = len(filtered_by_GENCODE_not_nan_100_perc_mTEC_false_length)
dp_mtc_2 = len(filtered_by_GENCODE_not_nan_more_95_perc_mTEC_false_length)
dp_mtc_3 = len(filtered_by_GENCODE_not_nan_less_95_perc_mTEC_false_length)

#%%
# create density plot 1 for GENCODE positive mTEC positive subset
plt.figure(figsize=(15,8))
sns.kdeplot(filtered_by_GENCODE_not_nan_100_perc_mTEC_true_length, bw_adjust = 0.5, label = '100 percent identity transcripts', fill = True)
sns.kdeplot(filtered_by_GENCODE_not_nan_more_95_perc_mTEC_true_length, bw_adjust = 0.5, label = 'More than 95 percent identity transcripts', fill = True)
sns.kdeplot(filtered_by_GENCODE_not_nan_less_95_perc_mTEC_true_length, bw_adjust = 0.5, label = 'Less than 95 percent identity transcripts', fill = True)
sns.kdeplot(trinity_CCLE_cancer_specific_non_gencode_length, bw_adjust = 0.5, label = 'Cancer-Specific Non-GENCODE, Non-mTEC BLAST transcripts', fill = True)
 
plt.xscale('log')
plt.title('Density plots for HCT116 TRINITY CCLE TRANSCRIPTS LENGTH: GENCODE MATCHED BY % IDENTITY (mTEC positive) + CANCER-SPECIFIC')
plt.xlabel('Transcripts length')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper left')
plt.show()
 
 
#%%
# create density plot 2 for GENCODE positive mTEC negative subset
 
plt.figure(figsize=(15,8))
sns.kdeplot(filtered_by_GENCODE_not_nan_100_perc_mTEC_false_length, bw_adjust = 0.5, label = '100 percent identity transcripts', fill = True)
sns.kdeplot(filtered_by_GENCODE_not_nan_more_95_perc_mTEC_false_length, bw_adjust = 0.5, label = 'More than 95 percent identity transcripts', fill = True)
sns.kdeplot(filtered_by_GENCODE_not_nan_less_95_perc_mTEC_false_length, bw_adjust = 0.5, label = 'Less than 95 percent identity transcripts', fill = True)
sns.kdeplot(trinity_CCLE_cancer_specific_non_gencode_length, bw_adjust = 0.5, label = 'Cancer-Specific Non-GENCODE, Non-mTEC BLAST transcripts', fill = True)
 
plt.xscale('log')
plt.title('Density plots for HCT116 TRINITY CCLE TRANSCRIPTS LENGTH: GENCODE MATCHED BY % IDENTITY (mTEC negative)+ CANCER-SPECIFIC')
plt.xlabel('Transcripts length')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper left')
plt.show()
#%%
# Additional statistics for transcripts length distribution

stats_tr_lengths = {}
dataframes_tr_lengths = [filtered_by_GENCODE_not_nan_100_perc_mTEC_true, filtered_by_GENCODE_not_nan_more_95_perc_mTEC_true, filtered_by_GENCODE_not_nan_less_95_perc_mTEC_true,
                          filtered_by_GENCODE_not_nan_100_perc_mTEC_false, filtered_by_GENCODE_not_nan_more_95_perc_mTEC_false, filtered_by_GENCODE_not_nan_less_95_perc_mTEC_false, 
                          trinity_CCLE_cancer_specific_non_gencode]
name_tr_lengths = ['GENCODE+mTEC+ 100%', 'GENCODE+mTEC+ >= 95%', 'GENCODE+mTEC+ <95%', 'GENCODE+mTEC- 100%', 'GENCODE+mTEC- >= 95%', 'GENCODE+mTEC- <95%', 'NON-MATCHED CANCER-SPECIFIC']

for name, df in zip(name_tr_lengths, dataframes_tr_lengths):
    stats_tr_lengths[name] = {
        'Mean': df['Transcript Length_x'].mean(),
        'Median': df['Transcript Length_x'].median(),
        'Max': df['Transcript Length_x'].max(),
        'Min': df['Transcript Length_x'].min(),
        'Std Deviation': df['Transcript Length_x'].std()
    }

# Convert the stats dictionary to a DataFrame for easier manipulation and visualization
stats_tr_lengths_final = pd.DataFrame(stats_tr_lengths)

#%%
#  Figure 3. Plot 5-6: RiboSeq and RNAseq Counts Distribution: GENCODE matches (reference) vs non-canonical

 # create Density Plot for RNAseq Reads Distribution


gencode_positive_mtec_positive_rnaseq_density =  final_dataframe.loc[(final_dataframe['perc_identical_x'] >= 95) & (final_dataframe['perc_identical_y'] >=95)]
gencode_positive_mtec_negative_rnaseq_density =  final_dataframe.loc[(final_dataframe['perc_identical_x'] >= 95) & ((final_dataframe['perc_identical_y'] < 95) | pd.isna(final_dataframe['perc_identical_y']))]
gencode_negative_mtec_positive_rnaseq_density =  final_dataframe.loc[((final_dataframe['perc_identical_x'] < 95) | pd.isna(final_dataframe['perc_identical_x'])) & (final_dataframe['perc_identical_y'] >=95)]

gencode_positive_mtec_positive_rnaseq_density_counts = gencode_positive_mtec_positive_rnaseq_density['est_counts_q_RNAseq'].tolist()
gencode_positive_mtec_negative_rnaseq_density_counts = gencode_positive_mtec_negative_rnaseq_density['est_counts_q_RNAseq'].tolist()
gencode_negative_mtec_positive_rnaseq_density_counts = gencode_negative_mtec_positive_rnaseq_density['est_counts_q_RNAseq'].tolist()
trinity_CCLE_cancer_specific_non_gencode_counts = trinity_CCLE_cancer_specific_non_gencode['est_counts_q_RNAseq'].tolist()

plt.figure(figsize=(15,8))
sns.kdeplot(gencode_positive_mtec_positive_rnaseq_density_counts, bw_adjust = 0.5, label = 'GENCODE positive, mTEC positive', fill = True)
sns.kdeplot(gencode_positive_mtec_negative_rnaseq_density_counts, bw_adjust = 0.5, label = 'GENCODE positive, mTEC negative', fill = True)
sns.kdeplot(gencode_negative_mtec_positive_rnaseq_density_counts, bw_adjust = 0.5, label = 'GENCODE negative, mTEC positive', fill = True)
sns.kdeplot(trinity_CCLE_cancer_specific_non_gencode_counts, bw_adjust = 0.5, label = 'Cancer-Specific Non-GENCODE, Non-mTEC positive transcripts', fill = True)
 
plt.xscale('log')
plt.title('RNAseq read counts distribution between GENCODE matched (mTEC negative and mTEC positive) and NON-CANONICAL CANCER SPECIFIC')
plt.xlabel('RNAseq counts')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper right')
plt.text(12000,0.001, '''GENCODE POSITIVE mTEC POSITIVE TRANSCRIPTS N = 10222\n\nGENCODE POSITIVE mTEC NEGATIVE TRANSCRIPTS N = 28298\n\n>GENCODE NEGATIVE mTEC POSITIVE TRANSCRIPTS N = 64243\n\nNON-MATCHED TRANSCRIPTS N = {}'''.format(dp_commmon) ,fontsize = 10, bbox = dict(facecolor = 'none', edgecolor = 'grey', boxstyle='round,pad=1', alpha = 0.5))
plt.show()

#%%
# Statistics RNASeq/RiboSeq

stats_rnaseq = {}
dataframes_rnaseq = [gencode_positive_mtec_positive_rnaseq_density, gencode_positive_mtec_negative_rnaseq_density, gencode_negative_mtec_positive_rnaseq_density, trinity_CCLE_cancer_specific_non_gencode]
name_stats_rnaseq = ['GENCODE+mTEC+', 'GENCODE+mTEC-', 'GENCODE-mTEC+', 'NON-CANONICAL CANCER-SPECIFIC']

for name, df in zip(name_stats_rnaseq, dataframes_rnaseq):
    stats_rnaseq[name] = {
        'Mean': df['est_counts_q_RNAseq'].mean(),
        'Median': df['est_counts_q_RNAseq'].median(),
        'Max': df['est_counts_q_RNAseq'].max(),
        'Min': df['est_counts_q_RNAseq'].min(),
        'Std Deviation': df['est_counts_q_RNAseq'].std()
    }

stats_stats_rnaseq_final = pd.DataFrame(stats_rnaseq)

stats_riboseq = {}
dataframes_riboseq = [gencode_positive_mtec_positive_rnaseq_density, gencode_positive_mtec_negative_rnaseq_density, gencode_negative_mtec_positive_rnaseq_density, trinity_CCLE_cancer_specific_non_gencode]
name_stats_riboseq = ['GENCODE+mTEC+', 'GENCODE+mTEC-', 'GENCODE-mTEC+', 'NON-CANONICAL CANCER-SPECIFIC']

for name, df in zip(name_stats_riboseq, dataframes_riboseq):
    stats_riboseq[name] = {
        'Mean': df['est_counts_q_RiboSeq'].mean(),
        'Median': df['est_counts_q_RiboSeq'].median(),
        'Max': df['est_counts_q_RiboSeq'].max(),
        'Min': df['est_counts_q_RiboSeq'].min(),
        'Std Deviation': df['est_counts_q_RiboSeq'].std()
    }

stats_stats_riboseq_final = pd.DataFrame(stats_riboseq)

#%%
# create violin plot for RNAseq estimated counts distribution

data_rnaseq_violin = [gencode_positive_mtec_positive_rnaseq_density_counts, gencode_positive_mtec_negative_rnaseq_density_counts, gencode_negative_mtec_positive_rnaseq_density_counts, 
        trinity_CCLE_cancer_specific_non_gencode_counts]
labels = ['GENCODE+mTEC+', 'GENCODE+mTEC-', 'GENCODE-mTEC+', 'NON-CANONICAL CANCER-SPECIFIC']

sns.violinplot(data=data_rnaseq_violin)

plt.xticks(range(len(labels)), labels, rotation = 30)

plt.show()

# create box plots for RiboSeq Reads Distribution
 
riboseq_box_plots = [gencode_positive_mtec_positive_rnaseq_density, gencode_positive_mtec_negative_rnaseq_density, gencode_negative_mtec_positive_rnaseq_density, trinity_CCLE_cancer_specific_non_gencode]
riboseq_box_plots_names = ['GENCODE + mTEC +', 'GENCODE + mTEC -', 'GENCODE - mTEC +', 'NON-CANONICAL CANCER SPECIFIC']

# Calculate '>0' and '=0' counts
summary_data_riboseq_plots = {}
for name, df in zip(riboseq_box_plots_names, riboseq_box_plots):
    greater_than_zero = (df['est_counts_q_RiboSeq'] > 0).sum()
    equal_zero = (df['est_counts_q_RiboSeq'] == 0).sum()
    summary_data_riboseq_plots[name] = {'RiboSeq counts > 0': greater_than_zero, 'RiboSeq counts = 0': equal_zero}

summary_df_riboseq_plots = pd.DataFrame(summary_data_riboseq_plots).T

plt.figure(figsize=(10, 6))

bar_width = 0.2

indices = np.arange(len(summary_df_riboseq_plots))

for i, col in enumerate(summary_df_riboseq_plots.columns):
    plt.bar(indices + i * bar_width, summary_df_riboseq_plots[col], width=bar_width, label=col)

plt.yscale('log')

plt.xlabel('Transcripts by subset type')
plt.ylabel('RiboSeq counts')
plt.title('RiboSeq counts in HCT116 De Novo Assembly trancsripts sorted by GENCODE/mTEC matching subset type')
plt.xticks(indices + bar_width, summary_df_riboseq_plots.index, rotation=30)
plt.legend()

plt.show()

#%%
# Figure 3. Plot 7 Transcripts Biotype

biotypes_transcripts = final_dataframe.copy()
biotypes_transcripts = biotypes_transcripts.reset_index()

biotypes_old = ['protein_coding_LoF', 'TR_V_gene','processed_transcript', 'TR_V_pseudogene']

biotypes_corrected = ['protein_coding', 'other', 'other', 'pseudogene']

biotypes_transcripts['Biotype'] =  biotypes_transcripts['Biotype'].replace(biotypes_old,biotypes_corrected)


biotypes_transcripts_100_perc =  biotypes_transcripts.loc[biotypes_transcripts['perc_identical_x'] == 100] 
biotypes_transcripts_more_95_perc=  biotypes_transcripts.loc[(biotypes_transcripts['perc_identical_x'] >= 95) & (biotypes_transcripts['perc_identical_x'] < 100)]
biotypes_transcripts_less_95_perc =  biotypes_transcripts.loc[biotypes_transcripts['perc_identical_x'] < 95]

biotypes_transcripts_100_perc_counts = biotypes_transcripts_100_perc['Biotype'].value_counts()
biotypes_transcripts_more_95_perc_counts = biotypes_transcripts_more_95_perc['Biotype'].value_counts()
biotypes_transcripts_less_95_perc_counts = biotypes_transcripts_less_95_perc['Biotype'].value_counts()

all_biotypes = ['lncRNA', 'protein_coding', 'other', 'decay', 'pseudogene']


biotypes_transcripts_100_perc_counts = biotypes_transcripts_100_perc_counts.reindex(all_biotypes,fill_value=0)
biotypes_transcripts_more_95_perc_counts = biotypes_transcripts_more_95_perc_counts.reindex(all_biotypes,fill_value=0)
biotypes_transcripts_less_95_perc_counts = biotypes_transcripts_less_95_perc_counts.reindex(all_biotypes,fill_value=0)

combined_counts = pd.DataFrame({'100% identity': biotypes_transcripts_100_perc_counts, '>=95% identity': biotypes_transcripts_more_95_perc_counts, '<95% identity': biotypes_transcripts_less_95_perc_counts})

combined_counts_t = combined_counts.T

ax_biotypes = combined_counts_t.plot(kind='bar', stacked=True, figsize=(10, 6))

plt.xlabel('GENCODE matched, by percentage of identity')
plt.ylabel('Transcripts count')
plt.title('Biotypes distribution in GENCODE matched HCT116 De Novo Assembly Transcripts')
plt.xticks(rotation=0)
plt.legend(title = 'Transcript Biotype')

for bar in ax_biotypes.patches:
    width, height = bar.get_width(), bar.get_height()
    x, y = bar.get_xy() 
    if height > 0:
        ax_biotypes.annotate(f'{height:.0f}', (x + width/2, y + height/2), ha='center', va='center')


plt.show()

#%%
# Figure 4. Plot 1 Immunopeptides Source distribution by DataBase (total number of peptides)

database_source_count_split = final_dataframe['File'].str.split(',').dropna()
 

database_source_count_list = [item.strip() for sublist in database_source_count_split for item in sublist]
 

database_source_count_series = pd.Series(database_source_count_list)
 

database_source_count_final = database_source_count_series.value_counts()
 
database_source_count_final_sns = database_source_count_final.reset_index()
database_source_count_final_sns.columns = ['Source DB', 'Peptides Counts']

plt.figure(figsize = (10,6))
database_source_count_final_plot = sns.barplot(x='Source DB', y = 'Peptides Counts', data = database_source_count_final_sns, palette = 'crest')

for p in database_source_count_final_plot.patches:
    database_source_count_final_plot.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.title('Total number of peptides counts by source experimental datasets')
plt.xlabel('Source Datasets')
plt.ylabel('Peptides Count')


#%%
# Figure 4. Plots 2-3 Immunopeptides Source distribution by DataBase (reference-matched and non-canonical)

peptides_by_source = final_dataframe.copy()

peptides_by_source['File'] = peptides_by_source['File'].str.split(',')

peptides_by_source_exploded = peptides_by_source.explode('File')

peptides_by_source_exploded_nc = peptides_by_source_exploded[peptides_by_source_exploded['Peptide_Type'] == 'UNKNOWN']
peptides_by_source_exploded_ref =  peptides_by_source_exploded[peptides_by_source_exploded['Peptide_Type'] == 'REF']

peptides_by_source_exploded_nc_counts = peptides_by_source_exploded_nc['File'].value_counts().reset_index()
peptides_by_source_exploded_ref_counts = peptides_by_source_exploded_ref['File'].value_counts().reset_index()

peptides_by_source_exploded_nc_counts.columns = ['Source DB', 'Peptides Counts']
peptides_by_source_exploded_ref_counts.columns = ['Source DB', 'Peptides Counts']

plt.figure(figsize = (10,6))
peptides_by_source_exploded_nc_counts_plot = sns.barplot(x='Source DB', y = 'Peptides Counts', data = peptides_by_source_exploded_nc_counts, palette = 'crest')

for p in peptides_by_source_exploded_nc_counts_plot.patches:
    peptides_by_source_exploded_nc_counts_plot.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.title('Non-canonical peptides counts by source experimental datasets')
plt.xlabel('Source datasets')
plt.ylabel('Peptides Count')

plt.figure(figsize = (10,6))
peptides_by_source_exploded_ref_counts_plot = sns.barplot(x='Source DB', y = 'Peptides Counts', data = peptides_by_source_exploded_ref_counts, palette = 'crest')

for p in peptides_by_source_exploded_ref_counts_plot.patches:
    peptides_by_source_exploded_ref_counts_plot.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.title('Reference-matched peptides counts by source experimental datasets')
plt.xlabel('Source datasets')
plt.ylabel('Peptides Count')


#%%
# Figure 4. Plots 5-6  bar plot with reproducible and not reproducibel peptides (reference-matched and non-canonical)

reproducible_peptides_count = final_dataframe['File'].str.split(',').str.len()



reproducible_peptides_count_positive = final_dataframe[reproducible_peptides_count >= 5]

reproducible_peptides_count_positive['Peptide_Type'] = reproducible_peptides_count_positive['Peptide_Type'].replace({'REF': 'Reference-matched', 'UNKNOWN': 'Non-canonical'})


reproducible_peptides_count_negative = final_dataframe[reproducible_peptides_count < 5]

reproducible_peptides_count_negative['Peptide_Type'] = reproducible_peptides_count_negative['Peptide_Type'].replace({'REF': 'Reference-matched', 'UNKNOWN': 'Non-canonical'})


reproducible_peptides_count_positive_count = reproducible_peptides_count_positive['Peptide_Type'].value_counts()

reproducible_peptides_count_negative_count = reproducible_peptides_count_negative['Peptide_Type'].value_counts()

ax_repoducible_positive = reproducible_peptides_count_positive_count.plot(kind = 'bar')

plt.title('Reprocudible Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)

for p in ax_repoducible_positive.patches:
    ax_repoducible_positive.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

ax_repoducible_negative = reproducible_peptides_count_negative_count.plot(kind = 'bar')

plt.title('Non-Reprocudible Peptide Counts by type')
plt.xlabel('Peptide Type')
plt.xlabel('Counts')
plt.xticks(rotation=0)

for p in ax_repoducible_negative.patches:
    ax_repoducible_negative.annotate(f'{int(p.get_height())}', 
                (p.get_x() + p.get_width() / 2., p.get_height()), 
                ha = 'center', va = 'center', 
                xytext = (0, 5), 
                textcoords = 'offset points')

plt.show()

#%%
# Figure 4. Plot 8  PSM_high density plot for non-canonical and reference-matched peptides

filtered_peptides_nc =  final_dataframe.loc[final_dataframe['Peptide_Type'] == 'UNKNOWN']
filtered_peptides_ref =  final_dataframe.loc[final_dataframe['Peptide_Type'] == 'REF']


filtered_peptides_nc_psm =  filtered_peptides_nc.loc[filtered_peptides_nc['PSMs_high'] > 0]
filtered_peptides_ref_psm =  filtered_peptides_ref.loc[filtered_peptides_ref['PSMs_high'] > 0]

density_plot_psm_high_nc_length = filtered_peptides_nc_psm['PSMs_high'].tolist()
density_plot_psm_high_ref_length = filtered_peptides_ref_psm['PSMs_high'].tolist()

dp_psm_nc = len(density_plot_psm_high_nc_length)
dp_psm_ref = len(density_plot_psm_high_ref_length)

plt.figure(figsize=(15,8))
sns.kdeplot(density_plot_psm_high_nc_length, bw_adjust = 0.5, label = 'Non-canonical peptides', fill = True)
sns.kdeplot(density_plot_psm_high_ref_length, bw_adjust = 0.5, label = 'Reference-matched peptides', fill = True)

plt.xscale('log')
plt.title('Density plots for HCT116 Immunopeptides distribution by PSMs_High filter')
plt.xlabel('PSMs_high value')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper right')
plt.show()


#%%
# Figure 4. Plot 9  PEP density plot for non-canonical and reference-matched peptides


density_plot_pep_nc_length = filtered_peptides_nc['BestPSM_PEP'].tolist()
density_plot_pep_ref_length = filtered_peptides_ref['BestPSM_PEP'].tolist()



plt.figure(figsize=(15,8))
sns.kdeplot(density_plot_pep_nc_length, bw_adjust = 0.5, label = 'Non-canonical peptides', fill = True)
sns.kdeplot(density_plot_pep_ref_length, bw_adjust = 0.5, label = 'Reference-matched peptides', fill = True)

plt.xscale('log')
plt.title('Density plots for HCT116 Immunopeptides distribution by BestPSM_PEP filter')
plt.xlabel('BestPSM_PEP value')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper right')
plt.show()

#%%
# Figure 4. Plot 10  RNAseq density plot for non-canonical and reference-matched peptides

density_plot_rnaseq_nc_length = filtered_peptides_nc['est_counts_q_RNAseq'].tolist()
density_plot_rnaseq_ref_length = filtered_peptides_ref['est_counts_q_RNAseq'].tolist()



plt.figure(figsize=(15,8))
sns.kdeplot(density_plot_rnaseq_nc_length, bw_adjust = 0.5, label = 'Non-canonical peptides', fill = True)
sns.kdeplot(density_plot_rnaseq_ref_length, bw_adjust = 0.5, label = 'Reference-matched peptides', fill = True)

plt.xscale('log')
plt.title('Density plots for HCT116 Immunopeptides distribution by RNAseq read counts filter')
plt.xlabel('RNAseq read counts')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper right')
plt.show()


#%%
# Figure 4. Plot 11 RiboSeq density plot for non-canonical and reference-matched peptides


density_plot_riboseq_nc_length = filtered_peptides_nc['est_counts_q_RiboSeq'].tolist()
density_plot_riboseq_ref_length = filtered_peptides_ref['est_counts_q_RiboSeq'].tolist()



plt.figure(figsize=(15,8))
sns.kdeplot(density_plot_riboseq_nc_length, bw_adjust = 0.5, label = 'Non-canonical peptides', fill = True)
sns.kdeplot(density_plot_riboseq_ref_length, bw_adjust = 0.5, label = 'Reference-matched peptides', fill = True)

plt.xscale('log')
plt.title('Density plots for HCT116 Immunopeptides distribution by RiboSeq read counts filter')
plt.xlabel('RiboSeq read counts')
plt.ylabel('Density')
plt.legend(fontsize = 12, loc = 'upper right')
plt.show()


#%%
# Figure 4. Plot 12 MHCBinder type stacked bar plot for non-canonical and reference-matched peptides

mhc_binder_dropna_peptides = final_dataframe.copy()
mhc_binder_dropna_peptides = mhc_binder_dropna_peptides.reset_index()
mhc_binder_dropna_peptides['BestBinder'] =  mhc_binder_dropna_peptides['BestBinder'].fillna("-")
mhc_binder_dropna_peptides['BestBinder'] =  mhc_binder_dropna_peptides['BestBinder'].replace("-",'Non-Binder')

mhc_binder_stacked_plot_data = pd.crosstab(mhc_binder_dropna_peptides['Peptide_Type'], mhc_binder_dropna_peptides['BestBinder'])

mhc_binder_stacked_plot_data_color = sns.color_palette('deep')

mhc_binder_stacked_plot_data_var = mhc_binder_stacked_plot_data.plot(kind = 'bar', stacked = True, figsize = (10,6), color=mhc_binder_stacked_plot_data_color)

plt.title('MHC Binding features distribution in non-canonical and reference matched peptides')
plt.xlabel('Peptide Type')
plt.ylabel('Counts')
plt.legend(title = 'MHC Binding type')

for bars in mhc_binder_stacked_plot_data_var.containers:
    mhc_binder_stacked_plot_data_var.bar_label(bars, label_type='center')


plt.show()

#%%
# Figure 4. Plot 13 best-candidate peptides upset plot (non-canonical, filtering at least 4/5, reproducible)


def apply_tests(data, report_all=False):
    
    data_subset = (data[["Peptide", "BestPSM_PEP", "PSMs_high",
                        "BestBinder", "est_counts_q_RNAseq", "est_counts_q_RiboSeq"]]
                    .dropna()           
                    .replace(np.nan, "-")

                    )

    data_tests = (pd
                .DataFrame({'Peptide': data_subset.Peptide,
                            'PEP <= 0.01' : (data_subset['BestPSM_PEP'] <= 0.01 ).astype('int'),
                            'PSM_high >= 2': (data_subset['PSMs_high'] >= 2 ).astype('int'),
                            'NetMHC Binder': (data_subset['BestBinder'].fillna("-") != "-" ).astype('int'),
                            'RNA Reads >= 5': ( data_subset['est_counts_q_RNAseq'].replace("-", 0).astype('float32') >= 5.0 ).astype('int'),
                            'RiboSeq Reads >= 1': ( data_subset['est_counts_q_RiboSeq'].replace("-", 0).astype('float32') >= 1 ).astype('int'),
                            })
                            .set_index("Peptide")
                            .assign(colsum = lambda df: df.apply(sum, axis = 1))
                    )

    if report_all:
        return data_tests
   
    ntests = data_tests.shape[1] - 2
    data_tests = data_tests.query("colsum >= @ntests").drop("colsum", axis=1)
    return data_tests

import upsetplot as upset
upset_df = apply_tests(reproducible_peptides_count_positive_drop.query("Peptide_Type == 'UNKNOWN'"))

fig = plt.figure(figsize=(5, 5))
upset.plot(upset_df.groupby(upset_df.columns.tolist()).size(), show_counts = "%d")

#%%
# TRINITY and ICR LAB INTERNAL DATABASE Peptides comparison for Vien Diagramm
upset_df.to_csv('daur_best_candidates_peptides.csv')


peptides_database = pd.read_csv('ABX_RUN2_PEPTIDES_dbs_sum_wfile.txt', sep = '\t')

novel_peptides = pd.merge(upset_df,peptides_database, left_on = 'Peptide', right_on = 'Peptide', how = 'left')
novel_peptides.to_csv('daur_best_candidates_peptides_novel.csv')
#%%
# create a dataframe of all unique peptides with 5% or less FDR rate

peptides_total_comparison = final_dataframe.copy()

columns_to_keep = ['Peptide','Peptide_FDR', 'BestPSM_PEP', 'PSMs_high', 'Peptide_Type']

peptides_total_comparison = peptides_total_comparison[columns_to_keep]

peptides_total_comparison_no_nan = peptides_total_comparison.dropna(subset = 'Peptide')

peptides_total_comparison_no_nan_fdr = peptides_total_comparison_no_nan.loc[peptides_total_comparison_no_nan['Peptide_FDR'] < 0.05]


peptides_total_comparison_no_nan_fdr.to_csv('daur_all_peptides_low_fdr.csv')

#%%
# just in case - drop duplicates

peptides_database_unique = peptides_database.drop_duplicates(subset = 'Peptide')
peptides_database_unique_fdr = peptides_database_unique.loc[peptides_database_unique['Peptide_FDR'] < 0.05]

peptides_total_comparison_no_nan_fdr_unique = peptides_total_comparison_no_nan_fdr.drop_duplicates(subset='Peptide')
#%%
# create Venn diagramm between TRINITY-based DB and ICR Peptides DB

import matplotlib_venn
from matplotlib_venn import venn2
 
trinity_peptides_set = set(peptides_total_comparison_no_nan_fdr_unique['Peptide'])
icr_peptides_set = set(peptides_database_unique_fdr['Peptide'])


venn2([trinity_peptides_set, icr_peptides_set], ('TRINITY-based immunopeptides database', 'ICR database'))

plt.show()

#%%
# # create Venn diagramm between TRINITY-based DB and ICR Peptides DB (best candidates by 2 filters: PSMs_high >= 2 and BestPSM_PEP < 0.1)


condition_1_trinity = peptides_total_comparison_no_nan_fdr_unique['PSMs_high'] >= 2
condition_2_trinity = peptides_total_comparison_no_nan_fdr_unique['BestPSM_PEP'] <= 0.01
condition_1_icrdb = peptides_database_unique_fdr['PSMs_high'] >= 2
condition_2_icrdb = peptides_database_unique_fdr['BestPSM_PEP'] <= 0.01

combinied_condition_trinity = condition_1_trinity & condition_2_trinity
combinied_condition_icrdb = condition_1_icrdb & condition_2_icrdb
trinity_peptides_high_confidence_two_filt = peptides_total_comparison_no_nan_fdr_unique[combinied_condition_trinity]
icr_peptides_high_confidence_two_filt = peptides_database_unique_fdr[combinied_condition_icrdb]


trinity_peptides_two_filt_set = set(trinity_peptides_high_confidence_two_filt['Peptide'])
icr_peptides_two_filt_set = set(icr_peptides_high_confidence_two_filt['Peptide'])


venn2([trinity_peptides_two_filt_set, icr_peptides_two_filt_set], ('TRINITY-based immunopeptides database', 'ICR database'))

plt.show()
#%%
# # create Venn diagramm between TRINITY-based DB and ICR Peptides DB (best candidates by at least 4/5 filters and reproducible in 5 or more datasets)

novel_peptides_high_confidence_five_filt_set = set(novel_peptides['Peptide'])


venn_final = venn2([novel_peptides_high_confidence_five_filt_set, icr_peptides_two_filt_set], ('TRINITY-based immunopeptides database', 'ICR database'))
for text in venn_final.subset_labels:
    if text:  # Check if the text object is not None
        text.set_text('')
plt.show()

# %%
