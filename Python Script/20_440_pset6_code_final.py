import os
import pandas as pd
import matplotlib.pyplot as plt

#modify ONLY this data_folder path to point to your path for this downloaded folder, after zipping
data_folder = "/Users/shashaankvenkatesh/Downloads/20.440-Pset-6-master"

#do not change these two paths
data_folder_a = data_folder + "/Dataset A files (merged expression data)"
data_folder_b = data_folder + "/Dataset B files (top 100 gene lists)"

#read each of our six gene list CSVs into a parseable DataFrame
dfs = {}
for filename in os.listdir(data_folder_b):
    if filename.__contains__("names"):
        filepath = os.path.join(data_folder_b, filename)
        try:
            df_name = filename[:-4]  # remove the .csv extension
            dfs[df_name] = pd.read_csv(filepath)
            print(f"Successfully read {filename} into DataFrame '{df_name}'")
        except pd.errors.ParserError:
            print(f"Error parsing {filename}. Skipping.")
        except Exception as e:
            print(f"An error occurred while reading {filename}: {e}")


#for each of these dataframes, generate a list (preserving the order) of GeneName
gene_lists = {}
for df_name, df in dfs.items():
  if 'GeneName' in df.columns:
    gene_lists[df_name] = df['GeneName'].tolist()
  else:
    print(f"DataFrame '{df_name}' does not have a 'GeneName' column.")

#enrichment dataset
enrichment_filename = "/Merged_Tumor_Normal_Expression_AllPatients.csv"
df = pd.read_csv(data_folder_a + enrichment_filename)


#parse through df, and eliminate the rows for which there is more than one NaN across all 4 "tumor" columns (or) more than one NaN across all 4 "normal" columns.
tumor_cols = [col for col in df.columns if 'tumor' in col.lower()]
normal_cols = [col for col in df.columns if 'normal' in col.lower()]

# Check if tumor_cols and normal_cols are not empty
if not tumor_cols or not normal_cols:
    print("Error: Could not find 'tumor' or 'normal' columns in the DataFrame.")
else:
    df = df[~((df[tumor_cols].isnull().sum(axis=1) > 1) | (df[normal_cols].isnull().sum(axis=1) > 1))]

# Calculate mean and standard deviation for tumor and normal samples
df['tumor_mean'] = df[tumor_cols].mean(axis=1)
df['tumor_std'] = df[tumor_cols].std(axis=1)
df['normal_mean'] = df[normal_cols].mean(axis=1)
df['normal_std'] = df[normal_cols].std(axis=1)

# Calculate signal-to-noise ratio for tumor and normal
df['tumor_snr'] = df['tumor_mean'] / df['tumor_std']
df['normal_snr'] = df['normal_mean'] / df['normal_std']

#function to calculate running enrichment score
def running_es_calculation(snr, genes, gene_set, p):
  N = len(genes)
  N_R = 0
  running_score = [0]

  for gene, snr_value in zip(genes, snr):
    if gene in gene_set:
      N_R += abs(snr_value)**p

  P_miss = 1/(N - len(gene_set))

  for gene, snr_value in zip(genes, snr):
    curr_score = running_score[-1]
    if gene in gene_set:
      P_hit = (abs(snr_value)**p)/(N_R)
      running_score.append(float(P_hit + curr_score))
    else:
      running_score.append(float(curr_score - P_miss))

  running_score.pop(0)
  max_ES = max(running_score, key=abs)

  return running_score, max_ES


#list of all geneset names
geneset_names = ['ER+_top100_upregulated_with_names', 'ER+_top100_downregulated_with_names', 'HER2+_top100_upregulated_with_names',
                 'HER2+_top100_downregulated_with_names', 'TNBC_top100_upregulated_with_names', 'TNBC_top100_downregulated_with_names']


# Create sorted_tumor_snr DataFrame
sorted_tumor_snr = df.sort_values(by='tumor_snr', ascending=False)

# Create sorted_normal_snr DataFrame
sorted_normal_snr = df.sort_values(by='normal_snr', ascending=False)

#create a plot for the GSEA visualization
fig, axes = plt.subplots(2, 6, figsize=(12,8)) # Create 2 rows and 6 columns of subplots

#perform running_es_calculation for the six genesets in geneset_names, with p = 1. 
# use df['Gene'] for genes and perform this twice, once with sorted_tumor_snr and once with sorted_normal_snr. 
for i, geneset_name in enumerate(geneset_names):
  geneset = dfs[geneset_name]['GeneName'].tolist()
  label = geneset_name[:-11]

  # Tumor
  running_score_tumor, max_es_tumor = running_es_calculation(sorted_tumor_snr['tumor_snr'], sorted_tumor_snr['Gene'], geneset, 1)
  axes[0, i].plot(running_score_tumor)
  axes[0, i].set_title(f'Tumor - \n{label}', fontsize = 8, fontweight = "bold")

  # Normal
  running_score_normal, max_es_normal = running_es_calculation(sorted_normal_snr['normal_snr'], sorted_normal_snr['Gene'], geneset, 1)
  axes[1, i].plot(running_score_normal)
  axes[1, i].set_title(f'Normal - \n{label}', fontsize = 8, fontweight = "bold")


for ax in axes.flat:
  ax.set(xlabel='Rank', ylabel='Running Enrichment Score')

plt.tight_layout()
plt.show()