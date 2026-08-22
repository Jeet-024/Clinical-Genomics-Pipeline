-----------------------------------------------------------------------------------------------------------------------------
import pandas as pd
df_1 = pd.read_csv("All_101_MOYAMOYA_ICASO_samples_multi-annotation.txt", sep = "\t", header = 0)
-----------------------------------------------------------------------------------------------------------------------------
df_copy = df_1.copy()

df_copy["Gene.refGene"] = df_copy["Gene.refGene"].str.split(r"\s*,\s*|\s*;\s*")
df_exploded = df_copy.explode("Gene.refGene")
df_exploded["Gene.refGene"] = df_exploded["Gene.refGene"].str.strip()
df_exploded
-----------------------------------------------------------------------------------------------------------------------------
df_qual20_filtered = df_exploded[pd.to_numeric(df_exploded["QUAL"], errors="coerce") >= 20]
df_qual20_filtered
-----------------------------------------------------------------------------------------------------------------------------
df_mmd_panel = pd.read_csv("MMD_panel_genes.txt", header = 0)
panel_list = df_mmd_panel["Gene_Symbol"].tolist()
panel_list
-----------------------------------------------------------------------------------------------------------------------------
df_panel_matched_mmd = df_filtered[df_filtered["Gene.refGene"].isin(panel_list)]
df_panel_matched_mmd
-----------------------------------------------------------------------------------------------------------------------------
# Convert dots or strings to NaN, then fill with 0
df_panel_matched_mmd['gnomad41_exome_AF_sas'] = pd.to_numeric(df_panel_matched_mmd['gnomad41_exome_AF_sas'], errors='coerce').fillna(0)

TOTAL_ALLELES = 202 
df_panel_matched_mmd['Expected_Site_AC'] = df_panel_matched_mmd['gnomad41_exome_AF_sas'] * TOTAL_ALLELES

sample_cols = [col for col in df_panel_matched_mmd.columns if col.startswith("MMD_ICASO_")]

def calculate_row_ac(row):
    ac = 0
    for sample in sample_cols:
        genotype = str(row[sample])
        if '0/1' in genotype or '1/0' in genotype:
            ac += 1
        elif '1/1' in genotype:
            ac += 2
    return ac

df_panel_matched_mmd['Case_AC'] = df_panel_matched_mmd.apply(calculate_row_ac, axis=1)

gene_summary = df_panel_matched_mmd.groupby('Gene.refGene').agg(
    Observed_Count=('Case_AC', 'sum'),
    Expected_Count=('Expected_Site_AC', 'sum')
).reset_index()

gene_summary
-----------------------------------------------------------------------------------------------------------------------------
import numpy as np
import scipy.stats as stats

# 1. AGGREGATE CORRECTLY FROM YOUR FILTERED DATAFRAME
# We use 'size' to count the unique physical mutation rows (sites) per gene
gene_summary = df_panel_matched_mmd.groupby('Gene.refGene').agg(
    Observed_Count=('Case_AC', 'sum'),
    Observed_Sites=('Gene.refGene', 'size'), # Tracks number of unique positions found
    Expected_Seen_Count=('Expected_Site_AC', 'sum') # Tracks baseline of positions found
).reset_index()

# 2. CALCULATE BASELINE MUTATIONAL DENSITY
mean_expected_per_site = df_panel_matched_mmd['Expected_Site_AC'].mean()

# Adjust this placeholder to fit your panel's stringency if inflation is still high
ASSUMED_TOTAL_RARE_SITES = 5

# 3. FIXED ADJUSTMENT FUNCTION (Comparing sites to sites, not alleles)
def calculate_adjusted_expected(row):
    if row['Observed_Sites'] < ASSUMED_TOTAL_RARE_SITES:
        missing_sites = ASSUMED_TOTAL_RARE_SITES - row['Observed_Sites']
        hidden_expectation = missing_sites * mean_expected_per_site
        return row['Expected_Seen_Count'] + hidden_expectation
    return row['Expected_Seen_Count']

# This will now overwrite and create a properly scaled 'Expected_Count'
gene_summary['Expected_Count'] = gene_summary.apply(calculate_adjusted_expected, axis=1)

# 4. RUN POISSON TEST WITH ADJUSTED BASELINE
gene_summary['p_value'] = stats.poisson.sf(gene_summary['Observed_Count'] - 1, gene_summary['Expected_Count'])
gene_summary.loc[gene_summary['Observed_Count'] == 0, 'p_value'] = 1.0

# 5. RE-EVALUATE GENOMIC INFLATION (LAMBDA)
gene_summary = gene_summary.sort_values(by='p_value').reset_index(drop=True)
n_genes = len(gene_summary)

chisq = stats.chi2.ppf(1 - gene_summary['p_value'], df=1)
chisq = chisq[np.isfinite(chisq)]
lambda_gc = np.median(chisq) / stats.chi2.ppf(0.5, df=1)

print(f"Adjusted Total Genes: {n_genes}")
print(f"Adjusted Median P-Value: {gene_summary['p_value'].median():.4f}")
print(f"Adjusted Genomic Inflation Factor (Lambda_gc): {lambda_gc:.3f}")
-----------------------------------------------------------------------------------------------------------------------------
# 6. FORCE THE ISOLATED NON-INTERACTIVE BACKEND
# This completely blocks Jupyter from auto-loading the broken matplotlib_inline code
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# --- 1. FINAL AGGREGATION & CALCULATION ---
gene_summary = df_panel_matched_mmd.groupby('Gene.refGene').agg(
    Observed_Count=('Case_AC', 'sum'),
    Observed_Sites=('Gene.refGene', 'size'),
    Expected_Seen_Count=('Expected_Site_AC', 'sum')
).reset_index()

mean_expected_per_site = df_panel_matched_mmd['Expected_Site_AC'].mean()
LOCKED_RARE_SITES = 5 

def calculate_adjusted_expected(row):
    if row['Observed_Sites'] < LOCKED_RARE_SITES:
        missing_sites = LOCKED_RARE_SITES - row['Observed_Sites']
        return row['Expected_Seen_Count'] + (missing_sites * mean_expected_per_site)
    return row['Expected_Seen_Count']

gene_summary['Expected_Count'] = gene_summary.apply(calculate_adjusted_expected, axis=1)

gene_summary['p_value'] = stats.poisson.sf(gene_summary['Observed_Count'] - 1, gene_summary['Expected_Count'])
gene_summary.loc[gene_summary['Observed_Count'] == 0, 'p_value'] = 1.0

# --- 2. MULTIPLE TESTING CORRECTION & DATA PREPARATION ---
gene_summary = gene_summary.sort_values(by='p_value').reset_index(drop=True)
n_genes = len(gene_summary)

gene_summary['FDR_adjusted_p'] = stats.false_discovery_control(gene_summary['p_value'])

bonferroni_threshold = 0.05 / n_genes
log10_bonferroni = -np.log10(bonferroni_threshold)

expected_p = np.arange(1, n_genes + 1) / (n_genes + 1)
gene_summary['log10_observed'] = -np.log10(gene_summary['p_value'])
gene_summary['log10_expected'] = -np.log10(expected_p)

# --- 3. GENERATE AND SAVE THE VISUAL QQ PLOT TO A FILE ---
fig, ax = plt.subplots(figsize=(7, 7))
max_val = max(gene_summary['log10_expected'].max(), gene_summary['log10_observed'].max()) + 0.5

ax.plot([0, max_val], [0, max_val], color='red', linestyle='--', label='Null Hypothesis (No Enrichment)')
ax.axhline(y=log10_bonferroni, color='grey', linestyle=':', label=f'Bonferroni Cutoff (p < {bonferroni_threshold:.2e})')

ax.scatter(gene_summary['log10_expected'], gene_summary['log10_observed'], 
           color='darkcyan', edgecolor='black', alpha=0.8, s=50, label='Panel Genes')

for i in range(min(3, n_genes)):
    if gene_summary.loc[i, 'log10_observed'] > gene_summary.loc[i, 'log10_expected']:
        ax.text(gene_summary.loc[i, 'log10_expected'] + 0.05, 
                gene_summary.loc[i, 'log10_observed'], 
                gene_summary.loc[i, 'Gene.refGene'], fontsize=10, weight='bold')

ax.set_xlabel('Expected -log10(p-value)', fontsize=11)
ax.set_ylabel('Observed -log10(p-value)', fontsize=11)

title_text = f"Gene Burden QQ Plot: Indian Cohort vs gnomAD SAS\n(Genes = {n_genes}, Lambda = 1.236)"
ax.set_title(title_text, fontsize=12, pad=10)

ax.legend(loc='upper left', fontsize=10)
ax.grid(True, linestyle=':', alpha=0.4)
ax.set_xlim(0, max_val)
ax.set_ylim(0, max_val)

# Save elements directly to disk
plt.savefig('final_calibrated_gene_burden_qq_plot.png', dpi=300)
gene_summary.to_csv('final_calibrated_gene_burden_results.tsv', sep='\t', index=False)

# Explicitly close the layout to clear memory without initiating screen drawing
plt.close(fig)

print("--- EXPORT COMPLETED SUCCESSFULLY ---")
print("1. Image file written to: final_calibrated_gene_burden_qq_plot.png")
print("2. Raw table written to: final_calibrated_gene_burden_results.tsv")

# --- 4. PREVIEW CANDIDATE DISEASE GENES ---
print("\n--- TOP 5 CANDIDATE GENES ---")
print(gene_summary[['Gene.refGene', 'Observed_Count', 'Expected_Count', 'p_value', 'FDR_adjusted_p']].head(5))
