import pandas as pd
import seaborn as sns
from scipy.stats import zscore
import matplotlib.pyplot as plt
import argparse as ap
import os

def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")
    
    # add non-optional argument for data
    parser.add_argument('-d', '--data', required = True, help = '.csv prs summary score files')
    
    # add non-optional argument for ancestry column
    parser.add_argument('-a', '--ancestry_list', required = True, help = 'ancestry list')
    
    # add non-optional argument for phenotype column
    parser.add_argument('-t', '--phenotype_list', required = True, help = 'phenotype list')
    
    # add non-optional argument for population
    parser.add_argument('-p', '--population', required = True, help = 'validation population')
    
    return parser

# parse args
args = make_arg_parser().parse_args()

data = args.data
ancestry_list = args.ancestry_list
pheno_list = args.phenotype_list
population = args.population

# parse ancestry and pheno lists
## ancestry
ancestry_string = ','.join(ancestry_list)
ancestry_string = ancestry_string.replace('[', '')
ancestry_string = ancestry_string.replace(']', '')
ancestry_string = ancestry_string.replace(',', '')
ancestry_new_list = list(ancestry_string.split(" "))

## pheno
pheno_string = ','.join(pheno_list)
pheno_string = pheno_string.replace('[', '')
pheno_string = pheno_string.replace(']', '')
pheno_string = pheno_string.replace(',', '')
pheno_new_list = list(pheno_string.split(" "))

# import Polygenic Risk Score Data
data = pd.read_table(data)

if '#FID' in data.columns:
    # drop columns
    to_drop = ['ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM', '#FID']
    data.drop(columns = to_drop, inplace = True)

    # set ID column as index
    data.set_index('IID', inplace = True)
elif '#IID' in data.columns:
    # drop columns
    to_drop = ['ALLELE_CT', 'NAMED_ALLELE_DOSAGE_SUM']
    data.drop(columns = to_drop, inplace = True)

    # set ID column as index
    data.rename(columns = {'#IID' : 'IID'}, inplace = True)
    data.set_index('IID', inplace = True)
else:
    raise ValueError('ID columns do not meet conditions')

# Z-transform the scores
data_zscore = data.apply(zscore)

# convert to long format
data_long = data_zscore.reset_index().melt(id_vars = 'IID')

# create ancestry column
ancestry_dfs = []
for ancestry in ancestry_new_list:
    match_cols = data.columns[data.columns.str.contains(ancestry)].tolist()
    ancestry_df = pd.DataFrame(data = {'variable' : match_cols, 'ancestry_col' : ancestry})
    ancestry_dfs.append(ancestry_df)
ancestry_df_all = pd.concat(ancestry_dfs, axis = 0)

# create pheno column
pheno_dfs = []
for pheno in pheno_new_list:
    match_cols = data.columns[data.columns.str.contains(pheno)].tolist()
    pheno_df = pd.DataFrame(data = {'variable' : match_cols, 'phenotype_col' : pheno})
    pheno_dfs.append(pheno_df)
pheno_df_all = pd.concat(pheno_dfs, axis = 0)

# merge dataframes
## ancestry and pheno
ancestry_pheno= ancestry_df_all.merge(pheno_df_all, on = 'variable')
## ancestry/pheno and data long
data_long = data_long.merge(ancestry_pheno, on = 'variable', how = 'left')

# Convert the 'ancestry_col' column to a categorical type with the custom sort order
data_long['ancestry_col'] = pd.Categorical(data_long['ancestry_col'], categories = ancestry_new_list, ordered = True)

# get color palette
pal = sns.color_palette('deep', data_long['ancestry_col'].nunique())

# BoxPlots
# ---------
for phenotype, subDF in data_long.groupby('phenotype_col'): # Iterate over phenos
    
    # Box Plots
    # ---------
    ax = sns.boxplot(data = subDF, y = 'ancestry_col', x = "value", hue = 'ancestry_col', palette = pal)
    ax.set_title(f"Z-Transformed PRS\nPopulation: {population}; Trait: {phenotype}")
    ax.set_xticklabels(ax.get_xticklabels(), rotation = 90)
    ax.set_ylabel("Cohort")
    ax.set_xlabel("PRS (z-score)")
    outfile = f"{population}.all_computed_PGS_scores.boxplots.{phenotype}.png"
    plt.savefig(outfile)
    plt.clf()
    
    # Density Plots
    # -------------
    # make sure axis background is transparent
    sns.set_theme(style = "white", rc = {"axes.facecolor": (0, 0, 0, 0)})

    # Initialize the FacetGrid object
    pal = sns.color_palette('deep', data_long['ancestry_col'].nunique())
    g = sns.FacetGrid(subDF, row = 'ancestry_col', hue = 'ancestry_col', height = 1.5, aspect = 7, palette = pal)

    # Draw the densities in a few steps
    g.map(sns.kdeplot, "value", bw_adjust = .5, clip_on = False, fill = True, alpha = 1, linewidth = 1.5)
    # draw another plot to outline the fill, either white or black
    g.map(sns.kdeplot, "value", clip_on = False, color = "w", lw = 2, bw_adjust = .5)

    # refline for each
    g.map(plt.axhline, y = 0, lw = 2, clip_on = False)


    # Define and use a simple function to label the plot in axes coordinates
    def label(x, color, label):
        ax = plt.gca()
        ax.text(0, .2, label, fontweight = "bold", color = color, ha = "left", va = "center", transform = ax.transAxes)

    g.map(label, "value")

    # Set the subplots to overlap
    g.fig.subplots_adjust(hspace = -.5)

    # Remove axes details that don't play well with overlap
    g.set_titles("")
    g.set(yticks = [], xlabel = "PRS (Z-Score)", ylabel = "")
    g.despine(bottom = True, left = True)

    # set labels
    # g.set_axis_labels("PRS (Z-Score)")

    # set big title
    plt.suptitle(f'PRS Scores by Cohort\nPopulation: {population}; Trait: {phenotype}', y = 0.98)

    # uncomment the following line if there's a tight layout warning
    # g.fig.tight_layout()
    
    outfile = f"{population}.all_computed_PGS_scores.densityplots.{phenotype}.png"
    plt.savefig(outfile)
    plt.clf()
    