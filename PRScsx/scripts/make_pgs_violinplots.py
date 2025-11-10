import pandas as pd
import seaborn as sns
import argparse as ap
import matplotlib.pyplot as plt
import re
from pathlib import Path

def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")

    parser.add_argument('-i', '--input_file_list', required = True, help = 'list of input files')

    parser.add_argument('-a', '--ancestry_list', required = True, help = 'list of ancestries')

    parser.add_argument('-p', '--pheno_score_group_list', required = True, help = 'list of phenotypes')

    parser.add_argument('--write_pst_samples', required = True, help = 'write posterior samples flag from PRScsx')

    parser.add_argument('--meta', required = True, help = 'meta flag from PRScsx')

    return parser

args = make_arg_parser().parse_args()

# parse arguments
input_file_list = args.input_file_list
ancestry_list = args.ancestry_list
pheno_score_group_list = args.pheno_score_group_list
write_pst_samples = args.write_pst_samples
meta = args.meta

# parse lists
## input files
input_file_string = ','.join(input_file_list)
input_file_string = input_file_string.replace('[', '')
input_file_string = input_file_string.replace(']', '')
input_file_string = input_file_string.replace(',', '')
input_file_new_list = list(input_file_string.split(" "))

## ancestry
ancestry_string = ','.join(ancestry_list)
ancestry_string = ancestry_string.replace('[', '')
ancestry_string = ancestry_string.replace(']', '')
ancestry_string = ancestry_string.replace(',', '')
ancestry_new_list = list(ancestry_string.split(" "))

## pheno score group
pheno_score_group_string = ','.join(pheno_score_group_list)
pheno_score_group_string = pheno_score_group_string.replace('[', '')
pheno_score_group_string = pheno_score_group_string.replace(']', '')
pheno_score_group_string = pheno_score_group_string.replace(',', '')
pheno_score_group_new_list = list(pheno_score_group_string.split(" "))

# read in pgs output files
pgs_output_dfs = []
for file in input_file_new_list:
    filename = Path(file).name
    df = pd.read_csv(file, sep = '\t')
    for ancestry in ancestry_new_list:
        if ancestry in filename:
            df['ancestry_group'] = ancestry
            if (ancestry == 'META') and (write_pst_samples.lower() == 'true') and (meta.lower() == 'true'):
                # subset columns to remove extra pst sample scores
                df = df[['CHR', 'RSID', 'POS', 'A1', 'A2', 'PGS.1', 'ancestry_group']]
                # rename pgs columns
                df = df.rename(columns = {'PGS.1' : 'PGS'})
    for pheno_score_group in pheno_score_group_new_list:
        if pheno_score_group in filename:
            df['phenotype_score_group'] = pheno_score_group
    print(df)
    pgs_output_dfs.append(df)
bigdf = pd.concat(pgs_output_dfs, axis = 0, ignore_index = True)
print(bigdf)

sns.catplot(kind = 'violin', data = bigdf, x = 'ancestry_group', y = 'PGS', col = 'phenotype_score_group', col_wrap = 2, sharey = False)

# save the plot as PNG file
plt.savefig("combined_PGS_violinplot.png", dpi = 300)