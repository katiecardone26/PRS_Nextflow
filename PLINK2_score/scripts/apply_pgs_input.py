#! /opt/conda/bin/python

# load packages
import pandas as pd
import argparse as ap
import sys
import os

# define arguments
def make_arg_parser():
    parser = ap.ArgumentParser(description = ".")
    # validation population
    parser.add_argument('--valPop', required = True, help = 'Validation population')

    # ancestry
    parser.add_argument('--ancestry', required = True, help = 'Ancestry')

    # phenotype
    parser.add_argument('--pheno', required = True, help = 'Phenotype')

    # plink flag (bfile or pfile)
    parser.add_argument('--plinkFlag', required = True, help = 'bfile or pfile')

    # pvar or bim files for all chromosomes
    parser.add_argument('--varFiles', required = True, help = '.pvar or .bim files for all chromosomes')

    # validation population variant ID formats
    parser.add_argument('--popsIds', required = True, nargs = '*', help = 'list of validation population variant ID formats')

    # score file variant ID format
    parser.add_argument('--scoreIdFormat', required = True, help = 'Format of the score file variant IDs')

    # score file name
    parser.add_argument('--scoreFile', required = True,  help = 'score file')

    # chromosome column name in score file
    parser.add_argument('--scoreChrCol', required = True, help = 'position/base pair column name in score files')

    # position column name in score file
    parser.add_argument('--scorePosCol', required = True, help = 'position/base pair column name in score files')

    # variant ID column name in score file
    parser.add_argument('--scoreIdCol', required = True, help = 'variant ID column name in score files')

    # A1 columns name in score file
    parser.add_argument('--scoreA1Col', required = True, help = 'A1 column name in score files')

    # A2 column name in score file
    parser.add_argument('--scoreA2Col', required = True, help = 'A2 column name in score files')
    
    # PGS column name in score file
    parser.add_argument('--scorePGSCol', required = True, help = 'PGS column name in score files')
    
    return parser


args = make_arg_parser().parse_args()

# make arguments into script variables
# reformat val_pop_variant_id_collect tuple into a dictionary
val_pop_variant_id_collect_list = args.popsIds
val_pop_list = val_pop_variant_id_collect_list[1::2]
var_id_type_list = val_pop_variant_id_collect_list[::2]

# variables
## plink flag
plink_flag = args.plinkFlag
print(plink_flag)
## validation population
validation_population = args.valPop
## ancestry
ancestry = args.ancestry
## phenotype
pheno = args.pheno

# column names
## chromosome
chr_colname = args.scoreChrCol
## position
pos_colname = args.scorePosCol
## variant ID
variant_id_colname = args.scoreIdCol
## A1
a1_colname = args.scoreA1Col
## A2
a2_colname = args.scoreA2Col
## PGS
pgs_colname = args.scorePGSCol

# input files
## score file name
score_file_name = args.scoreFile
## list of bim/pvar files
var_files_list = args.varFiles

# add score variant ID to list
var_id_type_list.append(args.scoreIdFormat)

# check if all variant ID formats are the same
if len(set(var_id_type_list)) == 1:
    print(f"{score_file_name} and {validation_population} have the same variant ID....making apply PGS input without reformatting variant IDs")

else:
    print(f"score files and {validation_population} have different variant IDs....converting bim/pvar file variant IDs to match score files")
    
    # reformat list of variant files
    var_files_string=','.join(var_files_list)
    var_files_string = var_files_string.replace('[', '')
    var_files_string = var_files_string.replace(']', '')
    var_files_string = var_files_string.replace(',', '')
    var_files_new_list = list(var_files_string.split(" "))

    # read in and concatenate variant files
    print('Reading in and concatenating chromosome separated plink variant files', flush = True)
    var_files_dfs = []
    for var_file in var_files_new_list:
        # bfile
        if plink_flag == 'bfile':
            plink = pd.read_csv(var_file, sep = None, engine = 'python', header = None, names = ['CHR', 'ID', 'CM', 'POS', 'A1', 'A2'])
            plink = plink.set_index(['POS', 'A1', 'A2', 'CHR'])
        # pfile
        elif plink_flag == 'pfile':
            plink = pd.read_csv(var_file, sep = None, engine = 'python', comment = '#', header = None)

            # get column names
            with open(var_file,"r") as fi:
                pvar_cols = []
                for ln in fi:
                    if ln.startswith("#CHROM"):
                        pvar_cols.append(ln[2:])
                        break
            pvar_cols = '\t'.join(pvar_cols)
            pvar_cols = pvar_cols.strip()
            pvar_cols = pvar_cols.replace('HROM', '#CHROM')
            pvar_cols = list(pvar_cols.split('\t'))
            plink.columns = pvar_cols
            plink.rename(columns={'REF' : 'A2', 'ALT' : 'A1', '#CHROM' : 'CHR'}, inplace = True)
            plink = plink.set_index(['POS', 'A1', 'A2', 'CHR'])
        else:
            raise ValueError("plink flag is not --bfile or --pfile")
        # append file to list
        var_files_dfs.append(plink)
    # concatenate chromosome separated files
    var_files_cat = pd.concat(var_files_dfs, axis = 0)
    
    # read in score file
    score_file = pd.read_table(score_file_name, sep = None, engine = 'python', dtype = {pos_colname: int, variant_id_colname: str, a1_colname: str, a2_colname: str})
    
    # get column order variable
    og_col_order = score_file.columns

    # check if plink chromosome column contains "chr" prefix (AOU has files formatted like this)
    # if so, reformat score file so they match
    if var_files_cat.index.get_level_values('CHR').astype(str).str.contains('chr').all():
        if not score_file[chr_colname].astype(str).str.contains('chr').all():
            score_file[chr_colname] = 'chr' + score_file[chr_colname].astype(str)
    else:
        if score_file[chr_colname].astype(str).str.contains('chr').all():
            score_file[chr_colname] = score_file[chr_colname].astype(str).str.replace('chr', '')

    # set score file index as position, A1, and A2
    ## A1- doesn't account for allele flips
    score_file_A1 = score_file.set_index([pos_colname, a1_colname, a2_colname, chr_colname])
    ## A2- accounts for allele flips
    score_file_A2 = score_file.set_index([pos_colname, a2_colname, a1_colname, chr_colname])
    
    # find intersections while account for allele flips
    forward_match = score_file_A1.index.intersection(var_files_cat.index)
    backward_match = score_file_A2.index.intersection(var_files_cat.index)

    print(f'Forward Match: {len(forward_match):,}')
    print(f'Backward Match: {len(backward_match):,}')

    # make copy of score file dataframe
    new_score_file=score_file.copy()
    
    # replace variant IDs that matched with plink file IDs
    ## forward match- replace IDs in new file
    new_score_file.set_index([pos_colname, a1_colname, a2_colname, chr_colname],inplace=True)
    new_score_file.loc[forward_match, variant_id_colname] = var_files_cat.loc[forward_match, 'ID']
    ## backward match- replace IDs but take inverse of score to account fot allele flips
    new_score_file=new_score_file.reset_index()[og_col_order]

    new_score_file.set_index([pos_colname, a2_colname, a1_colname, chr_colname], inplace = True)
    new_score_file=new_score_file[~new_score_file.index.duplicated()]
    new_score_file.loc[backward_match, variant_id_colname] = var_files_cat.loc[backward_match, 'ID']
    new_score_file.loc[backward_match, pgs_colname] = new_score_file.loc[backward_match, pgs_colname] * -1
    
    # reset index
    new_score_file = new_score_file.reset_index()[og_col_order]
    print(new_score_file)

# remove all columns except variant ID, A1, and score
new_score_file_sub = new_score_file[[variant_id_colname, a1_colname, pgs_colname]]

# rename columns
new_score_file_sub.rename(columns = {variant_id_colname : 'ID',
                                    a1_colname : 'ALLELE',
                                    pgs_colname : f'{ancestry}.{pheno}.SCORE'}, inplace = True)

# export file
new_score_file_sub.to_csv(f'{ancestry}.{pheno}.{validation_population}.Apply_PGS_Input.txt', sep = '\t', index = None)
    
    