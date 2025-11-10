#! /appl/nextflow-23.04.1.5866/bin/nextflow

// enable DSL2 syntax***
nextflow.enable.dsl = 2

// log info
log.info """\
NEXTFLOW - DSL2 - PRScsx - P I P E L I N E
============================================
run_as               :  ${workflow.commandLine}
run_location         :  ${launchDir}
started_at           :  ${workflow.start}
container            :  ${workflow.containerEngine}:${workflow.container}

Ancestries, Phenotypes and Chromosomes
============================================
chromosomes          :  ${params.chromosome_list}

Input Directories and LD Panel Type
============================================
genome_build         :  ${params.genome_build}
LD_directory         :  ${params.ld_directory}
LD_panel_type        :  ${params.ld_panel_type}
Sumstats colnames    :  ${params.sumstats_colnames}

PRScsx Parameters
============================================
prscsx_param_a       :  ${params.param_a}
prscsx_param_b       :  ${params.param_b}
prscsx_param_phi     :  ${params.param_phi}
prscsx_mcmc_iter     :  ${params.mcmc_iterations}
prscsx_mcmc_burn     :  ${params.mcmc_burnin}
prscsx_mcmc_thin     :  ${params.mcmc_thinning_factor}
prscsx_meta_flag     :  ${params.meta_flag}
prscsx_seed          :  ${params.seed}

"""
.stripIndent()

workflow {
    /*
    this section is called first if:
        PRScsx is the first workflow called
        PRScsx is called independently of other workflows
        Another GWAS workflow is not being stitched to PRScsx workflow
        Channels from another workflow are not inputs to the PRScsx workflow
    */
    prscsx_inputs = PRScsx_setup()
    prscsx_weights = PRScsx(prscsx_inputs)
}

workflow PRScsx_setup {
    /*
    sets up inputs if:
        PRScsx is the first workflow called
        PRScsx is called independently of other workflows
        Another GWAS workflow is not being stitched to PRScsx workflow
        Channels from another workflow are not inputs to the PRScsx workflow
    */
    main:
        // convert descriptor file to groovy file type
        descriptor_file = new File(params.input_descriptor_table_filename.toString())
        // read in descriptor file into a groovy tuple
        // descriptor_tuple = new Tuple(*descriptor_file.readLines().tail()*.split(','))
        // create channel from descriptor channel
        // descriptor_channel = Channel.fromList(descriptor_tuple)

        descriptor_lines = descriptor_file.readLines()
        if (descriptor_file.text.tokenize('\n').size() <= 2) {
            descriptor_lines[1..-1].each {
                line ->
                /* groovylint-disable-next-line GStringExpressionWithinString */
                lineG = line.replace('${launchDir}', "${launchDir}")
                line_parts = lineG.toString().trim().replace('[', '').replace(']', '').split(',') as List
                descriptor_tuple = line_parts
            }
        } else {
            descriptor_tuple = []
            descriptor_lines[1..-1].each {
                line ->
                /* groovylint-disable-next-line GStringExpressionWithinString */
                lineG = line.replace('${launchDir}', "${launchDir}")
                line_parts = lineG.toString().trim().replace('[', '').replace(']', '').split(',') as List
                descriptor_tuple.add(line_parts)
            }
        }

        descriptor_channel = Channel.fromList(descriptor_tuple)

        // parse header and retrieve indices of essential columns of descriptor file based on their column names
        descriptor_header = descriptor_file.readLines().get(0).split(',')
        descriptor_header.eachWithIndex { item, index ->
            if (item == params['input_descriptor_table_colnames']['cohort_colname']) {
                cohort_index = index
            }
            if (item == params['input_descriptor_table_colnames']['ld_ref_ancestry_colname']) {
                ancestry_index = index
            }
            if (item == params['input_descriptor_table_colnames']['phenotype_colname']) {
                pheno_index=index
            }
            if (item == params['input_descriptor_table_colnames']['sample_size_colname']) {
                sample_size_index = index
            }
            if (item == params['input_descriptor_table_colnames']['sumstats_full_filepath_colname']) {
                sumstats_index = index
            }
        }
        // makes each item (table row) in the list a separate tuple
        if (descriptor_file.text.tokenize('\n').size() <= 2) {
            descriptor_channel_collect = descriptor_channel.collect()
        } else {
            descriptor_channel_collect = descriptor_channel.map{ arr -> new Tuple(
                arr[cohort_index],
                arr[ancestry_index],
                arr[pheno_index],
                arr[sample_size_index],
                arr[sumstats_index]) }
        }
        // process prscsx score group to cohort map
        prscsx_group = Channel.fromList(params.prscsx_group_cohort_map.keySet().toList())
        prscsx_group_cohort=prscsx_group.map{ prscsx_group -> new Tuple (params['prscsx_group_cohort_map'][prscsx_group], prscsx_group) }.transpose()

        // add prscsx score group to tuple
        descriptor_score_group = descriptor_channel_collect.groupTuple(by: 0).join(prscsx_group_cohort, by: 0).transpose()

        // reorder tuple
        sumstats_summary = descriptor_score_group.map { cohort, ancestry, pheno, sample_size, sumstats, score_group -> new Tuple(
            score_group,
            cohort,
            ancestry,
            pheno,
            sample_size,
            sumstats
        ) }

        // create ancestry list
        ancestry_list = []
        if (descriptor_file.text.tokenize('\n').size() <= 2) {
            ancestry = descriptor_tuple[ancestry_index]
            ancestry_list.add(ancestry)
        } else {
            for (i in descriptor_tuple) {
                ancestry = i[ancestry_index]
                ancestry_list.add(ancestry)
            }
        }
        ancestry_list = ancestry_list.unique()
    emit:
        sumstats_summary
        ancestry_list
}

workflow PRScsx {
    /*
    PRScsx workflow includes:
        workflow code for when PRScsx is called independently or another workflow feeds into PRScsx
        process code
    If another workflow is stitched to PRScsx and happens BEFORE PRScsx, this section is called first
    If PRScsx is called independently, the workflow and PRScsx_setup sections are called first, then this section is called
    */
    take:
        cohort_pheno_sumstats
        pheno_table
    main:
        if (params.gwas_meta_workflow_stitching == true) {
            // parse pheno table from input workflows
            pheno_table_lines = pheno_table.map { file -> new Tuple(new File(file.toString()).readLines().tail()) }.transpose()
            pheno_table_lines_clean = pheno_table_lines.map { line -> new Tuple(line.toString().trim().replace('[', '').replace(']', '').split(',')) }
            // extract cohort, pheno and sample size from pheno table
            cohort_pheno_sample_size = pheno_table_lines_clean.map { arr -> new Tuple(arr[0], arr[1], arr[2]) }
            cohort_pheno_sample_size.view()
            // make ancestry list
            ancestry_list = params.ancestry_cohort_map.keySet().toList()
            // make ancestry channel
            ancestry = Channel.fromList(ancestry_list)
            // get cohort ancestry map
            cohort_ancestry = ancestry.map { ancestry -> new Tuple (params["ancestry_cohort_map"]["${ancestry}"], ancestry) }.transpose()

            // combine channels
            // cohort pheno sample size with cohort ancerstry
            cohort_pheno_sample_size_ancestry = cohort_pheno_sample_size.groupTuple(by: 0).join(cohort_ancestry, by: 0).transpose()
            // cohort pheno sample size ancestry with cohort pheno sumstats
            descriptor_channel = cohort_pheno_sample_size_ancestry.join(cohort_pheno_sumstats, by:[0,1])

            // process prscsx score group to cohort map
            prscsx_group = Channel.fromList(params.prscsx_group_cohort_map.keySet().toList())
            prscsx_group_cohort = prscsx_group.map{ prscsx_group -> new Tuple (params['prscsx_group_cohort_map'][prscsx_group], prscsx_group) }.transpose()

            // add prscsx score group to tuple
            descriptor_score_group = descriptor_channel.groupTuple(by: 0).join(prscsx_group_cohort, by: 0).transpose()
            // reorders tuple
            sumstats_summary = descriptor_score_group.map { cohort, pheno, sample_size, ancestry, sumstats, score_group -> new Tuple(
                score_group,
                cohort,
                ancestry,
                pheno,
                sample_size,
                sumstats
            ) }
            // make pgs input and bim
            // create sumstats input (doesn't include sumstats dir parameters since workflows are stitched)
            sumstats = sumstats_summary.map { score_group, cohort, ancestry, pheno, sample_size, sumstats -> new Tuple(
                score_group,
                cohort,
                ancestry,
                pheno,
                sumstats
            ) }
        } else {
            // rename take variables
            sumstats_summary = cohort_pheno_sumstats
            ancestry_list = pheno_table

            // make pgs input and bim
            // create sumstats input
            sumstats = sumstats_summary.map { score_group, cohort, ancestry, pheno, sample_size, sumstats -> new Tuple(
                score_group,
                cohort,
                ancestry,
                pheno,
                sumstats) }
        }
        // set up chromosome channel
        chromosome = Channel.fromList(params.chromosome_list)

        // make pgs input and bim
        // define sumstats column name channels
        chr_colname= "${params['sumstats_colnames']['chr_colname']}"
        pos_colname = "${params['sumstats_colnames']['pos_colname']}"
        a1_colname = "${params['sumstats_colnames']['a1_colname']}"
        a2_colname = "${params['sumstats_colnames']['a2_colname']}"
        quant_beta_or_colname = "${params['sumstats_colnames']['quant_beta_or_colname']}"
        bin_beta_or_colname = "${params['sumstats_colnames']['bin_beta_or_colname']}"
        pval_se_colname = "${params['sumstats_colnames']['pval_se_colname']}"
        // define genome build channel
        genome_build = "${params.genome_build}"
        // define LD panel type column name
        ld_panel_type = "${params.ld_panel_type}"
        // define my python
        my_python = "${params.my_python}"
        // call make pgs_input_and_bim process
        (pgs_input_tuple, bim_tuple) = make_pgs_input_and_bim(sumstats,
                                                               chr_colname,
                                                               pos_colname,
                                                               a1_colname,
                                                               a2_colname,
                                                               quant_beta_or_colname,
                                                               bin_beta_or_colname,
                                                               pval_se_colname,
                                                               genome_build,
                                                               ld_panel_type,
                                                               my_python)

        // cat_bim process
        // inputs
        // this command collects the bim files so they can be concatenated by phenotype in the process
        all_bim_files_tuple = bim_tuple.groupTuple(by: [0, 1])
        // call cat_bim process
        bim_final = cat_bim(all_bim_files_tuple, my_python)

        // prscsx
        /*
        map: extracts sample size & parallelization variables (pheno and score group) from sumstats summary channel
        first join: joins with PGS input channel
        groupTuple: groups sample sizes and PGS input files by score group and phenotype
        second join: joins with concatenated bim file channel
        combine: combines with chromosome channel
        */

        pgs_input_sample_size_pheno_bim_file_chr_group_tuple = sumstats_summary.map { 
            score_group, cohort, ancestry, pheno, sample_size, sumstats -> new Tuple (
            score_group, pheno, sample_size) 
            }.join(
                pgs_input_tuple, by: [0, 1]
            ).groupTuple(
                by: [0, 1]
            ).join(
                bim_final, by: [0, 1]
            ).combine(chromosome)
         // add meta to ancestry list if the META flag is true
        if (params.meta_flag == 'True') {
            ancestry_list.add('META')
        }

        // define LD directory channel
        ld_dir = "${params['ld_directory']}"
        // make PRScsx flags channels
        param_a = "${params.param_a}"
        param_b = "${params.param_b}"
        param_phi = "${params.param_phi}"
        mcmc_iterations = "${params.mcmc_iterations}"
        mcmc_burnin = "${params.mcmc_burnin}"
        mcmc_thinning_factor = "${params.mcmc_thinning_factor}"
        meta_flag = "${params.meta_flag}"
        write_posterior_samples = "${params.write_posterior_samples}"
        seed = "${params.seed}"
        // define PRScsx path
        my_prscsx = "${params.my_prscsx}"
        // call prscsx process
        prscsx_output = prscsx(pgs_input_sample_size_pheno_bim_file_chr_group_tuple,
                                ancestry_list,
                                ld_dir,
                                param_a,
                                param_b,
                                param_phi,
                                mcmc_iterations,
                                mcmc_burnin,
                                mcmc_thinning_factor,
                                meta_flag,
                                write_posterior_samples,
                                seed,
                                my_python,
                                my_prscsx)
        // cat pgs outputs process
        // tranpose PGS output tuple such that the format is [pheno,cohort,chr,file],[pheno,cohort,chr,file],etc.
        /*
        if META flag is true- create cohort name with META included (not an output of PRScsx process)
        create channel from sumstats summary with score group, ancestry, pheno and cohort
        */
        score_group_ancestry_pheno_cohort = sumstats_summary.map { score_group, cohort, ancestry, pheno, sample_size, filename -> new Tuple (score_group, ancestry, pheno, cohort) }
        // conditionally add META cohort if meta flag = true
        if (params.meta_flag == 'True') {
            meta = sumstats_summary.map{ score_group, cohort, ancestry, pheno, sample_size, filename -> 
            new Tuple (score_group, 'META', pheno, 'META_' + score_group) }
            meta = meta.unique()
            score_group_ancestry_pheno_cohort_meta = score_group_ancestry_pheno_cohort.concat(meta)
        } else {
            score_group_ancestry_pheno_cohort_meta = score_group_ancestry_pheno_cohort
        }
        /*
        tranpose: transpose channel so ancestry-specific outputs that were computed together in PRScsx are no longer grouped
        map: create a channel with score group, ancestry (extracted from output file name), pheno, and output file
        groupTuple: group channel by score group, ancestry, and phenotype
        join: join with channel with refactored cohort variables
        */
        prscsx_output_reformat = prscsx_output.transpose().map {
            score_group, pheno, chromosome, file ->
            new Tuple(
            score_group,
            file.toString().substring((file.toString().indexOf('.PRScsx_') + 8), (file.toString().indexOf('_pst_'))),
            pheno,
            file)
        }.groupTuple(
            by: [0, 1, 2]
        ).join(
            score_group_ancestry_pheno_cohort_meta, by: [0, 1, 2]
        )
        // true/false flag to convert output to b38
        b38_flag = "${params.b38_output_genome_build}"
        // call cat_pgs_outputs process
        prscsx_output_cat = cat_pgs_outputs(prscsx_output_reformat,
                                            my_python,
                                            b38_flag,
                                            ld_panel_type,
                                            meta_flag,
                                            write_posterior_samples,
                                            mcmc_burnin,
                                            mcmc_thinning_factor)
        
        // make_pgs_boxplots process
        // define plotting script path
        pgs_violinplots_script = "${moduleDir}/scripts/make_pgs_violinplots.py"
        /*
        map: extract concatenated output file paths from tuple
        collect: collect all output paths into one channel for plotting
        */
        
        all_summary_files = prscsx_output_cat.map {
            score_group, cohort, ancestry, pheno, path -> path
        }.collect()
        // make phenotype_score_group list
        pheno_score_group_list = sumstats_summary.map { score_group, cohort, ancestry, pheno, sample_size, sumstats -> pheno + '.' + score_group }.unique().collect()
        // call process
        violinplots = make_pgs_violinplots(all_summary_files,
                                            pgs_violinplots_script,
                                            ancestry_list,
                                            pheno_score_group_list,
                                            my_python,
                                            write_posterior_samples,
                                            meta_flag)

        // reformat prscsx output for emit
        prscsx_weights = prscsx_output_cat.map { score_group, cohort, ancestry, pheno, file -> new Tuple(cohort, ancestry, pheno, file) }
        
        json_params = dump_params_to_json(params)

    emit:
        prscsx_weights

}

process make_pgs_input_and_bim {
    publishDir "${launchDir}/${cohort}/"
    memory = '65 GB'
    machineType 'n2-standard-16'
    input:
        tuple val(score_group), val(cohort), val(ancestry), val(pheno), path(sumstats)
        val(chr_colname)
        val(pos_colname)
        val(a1_colname)
        val(a2_colname)
        val(quant_beta_or_colname)
        val(bin_beta_or_colname)
        val(pval_se_colname)
        val(genome_build)
        val(ld_panel_type)
        val(my_python)
    output:
        tuple val(score_group), val(pheno), path("${ancestry}.${pheno}.${score_group}.PRScsx_input.txt")
        tuple val(score_group), val(pheno), path("${ancestry}.${pheno}.${score_group}.bim_file.bim")
    script:
        """
        #! ${my_python}

        # import pandas packages
        import pandas as pd
        import sys

        # import sumstats and snp info file
        sumstats = pd.read_csv("${sumstats}", sep = None, engine = 'python', dtype = {'${pval_se_colname}' : str})
        snp = pd.read_csv("/app/snpinfo_mult_${ld_panel_type}_hm3_map_b37_b38", sep = '\t')

        # clean sumstats
        ## remove insertions and deletions- I need to add other annotations for insertions and deletions as well
        sumstats = sumstats[sumstats['${a1_colname}'] != 'I']
        sumstats = sumstats[sumstats['${a1_colname}'] != 'D']
        sumstats = sumstats[sumstats['${a1_colname}'] != '-']
        sumstats = sumstats[sumstats['${a1_colname}'] != 'i']
        sumstats = sumstats[sumstats['${a1_colname}'] != 'd']

        ## capitalize alleles
        sumstats['${a1_colname}'] = sumstats['${a1_colname}'].str.upper()
        sumstats['${a2_colname}'] = sumstats['${a2_colname}'].str.upper()

        ## remove missing values from all important columns except beta/OR (that will be done later)
        sumstats.dropna(subset = ['${chr_colname}',
                                    '${pos_colname}',
                                    '${a1_colname}',
                                    '${a2_colname}',
                                    '${pval_se_colname}'], inplace = True)

        ## remove 'chr' from chromosome column in sumstats if it is there
        if sumstats['${chr_colname}'].dtype == 'object':
            print('chr character is present')
            sumstats['${chr_colname}'] = sumstats['${chr_colname}'].str.replace('chr', '')
        else:
            print('chr character is not present')

        ## convert chromosome and position columns to int
        sumstats['${chr_colname}'] = sumstats['${chr_colname}'].astype(int)
        sumstats['${pos_colname}'] = sumstats['${pos_colname}'].astype(int)

        # add CHR_BP column to both files
        ## sumstats
        sumstats['${genome_build}_CHR_BP'] = sumstats['${chr_colname}'].astype(str) + ':' + sumstats['${pos_colname}'].astype(str)

        # subset snp file
        snp = snp[['SNPINFO_ID', '${genome_build}_CHR_BP']]

        # merge sumstats and SNP file
        merge = sumstats.merge(snp, how = 'inner', on = '${genome_build}_CHR_BP')

        # create bim file
        bim = merge.copy()

        # add zero column to bim file
        bim['zero'] = 0

        # drop extra columns in bim file
        bim = bim[['${chr_colname}',
                    'SNPINFO_ID',
                    'zero',
                    '${pos_colname}',
                    '${a1_colname}',
                    '${a2_colname}']]

        # sort bim file
        bim.sort_values(by = ['${chr_colname}', '${pos_colname}'], inplace = True)

        # export bim file
        bim.to_csv("${ancestry}.${pheno}.${score_group}.bim_file.bim", sep = '\t', index = False, header = False)

        # check for beta vs. odds ratio and do some checks
        double_checked_or_beta = False

        if '${bin_beta_or_colname}' in merge.columns:
            print('using binary column')
            beta_vs_or = 'OR' if merge['${bin_beta_or_colname}'].min() >= 0 else 'BETA'

            if beta_vs_or == 'OR':
                if 'o' in '${bin_beta_or_colname}'.lower():
                    double_checked_or_beta = True
            elif beta_vs_or == 'BETA':
                if 'b' in '${bin_beta_or_colname}'.lower():
                    double_checked_or_beta = True

            merge.rename(columns={'${bin_beta_or_colname}' : beta_vs_or}, inplace = True)

        elif '${quant_beta_or_colname}' in merge.columns:
            print('using quantitative column')
            beta_vs_or = 'OR' if merge['${quant_beta_or_colname}'].min() >= 0 else 'BETA'

            if beta_vs_or == 'OR':
                if 'o' in '${quant_beta_or_colname}'.lower():
                    double_checked_or_beta = True
            elif beta_vs_or == 'BETA':
                if 'b' in '${quant_beta_or_colname}'.lower():
                    double_checked_or_beta = True

            merge.rename(columns = {'${quant_beta_or_colname}' : beta_vs_or}, inplace = True)

        else:
            sys.exit("BETA or OR columns specified were not found in file")
        print(beta_vs_or)

        if not double_checked_or_beta:
            sys.exit("Column renaming - could not detect OR vs BETA")

        # check for p-value vs. standard error
        merge['P_SE_CHECK'] = merge['${pval_se_colname}'].astype(float)
        p_vs_se = 'P' if merge['P_SE_CHECK'].max() <= 1 else 'SE'
        merge.drop(columns = ['P_SE_CHECK'], inplace = True)

        # also do a little bit of regex for protection
        double_checked_p_se = False
        if p_vs_se == 'P':
            if 'p' in '${pval_se_colname}'.lower():
                double_checked_p_se = True
        elif p_vs_se == 'SE':
            if 's' in '${pval_se_colname}'.lower():
                double_checked_p_se = True
        else:
            sys.exit("P or SE columns specified were not found in file")
        print(p_vs_se)

        if not double_checked_p_se:
            sys.exit("Column renaming - could not detect P vs SE")

        # drop extra columns in pgs input
        merge = merge[['SNPINFO_ID',
                        '${a1_colname}',
                        '${a2_colname}',
                        beta_vs_or,
                        '${pval_se_colname}']]

        # drop rows with missing data in beta/OR column
        merge.dropna(subset = [beta_vs_or], inplace = True)

        # rename columns in pgs input
        merge.rename(columns = {'SNPINFO_ID' : 'SNP',
                                '${a1_colname}' : 'A1',
                                '${a2_colname}' : 'A2',
                                '${pval_se_colname}' : p_vs_se}, inplace = True)

        # export pgs input file
        merge.to_csv("${ancestry}.${pheno}.${score_group}.PRScsx_input.txt", sep = '\t', index = False)
        """

    stub:
        """
        touch ${ancestry}.${pheno}.${score_group}.PRScsx_input.txt
        touch ${ancestry}.${pheno}.${score_group}.bim_file.bim
        """
}

process cat_bim {
    publishDir "${launchDir}/bim_files/"
    machineType 'n2-standard-4'
    input:
        tuple val(score_group), val(pheno), path(all_bim_files)
        val(my_python)
    output:
        tuple val(score_group), val(pheno), path("${score_group}.${pheno}.bim_file.bim")
    script:
        """
        #! ${my_python}

        # concatenate all bim files into one
        import pandas as pd

        # reformat all_bim_files channel object into a list that python can loop through
        string = "${all_bim_files}"
        print(string)
        list = string.split()
        print(list)

        # create empty list in which files that need to be concatenated can be added to
        dfs = []

        # in a loop, read in each bim file and add to the empty list
        for f in list:
            print(f)
            dfs.append(pd.read_table(f, sep = '\t', header = None))

        # concatenate all the dataframes in the list (all bim files)
        df = pd.concat(dfs)
        print(df)

        # drop duplicate variant IDs
        df.drop_duplicates(keep = 'first',inplace = True)

        # sort concatenated dataframe by chromosome and then position
        df.sort_values(by = [0, 3], inplace = True)

        # export concatenated dataframe
        df.to_csv("${score_group}.${pheno}.bim_file.bim", sep = '\t', index = False, header = False)
        """

    stub:
        """
        touch ${score_group}.${pheno}.bim_file.bim
        """
}
process prscsx {
    publishDir "${launchDir}/PRScsx_output/chromosome_separated_outputs/"
    machineType 'n2-standard-4'

    input:
        tuple val(score_group), val(pheno), val(sample_sizes), path(pgs_input_files), path(bim_output_file), val(chromosome)
        val(ancestry_list)
        path(ld_directory)
        val(param_a)
        val(param_b)
        val(param_phi)
        val(mcmc_iterations)
        val(mcmc_burnin)
        val(mcmc_thinning_factor)
        val(meta_flag)
        val(write_posterior_samples)
        val(seed)
        val(my_python)
        val(my_prscsx)
    output:
        tuple val(score_group), val(pheno), val(chromosome), path("${score_group}.${pheno}.PRScsx_${ancestry_list.size() < 2 ? ancestry_list.join(',') : '{' + ancestry_list.join(',') + '}'}_pst_eff_a${params.param_a}_b${params.param_b}_phi${params.param_phi}_chr${chromosome}.txt")
    shell:
    """
    # makes ancestry list for PRScsx from tuple so the order is correct
    touch temp1
    echo ${pgs_input_files} | sed 's/[][]//g' | tr , '\n' | grep -v 'META' | sed 's/[.].*//g' | tr '\n' , | sed 's/ //g' | sed 's/,\$//g' >> temp1
    ancestry_list=\$(cat "temp1")
    echo \$ancestry_list
    rm temp1

    # makes prscsx input list for PRScsx from tuple so the order is correct
    touch temp2
    echo ${pgs_input_files} | sed 's/[][]//g' | sed 's/ //g' >> temp2
    pgs_input=\$(cat "temp2")
    echo \$pgs_input
    rm temp2

    # makes sample size list for PRScsx from tuple so the order is correct
    touch temp3
    echo ${sample_sizes} | sed 's/[][]//g' | sed 's/ //g'  >> temp3
    sample_size=\$(cat "temp3")
    echo \$sample_size
    rm temp3

    # PRScsx command if phi parameter is auto (not specified)
    if [[ '${param_phi}' == 'auto' ]]; then
        ${my_python} ${my_prscsx} \
        --ref_dir=${ld_directory} \
        --bim_prefix=${score_group}.${pheno}.bim_file \
        --sst_file=\$pgs_input \
        --n_gwas=\$sample_size \
        --pop=\$ancestry_list \
        --out_dir=. \
        --out_name=${score_group}.${pheno}.PRScsx \
        --a=${param_a} \
        --b=${param_b} \
        --n_iter=${mcmc_iterations} \
        --n_burnin=${mcmc_burnin} \
        --thin=${mcmc_thinning_factor} \
        --chrom=${chromosome} \
        --meta=${meta_flag} \
        --write_pst=${write_posterior_samples} \
        --seed=${seed}
    # PRScsx command if phi parameter is specified
    else
        ${my_python} ${my_prscsx} \
        --ref_dir=${ld_directory} \
        --bim_prefix=${score_group}.${pheno}.bim_file \
        --sst_file=\$pgs_input \
        --n_gwas=\$sample_size \
        --pop=\$ancestry_list \
        --out_dir=. \
        --out_name=${score_group}.${pheno}.PRScsx \
        --a=${param_a} \
        --b=${param_b} \
        --phi=${param_phi} \
        --n_iter=${mcmc_iterations} \
        --n_burnin=${mcmc_burnin} \
        --thin=${mcmc_thinning_factor} \
        --chrom=${chromosome} \
        --meta=${meta_flag} \
        --write_pst=${write_posterior_samples} \
        --seed=${seed}
    fi
    """

    stub:
        """
        # reformat ancestry list
        echo ${ancestry_list} | sed 's/[][]//g' | sed 's/,//g' > temp
        ancestry_list=\$(cat "temp")
        echo \$ancestry_list
        echo \$ancestry_list | tr ' ' '\n' | wc -l > num_ancestries_file

        # output files
        for ancestry in \$ancestry_list
        do
        filename='${score_group}.${pheno}.PRScsx_'\$ancestry'_pst_eff_a${param_a}_b${param_b}_phi${param_phi}_chr${chromosome}.txt'
        touch \$filename
        done
        """
}

process cat_pgs_outputs {
    publishDir "${launchDir}/PRScsx_output/concatenated_outputs/",mode:'copy'
    machineType 'n2-standard-4'
    input:
        tuple val(score_group), val(ancestry), val(pheno), file(prscsx_output_files), val(cohort)
        val(my_python)
        val(b38_flag)
        val(ld_panel_type)
        val(meta_flag)
        val(write_posterior_samples)
        val(mcmc_burnin)
        val(mcmc_thinning_factor)
    output:
        tuple val(score_group), val(cohort), val(ancestry), val(pheno), path("${ancestry}.${pheno}.${score_group}.PRScsx_output.txt")
    script:
    """
    #! ${my_python}
    # import pandas
    import pandas as pd

    # reformat prscsx_output_files channel object into a list that python can loop through
    # save prscsx output files into a string
    string = "${prscsx_output_files}"
    # create a list of prscsx output files from that string
    list = string.split()
    print(list)

    # create empty list in which files that need to be concatenated can be added to
    dfs = []

    # loop through prscsx output files, read them in and append them to the empty list
    for f in list:
        dfs.append(pd.read_table(f, sep = '\t', header = None))

    # concatenate all prscsx output files for each ancestry_pheno pair (output files are chromosome separated)
    df = pd.concat(dfs)

    # sort concatenated file by chromosome and then position
    df.sort_values(by = [0, 2], inplace = True)
    print(df)

    # check if write pst samples == True
    if ("${write_posterior_samples}".lower() == "true") and ("${meta_flag}".lower() == "true") and ("${ancestry}" == "META"):
        # add column names except PGS
        df = df.rename(columns = {0 : 'CHR', 1 : 'RSID', 2 : 'POS', 3 : 'A1', 4 : 'A2'})

        # identify number of PGS columns
        num_pgs_cols = int("${mcmc_burnin}") / int("${mcmc_thinning_factor}") + 5

        # loop through PGS columns and rename by index
        pgs_col_list = []
        for index, col in enumerate(range(5, int(num_pgs_cols))):
            colnum = index + 1
            colname = 'PGS.' + str(colnum)
            pgs_col_list.append(colname)
            df = df.rename(columns = {col : colname})
    else:
        # add column names with one score column name
        df = df.rename(columns = {0 : 'CHR', 1 : 'RSID', 2 : 'POS', 3 : 'A1', 4 : 'A2', 5 : 'PGS'})

        # subset columns
        df = df[['CHR', 'RSID', 'POS', 'A1', 'A2', 'PGS']]

    # format to b38 if flag is true
    if "${b38_flag}" == "true":
        # read in snp info file
        snp = pd.read_csv("/app/snpinfo_mult_${ld_panel_type}_hm3_map_b37_b38", sep = '\t')

        # rename columns in snp file
        snp.rename(columns = {'SNPINFO_ID' : 'RSID'}, inplace = True)

        # split b38 chr:pos column in snp file
        snp[['CHR', 'POS']] = snp['b38_CHR_BP'].str.split(':', expand = True)

        # drop b37 id in snp file
        snp.drop(columns = ['b37_CHR_BP'], inplace = True)

        # drop chr and pos columns in pgs df
        df.drop(columns = ['CHR', 'POS'], inplace = True)

        # merge snp file and pgs output
        df = df.merge(snp, on = 'RSID', how = 'inner')

        if ("${write_posterior_samples}".lower() == "true") and ("${meta_flag}".lower() == "true") and ("${ancestry}" == "META"):
            # create column list
            all_col_list = ['CHR', 'RSID', 'POS', 'A1', 'A2'] + pgs_col_list

            # reorder columns with multiple PGS
            df = df[all_col_list]
        else:
            # reorder columns with one PGS
            df = df[['CHR', 'RSID', 'POS', 'A1', 'A2', 'PGS']]

    # export concatenated file
    df.to_csv("${ancestry}.${pheno}.${score_group}.PRScsx_output.txt", sep = '\t', index = False)
    """

    stub:
    """
    touch ${ancestry}.${pheno}.${score_group}.PRScsx_output.txt
    """
}

process make_pgs_violinplots {
    publishDir "${launchDir}/Plots/",mode:'copy'
    machineType 'n2-standard-4'
    input:
        path(summary_files)
        path(pgs_violinplots_script)
        val(ancestry_list)
        val(pheno_score_group_list)
        val(my_python)
        val(meta_flag)
        val(write_posterior_samples)

    output:
        path('combined_PGS_violinplot.png')

    shell:
        """
        echo "Triggered once after all files are concatenated"
        # call plotting script
        ${my_python} ${pgs_violinplots_script} \
        --input_file_list '${summary_files}' \
        --ancestry_list '${ancestry_list}' \
        --pheno_score_group_list '${pheno_score_group_list}' \
        --write_pst_samples '${write_posterior_samples}' \
        --meta '${meta_flag}'
        """

    stub:
        '''
        touch 'combined_PGS_violinplot.png'
        '''
}

import groovy.json.JsonBuilder
process dump_params_to_json {
    publishDir "${launchDir}/Summary", mode: 'copy'
    machineType 'n2-standard-2'

    input:
        val params_dict
    output:
        path('prscsx_params.json')
    shell:
        """
        echo '${new JsonBuilder(params_dict).toPrettyString().replace(';', '|')}' > prscsx_params.json
        """
}

