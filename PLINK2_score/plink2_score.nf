#! /appl/nextflow-23.04.1.5866/bin/nextflow

// enable DSL2 syntax***
nextflow.enable.dsl = 2

// log info
log.info """\
NEXTFLOW - DSL2 - Apply PGS - P I P E L I N E
============================================
run_as                       :  ${workflow.commandLine}
run_location                 :  ${launchDir}
started_at                   :  ${workflow.start}
container                    :  ${workflow.containerEngine}:${workflow.container}

Chromosomes
============================================
chromosomes                  :  ${params.chromosome_list}

Polygenic Score File Information
============================================
score variant ID format      :  ${params.score_variant_id_format}
chromosome column name       :  ${params.score_file_colnames['chr_colname']}
position column name         :  ${params.score_file_colnames['pos_colname']}
A1 column name               :  ${params.score_file_colnames['a1_colname']}
A2 column name               :  ${params.score_file_colnames['a2_colname']}
variant ID column name       :  ${params.score_file_colnames['variant_id_colname']}
polygenic score column name  :  ${params.score_file_colnames['pgs_colname']}

PLINK 2 Score Parameters
============================================
dosage transformation        :  ${params.dosage_transformation}
chrX model                   :  ${params.xchr_model}
allele frequency file        :  ${params.read_freq}
no mean imputatation         :  ${params.no_mean_imputation}
independent standard error   :  ${params.independent_se}

"""
.stripIndent()

workflow {
    /*
    this section is called first if:
        Apply PGS is the first workflow called
        Apply PGS is called independently of other workflows
        Another PRScsx workflow is not being stitched to Apply PGS workflow
        Channels from another workflow are not inputs to the Apply PGS workflow
    */
    score_weights = PLINK2_SCORE_setup()
    indiv_polygenic_scores = PLINK2_SCORE(score_weights)
}

workflow PLINK2_SCORE_setup {
    main:
        // make summary table file name into a string type
        descriptor_file = new File(params.input_descriptor_table_filename.toString())
        // reads summary table file into groovy making it into a list
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
        // creates a channel out of the summary table list
        descriptor_channel = Channel.fromList(descriptor_tuple)
        // parse header and retrieve indices of essential columns based on their column names
        descriptor_header = descriptor_file.readLines().get(0).split(',')
        descriptor_header.eachWithIndex { item, index ->
            if (item == params['input_descriptor_table_colnames']['cohort_colname']) {
                cohort_index = index
            }
            if (item == params['input_descriptor_table_colnames']['ancestry_colname']) {
                ancestry_index = index
            }
            if (item == params['input_descriptor_table_colnames']['phenotype_colname']) {
                pheno_index = index
            }
            if (item == params['input_descriptor_table_colnames']['pgs_weights_full_filpath_colname']) {
                score_file_index = index
            }
        }
        // find out indices
        // makes each item (table row) in the list a separate tuple
        if (descriptor_file.text.tokenize('\n').size() <= 2) {
            descriptor_channel_collect = descriptor_channel.collect()
            score_weights = descriptor_channel_collect.map { arr -> new Tuple(
                arr[cohort_index],
                arr[ancestry_index],
                arr[pheno_index],
                arr[score_file_index]) }
        } else {
            score_weights = descriptor_channel.map { arr -> new Tuple(
                arr[cohort_index],
                arr[ancestry_index],
                arr[pheno_index],
                arr[score_file_index]) }
        }
    emit:
        score_weights
}
workflow PLINK2_SCORE {
    take:
        score_weights
    main:
        // create validation cohort list
        validation_population_list = params.validation_populations.keySet().toList()

        // create validation_population, and chromomosome channels to parallelize by
        // could we pass chromosome from prscsx.nf??
        validation_population = Channel.fromList(validation_population_list)
        chromosome = Channel.fromList(params.chromosome_list)

        // apply PGS input process
        /*
        combine: combine val pop with chromosome
        map: add plink flag and bim/pvar file
        */
        val_pop_variant_files = validation_population.combine(chromosome).map { val_pop, chr ->
        new Tuple(
            val_pop,
            params.validation_populations[val_pop]['plink_file_flag'],
            params.validation_populations[val_pop]['plink_file_flag'] == '--bfile' ? "${params.validation_populations[val_pop]['plink_prefix']}${chr}${params.validation_populations[val_pop]['plink_suffix']}.bim" : "${params.validation_populations[val_pop]['plink_prefix']}${chr}${params.validation_populations[val_pop]['plink_suffix']}.pvar"
        ) }.groupTuple(by:[0, 1])
        // combine score weights and val pop variant files tuples for parallelization
        score_weights_val_pop_variant_files = score_weights.combine(val_pop_variant_files)
        /*
        map: create variant_id_format tuple WITH validation population variable
        collect: collect val_pop_variant_id tuples into one tuple (so the process doesn't parallelize by validation population)
        */
        val_pop_variant_ids = validation_population.map { val_pop ->
        new Tuple(val_pop, params.validation_populations[val_pop]['variant_id_format']) }.collect()
        // make column names channels
        score_chr_colname="${params['score_file_colnames']['chr_colname']}"
        score_pos_colname="${params['score_file_colnames']['pos_colname']}"
        score_a1_colname="${params['score_file_colnames']['a1_colname']}"
        score_a2_colname="${params['score_file_colnames']['a2_colname']}"
        score_id_colname="${params['score_file_colnames']['variant_id_colname']}"
        score_pgs_colname="${params['score_file_colnames']['pgs_colname']}"
        // make script a channel
        apply_pgs_input_script = "${moduleDir}/scripts/apply_pgs_input.py"
        // define python path channel
        my_python = "${params.my_python}"
        // define score variant ID type channel
        score_variant_id_format = "${params.score_variant_id_format}"
        // call process
        apply_pgs_input = make_apply_pgs_input(score_weights_val_pop_variant_files,
                                                val_pop_variant_ids,
                                                validation_population_list,
                                                apply_pgs_input_script,
                                                score_chr_colname,
                                                score_pos_colname,
                                                score_a1_colname,
                                                score_a2_colname,
                                                score_id_colname,
                                                score_pgs_colname,
                                                my_python,
                                                score_variant_id_format)

        // join score files process
        // collect separate apply PGS input files into one channel
        all_apply_pgs_inputs = apply_pgs_input.collect()
        // call process
        combine_apply_pgs_input = join_score_files(all_apply_pgs_inputs, my_python)

        // clean_sample_list process
        // create tuple with validation population, sample list, and fam/pvar file
        val_pop_sample_list = validation_population.combine(chromosome).map { val_pop, chr ->
        new Tuple(
            val_pop,
            params.validation_populations[val_pop]['population_subset_file'],
            params.validation_populations[val_pop]['population_subset_file_id_col'],
            params.validation_populations[val_pop]['population_subset_file_delim'],
            params.validation_populations[val_pop]['plink_file_flag'] == '--bfile' ? "${params.validation_populations[val_pop]['plink_prefix']}${chr}${params.validation_populations[val_pop]['plink_suffix']}.fam" : "${params.validation_populations[val_pop]['plink_prefix']}${chr}${params.validation_populations[val_pop]['plink_suffix']}.psam"
        ) }.groupTuple(by: [0, 1, 2])

        // call process
        cleaned_sample_list = clean_sample_list(val_pop_sample_list, my_python)

        // compute_scores process
        // make tuple with validation population, chromosome, plink flag, plink files, plink file prefix, and sample list
        plink_files = cleaned_sample_list.combine(chromosome).map { val_pop, sample_list, chr ->
        new Tuple(val_pop, chr,
        params.validation_populations[val_pop]['plink_file_flag'],
        new Tuple(
            *(params.validation_populations[val_pop]['plink_file_flag'] == '--bfile' ? ['.bim', '.bed', '.fam'] : ['.pvar', '.pgen', '.psam']).collect {
            ext -> "${params.validation_populations[val_pop]['plink_prefix']}${chr}${params.validation_populations[val_pop]['plink_suffix']}${ext}"
            } ),
        sample_list
        ) }

        // create --read-freq plink line
        plink_read_freq_line = params['read_freq'] == 'auto' ? '' : "--read-freq ${params.read_freq}"
        // make channels from other PLINK score parameters
        dosage_transformation="${params.dosage_transformation}"
        xchr_model="${params.xchr_model}"
        no_mean_imputation="${params.no_mean_imputation}"
        independent_se="${params.independent_se}"
        // define my_plink2 executable channel
        my_plink2 = "${params.my_plink2}"
        // call process
        plink_score_output = compute_scores(plink_files,
                                            combine_apply_pgs_input,
                                            plink_read_freq_line,
                                            dosage_transformation,
                                            xchr_model,
                                            no_mean_imputation,
                                            independent_se,
                                            my_plink2)

        // concatenate_plink_score_outputs
        // group chromosome separated files together
        plink_score_output_grouped = plink_score_output[0].groupTuple(by: 0, size: params.chromosome_list.size())
        // make cat outputs script channel
        cat_outputs_script = "${moduleDir}/scripts/cat_plink_score_outputs.py"
        // call process
        cat_scores = concatenate_plink_score_outputs(plink_score_output_grouped,
                                                        cat_outputs_script,
                                                        my_python)

        /*
        make distribution plots process
        output from concatenate_plink_score_outputs
        Need: combined scores, population, summary weights file
        */

        // define script
        summary_plots_script = "${moduleDir}/scripts/pgs_plotting.py"
        // make ancestry and phenotype lists
        ancestry_list = score_weights.map { cohort, ancestry, pheno, score_file -> ancestry }.unique().collect()
        pheno_list = score_weights.map { cohort, ancestry, pheno, score_file -> pheno }.unique().collect()
        // call process
        plots = make_summary_plots(summary_plots_script,
                                    cat_scores,
                                    ancestry_list,
                                    pheno_list,
                                    my_python)

        // reformat output for emit
        indiv_polygenic_scores = cat_scores

        json_params = dump_params_to_json(params)

    emit:
       indiv_polygenic_scores
}

process make_apply_pgs_input {
    publishDir "${launchDir}/${cohort}/"
    memory = '64 GB'
    machineType 'n2-standard-16'
    input:
        tuple val(cohort), val(ancestry), val(pheno), path(score_weights_file), val(validation_population), val(plink_flag), path(bim_pvar_files)
        val(val_pop_variant_id_collect)
        val(validation_population_list)
        path(apply_pgs_input_script)
        val(score_chr_colname)
        val(score_pos_colname)
        val(score_a1_colname)
        val(score_a2_colname)
        val(score_id_colname)
        val(score_pgs_colname)
        val(my_python)
        val(score_variant_id_format)
    output:
        path("${ancestry}.${pheno}.${validation_population}.Apply_PGS_Input.txt")
    shell:
        """
        # reformat plink flag
        echo '${plink_flag}' | sed 's/--//g' > temp
        plink_flag_reformatted=\$(cat temp)
        echo \$plink_flag_reformatted
        rm temp

        ${my_python} ${apply_pgs_input_script} \
        --valPop '${validation_population}' \
        --ancestry '${ancestry}' \
        --pheno '${pheno}' \
        --plinkFlag \$plink_flag_reformatted \
        --varFiles '${bim_pvar_files}' \
        --popsIds '${val_pop_variant_id_collect}' \
        --scoreIdFormat '${score_variant_id_format}' \
        --scoreFile '${score_weights_file}' \
        --scoreChrCol '${score_chr_colname}' \
        --scorePosCol '${score_pos_colname}' \
        --scoreIdCol '${score_id_colname}' \
        --scoreA1Col '${score_a1_colname}' \
        --scoreA2Col '${score_a2_colname}' \
        --scorePGSCol '${score_pgs_colname}' \
        """
    stub:
        """
        touch ${ancestry}.${pheno}.${validation_population}.Apply_PGS_Input.txt
        """
}

process join_score_files {
    publishDir "${launchDir}/"
    machineType 'n2-standard-4'
    input:
        path(all_apply_pgs_input_files)
        val(my_python)
    output:
        path('combined_score_file_for_plink.txt')
    script:
        """
        #! ${my_python}

        # load packages
        import pandas as pd

        # create an empty list for dataframes to be added to
        dfs = []

        # reformat list of files so python can loop through it
        input_file_string = '${all_apply_pgs_input_files}'
        input_file_list = list(input_file_string.split(" "))
        print(input_file_list)

        # in a loop, read in all of the PGS input files and add them to a list of dataframes
        for f in input_file_list:
            print(f)
            temp_df = pd.read_table(f, sep='\t', index_col=['ID', 'ALLELE'])
            temp_df = temp_df[~temp_df.index.duplicated()]
            dfs.append(temp_df)

        # concatenate the COLUMNS dataframes in the list (combining separate apply PGS input files)
        # this makes a score file with one ID column, one allele column, and one score column for each training cohort
        df = pd.concat(dfs, axis=1)

        # fill missing values with zeros (if a cohort doesn't have a SNP, we can assume the score for that SNP in that cohort is zero)
        df = df.fillna(0)

        # export dataframe
        df.to_csv('combined_score_file_for_plink.txt', sep='\t')
        """
    stub:
        '''
        touch combined_score_file_for_plink.txt
        '''
}

process clean_sample_list {
    publishDir "${launchDir}/sample_lists/"
    machineType 'n2-standard-4'
    input:
        tuple val(validation_population), path(sample_list), val(id_col), val(sample_delim), path(fam_psam_file_list)
        val(my_python)
    output:
        tuple val(validation_population), path("${validation_population}.fam_format_sample_list.txt")
    script:
        """
        #! ${my_python}

        # load packages
        import pandas as pd

        # reformat list of fam/psam files
        fam_psam_file_string = "${fam_psam_file_list}"
        fam_psam_file_list_reformatted = fam_psam_file_string.split(' ')

        # extract first item in list
        first_fam_psam_file = fam_psam_file_list_reformatted[0]

        # import fam/psam file
        fam_psam = pd.read_csv(first_fam_psam_file, sep = None, engine = 'python', header = None, comment = '#', usecols = [0, 1], names = ['FID', 'IID'], dtype = str)
        print(fam_psam)

        if fam_psam['IID'].isnull().all():
            fam_psam.drop(columns = ['IID'], inplace = True)
            fam_psam['IID'] = fam_psam['FID']
            fam_psam.drop(columns = ['FID'], inplace = True)
            print(fam_psam)

        # check if sample file has a header
        if "${id_col}".isdigit():
            # sample file does not have a header
            ## read in file
            sample=pd.read_csv("${sample_list}", sep = "${sample_delim}", engine = 'python', header = None, dtype = str)
            ## filter fam/psam file to IDs in sample file
            fam_psam_sample=fam_psam[fam_psam['IID'].isin(sample[${id_col}])]
        else:
            # sample file has a header
            ## read in file
            sample=pd.read_csv("${sample_list}", sep = "${sample_delim}", engine = 'python', dtype = str)
            ## filter fam/psam file to IDs in sample file
            fam_psam_sample = fam_psam[fam_psam['IID'].isin(sample["${id_col}"])]

        print(fam_psam_sample)

        # export sample list
        fam_psam_sample.to_csv("${validation_population}.fam_format_sample_list.txt", sep = '\t', index = None, header = None)
        """
    stub:
        """
        touch ${validation_population}.fam_format_sample_list.txt
        """
}

process compute_scores {
    publishDir "${launchDir}/plink_score_output/chromosome_separated_outputs/"
    disk params.max_chr_bed_size
    machineType 'n2-standard-4'
    input:
        tuple val(validation_population), val(chromosome), val(plink_flag), path(plink_files), path(sample_list_file)
        path(combined_apply_pgs_input)
        val(plink_read_freq_line)
        val(dosage_transformation)
        val(xchr_model)
        val(no_mean_imputation)
        val(independent_se)
        val(my_plink2)
    output:
        tuple val(validation_population), path("${validation_population}.Apply_PGS_Output.chr${chromosome}.sscore"), path("${validation_population}.Apply_PGS_Output.chr${chromosome}.sscore.vars")
        path("${validation_population}.Apply_PGS_Output.chr${chromosome}.log")
    shell:
        """
        # get column numbers of scores in score file
        sed 's/\t/,/g' ${combined_apply_pgs_input} | head -n1 | sed 's/[^,]//g' | wc -c > temp1
        total_colnums=\$(cat "temp1")
        echo \$total_colnums
        seq 3 \$total_colnums | tr '\n' , | sed 's/,\$//g' > temp2
        plink_score_colnums=\$(cat "temp2")
        echo \$plink_score_colnums
        rm temp1
        rm temp2

        # get plink prefix
        echo ${plink_files} | sed 's/[][]//g' | cut -d',' -f1 | sed -e 's/.bed//g' -e 's/.bim//g' -e 's/.fam//g' -e 's/.psam//g' -e 's/.pvar//g' -e 's/.pgen//g' | sed 's/,//g' | sed 's/ //g' > temp
        plink_prefix=\$(cat "temp")

        # plink score command
            # first line includes score input, and opts to read the score input header
            # first line also includes dosage transformation, no mean imputation, and independent SE parameters
            # first line also specifies the output files to include score sums, not score averages
            # finally, first line directs creation of an output file with the variants plink score was computed with
            # second line specifies the X chromosome model parameter
            # third line specifies the column numbers of the score columns
            # fourth line only keeps individuals specified in the sample list
            # fourth line include --read-freq flag is user opted to use it
            # fifth line specifies the input plink files
            # sixth line removes unequal duplicates- duplicate variant IDs were causing an error
            # seventh line specifies the output file format
        ${my_plink2} --score ${combined_apply_pgs_input} header-read ${dosage_transformation} ${no_mean_imputation} ${independent_se} cols=+scoresums,-scoreavgs list-variants \
        --xchr-model ${xchr_model} \
        --score-col-nums \$plink_score_colnums \
        --keep ${sample_list_file} ${plink_read_freq_line} \
        ${plink_flag} \$plink_prefix \
        --rm-dup exclude-mismatch \
        --out ${validation_population}.Apply_PGS_Output.chr${chromosome} > ${validation_population}.Apply_PGS_Output.chr${chromosome}.log
        """
    stub:
        """
        touch ${validation_population}.Apply_PGS_Output.chr${chromosome}.sscore
        touch ${validation_population}.Apply_PGS_Output.chr${chromosome}.sscore.vars
        touch ${validation_population}.Apply_PGS_Output.chr${chromosome}.log
        """
}

process concatenate_plink_score_outputs {
    publishDir "${launchDir}/plink_score_output/concatenated_outputs/", mode:'copy'
    machineType 'n2-standard-4'
    input:
        tuple val(validation_population), path(plink_score_output_files), path(plink_score_var_lists_files)
        path(cat_outputs_script)
        val(my_python)
    output:
        tuple val(validation_population), path("${validation_population}.all_computed_PGS_scores.txt"), path("${validation_population}.all_computed_PGS_scores.variant_list.txt")
    shell:
        """
        ${my_python} ${cat_outputs_script} \
        --valPop '${validation_population}' \
        --score_file_list '${plink_score_output_files}' \
        --variant_list '${plink_score_var_lists_files}'
        """
    stub:
        """
        touch ${validation_population}.all_computed_PGS_scores.txt
        touch ${validation_population}.all_computed_PGS_scores.variant_list.txt
        """
}

process make_summary_plots {
    publishDir "${launchDir}/Plots/",mode:'copy'
    machineType 'n2-standard-4'

    input:
        path(summary_plots_script)
        tuple val(validation_population), path(concatenated_scores), path(variant_list)
        val(ancestry_list)
        val(pheno_list)
        val(my_python)
    output:
        path("${validation_population}.all_computed_PGS_scores.boxplots.*.png")
        path("${validation_population}.all_computed_PGS_scores.densityplots.*.png")
    shell:
        """
        ${my_python} ${summary_plots_script} \
        --data ${concatenated_scores} \
        --population ${validation_population} \
        --ancestry_list '${ancestry_list}' \
        --phenotype_list '${pheno_list}'
        """
    stub:
        """
        touch ${validation_population}.all_computed_PGS_scores.boxplots.stub.png
        touch ${validation_population}.all_computed_PGS_scores.densityplots.stub.png
        """
}

import groovy.json.JsonBuilder
process dump_params_to_json {
    publishDir "${launchDir}/Summary", mode: 'copy'
    machineType 'n2-standard-2'

    input:
        val params_dict
    output:
        path('plink_score_params.json')
    shell:
        """
        echo '${new JsonBuilder(params_dict).toPrettyString().replace(';', '|')}' > plink_score_params.json
        """
}
