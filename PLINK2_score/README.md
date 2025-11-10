# Apply PGS Nextflow Pipeline
The Apply PGS Nextflow pipeline utilizes PLINK2 <code>--score</code> function (linear scoring) to generate individual-level polygenic scores (PGS). This works by multiplying polygenic score weights, or variant-level scores, by individual genotypes, and then sums and averages these values for each individual. This pipeline makes inputs for PLINK2 <code>--score</code>, runs the tool, and concatenates chromosome-separated outputs.

## Files in this Directory
* <code>README.md</code>: User-side documentation
* <code>Docker_Singularity_README.md</code>: Information on docker and singularity containers
* <code>Dockerfile</code>: Dockerfile used to build the docker image
* <code>apply_pgs.sif</code>: Singularity image
  * NOTE: Singularity image might need to be rebuilt after cloning the github. Errors have arisen from pushing and pulling the image to/from github
* <code>apply_pgs.nf</code>: Workflow code
* <code>apply_pgs.config</code>: Configuration file to modify with personalized filepaths and parameters
* <code>nextflow.config</code>: Configuration file for compute profiles

## Software Requirements
* Nextflow version 23.04.1.5866
* Singularity version 3.8.3 OR Docker version 4.30.0

## Basic Commands to Run This Pipeline
* Build a singularity container from the docker image: <code>singularity build apply_pgs.sif docker://pennbiobank/apply_pgs:latest</code>
  * A docker or singualrity container is used to provide the necessarily software for this pipeline
  * The docker image will be directly pulled from dockerhub automatically whereas the singularity container needs to be built manually before running the pipeline
  * See [Docker_Singularity_README.md](https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench/blob/main/Apply_PGS/Docker_Singularity_README.md) for more information
* Run nextflow pipeline: <code>nextflow run apply_pgs.nf -resume -profile <standard/cluster/all_of_us></code>
  * <code>-resume</code> flag picks up the workflow where it left off, otherwise, the workflow will rerun from the beginning
  * <code>-profile</code> selects the compute profiles we set up in <code>nextflow.config</code>
  * <code>-profile standard</code> uses the docker image and executes the processes locally
  * <code>-profile cluster</code> uses the singularity container and submits processes to a queue- optimal for HPC or LPC computing systems
  * <code>-profile all_of_us</code> uses the docker image and executes pipelines on the All of Us Researcher Workbench
  * For more information, visit [Nextflow's Documentation](https://www.nextflow.io/docs/latest/cli.html)

## Inputs
* PGS weights
* Validation population <code>plink1</code> or <code>plink2</code> files

## Outputs
* Chromosome/validation population separated output files (raw output files) will be symlinked in <code>Apply_PGS_Output/</code> (each score is a column in the file)
* Validation population separated output files (chromosome separated outputs concatenated and sorted) will be symlinked in in <code>Summary/</code> (each score is a column in the file)

## Configuration
* This pipeline requires specification of the following parameters in <code>apply_pgs.config</code>:
  * List of chromosomes
  * Full file path to directory containing PGS weights
  * Path to descriptor file containing columns with cohort nickname, ancestry, phenotype, GWAS sample size, and summary statistics filename
    * Each row represents a different GWAS for PRScsx to compute a score based on
    * The file must contain the following columns:
      * Cohort: a nickname for the cohort of interest
      * Ancestry: The PGS genetically inferred ancestry group
      * Phenotype: The PGS phenotype
      * Filename: Name of the PGS weights filename
    * Example File:
      * ![Title](images/Apply_PGS_Descriptor_File_Example.png) 
  * Groovy dictornary specifying column names of of the descriptor file
    * Necessary columns:
      * Cohort
      * Ancestry
      * Phenotype
      * PGS weights filename
  * Delimiter of the descriptor file
  * Variant ID format in PGS weights files
    * This method assumes that all PGS weights files have the same variant ID format
    * Modeled based on PLINK formats
    * Example formats:     
      * <code>RSID</code> = RSID (ex: rs6893237)
      * <code>@:#:$r:$a</code> = chr:pos:ref:alt (ex: 1:100000:G:A)
      * <code>chr@:#:$r:$a</code>  = chr:pos:ref:alt (ex: chr1:100000:G:A)
      * <code>@:#</code> = chr:pos (ex: 1:100000)
      * <code>chr@:#</code> = chr:pos (ex: chr1:100000)
  * Delimeter in score file
    * This method assumes that all PGS weights files have the same delimiter
    * Add delimiter descriptions !!
  * Groovy dictionary specifying validation population information, including:
    * Validation population nickname
    * PLINK file prefix (before chromosome number)
      * PLINK files MUST be chromosome separated
    * PLINK file suffix (after chromosome number and before plink suffix (ex: <code>.bim</code>)
    * Validation population sample list
    * Variant ID format
      * Example formats described above
  * <code>plink --score</code> parameters
    * Parameters are described in <code>apply_pgs.config</code>. For more detailed descriptions, visit [PLINK website](https://www.cog-genomics.org/plink/2.0/score)
      * <code>dosage_transformation</code>
      * <code>xchr_model</code>
      * <code>read_freq</code>
      * <code>no_mean_imputation</code>
      * <code>independent_se</code>


## Workflow Processes
1. <code>variant_id_map</code>
  * Checks that PGS weights and validation population bim or pvar files have the same variant ID
  * Parallelizes by <code>validation_population</code>
2. <code>make_apply_pgs_input</code>
  * Reformats PGS weights into format that <code>plink --score</code> requires 
  * Parallelizes by <code>score</code>
3. <code>join_score_files</code>
  * Combines score-separated input files into one file (<code>plink --score</code> can compute many scores at once)
  * No parallelization
4. <code>compute_scores</code>
  * Uses <code>plink --score</code> to multiply variant-level scores by individual genotypes and sums all products for a given individual
  * Parallelizes by <code>validation_population</code> and <code>chromosome</code>
5. <code>concatenate_plink_score_outputs</code>
  * Concatenates chromosome separated <code>plink --score</code> outputs by <code>validation_population</code>
  * Parallelizes by <code>validation_population</code>

## Support
Please direct any problems or questions to Katie Cardone ([Katie.Cardone@pennmedicine.upenn.edu](mailto:Katie.Cardone@pennmedicine.upenn.edu))
