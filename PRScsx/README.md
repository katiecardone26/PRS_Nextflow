# PRScsx Nextflow Pipeline
PRScsx is a tool used to calculate polygenic scores. It integrates GWAS summary statistics and LD panels across multiple ancestry groups to improve global polygenic prediction. This nextflow pipeline seemlessly takes GWAS summary statistics, makes the input files for PRScsx, runs PRScsx, concatenates the chromosome-separated outcomes, and produces a summary violin plot, all in a parallelized fashion.

## Files in this Directory
* <code>README.md</code>: User-side documentation
* <code>Docker_Singularity_README.md</code>: Information on docker and singularity containers
* <code>Dockerfile</code>: Dockerfile used to build the docker image
* <code>prscsx.sif</code>: Singularity image
  * NOTE: Singularity image might need to be rebuilt after cloning the github. Errors have arisen from pushing and pulling the image to/from github
* <code>prscsx.nf</code>: Workflow code
* <code>prscsx.config</code>: Configuration file to modify with personalized filepaths and parameters
* <code>nextflow.config</code>: Configuration file for compute profiles

## Software Requirements
* Nextflow version 23.04.1.5866
* Singularity version 3.8.3 OR Docker version 4.30.0

## Basic Commands to Run This Pipeline
* Build a singularity container from the docker image: <code>singularity build prscsx.sif docker://pennbiobank/prscsx:latest</code>
  * A docker or singualrity container is used to provide the necessarily software for this pipeline
  * The docker image will be directly pulled from dockerhub automatically whereas the singularity container needs to be built manually before running the pipeline
  * See [Docker_Singularity_README.md](https://github.com/PMBB-Informatics-and-Genomics/geno_pheno_workbench/blob/main/PRScsx/Docker_Singularity_README.md) for more information
* Run nextflow pipeline: <code>nextflow run prscsx.nf -resume -profile <standard/cluster/all_of_us></code>
  * <code>-resume</code> flag picks up the workflow where it left off, otherwise, the workflow will rerun from the beginning
  * <code>-profile</code> selects the compute profiles we set up in <code>nextflow.config</code>
  * <code>-profile standard</code> uses the docker image and executes the processes locally
  * <code>-profile cluster</code> uses the singularity container and submits processes to a queue- optimal for HPC or LPC computing systems
  * <code>-profile all_of_us</code> uses the docker image and executes pipelines on the All of Us Researcher Workbench
  * For more information, visit [Nextflow's Documentation](https://www.nextflow.io/docs/latest/cli.html)

## Inputs
* GWAS summary statistics
  * Must have the following columns (specify their column names in the config file)
    * Chromosome
    * Position/Base Pair
    * Allele 1/A1
    * Allele 2/A2
    * Beta OR Odds Ratio
    * P-Value OR Standard Error
* LD Reference Panel and Multi-Ancestry SNP Info File
  * Available at [PRScsx github](https://github.com/getian107/PRScsx/tree/master)
  * All data provided on github is in build 37 (b37/hg19)
  * Follow the instructions on the github to download and decompress the files
  * LD panels are ancestry-specific
  * SNP info file contains SNPs in all ancestry-specific LD panels and their population specific allele frequencies and allele flips

## Outputs
* Chromosome/cohort separated output files (raw output files) will be symlinked in <code>PRScsx_output/</code>
* Cohort separated output files (chromosome separated outputs concatenated and sorted) will be symlinked in in <code>Summary/</code>
* A violin plot comparing the variant level scores across cohorts will be symlinked as <code>plots/combined_PGS_violinplot.png</code>

## Configuration
* This pipeline requires specification of the following parameters in <code>prscsx.config</code>
  * List of chromosomes
  * Full file path to directory containing GWAS summary statistics
  * Path to summary/descriptor file containing columns with cohort nickname, ancestry, phenotype, GWAS sample size, and summary statistics filename
    * Each row represents a different GWAS for PRScsx to compute a score based on
    * The file must contain the following columns:
      * Cohort: a nickname for the cohort of interest
      * Ancestry: The GWAS genetically inferred ancestry group
      * Phenotype: The GWAS phenotype
      * Sample size: The GWAS sample size
      * Filename: Name of the GWAS summary statistics filename
    * Example File:    
   ![Title](images/Summary_Descriptor_File_Example.png)
  * Groovy dictornary specifying column names of of the summary/descriptor file
    * Necessary columns:
      * Cohort
      * Ancestry
      * Phenotype
      * Sample size
      * Summary statistics filename
  * Delimiter of the summary/descriptor file
  * Genome build of the summary statistics
    * Accepts either "b38" (build 38) or "b37" (build 37)
    * This method assumes that all summary statistic files are in the same genome build
  * Groovy dictionary specifying column names in GWAS summary statistics files
    * This method assumes that all summary statistic files have the same column names
    * Necessarily input columns specified above
  * Full file path to directory containing LD reference panels AND multi-ancestry SNP info file
    * The multi-ancestry SNP info file MUST be in the same directory as the LD panel directories, or else PRScsx will error out
  * Type of LD panel being used
    * Either 1000 genomes (label as 1kg) or UK Biobank (label as ukbb) LD panels will be accepted by PRScsx
  * PRScsx parameters
    * Parameters are descibed in <code>prscsx.config</code>. Full descriptions of these parameters are on the [PRScsx github](https://github.com/getian107/PRScsx/tree/master)
    * <code>param_a</code>
    * <code>param_b</code>
    * <code>param_phi</code>
    * <code>mcmc_iterations</code>
    * <code>mcmc_burnin</code>
    * <code>mcmc_thinning_factor</code>
    * <code>meta_flag</code>
    * <code>seed</code>

## Workflow Processeses
1. <code>make_pgs_input_and_bim</code>
  * Clean sumamry statistics file for PRScsx by removing insertions and deletions, rows with missing values, and converting numbers in scientific format to integers
  * Adds RSIDs to summary statistics to ensure the variant IDs match the LD panels provided by PRScsx using multi-ancestry SNP info files built into the docker/singularity image (identical to the ones on the PRScsx github (both UKBB and 1KG files) but the coordinates are in either b37 or b38 depending on the genome build specified)
  * Reformats summary statistics to make PRScsx input file and bim file
  * Parallelizes by <code>ancestry_phenotype</code>
2. <code>cat_bim</code>
  * Concatenates <code>ancestry_phenotype</code> separated bim files by <code>phenotype</code>
  * PRScsx only takes one bim file as input and it computes scores for multiple ancestry groups simultaneously
  * Parallelizes by <code>phenotype</code>
3. <code>prscsx</code>
  * Runs PRScsx
  * Parallelizes by <code>phenotype</code> and <code>chromosome</code>
4. <code>cat_pgs_outputs</code>
  * Concatenates PRScsx outputs by <code>ancestry_phenotype</code>
  * Parallelizes by <code>ancestry_phenotype</code>
5. <code>make_pgs_violinplots</code>
  * Creates violin plots of variant-level scores for each <code>ancestry_phenotype</code> group to compare distributions.

## Support
Please direct any problems or questions to Katie Cardone ([Katie.Cardone@pennmedicine.upenn.edu](mailto:Katie.Cardone@pennmedicine.upenn.edu))
