// include both PRScsx and Apply PGS workflows
include {PRScsx_setup} from './prscsx.nf'
include {PRScsx} from './prscsx.nf'
include {PLINK2_SCORE} from './plink2_score.nf'

// run PRScsx workflow then Apply PGS workflow
workflow {
    prscsx_input=PRScsx_setup()
    prscsx_weights = PRScsx(prscsx_input)

    indiv_polygenic_scores = PLINK2_SCORE(prscsx_weights)
}