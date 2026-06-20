# R Run Log
Started:  2026-06-20T13:58:17+0700
Finished: 2026-06-20T13:58:29+0700
R:        R version 4.3.3 (2024-02-29 ucrt)

## Benchmark checks
                                   check passed                    value
         A_additive_is_additive_dominant   TRUE PCV=100.0; VPC_main=0.00
             B_interaction_increases_vpc   TRUE          A=4.50; B=17.20
  B_interaction_leaves_residual_variance   TRUE                 PCV=51.7
 C_detection_reduces_observed_prevalence   TRUE         A=23.6%; C=11.6%
    D_detection_can_mask_interaction_vpc   TRUE         B=17.20; D=16.06
             E_sparse_strata_are_flagged   TRUE   B min_n=130; E min_n=2
                                                    criterion
                                   PCV >= 80 and VPC_main < 1
                B VPC_null > A VPC_null + 5 percentage points
                                                     PCV < 70
                C observed prevalence < A observed prevalence
 D VPC_null < B VPC_null despite same residual-interaction SD
              E minimum stratum size < B minimum stratum size

## Scenario results
 scenario overall_prevalence   vpc_null  vpc_main       pcv min_stratum_n
        A          0.2365000  4.4953400 0.0000000 100.00000           130
        B          0.2633333 17.1978511 9.1215635  51.67455           130
        C          0.1165000  0.8175855 0.1806489  78.04557           130
        D          0.1303333 16.0581117 8.0152621  54.45023           130
        E          0.1162857 19.4709700 2.9820912  87.28741             2

## Note on cross-language comparison
R uses Mersenne Twister (set.seed); Python uses PCG64 (numpy.random.default_rng).
Exact numerical values differ, but statistical patterns and benchmark pass/fail
should be equivalent. If benchmark pass/fail disagrees, investigate seed sensitivity.
