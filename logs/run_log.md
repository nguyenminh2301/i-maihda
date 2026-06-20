# Run log
Started: 2026-06-20T00:10:42.530650+00:00
Finished: 2026-06-20T00:10:45.063135+00:00
Python: 3.13.6

## Benchmark checks
| check                                   | passed   | value                    | criterion                                                    |
|:----------------------------------------|:---------|:-------------------------|:-------------------------------------------------------------|
| A_additive_is_additive_dominant         | True     | PCV=100.0; VPC_main=0.00 | PCV >= 80 and VPC_main < 1                                   |
| B_interaction_increases_vpc             | True     | A=4.32; B=22.58          | B VPC_null > A VPC_null + 5 percentage points                |
| B_interaction_leaves_residual_variance  | True     | PCV=35.8                 | PCV < 70                                                     |
| C_detection_reduces_observed_prevalence | True     | A=23.3%; C=11.3%         | C observed prevalence < A observed prevalence                |
| D_detection_can_mask_interaction_vpc    | True     | B=22.58; D=13.68         | D VPC_null < B VPC_null despite same residual-interaction SD |
| E_sparse_strata_are_flagged             | True     | B min_n=144; E min_n=1   | E minimum stratum size < B minimum stratum size              |

## Scenario results
| scenario   |   overall_prevalence |   vpc_null |   vpc_main |      pcv |   min_stratum_n |
|:-----------|---------------------:|-----------:|-----------:|---------:|----------------:|
| A          |            0.233333  |    4.31889 |    0       | 100      |             144 |
| B          |            0.271     |   22.5795  |   15.7797  |  35.7573 |             144 |
| C          |            0.112833  |    0       |    0       | nan      |             144 |
| D          |            0.136833  |   13.6776  |    8.79737 |  39.1222 |             144 |
| E          |            0.0905714 |   14.6996  |    9.44231 |  39.4942 |               1 |