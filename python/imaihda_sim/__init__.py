"""Synthetic-data I-MAIHDA sensitivity demo for HIC-MIC comparisons."""
from .simulate import simulate_intersectional_data, scenario_grid
from .fit import fit_imaihda, vpc_latent, pcv, LOGISTIC_L1_VARIANCE
from .detection import (
    correct_detection_bias,
    vpc_detection_bounds,
    detection_tipping_point,
)
from .uncertainty import sparse_strata_vpc
