from pathlib import Path
import pytest

import src.integrator 
import src.CH_file_parser
import numpy as np
from src.plot_helpers import build_time_axis
import write_ch as wch

# ------------------------- CREATING DUMMY .CH FILES ------------------------- # 

t1 = np.linspace(0, 12, 2400)  # times in minutes
values1  = np.sin(2 * np.pi * t1) * 500   # mAU
true_integral_1 = 0


params1 = wch.write_ch(
    path        = "DAD1A_test_1.ch",
    times_min   = t1,
    values      = values1,
    units       = "mAU",
    signal      = "DAD1A, Sig=215.0,16.0 Ref=off",
    method      = "DEMO_METHOD.M",
    instrument  = "Agilent 1260 Infinity II",
    date        = "28-Jun-26, 12:00:00",
    notebook    = "DEMO_NOTEBOOK",
    parent_dir  = "demo"
)


t2 = np.linspace(0, 11.5, 2400)   # times in minutes
values2  = np.sin(2 * np.pi * t2) * 500   # mAU
true_integral_2 = - 1/(2 * np.pi) * 500 * (np.cos(2 * np.pi * t2[-1]) - np.cos(2 * np.pi *t2[0]))


params2 = wch.write_ch(
    path        = "DAD1A_test_2.ch",
    times_min   = t2,
    values      = values2,
    units       = "mAU",
    signal      = "DAD1A, Sig=215.0,16.0 Ref=off",
    method      = "DEMO_METHOD.M",
    instrument  = "Agilent 1260 Infinity II",
    date        = "28-Jun-26, 12:00:00",
    notebook    = "DEMO_NOTEBOOK",
    parent_dir  = "demo",
)


# ------------------------- TEST ------------------------- # 




@pytest.mark.parametrize("time, signal, true_integral", [
    (t1, values1, true_integral_1 ),
    (t2, values2, true_integral_2 )
])
def test_calculate_integral_absolute_precision(time, signal, true_integral):
    # computing the integral
    computed_integral = src.integrator.calculate_integral(time, signal)
    assert np.allclose(computed_integral, true_integral, atol=1e-5)
    





@pytest.mark.parametrize("filename, time_start_s, use_header_time, true_integral", [
    ("DAD1A_test_1.ch", 0, True, true_integral_1 ),
    ("DAD1A_test_2.ch", 0, True, true_integral_2 )
])
def test_calculate_integral_in_context(filename, time_start_s, use_header_time, true_integral):
    ch_file = Path(filename)

    # getting the scaling factor
    data = ch_file.read_bytes()
    retrieved_header =  src.CH_file_parser.extract_header(data) 
    scaling_factor = retrieved_header["scaling_factor"]

    # retrieving the raw data and the time
    _, t, raw, signal = src.CH_file_parser.read_ch_data(ch_file, None, None, None, time_start_s, use_header_time)

    # computing the integral
    computed_raw_integral = src.integrator.calculate_integral(t, raw)
    computed_integral = computed_raw_integral * scaling_factor

    assert np.allclose(computed_integral, true_integral*60, atol=1e-3)
    

