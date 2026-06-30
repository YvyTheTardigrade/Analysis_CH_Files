from pathlib import Path
import pytest

import src.CH_file_parser 
import numpy as np
from src.plot_helpers import build_time_axis
import write_ch as wch

# ------------------------- CREATING DUMMY .CH FILES ------------------------- # 

t1 = np.linspace(0, 12, 2400) # times in minutes 
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


t2 = np.linspace(0, 11.5, 2400)  # times in minutes
values2  = np.sin(2 * np.pi * t2) * 500   # mAU
true_integral_2 = - 1/(2 * np.pi) * 500 * (np.cos(2 * np.pi * t2[-1]) - np.cos(2 * np.pi * t2[0]))


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


@pytest.mark.parametrize("filename,params", [
    ("DAD1A_test_1.ch", params1),
    ("DAD1A_test_2.ch", params2)
])
def test_extract_header(filename, params):
    ch_file = Path(filename)
    data = ch_file.read_bytes()
    retrieved_header =  src.CH_file_parser.extract_header(data) 
    assert retrieved_header["file_type_number"] == "130"
    assert retrieved_header["file_type_name"] == "LC DATA FILE"
    assert retrieved_header["notebook_name"] == params["notebook"]
    assert retrieved_header["parent_directory"] == params["parent_dir"]
    assert retrieved_header["date"] == params["date"]
    assert retrieved_header["method"] == params["method"]
    assert retrieved_header["instrument"] ==  params["instrument"]
    assert retrieved_header["units"] ==  params["units"]
    assert retrieved_header["signal_name"] == params["signal"]
    assert retrieved_header["first_time_ms"] == params["t_start_ms"]
    assert retrieved_header["last_time_ms"] == params["t_end_ms"]
    assert retrieved_header["scaling_factor"] == params["scaling_factor"]
    assert retrieved_header["header_size_bytes"] == wch.HEADER_SIZE


@pytest.mark.parametrize("filename,values", [
    ("DAD1A_test_1.ch", values1 ),
    ("DAD1A_test_2.ch", values2 )
])
def test_decode_signal_from_offset(filename, values):
    ch_file = Path(filename)
    data = ch_file.read_bytes()
    retrieved_header = src.CH_file_parser.extract_header(data)
    HEADER_SIZE = int(retrieved_header["header_size_bytes"])
    scaling_factor = retrieved_header["scaling_factor"]
    retrieved_raw = src.CH_file_parser.decode_signal_from_offset(data, HEADER_SIZE)
    effective_target_points = src.CH_file_parser.deduce_output_point_count(
            retrieved_header.get("first_time_ms"),
            retrieved_header.get("last_time_ms"),
            len(retrieved_raw),
        )

    signal_raw = src.CH_file_parser.align_signal(retrieved_raw, drop=0, target_points=effective_target_points)
    scaling = float(retrieved_header.get("scaling_factor", 1.0))
    signal_scaled = [x * scaling for x in signal_raw]
    #import matplotlib.pyplot as plt 
    #plt.plot(values)
    #plt.plot(signal_scaled)
    #plt.show()
    
    assert np.allclose(values, signal_scaled)


@pytest.mark.parametrize("signal,drop,target_pts,expected_res", [
    (np.arange(10), 2, None, np.arange(2,10)),
    (np.arange(10), 2, 4, np.arange(2,6)),
    (np.arange(10), 2, 10, np.arange(2,12))
])
def test_align_signal(signal, drop, target_pts, expected_res):
    aligned = np.array(src.CH_file_parser.align_signal(signal.tolist(), drop, target_pts))
    assert np.all(aligned == expected_res)

#def test_deduce_output_point_count():
#    pass 

@pytest.mark.parametrize("n_points, total_time_s, time_start_s, first_time_ms, last_time_ms, use_header_time, expected", [
    (2400, (t1[-1] - t1[0])/1000., t1[0]/1000., t1[0], t1[-1], True, t1),
    (2400, (t1[-1] - t1[0])/1000., t1[0]/1000., None, None, False, t1),
    (2400, (t2[-1] - t2[0])/1000., t2[0]/1000., t2[0], t2[-1], True, t2),
    (2400, (t2[-1] - t2[0])/1000., t2[0]/1000., None, None, False, t2),
])
def test_build_time_axis(n_points, total_time_s, time_start_s, first_time_ms, last_time_ms, use_header_time, expected):
    time_axis = src.CH_file_parser.build_time_axis(n_points, total_time_s, time_start_s, first_time_ms, last_time_ms, use_header_time)
    assert np.allclose(time_axis, expected/1000.)


@pytest.mark.parametrize("filename, off, nb_pts, duration, time_start_s, use_header_time, expected_time, expected_signal", [
    ("DAD1A_test_1.ch", None, None, None, 0, True, t1, values1),
    ("DAD1A_test_2.ch", None, None, None, 0, True, t2, values2)
])
def test_read_ch_data(filename, off, nb_pts, duration, time_start_s, use_header_time, expected_time, expected_signal):
    ch_file = Path(filename)
    _, t, _, signal = src.CH_file_parser.read_ch_data(ch_file, off, nb_pts, duration, time_start_s, use_header_time)
    assert np.allclose(t, expected_time*60.)  # expected times is given in minutes !!!
    assert np.allclose(signal, expected_signal)



