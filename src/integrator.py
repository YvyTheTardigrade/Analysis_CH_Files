from typing import Any, Dict, List, Optional, Set, Tuple


def calculate_integral(times: List[float], signal: List[int]) -> float:
    if len(times) != len(signal):
        raise ValueError("times and signal must have the same length")

    n = len(signal)
    if n < 2:
        return 0.0

    dt = times[1] - times[0]
    if n == 2:
        return dt * (signal[0] + signal[1]) / 2.0

    intervals = n - 1
    if intervals % 2 == 0:
        odd_sum = sum(signal[i] for i in range(1, n - 1, 2))
        even_sum = sum(signal[i] for i in range(2, n - 1, 2))
        return (dt / 3.0) * (signal[0] + signal[-1] + 4.0 * odd_sum + 2.0 * even_sum)

    usable_n = n - 1
    odd_sum = sum(signal[i] for i in range(1, usable_n - 1, 2))
    even_sum = sum(signal[i] for i in range(2, usable_n - 1, 2))
    simpson_area = (dt / 3.0) * (signal[0] + signal[usable_n - 1] + 4.0 * odd_sum + 2.0 * even_sum)
    trapezoid_tail = dt * (signal[usable_n - 1] + signal[usable_n]) / 2.0
    return simpson_area + trapezoid_tail


