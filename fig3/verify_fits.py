"""Re-evaluate saved MATLAB predictions with an independent split LSODA solve."""
import csv
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp


def replay(times, dose, wash, initial_count, parameters, kind, pk_mode, panel):
    k_fe, k_fb = 1.25, 0.21
    input_gain = 3.5e-6 * 250 / 0.004629
    kp, theta = parameters['kp'], parameters[f'theta_{panel}']
    k_bf = parameters.get('kBF', 0)
    exponent = parameters.get('n', 1)
    plateau = input_gain * dose * wash * k_fb / (k_fe + k_fb)

    def kill(cb):
        if kind == 'linear':
            return parameters['beta'] * cb
        if cb == 0:
            return 0
        return parameters['Kmax'] / (1 + (parameters['Khalf'] / cb) ** exponent)

    def pk(state, ce):
        cf, cb = state[:2]
        return [input_gain * ce - (k_fe + k_fb) * cf + k_bf * cb, k_fb * cf - k_bf * cb]

    def rhs(_, state, ce):
        # Guard tiny negative numerical CB in the independent ODE reference.
        cb = plateau if pk_mode == 'fixed' else max(0, state[1])
        return [*pk(state, ce), kp * (1 - np.exp(state[2]) / theta) - kill(cb)]

    options = dict(method='LSODA', rtol=2e-11, atol=1e-13, dense_output=True)

    def solve(rhs, interval, state, ce):
        result = solve_ivp(rhs, interval, state, args=(ce,), **options)
        assert result.success, result.message
        return result

    first_time = times[0]
    pk_state = [0, 0]
    if first_time > 0:
        pk_rhs = lambda _, state, ce: pk(state, ce)
        pk_state = solve(pk_rhs, (0, min(wash, first_time)), pk_state, dose).y[:, -1]
        if first_time > wash:
            pk_state = solve(pk_rhs, (wash, first_time), pk_state, 0).y[:, -1]
    state = [*pk_state, np.log(initial_count)]
    before = times <= wash
    predicted = np.empty_like(times)
    if first_time < wash:
        pre = solve(rhs, (first_time, min(times[-1], wash)), state, dose)
        predicted[before] = np.exp(pre.sol(times[before])[2])
        state = pre.y[:, -1]
    elif first_time == wash:
        predicted[before] = initial_count
    if np.any(~before):
        post = solve(rhs, (max(first_time, wash), times[-1]), state, 0)
        predicted[~before] = np.exp(post.sol(times[~before])[2])
    return predicted


def rows(path):
    with path.open(newline='') as stream:
        return list(csv.DictReader(stream))


def main(output):
    observations = defaultdict(list)
    for row in rows(Path(__file__).with_name('observations.csv')):
        observations[row['panel']].append(row)

    output = Path(output)
    estimates = defaultdict(dict)
    for row in rows(output / 'parameters.csv'):
        estimates[row['model']][row['parameter']] = float(row['estimate'])
    residuals = defaultdict(list)
    for row in rows(output / 'residuals.csv'):
        residuals[int(row['model_id']), int(row['panel_id'])].append(row)

    maximum_error = 0
    pair_count = 0
    for summary in rows(output / 'summary.csv'):
        model = summary['model']
        mode, kind = model.split('_', 1)
        for panel, group in observations.items():
            times = np.array([float(row['time_h']) for row in group])
            counts = np.array([float(row['cells']) for row in group])
            dose = float(group[0]['dose_nM'])
            wash = float(group[0]['wash_h'])
            reference = replay(times, dose, wash, counts[0], estimates[model], kind, mode, panel)
            panel_id = ord(panel) - ord('a') + 1
            fitted_rows = residuals[int(summary['model_id']), panel_id]
            assert len(times) == len(fitted_rows)
            np.testing.assert_allclose(times, [float(row['time_h']) for row in fitted_rows], atol=1e-10)
            np.testing.assert_allclose(counts, [float(row['observed_N']) for row in fitted_rows], atol=1e-10)
            fitted = np.array([float(row['predicted_N']) for row in fitted_rows])
            maximum_error = max(maximum_error, np.max(np.abs(reference - fitted)) / max(counts))
            pair_count += 1
    assert pair_count == 81 and maximum_error < 2e-7
    print(f'{pair_count} model/condition pairs, maximum scaled difference vs LSODA: {maximum_error:.6g}')


if __name__ == '__main__':
    if len(sys.argv) != 2:
        raise SystemExit('Usage: python fig3/verify_fits.py RESULTS_DIRECTORY')
    main(sys.argv[1])
