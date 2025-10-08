"""
@author: Nikolay S. Markov
"""

import functools
import os
import numpy as np
import pandas as pd

from . import tlc, fev1, fvc, fev1_fvc, dlco


PARAMS = {
    'tlc': tlc,
    'fev1': fev1,
    'fvc': fvc,
    'fev1_fvc': fev1_fvc,
    'dlco': dlco
}


def load_splines(param):
    data = []
    for sex in ('female', 'male'):
        table = pd.read_csv(os.path.join(
            os.path.dirname(__file__),
            param,
            f'{sex}.csv'
        ))
        table = table.iloc[:, :4] # keep only first 4 columns
        table.columns = ('age', 'Mspline', 'Sspline', 'Lspline')
        data.append(table)
        data[-1]['sex'] = sex
    data = pd.concat(data)
    return data.set_index(['sex', 'age'])


def compute_spline(param, formula, sex, age, interpolate=False):
    SPLINES = {
        'median': 'Mspline'
    }
    formula = SPLINES[formula]
    lookup = load_splines(param)
    if interpolate:
        pass
    else:
        if isinstance(age, pd.Series):
            age = (age * 4).round() / 4
            spline = pd.Series(np.nan, index=age.index)
            idx = pd.MultiIndex.from_arrays([sex, age])
            spline[age.notna()] = lookup[formula][idx[age.notna()]].values
        else:
            age = round(age * 4) / 4
            spline = lookup[formula][(sex, age)]
    return spline


def get_reference(param, formula, sex, age, height):
    param_impl = PARAMS[param]
    formula_impl = {}
    for s in ('female', 'male'):
        formula_impl[s] = getattr(getattr(param_impl, s), formula)
    spline = compute_spline(param, formula, sex, age)
    if isinstance(sex, pd.Series):
        result = pd.Series(np.nan, index=sex.index)
        for s in ('female', 'male'):
            idx = sex.eq(s)
            result[idx] = formula_impl[s](np, age[idx], height[idx], spline[idx])
        return result
    return formula_impl[sex](np, age, height, spline)


def percent_predicted(param, sex, age, height, value):
    return 100 * value / get_reference(param, 'median', sex, age, height)


for param in PARAMS.keys():
    locals()[f'{param}_percent_predicted'] = functools.partial(percent_predicted, param)
