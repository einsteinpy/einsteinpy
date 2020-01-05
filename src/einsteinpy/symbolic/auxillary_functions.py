import warnings

import sympy


def raise_warning(WarningType, message):
    warnings.warn(message, WarningType)


def _flatten_list(seq):
    # flatten an arbitarily nested list
    if not seq:
        return []
    if not isinstance(seq[0], list):
        return [seq[0]] + _flatten_list(seq[1:])
    return _flatten_list(seq[0]) + _flatten_list(seq[1:])


def simplify_sympy_array(arr):
    flattened_list = _flatten_list(arr.tolist())
    simplified_flattened_list = [sympy.simplify(e) for e in flattened_list]
    return sympy.Array(simplified_flattened_list, arr.shape)
