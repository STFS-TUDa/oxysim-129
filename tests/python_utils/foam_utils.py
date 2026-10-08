"""
Minimal OpenFOAM-to-pandas helpers used by the integration tests.

Vendored from the (now-removed) python_toolbox submodule's
python_utilities.pyvistaTools module — only the call chain actually used
by the test suite (foam_to_df) is kept.
"""

import numpy as np
import pandas as pd
import pyvista as pv


def foam_to_df(path, time, index=True, keep_boundaries=[], case_type="reconstructed", variables=None, point_or_cell="point"):
    """
    Reads an OpenFOAM case with pyvista and converts it to a pandas DataFrame.

    Parameters
    ----------
    path: str
        Path to the OpenFOAM case file (an empty file with a .foam extension).
    time: int or float
        Time index or time value to read.
    index: bool
        Use time index (True) or time value (False).
    keep_boundaries : list, optional
        Boundaries to include in the DataFrame, by default [].
    case_type: str
        Case type to read: "reconstructed" or "decomposed".
    variables: list(str)
        Names of variables to include.
    point_or_cell: str
        Use point data ("point") or cell data ("cell").

    Returns
    -------
    pd.DataFrame
    """
    multiblock = _foam_to_pyvista(path, time=time, index=index, case_type=case_type)
    return _foam_multiblock_to_df(multiblock, keep_boundaries, variables=variables, point_or_cell=point_or_cell)


def _foam_to_pyvista(path, time, index=True, case_type="reconstructed"):
    reader = pv.POpenFOAMReader(path)
    reader.cell_to_point_creation = False
    reader.case_type = case_type

    if index:
        time = reader.time_values[time]
    reader.set_active_time_value(time)

    return reader.read()


def _foam_multiblock_to_df(multiblock, keep_boundaries=[], variables=None, point_or_cell="point"):
    for boundary_name in multiblock[1].keys():
        if boundary_name not in keep_boundaries:
            multiblock[1].pop(boundary_name)

    grid = multiblock.combine()
    grid = grid.cell_centers()
    grid.clear_cell_data()
    return _convert_to_df(grid, variables=variables, point_or_cell=point_or_cell)


def _convert_to_df(grid, variables=None, point_or_cell="point"):
    df_columns = ["x", "y", "z"]
    if point_or_cell == "point":
        grid_data_fn = grid.point_data
        var_array = grid.points
    elif point_or_cell == "cell":
        grid_data_fn = grid.cell_data
        var_array = grid.cell_centers().points

    if not variables:
        variables = grid_data_fn.keys()
    variables = [v for v in variables if v not in df_columns]

    for v in variables:
        data = grid_data_fn[v]
        try:
            number_of_vars = np.shape(data)[1]
            var_array = np.append(var_array, data, axis=-1)
            df_columns += ["{}_{}".format(v, i) for i in range(number_of_vars)]
        except IndexError:
            df_columns.append(v)
            data = np.expand_dims(data, axis=-1)
            var_array = np.append(var_array, data, axis=-1)

    return pd.DataFrame(var_array, columns=df_columns)
