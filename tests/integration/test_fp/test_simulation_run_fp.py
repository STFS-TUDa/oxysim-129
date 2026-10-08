import subprocess
import pytest
import os
import pyvista as pv
import pandas as pd
import numpy as np
from scipy.interpolate import interp1d
from scipy.integrate import trapezoid
from scipy.stats import pearsonr
import pathlib
import shutil
from matplotlib import pyplot as plt
from conftest import plot_line_profiles
from python_utils.foam_utils import foam_to_df

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(TEST_DIR, "report/figures")

# Create FIGURES_DIR if it does not exists
pathlib.Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)


def calc_flame_position(df):
    f = interp1d(df["yc"], df["x"])
    return f((df["yc"].max() + df["yc"].min()) / 2)


def calc_consumption_speed(df):
    integral = trapezoid(df["omega_yc"], x=df["x"])
    return 1.0 / (df["rho"][0] * (df["yc"][len(df) - 1] - df["yc"][0])) * integral


def calc_displacement_speed(df):
    return (df["U_0"][len(df) - 1] - df["U_0"][0]) / (
        df["rho"][0] / df["rho"][len(df) - 1] - 1
    )


def calc_flame_thickness(df):
    """
    calculate the flame thickness based on the normalized maximum gradient of the progress variable
    """
    # calc flame thickness
    PV = df["yc"]
    dPVdx = np.gradient(PV, df["x"])
    ft = (PV.max() - PV.min()) / dPVdx.max()
    return ft


def validate_freelyPropagatingFlame(result, reference, yvs=[]):

    # Compare global flame properties
    np.testing.assert_allclose(
        calc_flame_thickness(result),
        calc_flame_thickness(reference),
        rtol=5e-2,
        err_msg=f"Mismatch in the thermal flame thickness.",
    )
    np.testing.assert_allclose(
        calc_displacement_speed(result),
        calc_displacement_speed(reference),
        rtol=5e-2,
        err_msg="Missmatch in the displacement speed.",
    )
    np.testing.assert_allclose(
        calc_consumption_speed(result),
        calc_consumption_speed(reference),
        rtol=5e-2,
        err_msg="Missmatch in the consumption speed.",
    )

    # Compare flame position
    dx = np.round(np.max(np.diff(result["x"])), 9)
    np.testing.assert_allclose(
        calc_flame_position(result),
        calc_flame_position(reference),
        rtol=0,
        atol=dx * 2,
        err_msg=f"Mismatch in the flame position (defined by the Temperature isoline) by more than 2 cells (dx={dx}).",
    )
    # Define shifted x-position for integral value comparison
    result["x_rel"] = result["x"] - calc_flame_position(result)
    reference["x_rel"] = reference["x"] - calc_flame_position(reference)

    # Compare the min, max and integral value of quantities of intere
    for yv in yvs:
        # Difference of single point values
        rtol = 2e-2
        atol = np.max(np.abs(result[yv])) * 1e-5
        np.testing.assert_allclose(
            result[yv][0],
            reference[yv][0],
            rtol=rtol,
            atol=atol,
            err_msg=f"Mismatch in fresh gas value of {yv}.",
        )  # fresh gas state
        np.testing.assert_allclose(
            result[yv][len(result) - 1],
            reference[yv][len(reference) - 1],
            rtol=rtol,
            atol=atol,
            err_msg=f"Mismatch in burnt gas value of {yv}.",
        )  # burnt gas state
        np.testing.assert_allclose(
            np.min(result[yv]),
            np.min(reference[yv]),
            rtol=rtol,
            atol=atol,
            err_msg=f"Mismatch in minimum value of {yv}.",
        )  # minimum value
        np.testing.assert_allclose(
            np.max(result[yv]),
            np.max(reference[yv]),
            rtol=rtol,
            atol=atol,
            err_msg=f"Mismatch in maximum value of {yv}.",
        )  # maximum value

        int_ref = trapezoid(reference[yv], reference["x"])
        int_result = trapezoid(result[yv], result["x"])

        # Maximum difference of 1% in the integral values of the curves
        np.testing.assert_allclose(
            int_result,
            int_ref,
            rtol=rtol,
            atol=atol,
            err_msg=f"Mismatch in integral value of {yv}.",
        )

        # Shift flame position and interpolate onto new grid
        x_min = np.max([reference["x_rel"].min(), result["x_rel"].min()])
        x_max = np.min([reference["x_rel"].max(), result["x_rel"].max()])
        x_grid = np.linspace(x_min, x_max, len(reference))
        y_reference = interp1d(reference["x_rel"], reference[yv])(x_grid)
        y_result = interp1d(result["x_rel"], result[yv])(x_grid)

        # Area between the shifted curves normalized by the integral
        norm_area_between_curves = (
            trapezoid(np.abs(y_result - y_reference), x_grid) / int_ref
        )
        assert (
            norm_area_between_curves < rtol
        ), f"Normalized area for variable {yv} between shifted curve is {norm_area_between_curves}, which is greater than the threshold of {rtol}."

        # Pearson correlation coefficient (for variables that are not approximately constant)
        if not np.allclose(reference[f"{yv}"].max(), reference[f"{yv}"].min()):
            pearson_correlation, _ = pearsonr(y_reference, y_result)
            assert (
                pearson_correlation >= 0.999
            ), f"Pearson correlation coeffiecient for shifted variable {yv} is {pearson_correlation}, which is greater than the threshold of 0.99."

    # Check pressure fluctuations
    np.testing.assert_allclose(
        result["p"],
        101325,
        rtol=2e-3,
        err_msg="Pressure fluctuations are higher than the tolerance.",
    )


def test_fpFlame(testdir, write_report):
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun", "-p"], cwd=testdir, check=True)

    # load and compare results
    time = 0.01
    result = foam_to_df(
        path=os.path.join(testdir, "case.foam"),
        time=time,
        index=False,
        keep_boundaries=["inlet", "outlet"],
    )
    result.sort_values("x", inplace=True, ignore_index=True)
    reference = pd.read_parquet(
        os.path.join(TEST_DIR, "reference/foam/t_{:.3f}.parq".format(time))
    )
    reference.sort_values("x", inplace=True, ignore_index=True)

    quantities_of_interest = ["Z", "yc", "omega_yc", "U_0", "rho"]

    # plot results
    if write_report:
        plot_line_profiles(
            reference=reference,
            result=result,
            xv="x",
            yvs=[*quantities_of_interest, "p"],
            filename=f"{FIGURES_DIR}/t_{time}_physicalSpace.png",
        )
        plot_line_profiles(
            reference=reference,
            result=result,
            xv="yc",
            yvs=[*quantities_of_interest, "p"],
            filename=f"{FIGURES_DIR}/t_{time}_stateSpace.png",
        )

    # test
    validate_freelyPropagatingFlame(
        reference=reference, result=result, yvs=quantities_of_interest
    )
