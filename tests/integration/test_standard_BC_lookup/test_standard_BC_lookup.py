import subprocess
import pytest
import os
import pandas as pd
import numpy as np
import pyvista as pv

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

def prepare_testcase(testdir):
    # There were issues with copying the FLUT -> therefore it's linked 
    subprocess.run('ln -s ' + os.path.join(TEST_DIR,'FLUT.h5') + ' ' + testdir, shell=True)
    subprocess.run('ln -s ' + os.path.join(TEST_DIR,'FLUT_BC.h5') + ' ' + testdir, shell=True)

def read_foam_case(testdir):
    reader = pv.POpenFOAMReader(os.path.join(testdir,'case.foam'))
    reader.case_type = 'reconstructed'
    reader.set_active_time_value(0.001)
    vtk = reader.read()

    wall = vtk['boundary']['WALL']

    T_BC = wall["T_BC"][:100]
    T = wall["T"][:100]
    coords = wall.cell_centers().points[:100,2]

    return T_BC, T, coords

def comparison(T_BC, T, coords, write_report=False):
    if write_report:
        if not os.path.isdir(os.path.join(TEST_DIR, "report")):
            os.mkdir(os.path.join(TEST_DIR, "report"))
        d = {"T_BC": T_BC, "T": T, "z": coords, "err": T_BC - T}
        pd.DataFrame(d).to_csv(os.path.join(TEST_DIR,'report/report.csv'),index=False)
    np.testing.assert_allclose(T_BC,T,atol=10)


def test_standard_BC_lookup(testdir, write_report):
    # run OpenFOAM
    prepare_testcase(testdir)
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    T_BC, T, coords = read_foam_case(testdir)
    comparison(T_BC, T, coords, write_report)
