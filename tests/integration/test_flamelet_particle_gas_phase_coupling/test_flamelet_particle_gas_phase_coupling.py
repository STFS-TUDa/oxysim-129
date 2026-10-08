import subprocess
import pytest
import os 
import pandas as pd
import numpy as np
import pathlib
import os 
import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv
import subprocess

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

FIGURES_DIR = os.path.join(TEST_DIR, "report/figures")

# Create FIGURES_DIR if it does not exists
pathlib.Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)

OF_settings = pd.read_csv(os.path.join(TEST_DIR,'OF_settings.csv'),index_col=0)
reference_df = pd.read_csv(os.path.join(TEST_DIR,'reference/data/reference_data.csv'),index_col=0)

def prepare_testcase(index, testdir):
    for c in OF_settings.columns:
        command = 'grep -rl '+c+ ' ' + testdir + ' | xargs sed -i "s/'+c+'/'+str(OF_settings[c].iloc[index]).replace('nan','')+'/g"'
        subprocess.run(command,shell=True)
    # There were issues with copying the FLUT -> therefore it's linked 
    subprocess.run('ln -s ' + os.path.join(TEST_DIR,'Z_Y_YY_hanorm_non_reactive.h5') + ' ' + testdir, shell=True)

def read_foam_case(testdir):
    reader = pv.POpenFOAMReader(os.path.join(testdir,'case.foam'))
    reader.case_type = 'reconstructed'
    reader.set_active_time_value(0.3)
    vtk = reader.read()

    internal = vtk['internalMesh']
    x, y = 0.e-3, 0.e-3
    line = internal.sample_over_line((x,y,0),(x,y,np.max(internal.points[:,2])),resolution=100)

    df_in = pd.read_csv(os.path.join(testdir,"postProcessing/inletValue/0/surfaceFieldValue.dat"), skiprows=4,sep=r'\s+',usecols=[0,1])
    df_out = pd.read_csv(os.path.join(testdir,"postProcessing/outletValue/0/surfaceFieldValue.dat"), skiprows=4,sep=r'\s+',usecols=[0,1])
    phi_in = -df_in['Time'].iloc[-1]
    phi_out = df_out['Time'].iloc[-1]

    d = {var:[line[var][-1]] for var in reference_df.columns if var in line.array_names}
    d.update({'massflow_inlet':[phi_in],'massflow_outlet':[phi_out]})
    result_df = pd.DataFrame(d)
    return result_df

def comparison(result_df,reference_df,i,write_report=False):
    if write_report:
        d = {'':['result','reference']}
        d.update({col:[result_df[col].iloc[0],reference_df[col].iloc[i]] for col in result_df.columns})
        pd.DataFrame(d).to_csv(os.path.join(TEST_DIR,'report/report'+str(i)+'.csv'),index=False)
    for col in result_df.columns:
        np.testing.assert_allclose(result_df[col].iloc[0],reference_df[col].iloc[i],rtol=5e-3 if col != 'CH3O' else 20e-3)


def test_volatile(testdir, write_report):
    i = 0
    # Prepare testcase
    prepare_testcase(i, testdir)
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    result_df = read_foam_case(testdir)
    comparison(result_df,reference_df,i,write_report)

def test_volatile_plus_gas(testdir, write_report):
    i = 1
    # Prepare testcase
    prepare_testcase(i, testdir)
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    result_df = read_foam_case(testdir)
    comparison(result_df,reference_df,i,write_report)

def test_volatile_plus_gas_plus_char_off(testdir, write_report):
    i = 2
    # Prepare testcase
    prepare_testcase(i, testdir)
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    result_df = read_foam_case(testdir)
    comparison(result_df,reference_df,i,write_report)