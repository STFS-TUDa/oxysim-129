import subprocess
import pytest
import os 
import pandas as pd
import numpy as np
import pathlib
import json
import os 
import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

FIGURES_DIR = os.path.join(TEST_DIR, "report/figures")

# Create FIGURES_DIR if it does not exists
pathlib.Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)

def read_time_resolved_particle_data(pp_path):
    prt_dfs = {}
    for cloud_path, cloud_name in zip(['writeCarbCloud'],['carbCloud']):
        path = os.path.join(pp_path, cloud_path)
        with open(os.path.join(path,cloud_name + '.vtp.series')) as f:
            d = json.load(f)

        time_file_dict = {}
        times = []
        files = []
        for di in d['files']:
            times.append(di['time'])
            files.append(os.path.join(path,di['name']))

        sorting = np.argsort(np.array(times))
        files = list(np.array(files)[sorting])

        data_list = []
        for file in files:
            vtk = pv.read(file)
            data_list.append(np.concatenate((vtk['T'][:,np.newaxis],vtk['origId'][:,np.newaxis],vtk['d'][:,np.newaxis],vtk['YC(s)'][:,np.newaxis],vtk['TimeValue']*np.ones(vtk.n_points)[:,np.newaxis]),axis=1))
        data = np.concatenate(tuple(data_list))
        prt_dfs[cloud_name] = pd.DataFrame({'T':data[:,0],'origId':data[:,1],'d':data[:,2],'YC(s)':data[:,3],'Time':data[:,4]})

    rho_prt = 707
    mass0 = 1e-11
    prt_df = prt_dfs['carbCloud']
    prt_df['mass'] = prt_df['d']**3/6*np.pi*rho_prt
    prt_df['mass/mass0'] = prt_df['mass']/mass0
    return prt_df

def calculate_volatiles(testdir):

    df = pd.read_csv(os.path.join(testdir,'postProcessing/sumZ/0/surfaceFieldValue.dat'),header=5,sep=r'\t',engine='python')
    return np.sum(df['weightedSum(Z)']*5e-6)

def test_C2SM(testdir, write_report):
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # Load particles
    prt_df = read_time_resolved_particle_data(os.path.join(testdir,'postProcessing'))

    # Calculate volatile mass passing the outlet
    mass_gas_volatiles = calculate_volatiles(testdir)

    # Load particle reference data
    particle_reference_df = pd.read_csv(os.path.join(TEST_DIR,'reference/foam/profiles.csv'),index_col=0)

    # Load final reference
    final_reference_df = pd.read_csv(os.path.join(TEST_DIR,'reference/foam/final.csv'),index_col=0)

    # plot results
    if write_report:
        fig, ax = plt.subplots(1,1)
        ax.plot(prt_df['Time'],prt_df['mass/mass0'],label='mass/mass0 OpenFOAM')
        ax.plot(particle_reference_df['Time'],particle_reference_df['mass/mass0'],label='Mass/Mass0 Python',linestyle='--')

        ax.plot(prt_df['Time'],prt_df['YC(s)'],label='Ychar OpenFOAM')
        ax.plot(particle_reference_df['Time'],particle_reference_df['YC(s)'],label='YChar Python',linestyle='--')

        ax.legend()
        ax.set_ylabel('Y (-)')
        ax.set_xlabel('t (s)')
        fig.savefig(os.path.join(FIGURES_DIR,'comparison.png'))

        if write_report:
            # particle profiles
            prt_df.to_csv(os.path.join(TEST_DIR,'report/profile.csv'))

            # final values
            d = {'':['result','reference']}
            d.update({'Volatiles':[mass_gas_volatiles,final_reference_df['mass_volatiles_final'].iloc[0]],'Final particle mass':[prt_df['mass'].iloc[-1],final_reference_df['mass_final_particle'].iloc[0]]})
            pd.DataFrame(d).to_csv(os.path.join(TEST_DIR,'report/final.csv'),index=False)


    # test
    np.testing.assert_allclose(prt_df['mass/mass0'],particle_reference_df['mass/mass0'],atol=1.1e-2)
    np.testing.assert_allclose(prt_df['YC(s)'],particle_reference_df['YC(s)'],atol=0.016)
    np.testing.assert_allclose(np.array(mass_gas_volatiles), final_reference_df['mass_volatiles_final'],rtol=1e-3)
    np.testing.assert_allclose(np.array(prt_df['mass'].iloc[-1]), final_reference_df['mass_final_particle'].iloc[0],rtol=1e-3)


if __name__ == "__main__":
    test_C2SM('foam_template', True)
