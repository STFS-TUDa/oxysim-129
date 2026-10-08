import subprocess
import pytest
import os
import pandas as pd
import numpy as np
import pathlib
import json
import numpy as np
from matplotlib import pyplot as plt
import pyvista as pv

TEST_DIR = os.path.dirname(os.path.abspath(__file__))

FIGURES_DIR = os.path.join(TEST_DIR, "report/figures")
pathlib.Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)

def read_time_resolved_particle_data(pp_path):
    prt_dfs = {}
    for cloud_path, cloud_name in zip(['writeCarbCloud'], ['carbCloud']):
        path = os.path.join(pp_path, cloud_path)
        with open(os.path.join(path, cloud_name + '.vtp.series')) as f:
            d = json.load(f)

        times = []
        files = []
        for di in d['files']:
            times.append(di['time'])
            files.append(os.path.join(path, di['name']))

        sorting = np.argsort(np.array(times))
        files = list(np.array(files)[sorting])

        data_list = []
        for file in files:
            vtk = pv.read(file)
            data_list.append(np.concatenate((
                vtk['T'][:, np.newaxis],
                vtk['origId'][:, np.newaxis],
                vtk['TimeValue'] * np.ones(vtk.n_points)[:, np.newaxis]
            ), axis=1))
        data = np.concatenate(tuple(data_list))
        prt_dfs[cloud_name] = pd.DataFrame({'T': data[:, 0], 'origId': data[:, 1], 'Time': data[:, 2]})
    return prt_dfs

def test_heat_transfer_particle(testdir, write_report):
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # Load probe data
    T_probes = pd.read_csv(os.path.join(testdir, 'postProcessing/probes/0/T'),
                           skiprows=4, sep=r'\s+', usecols=[0, 1, 2, 3])
    T_probes = T_probes.set_axis(['Time', 'Position 1', 'Position 2', 'Position 3'], axis=1)

    # Load carbonaceous particle data
    prt_dfs = read_time_resolved_particle_data(os.path.join(testdir, 'postProcessing'))
    prt_dfs['carbCloud'] = prt_dfs['carbCloud'].groupby('origId')

    result_df = pd.DataFrame({
        'Time': T_probes['Time'],
        'Probe 1': T_probes['Position 1'],
        'Probe 2': T_probes['Position 2'],
        'Probe 3': T_probes['Position 3'],
    })
    grouped_df = prt_dfs['carbCloud']
    for i, pid in enumerate(grouped_df.groups.keys()):
        col = 'carbCloud Position ' + str(i + 1)
        tmp_df = grouped_df.get_group(pid)[['T']].rename(columns={'T': col}).iloc[1:].reset_index(drop=True)
        result_df = pd.concat([result_df, tmp_df], axis=1)

    reference_data = pd.read_csv(os.path.join(TEST_DIR, 'reference/foam/reference_data.csv'), index_col=0)

    if write_report:
        fig, ax = plt.subplots(1, 1)
        for carb_id, pos in zip(prt_dfs['carbCloud'].groups.keys(),
                                ['Position 1', 'Position 2', 'Position 3']):
            p = ax.plot(T_probes['Time'], T_probes[pos])
            prt = prt_dfs['carbCloud'].get_group(carb_id)
            ax.plot(prt['Time'], prt['T'], c=p[0].get_color(), linestyle='--')
        ax.set_ylabel('$T$ (K)')
        ax.set_xlabel('$t$ (s)')
        ax.grid()
        fig.savefig(os.path.join(FIGURES_DIR, 'comparison.png'))

    report_df = pd.DataFrame()
    for col in reference_data.columns:
        report_df[col] = result_df[col]
        report_df[col + '_ref'] = reference_data[col]
    report_df.to_csv(os.path.join(TEST_DIR, 'report/report.csv'))

    for col in reference_data.columns:
        reference = reference_data[col].to_numpy()
        test_result = result_df[col].to_numpy()
        np.testing.assert_allclose(reference, test_result, atol=5e-1)


if __name__ == "__main__":
    test_heat_transfer_particle('foam_template', True)
