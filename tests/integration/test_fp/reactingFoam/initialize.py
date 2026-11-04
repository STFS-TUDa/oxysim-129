import os
import numpy as np
import pandas as pd
import python_utilities as pytb

os.makedirs("0.org", exist_ok=True)

n_cells = 2500
length = 0.025
x = (np.arange(n_cells) + 0.5) * length / n_cells

flame = pd.read_parquet("../reference/cantera/fp_Z0.05000.parq")
species = ["H2", "H", "O2", "O", "OH", "HO2", "H2O2", "H2O", "CH4", "CH3", "CH2O", "CH3O", "HCO", "CO", "CO2", "N2"]

x_coarse = np.linspace(0, length, 10)
df = pd.DataFrame({"x": x})
for var in species + ["T", "sl"]:
    df[var] = np.interp(x, x_coarse, np.interp(x_coarse, flame["X"], flame[var]))

for s in species:
    field_info = {
        'class': 'volScalarField',
        'location': '0.org',
        'object': s,
        'dimension': '[0 0 0 0 0 0 0]',
        'values': df[s].values,
        'boundary': {
            'patch': ['inlet', 'outlet', 'frontAndBack'],
            'type': ['fixedValue', 'zeroGradient', 'empty'],
            'value': [df[s][0], None, None]}
    }
    pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volScalarField',
    'location': '0.org',
    'object': 'T',
    'dimension': '[0 0 0 1 0 0 0]',
    'values': df['T'].values,
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['fixedValue', 'zeroGradient', 'empty'],
        'value': [df["T"][0], None, None]}
}
pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volVectorField',
    'location': '0.org',
    'object': 'U',
    'dimension': '[0 1 -1 0 0 0 0]',
    'values': np.column_stack([df["sl"].values, np.zeros(n_cells), np.zeros(n_cells)]),
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['fixedValue', 'zeroGradient', 'empty'],
        'value': [np.array([df["sl"][0], 0, 0]), None, None]}
}
pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volScalarField',
    'location': '0.org',
    'object': 'p',
    'dimension': '[1 -1 -2 0 0 0 0]',
    'values': 101325,
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['zeroGradient', 'waveTransmissive', 'empty'],
        'value': [None, {'field': 'p', 'psi': 'thermo:psi', 'gamma': 1.4, 'fieldInf': 101325, 'lInf': 0.025, 'value': 101325}, None]}
}
pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volScalarField',
    'location': '0.org',
    'object': 'Ydefault',
    'dimension': '[0 0 0 0 0 0 0]',
    'values': 0,
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['fixedValue', 'zeroGradient', 'empty'],
        'value': [0, None, None]}
}
pytb.foamTools.write_foam_file(field_info)
