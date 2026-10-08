import pyFLUT
import pandas as pd
import scipy
import cantera as ct
import os
import numpy as np
import python_utilities as pytb
import sys

def get_coordinates(lines):
    firstLine = 0
    for i, line in enumerate(lines):
        if ("(" in line) and (firstLine==0):
            firstLine = i+1
            break

    nPoints = int(lines[firstLine-2])

    coordinates = np.empty((3, nPoints))
    for i, line in enumerate(lines[firstLine:firstLine+nPoints]):
        line = line.strip(")\n")
        line = line.strip("(")
        coordinates[:, i] = [float(coord) for coord in line.split(" ")]

    return coordinates

# Write out cell centres in OF
os.system("cp -r 0.org 0; blockMesh; postProcess -func writeCellCentres;")

# Define initial flame profile
# Read in grid coordinates
with open(os.path.join("./", "0/C")) as f:
    coordinate_lines = f.readlines()

# Field interpolation
coordinates = get_coordinates(coordinate_lines)
df = pd.DataFrame(data=coordinates[0,:].T, columns=["x"])

flame = pyFLUT.read_flame("../reference/cantera/fp_Z0.05000.ct")
flame = flame.map_variables(from_inp="X", to_inp="X", n_points=10)
flame = flame.map_variables(from_inp="X", to_inp="X", n_points=df["x"].values)
flame["Z"] = 0.05
flame["PV"] = flame["CO"] + flame["CO2"]
flame["p"] = 101325
flame["Ux"] = flame["sl"]
flame["Uy"] = flame["Uz"] = 0

df = flame.convert_to_DataFrame()

for var in ["PV", "Z"]:
    field_info = {
        'class': 'volScalarField', 
        'location': '0', 
        'object': var,
        'dimension': '[0 0 0 0 0 0 0]',
        'values': df[var].values, 
        'boundary': {
            'patch': ['inlet', 'outlet', 'frontAndBack'],
            'type': ['fixedValue', 'zeroGradient', 'empty', 'empty'], 
            'value': [df[var][0], None, None]}
    }
    pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volVectorField', 
    'location': '0', 
    'object': 'U', 
    'dimension': '[0 1 -1 0 0 0 0]',
    'values': df[["Ux", "Uy", "Uz"]].values, 
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['fixedValue', 'zeroGradient', 'empty'], 
        'value': [np.array([df['Ux'].values[0], 0, 0]), None, None]}
}
pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volScalarField', 
    'location': '0', 
    'object': 'p', 
    'dimension': '[1 -1 -2 0 0 0 0]',
    'values': 101325, 
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['zeroGradient', 'fixedValue', 'empty'], 
        'value': [None, 101325, None]}
}
pytb.foamTools.write_foam_file(field_info)

field_info = {
    'class': 'volScalarField', 
    'location': '0', 
    'object': 'omegaPV', 
    'dimension': '[1 -3 -1 0 0 0 0]',
    'values': 0, 
    'boundary': {
        'patch': ['inlet', 'outlet', 'frontAndBack'],
        'type': ['calculated', 'calculated', 'empty'], 
        'value': [0, 0, None]}
}
pytb.foamTools.write_foam_file(field_info)
