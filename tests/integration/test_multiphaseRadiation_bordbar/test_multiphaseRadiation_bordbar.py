import subprocess
import os
import pandas as pd
import numpy as np
import pathlib
import matplotlib.pyplot as plt 
import re
import sys
import pyvista as pv

TEST_DIR = os.path.dirname(os.path.abspath(__file__))
FIGURES_DIR = os.path.join(TEST_DIR, "report/figures")

# Create FIGURES_DIR if it does not exists
pathlib.Path(FIGURES_DIR).mkdir(parents=True, exist_ok=True)

def particle_emissive_power(
    radAreaPT4,
    dt,
    cell_volume,
    emissivity
):
    """
    Total particle emissive power density E_p [W/m³].
    """
    SIGMA = 5.670374419e-8  # W/(m² K⁴)
    radAreaPT4 = np.asarray(radAreaPT4, dtype=float)
    cell_volume = np.asarray(cell_volume, dtype=float)

    return (4.0 * emissivity * SIGMA * radAreaPT4 / (dt * cell_volume))

def particle_absorption(
    radAreaP,
    dt,
    cell_volume,
    emissivity,
    G
):
    """
    Particle radiative absorption [W/m³].
    """
    # Projected particle area concentration Ap [1/m].
    radAreaP = np.asarray(radAreaP, dtype=float)
    cell_volume = np.asarray(cell_volume, dtype=float)
    Ap = radAreaP / (dt * cell_volume)
    
    ap = emissivity * np.asarray(Ap, dtype=float)
    G = np.asarray(G, dtype=float)

    return ap * G

def get_rad_source(internalMesh, emissivity, dt):
    internalMesh = internalMesh.compute_cell_sizes(length=False, area=False, volume=True)
    internalMesh['particle_emission_coal'] = particle_emissive_power(radAreaPT4=internalMesh['coalCloud:radAreaPT4'],dt=dt, cell_volume=internalMesh.cell_data["Volume"], emissivity=emissivity)
    internalMesh['particle_absorption_coal'] = particle_absorption(radAreaP=internalMesh['coalCloud:radAreaP'], dt=dt, cell_volume=internalMesh.cell_data["Volume"], emissivity=emissivity, G=internalMesh['G'])
    internalMesh['particle_emission_ash'] = particle_emissive_power(radAreaPT4=internalMesh['ashCloud:radAreaPT4'],dt=dt, cell_volume=internalMesh.cell_data["Volume"], emissivity=0.7)
    internalMesh['particle_absorption_ash'] = particle_absorption(radAreaP=internalMesh['ashCloud:radAreaP'], dt=dt, cell_volume=internalMesh.cell_data["Volume"], emissivity=0.7, G=internalMesh['G'])
    internalMesh['particle_emission'] = internalMesh['particle_emission_coal'] + internalMesh['particle_emission_ash']
    internalMesh['particle_absorption'] = internalMesh['particle_absorption_coal'] + internalMesh['particle_absorption_ash']
    return internalMesh['ShGas'] + internalMesh['particle_absorption'] - internalMesh['particle_emission'] # ShGas = gas_absorption - gas_emission

def plot_report(result, reference, output_path: str, label: str):
    fig, ax = plt.subplots(nrows=1, ncols=1, tight_layout=True)
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    x = pd.Series(np.asarray(line1.points[:,2]))
    y = pd.Series(np.asarray(line1["rad_source"]/1000))
    ax.plot(x, y, "-",label=label)
    ax.plot(reference["z"], reference["Sh"], '--', label=f"Referenz Gronarz 2017")
    ax.set_xlabel(r'$z$ in m')
    ax.legend()
    fig.set_size_inches(16, 4)
    fig.savefig(output_path)

def set_epsilon0(epsilon0,testdir):
    path = pathlib.Path(f"{testdir}/constant/coalCloudProperties")
    text = path.read_text()

    pattern = r'^(\s*epsilon0\s+)[0-9.+-eE]+(\s*;?)$'
    if not re.search(pattern, text, flags=re.MULTILINE):
        sys.exit("epsilon0 not found")

    text_new = re.sub(
        pattern,
        rf'\g<1>{epsilon0}\2;',
        text,
        flags=re.MULTILINE
    )
    path.write_text(text_new)

def set_dParticle(dParticle,testdir):
    path = pathlib.Path(f"{testdir}/0.org/lagrangian/coalCloud/d")
    text = path.read_text()

    pattern = r'\{2e-2\}'

    if not re.search(pattern, text):
        sys.exit("'{'" + "2e-2" +"}' not found")

    text_new = re.sub(
        pattern,
        dParticle,
        text
    )
    
    path.write_text(text_new)

def set_deltaT(deltaT, testdir):
    path = pathlib.Path(f"{testdir}/0.org/lagrangian/coalCloud/T")
    lines = path.read_text().splitlines(keepends=True)
    
    inside = False
    found_block = False
    new_lines = []

    for line in lines:

        stripped = line.strip()

        if stripped == "(":
            inside = True
            found_block = True
            new_lines.append(line)
            continue

        if stripped == ")":
            inside = False
            new_lines.append(line)
            continue

        if inside and stripped:
            value = float(stripped)
            value += deltaT - 100
            newline = "\n" if line.endswith("\n") else ""
            new_lines.append(f"{value:g}{newline}")
        else:
            new_lines.append(line)

    if not found_block:
        sys.exit("Opening '(' not found")

    if inside:
        sys.exit("Opening '(' found, but closing ')' not found")

    path.write_text("".join(new_lines))

def test_standard(testdir, write_report):
    # standard T particle = 100, epsilon0 = 0.29, d=2e-2    
    epsilon0 = 0.29
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=epsilon0, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/epsi_coal029.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/eps029.png", label="eps = 0.29")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/epsi_coal029.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")

def test_eps0(testdir, write_report):
    # set epsilon0
    epsilon0 = 0
    set_epsilon0(epsilon0,testdir)
    
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=epsilon0, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/epsi_coal0.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/eps0.png", label="eps = 0")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/epsi_coal0.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")

def test_eps1(testdir, write_report):
    # set epsilon0
    epsilon0 = 1
    set_epsilon0(epsilon0,testdir)
    
    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=epsilon0, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/epsi_coal1.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/eps1.png", label="eps = 1")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/epsi_coal1.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")

def test_Ap025(testdir, write_report):
    # set d
    dParticle='{1e-2}'
    set_dParticle(dParticle,testdir)    

    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=0.29, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/Ap025.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/Ap025.png", label="Ap = 0.25")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/Ap025.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")

def test_Ap4(testdir, write_report):
    # set d
    dParticle='{4e-2}'
    set_dParticle(dParticle,testdir)    

    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=0.29, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/Ap4.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/Ap4.png", label="Ap = 4")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/Ap4.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")

def test_deltaT0(testdir, write_report):
    # set T particle
    deltaT = 0
    set_deltaT(deltaT=deltaT, testdir=testdir) 

    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=0.29, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/deltaT0.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/deltaT0.png", label="deltaT = 0 K")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/deltaT0.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")

def test_deltaT200(testdir, write_report):
    # set T particle
    deltaT = 200
    set_deltaT(deltaT=deltaT, testdir=testdir) 

    # run OpenFOAM
    subprocess.run(["bash", "./Allrun"], cwd=testdir, check=True)

    # load results
    time = 0.002
    dt = 0.001
    reader = pv.POpenFOAMReader(pathlib.Path(f"{testdir}/case.foam"))
    reader.set_active_time_value(time)
    mesh = reader.read()
    result = mesh["internalMesh"]
    result['rad_source'] = get_rad_source(internalMesh=result, emissivity=0.29, dt=dt)
    reference = pd.read_parquet(os.path.join(TEST_DIR, "reference/Gronarz/deltaT200.parq"))
    line1 = result.sample_over_line((0, 10, 0), (0, 10, 40),resolution=1000)
    
    # plot results
    if write_report:
        plot_report(result=result, reference=reference, output_path=f"{FIGURES_DIR}/deltaT200.png", label="deltaT = 200 K")

    # test
    reference_line1 = pd.read_parquet(os.path.join(TEST_DIR, "reference/foam/deltaT200.parq"))
    rtol=1e-02
    atol=0.01
    np.testing.assert_allclose(line1["rad_source"]/1000, reference_line1["Sh"], atol=atol, rtol=rtol, err_msg=f"Test failed for line1 field Sh")