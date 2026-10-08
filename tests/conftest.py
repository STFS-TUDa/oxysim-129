"""!
In the conftest.py, custom fixtures for the usage in pytest testcases are defined. This part of the documentation handles these custom definitions.
For general information on how to write testcases, see the [testing documentation](@ref testing).
For general information on how fixtures work, see the [pytest fixtures documentation](https://docs.pytest.org/en/stable/how-to/fixtures.html)
"""
import pytest
import pathlib
import shutil
import os
import matplotlib.pyplot as plt
import numpy as np

label_dict = {"x": "$x$ (m)", 
              "y": "$y$ (m)", 
              "z": "$z$ (m)", 
              "T": "$T$ (K)", 
              "U_0": "$u$ (m/s)", 
              "p": "$p$ (Pa)", 
              "pd": "$p_\\mathrm{dynamic}$ (Pa)"
              }

def pytest_addoption(parser):
    parser.addoption("--write_report", action="store", default="True")
    parser.addoption("--clean", action="store", default="False")

@pytest.fixture(scope="session")
def write_report(pytestconfig):
    return pytestconfig.getoption("--write_report")=="True"

@pytest.fixture(scope="session")
def clean(pytestconfig):
    """!
    The clean fixture provides an easy access point for the `--clean` flag. If pytest is called with this flag (e.g. `pytest test_mystuff.py --clean=True` in the terminal), the value for `clean` will evaluate to true.

    EXAMPLE:

    ```python
    import pytest

    def test_something(testdir, clean):
        if clean==True:
            print("--clean=True has been used")
        elif clean==False:
            print("--clean=True has not been used")
    ```
    """
    return pytestconfig.getoption("--clean")=="True"

@pytest.fixture
def testdir(request, clean):
    """!
    This fixture creates a temporary copy of the `foam_template` subdirectory next to the pytest. The copy is placed inside a `temp` directory, which is created as necessary. The fixture makes sure that
    - Any remnant data from previous runs is deleted
    - A fresh copy of `foam_template` is created inside `temp`, named after the function calling the fixture.
    - The OpenFOAM run is executed inside this subdirectory
    - The data from this run is deleted, except for cases with very hard errors
    - the empty `temp` directory will remain, but empty directories are ignored by git anyways

    NOTE: Using this function, the OpenFOAM run is executed **once for every test function**. All tests which request the result data will obtain a distinctly generated object. Changes made to the data by one test can therefore never influence the results of the other tests. This comes at the cost of running a simulation several times.

    The `testdir_module` fixture provides similar functionality with a module scope, creating and running only one simulation for all tests in one test module.

    EXAMPLE:

    Calling the fixture like

    ```python
    def test_heat_conservation(testdir, write_report):
        pass

    def test_mass_conservation(testdir, write_report):
        pass
    ```
    in a test directory such as

    ```
    ├── test_mystuff.py
    └── foam_template
        ├── input_file_1
        └── input_file_2
    ```
    
    will create new temporary working directories

    ```
    ├── test_mystuff.py
    ├── foam_template
    |   ├── input_file_1
    |   └── input_file_2
    └── temp
        └── test_heat_conservation
        |   ├── input_file_1
        |   └── input_file_2
        └── test_mass_conservation
            ├── input_file_1
            └── input_file_2
    ```
    """
    # Get the folder where the test_*.py file is located
    TEST_FILE_DIR = os.path.dirname(request.fspath)
    # Create testfolder root, if it does not exist
    TEST_FOLDER_ROOT = f"{TEST_FILE_DIR}/temp"
    pathlib.Path(f"{TEST_FOLDER_ROOT}").mkdir(parents=True, exist_ok=True)

    # Get the name of the function being executed
    FOLDER_NAME = request.node.originalname
    TEST_DIR = os.path.join(f"{TEST_FOLDER_ROOT}", FOLDER_NAME)

    # Clean up, if folder already exists
    if os.path.exists(TEST_DIR):
        shutil.rmtree(TEST_DIR)
    shutil.copytree(f'{TEST_FILE_DIR}/foam_template', TEST_DIR)

    # Return value of fixture
    yield TEST_DIR

    # Cleanup after test
    if clean:
        shutil.rmtree(TEST_DIR)

@pytest.fixture(scope="module")
def testdir_module(request, clean):
    """!
    The `testdir_module` fixture is similar to `testdir`, but has a module scope. It creates a temporary copy of the `foam_template` subdirectory next to the pytest. The copy is placed inside a `temp` directory, which is created as necessary. The fixture makes sure that
    - Any remnant data from previous runs is deleted
    - A fresh copy of `foam_template` is created inside `temp`, named after the module calling the fixture.
    - The OpenFOAM run is executed inside this subdirectory
    - The data from this run is deleted, except for cases with very hard errors
    - the empty `temp` directory will remain, but empty directories are ignored by git anyways

    NOTE: Using this function, the OpenFOAM run is executed **only once for the entire module**. All tests which request the result data will obtain a reference on the exact same object. Changes made to the data by one test can therefore influence the results of the other tests. 

    EXAMPLE:

    Calling the fixture like

    ```python
    def test_heat_conservation(testdir, write_report):
        pass

    def test_mass_conservation(testdir, write_report):
        pass
    ```

    in a test directory such as

    ```
    ├── test_mystuff.py
    └── foam_template
        ├── input_file_1
        └── input_file_2
    ```

    will create a new temporary working directory

    ```
    ├── test_mystuff.py
    ├── foam_template
    |   ├── input_file_1
    |   └── input_file_2
    └── temp
        └── test_mystuff
            ├── input_file_1
            └── input_file_2
    ```
    """
    # Get the folder where the test_*.py file is located
    TEST_FILE_DIR = os.path.dirname(request.fspath)
    # Create testfolder root, if it does not exist
    TEST_FOLDER_ROOT = f"{TEST_FILE_DIR}/temp"
    pathlib.Path(f"{TEST_FOLDER_ROOT}").mkdir(parents=True, exist_ok=True)

    # Get the name of the function being executed
    FOLDER_NAME = request.module.__name__
    TEST_DIR = os.path.join(f"{TEST_FOLDER_ROOT}", FOLDER_NAME)

    # Clean up, if folder already exists
    if os.path.exists(TEST_DIR):
        shutil.rmtree(TEST_DIR)
    shutil.copytree(f'{TEST_FILE_DIR}/foam_template', TEST_DIR)

    # Return value of fixture
    yield TEST_DIR

    # Cleanup after test
    if clean:
        shutil.rmtree(TEST_DIR)

def plot_line_profiles(result, reference, xv, yvs, filename, title="", ncols_max=4, labels=["reactiveFoam", "Reference"], show=False):
    yvs = [yv for yv in yvs if xv != yv]

    nvars = len(yvs)
    if nvars > ncols_max:
        ncols = ncols_max
        nrows = int((nvars - 1) / ncols) + 1
    else:
        ncols = nvars
        nrows = 1
    
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, sharex=True)
    axs = axs.flatten()

    for yv, ax in zip(yvs, axs):
        reference.plot(xv, yv, ax=ax, legend=False, color="k", label=labels[1])
        result.plot(xv, yv, ax=ax, legend=False, linestyle="-.", color="C1", label=labels[0])
        ax.set_xlabel(label_dict.get(xv, xv))
        ax.set_ylabel(label_dict.get(yv, yv))
    axs[0].legend()

    fig.suptitle(title)
    fig.set_size_inches(ncols*4,nrows*3)
    fig.tight_layout()
    fig.savefig(f"{filename}")
    if not show:
        plt.close()

def plot_line_profiles_time_series(result, reference, xv, yvs, filename, labels=["reactiveFoam", "Reference"], title="", ncols_max=4, time_name="t"):
    times = result[f"{time_name}"].unique()
    yvs = [yv for yv in yvs if xv != yv]

    nvars = len(yvs)
    if nvars > ncols_max:
        ncols = ncols_max
        nrows = int((nvars - 1) / ncols) + 1
    else:
        ncols = nvars
        nrows = 1
    
    fig, axs = plt.subplots(ncols=ncols, nrows=nrows, sharex=True)
    axs = axs.flatten()

    for t in times:
        tmp_reference = reference[np.isclose(reference["t"],t)].sort_values("x", ignore_index=True)
        tmp_result = result[np.isclose(result["t"],t)].sort_values("x", ignore_index=True)
        for yv, ax in zip(yvs, axs):
            tmp_reference.plot(xv, yv, ax=ax, legend=False, color="k", label=labels[0])
            tmp_result.plot(xv, yv, ax=ax, legend=False, linestyle="-.", color="C1", label=labels[1])
            ax.set_xlabel(label_dict.get(xv, xv))
            ax.set_ylabel(label_dict.get(yv, yv))
        axs[0].legend(labels)

    fig.suptitle(title)
    fig.set_size_inches(ncols*4,nrows*3)
    fig.tight_layout()
    fig.savefig(f"{filename}")
    plt.close()
