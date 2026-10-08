/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    flameletFoam

Group
    grpCombustionSolvers

Description
    Solver for reactive gas-phase simulations with tabulated chemistry.

\*---------------------------------------------------------------------------*/

#include "fvCFD.H"

#include "pimpleControl.H"
#include "fvOptions.H"
#include "fvcSmooth.H"

// Turbulence model
#include "turbulentFluidThermoModel.H"

#include "multivariateScheme.H"
#include "pressureControl.H"
#include "localEulerDdtScheme.H"

// STFS inhouse libraries
#include "flameletThermo.H"
#include "EquationHandler.H"
#include "CloudHandler.H"
#include "canteraThermo.H"

// Radiation model
#include "radiationModel.H"

#define FLOW
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{
    argList::addNote(
        "Solver for combustion using tabulated chemistry.");

#include "postProcess.H"

#include "addCheckCaseOptions.H"
#include "setRootCaseLists.H"
#include "createTime.H"
#include "createMesh.H"
#include "createControl.H"
#include "initContinuityErrs.H"
#include "createFields.H"
#include "createInputVariableEquations.H"
#include "createRhoUfIfPresent.H"
#include "createTimeControls.H"

    turbulence->validate();

    // if (!LTS)
    //{
#include "compressibleCourantNo.H"
#include "setInitialDeltaT.H"
    //}


    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    if (pimple.consistent())
    {
        Info << "Using pimple consistent" << endl;
    }
    Info << "\nStarting time loop\n"
         << endl;

    scalar executeInterval =
        controlDict.lookupOrDefault(
            "deleteParticlesInterval",
            readScalar(controlDict.lookup("writeInterval")));
    scalar nextExecuteTime = runTime.value() + executeInterval;

    while (runTime.run())
    {
#include "readTimeControls.H"

        // if (LTS)
        //{
        //     #include "setRDeltaT.H"
        // }
        // else
        //{
#include "compressibleCourantNo.H"
#include "setDeltaT.H"
        //}

        ++runTime;

        Info << "Time = " << runTime.timeName() << nl << endl;

        cloudHandler.evolveClouds(); // update particle cloud

#include "rhoEqn.H"

        // --- Pressure-velocity PIMPLE corrector loop
        while (pimple.loop())
        {
#ifdef FLOW
#include "UEqn.H"
#endif

            // Call update functions for flamelet
            equationHandler.solve();
            Info << "Equations solved" << endl;
            flameletThermo.update();
            thermo.alpha() = flameletThermo.alpha();

            Info << "Flamelet update successful" << endl;
            radiation->correct();
            // --- Pressure corrector loop
            while (pimple.correct())
            {
                if (lowMachPressure)
                {
                    // Default option.
                    // Recommended for tables with constant pressure.
#include "pdEqn.H"
                }
                else
                {
                    // Necessary for tables with pressure
                    // Needs to be activated by setting the keyword "lowMach" to false in the pimple dict (fvSolution).
                    // lowMach false;
                    if (pimple.consistent())
                    {
#include "pcEqn.H"
                    }
                    else
                    {
#include "pEqn.H"
                    }
                }
            }

            if (pimple.turbCorr())
            {
                // rho = flameletThermo.rho();
                turbulence->correct();
            }
        }

        rho = psi * p;
        // runTime.write();
        if (runTime.value() >= nextExecuteTime)
        {
            if (controlDict.lookupOrDefault("deleteParticles", false))
            {
                cloudHandler.deleteParticles(controlDict.subDict("deleteParticlesDict"));
            }
            nextExecuteTime += executeInterval;
        }
        if (runTime.writeTime())
        {
            cloudHandler.particleProcessorDistributionInfo();

            runTime.write();
            flameletThermo.getField("ha").write();
        };

        Info << "T gas min/max   = " << min(T).value() << ", " << max(T).value() << endl;
        Info << "rho min/max : " << min(rho).value() << " " << max(rho).value() << endl;
        // Info << "p gas min/max   = " << min(p).value() << ", " << max(p).value() << endl;
        Info << "p or pd min/max   = " << min(*pd).value() << ", " << max(*pd).value() << endl;
        Info << "mag U min/max   = " << min(mag(U)).value() << ", " << max(mag(U)).value() << endl;

        runTime.printExecutionTime(Info);
    }

    Info << "End\n"
         << endl;

    return 0;
}


// ************************************************************************* //
