/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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

\*---------------------------------------------------------------------------*/

#include "basicReactingCloudSTFS.H" // typedef of basicReactingCloudSTFS

// // Definition of necessary Macros (for buidling the C files)
#include "makeReactingParcelCloudFunctionObjects.H"     // Reacting variant
#include "makeReactingParcelSTFSCloudFunctionObjects.H" // Reacting variant

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingParcelInjectionModels.H" // Reacting variant
#include "makeParcelPatchInteractionModels.H"
#include "makeParcelStochasticCollisionModels.H"
#include "makeReactingParcelSurfaceFilmModels.H" // Reacting variant

// Thermodynamic
#include "makeParcelHeatTransferModels.H"
#include "makeParcelSTFSHeatTransferModels.H"

// Reacting
#include "makeReactingParcelCompositionModels.H"
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSTFSPhaseChangeModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactingParcelCloudFunctionObjects(basicReactingCloudSTFS);
makeSTFSParcelCloudFunctionObjects(basicReactingCloudSTFS);

// Kinematic sub-models
// Note:
// We do not have to build this, since we are using the
// standard KinamticParcel Template class
// for which the models where still build during OF source Code build
makeThermoParcelForces(basicReactingCloudSTFS);
makeParcelDispersionModels(basicReactingCloudSTFS);
makeReactingParcelInjectionModels(basicReactingCloudSTFS);
makeParcelPatchInteractionModels(basicReactingCloudSTFS);
makeParcelStochasticCollisionModels(basicReactingCloudSTFS);
makeReactingParcelSurfaceFilmModels(basicReactingCloudSTFS);

// Thermo sub-models
makeParcelHeatTransferModels(basicReactingCloudSTFS);
makeParcelSTFSHeatTransferModels(basicReactingCloudSTFS);

// Reacting sub-models
// Note:
// We need to build all here, since we are using a new ReactingParcel Template class
// (ReactingParcelSTFS)
makeReactingParcelCompositionModels(basicReactingCloudSTFS);
makeReactingParcelPhaseChangeModels(basicReactingCloudSTFS);
makeReactingParcelSTFSPhaseChangeModels(basicReactingCloudSTFS);

// MPPIC sub-models
makeMPPICParcelDampingModels(basicReactingCloudSTFS);
makeMPPICParcelIsotropyModels(basicReactingCloudSTFS);
makeMPPICParcelPackingModels(basicReactingCloudSTFS);


// ************************************************************************* //
