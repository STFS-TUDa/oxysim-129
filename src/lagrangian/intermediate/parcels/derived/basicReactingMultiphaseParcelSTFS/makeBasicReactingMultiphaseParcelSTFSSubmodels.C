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

#include "basicReactingMultiphaseCloud.H"     // typedef of basicReactingMultiphaseCloud
#include "basicReactingMultiphaseCloudSTFS.H" // typedef of basicReactingMultiphaseCloudSTFS

// // Definition of necessary Macros (for buidling the C files)
#include "makeReactingParcelCloudFunctionObjects.H"     // Reacting variant
#include "makeReactingParcelSTFSCloudFunctionObjects.H" // Reacting variant

// Kinematic
#include "makeThermoParcelForces.H" // thermo variant
#include "makeParcelDispersionModels.H"
#include "makeReactingMultiphaseParcelInjectionModels.H" // MP variant
#include "makeParcelPatchInteractionModels.H"
#include "makeReactingMultiphaseParcelStochasticCollisionModels.H" // MP variant
#include "makeReactingParcelSurfaceFilmModels.H"                   // Reacting variant
#include "makeReactingMultiphaseParcelSTFSInjectionModels.H"       //Custom Injection Models

// Thermodynamic
#include "makeParcelHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSTFSPhaseChangeModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactingParcelCloudFunctionObjects(basicReactingMultiphaseCloudSTFS);
makeSTFSParcelCloudFunctionObjects(basicReactingMultiphaseCloudSTFS);

// Kinematic sub-models
makeThermoParcelForces(basicReactingMultiphaseCloudSTFS);
makeParcelDispersionModels(basicReactingMultiphaseCloudSTFS);
makeReactingMultiphaseParcelInjectionModels(basicReactingMultiphaseCloudSTFS);
makeReactingMultiphaseParcelSTFSInjectionModels(basicReactingMultiphaseCloudSTFS);
makeParcelPatchInteractionModels(basicReactingMultiphaseCloudSTFS);
makeReactingMultiphaseParcelStochasticCollisionModels(basicReactingMultiphaseCloudSTFS);
makeReactingParcelSurfaceFilmModels(basicReactingMultiphaseCloudSTFS);

// Thermo sub-models
makeParcelHeatTransferModels(basicReactingMultiphaseCloudSTFS);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels(basicReactingMultiphaseCloudSTFS);
makeReactingParcelPhaseChangeModels(basicReactingMultiphaseCloudSTFS);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelDevolatilisationModels(basicReactingMultiphaseCloudSTFS);
makeReactingMultiphaseParcelSurfaceReactionModels(basicReactingMultiphaseCloudSTFS);

// MPPIC sub-models
makeMPPICParcelDampingModels(basicReactingMultiphaseCloudSTFS);
makeMPPICParcelIsotropyModels(basicReactingMultiphaseCloudSTFS);
makeMPPICParcelPackingModels(basicReactingMultiphaseCloudSTFS);

// ************************************************************************* //
