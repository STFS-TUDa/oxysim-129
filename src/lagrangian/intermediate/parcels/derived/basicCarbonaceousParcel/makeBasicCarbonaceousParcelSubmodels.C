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


#include "basicCarbonaceousParcel.H"
#include "PtrDynList.H"
#include "basicCarbonaceousCloud.H" // typedef of basicCarbonaceousCloud


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
#include "makeParcelSTFSHeatTransferModels.H"

// Reacting
#include "makeReactingMultiphaseParcelCompositionModels.H" // MP Variant
#include "makeReactingParcelPhaseChangeModels.H"
#include "makeReactingParcelSTFSPhaseChangeModels.H"

// Reacting multiphase
#include "makeReactingMultiphaseParcelDevolatilisationModels.H"
#include "makeReactingMultiphaseParcelSurfaceReactionModels.H"

// Carbonaceous
#include "makeCarbonaceousParcelDevolatilisationModels.H"
#include "makeCarbonaceousParcelHeatCapacityModels.H"
#include "makeCarbonaceousParcelSurfaceReactionModels.H"

// MPPIC sub-models
#include "makeMPPICParcelDampingModels.H"
#include "makeMPPICParcelIsotropyModels.H"
#include "makeMPPICParcelPackingModels.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makeReactingParcelCloudFunctionObjects(basicCarbonaceousCloud);
makeSTFSParcelCloudFunctionObjects(basicCarbonaceousCloud);

// Kinematic sub-models
makeThermoParcelForces(basicCarbonaceousCloud);
makeParcelDispersionModels(basicCarbonaceousCloud);
makeReactingMultiphaseParcelInjectionModels(basicCarbonaceousCloud);
makeReactingMultiphaseParcelSTFSInjectionModels(basicCarbonaceousCloud);
makeParcelPatchInteractionModels(basicCarbonaceousCloud);
makeReactingMultiphaseParcelStochasticCollisionModels(
    basicCarbonaceousCloud);
makeReactingParcelSurfaceFilmModels(basicCarbonaceousCloud);

// Thermo sub-models
makeParcelHeatTransferModels(basicCarbonaceousCloud);

// Reacting sub-models
makeReactingMultiphaseParcelCompositionModels(
    basicCarbonaceousCloud);
makeReactingParcelPhaseChangeModels(basicCarbonaceousCloud);

// Reacting multiphase sub-models
makeReactingMultiphaseParcelDevolatilisationModels(
    basicCarbonaceousCloud);

makeReactingMultiphaseParcelSurfaceReactionModels(
    basicCarbonaceousCloud);

// Carbonaceous sub-models
makeCarbonaceousDevolatilisationModels(basicCarbonaceousCloud);
makeCarbonaceousParcelHeatCapacityModels(basicCarbonaceousCloud);
makeCarbonaceousParcelSurfaceReactionModels(basicCarbonaceousCloud);

// MPPIC sub-models
makeMPPICParcelDampingModels(basicCarbonaceousCloud);
makeMPPICParcelIsotropyModels(basicCarbonaceousCloud);
makeMPPICParcelPackingModels(basicCarbonaceousCloud);

// ************************************************************************* //
