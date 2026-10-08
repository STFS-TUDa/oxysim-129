/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  |
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

#include "wsggAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "unitConversion.H"
#include "basicSpecieMixture.H"
#include "thermoPhysicsTypes.H"
#include "reactingMixture.H"
#include "surfaceFields.H"
#include "symmetryFvPatch.H"
#include "cyclicFvPatch.H"
#include "processorFvPatch.H"
#include "zeroGradientFvPatchFields.H"
#include "fixedInternalValueFvPatchFields.H"


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
defineTypeNameAndDebug(wsggAbsorptionEmission, 0);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggAbsorptionEmission::wsggAbsorptionEmission(
    const dictionary& dict,
    const fvMesh& mesh)
: absorptionEmissionModel(dict, mesh),
  thermoPtr_(nullptr)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggAbsorptionEmission::~wsggAbsorptionEmission()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::wsggAbsorptionEmission::setThermoPtr(
    const basicThermo& thermo)
{
    // make absorptionEmission model thermo aware
    thermoPtr_ = &thermo;
}

// ************************************************************************* //
