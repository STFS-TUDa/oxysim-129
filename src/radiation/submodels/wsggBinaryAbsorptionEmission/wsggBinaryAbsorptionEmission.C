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

\*---------------------------------------------------------------------------*/

#include "wsggBinaryAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "wsggAbsorptionEmission.H"
// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
defineTypeNameAndDebug(wsggBinaryAbsorptionEmission, 0);

addToRunTimeSelectionTable(
    absorptionEmissionModel,
    wsggBinaryAbsorptionEmission,
    dictionary);

// Construct model1 and verify that it is a wsggAbsorptionEmission
static autoPtr<wsggAbsorptionEmission> newWsggAbsorptionEmissionModel(
    const dictionary& dict,
    const fvMesh& mesh)
{
    autoPtr<absorptionEmissionModel> absorptionEmissionModelPtr(
        absorptionEmissionModel::New(dict, mesh));

    wsggAbsorptionEmission* wsggAbsorptionEmissionModelPtr =
        dynamic_cast<wsggAbsorptionEmission*>(absorptionEmissionModelPtr.get());

    if (!wsggAbsorptionEmissionModelPtr)
    {
        FatalIOErrorInFunction(dict)
            << "model1 of " << wsggBinaryAbsorptionEmission::typeName
            << " must be of type " << wsggAbsorptionEmission::typeName
            << "but is of type " << absorptionEmissionModelPtr->type() << nl
            << "This is necessary for combining it with particle absorption emission models." << nl
            << exit(FatalIOError);
    }

    // Cast succeeded: transfer ownership
    absorptionEmissionModelPtr.release();
    return autoPtr<wsggAbsorptionEmission>(wsggAbsorptionEmissionModelPtr);
}

} // End namespace radiation
} // End namespace Foam


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggBinaryAbsorptionEmission::wsggBinaryAbsorptionEmission(
    const dictionary& dict,
    const fvMesh& mesh)
: absorptionEmissionModel(dict, mesh),
  coeffsDict_(dict.optionalSubDict(typeName + "Coeffs")),
  model1_(newWsggAbsorptionEmissionModel(coeffsDict_.subDict("model1"), mesh)),
  model2_(
      absorptionEmissionModel::New(coeffsDict_.subDict("model2"), mesh)),
  nBands_(model1_->nBands()),
  thermoPtr_(nullptr)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggBinaryAbsorptionEmission::~wsggBinaryAbsorptionEmission()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBinaryAbsorptionEmission::aCont(const label bandI) const
{
    return model1_->aCont(bandI) + model2_->aCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBinaryAbsorptionEmission::aDisp(const label bandI) const
{
    return model1_->aDisp(bandI) + model2_->aDisp(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBinaryAbsorptionEmission::eCont(const label bandI) const
{
    return model1_->eCont(bandI) + model2_->eCont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBinaryAbsorptionEmission::eDisp(const label bandI) const
{
    return model1_->eDisp(bandI) + model2_->eDisp(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBinaryAbsorptionEmission::ECont(const label bandI) const
{
    return model1_->ECont(bandI) + model2_->ECont(bandI);
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBinaryAbsorptionEmission::EDisp(const label bandI) const
{
    return model1_->EDisp(bandI) + model2_->EDisp(bandI);
}

void Foam::radiation::wsggBinaryAbsorptionEmission::correct(
    Foam::volScalarField& aa,
    PtrList<Foam::volScalarField>& aLambda

) const
{
    model1_->correct(aa, aLambda);
    model2_->correct(aa, aLambda);
}

void Foam::radiation::wsggBinaryAbsorptionEmission::setThermoPtr(
    const basicThermo& thermo)
{
    // make absorptionEmission model thermo aware
    thermoPtr_ = &thermo;
    if (isA<wsggAbsorptionEmission>(model1_()))
    {
        model1_->setThermoPtr(*thermoPtr());
    }
    if (isA<wsggAbsorptionEmission>(model2_()))
    {
        dynamic_cast<wsggAbsorptionEmission&>(model2_())
            .setThermoPtr(*thermoPtr());
    }
}

// ************************************************************************* //
