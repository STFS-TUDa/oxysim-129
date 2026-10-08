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

#include "wsggDorigonAbsorptionEmission.H"
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
defineTypeNameAndDebug(wsggDorigonAbsorptionEmission, 0);

addToRunTimeSelectionTable(
    absorptionEmissionModel,
    wsggDorigonAbsorptionEmission,
    dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggDorigonAbsorptionEmission::wsggDorigonAbsorptionEmission(
    const dictionary& dict,
    const fvMesh& mesh)
: wsggAbsorptionEmission(dict, mesh),
  coeffsDict_(dict.subDict(typeName + "Coeffs")),
  bi_(coeffsDict_.lookup("biCoeffs")),
  ki_(coeffsDict_.lookup("kiCoeffs")),
  // These variables are only required for grey implementation. This is the non-gray implementation, so we don't use these.
  //     pathLength_(coeffsDict_.lookup("pathLength")),
  //     meanBeamPathAutoCalcMode_(coeffsDict_.lookupOrDefault<bool>("meanBeamPathAutoCalcMode", false)),
  //     sector_(coeffsDict_.lookup("sector")),
  nBands_(ki_.size()),
  efields(nBands_),
  afields(nBands_)
{
    Info << "Number of grey gasses available (Excluding the clear gas): " << nBands_ - 1 << "\n";

    forAll(efields, i)
    {
        efields.set(
            i,
            new volScalarField(
                IOobject(
                    "eCont" + name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("e", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)
                // zeroGradientFvPatchScalarField::typeName
                ));

        afields.set(
            i,
            new volScalarField(

                IOobject(
                    "aCont" + name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("a", dimensionSet(0, -1, 0, 0, 0, 0, 0), 0.0),
                zeroGradientFvPatchScalarField::typeName));
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggDorigonAbsorptionEmission::~wsggDorigonAbsorptionEmission()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggDorigonAbsorptionEmission::aCont(const label bandI) const
{
    return afields[bandI];
}

Foam::tmp<Foam::volScalarField> // This is actually the weight, which is required for the emissitivy. Emissitivy = Weight x Absorptivity --->>>  eCont x aCont
Foam::radiation::wsggDorigonAbsorptionEmission::eCont(const label bandI) const
{
    return efields[bandI];
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggDorigonAbsorptionEmission::ECont(const label bandI) const
{
    tmp<volScalarField> tE(
        new volScalarField(
            IOobject(
                "ECont" + name(bandI),
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("E", dimMass / dimLength / pow3(dimTime), 0.0)));

    return tE;
}

void Foam::radiation::wsggDorigonAbsorptionEmission::correct // Correct updates the a and e fields. They are kept in memmory and accessed with aCont and eCont
    (
        Foam::volScalarField& aa,
        PtrList<Foam::volScalarField>& aLambda) const
{
    const basicThermo* thermoPtr = this->thermoPtr();
    if (!thermoPtr)
    {
        FatalErrorInFunction
            << "WsggAbsorptionEmission model is not thermo aware and has to be initialised first."
            << exit(FatalError);
    }

    const basicThermo& thermo_ = *thermoPtr;
    // Access Mr directly from table
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();
    const volScalarField& CO2 = mesh_.lookupObject<volScalarField>("CO2");
    const volScalarField& H2O = mesh_.lookupObject<volScalarField>("H2O");
    const volScalarField& Wt = thermo_.W(); // mean molecular weight

    // Loop all cells
    forAll(afields[0], celli)
    {
        scalar XkCO2 = CO2[celli] * Wt[celli] / WCO2_;
        scalar XkH2O = H2O[celli] * Wt[celli] / WH2O_;
        scalar Mr = XkH2O / (XkCO2 + 1e-18);
        //--------      Calculate a (Absorption coeff)
        afields[0][celli] = 0.0;
        scalarList weightList = getWeights(T[celli], Mr);
        for (label bandI = 1; bandI < nBands_; bandI++)
        {
            afields[bandI][celli] = 0;
            efields[bandI][celli] = weightList[bandI];
            // equation 11 in Dorigon
            afields[bandI][celli] += p[celli] / 101325 * (XkCO2 + XkH2O) * (ki_[bandI]);
        }
        efields[0][celli] = weightList[0];


    } // Finished looping all cells


    // Loop over all boundary faces, and correct boundary weights, to be used in BC
    forAll(mesh().boundaryMesh(), patchI)
    {
        forAll(efields[0].boundaryFieldRef()[patchI], faceI)
        {
            scalarList weightList = getWeights(T.boundaryField()[patchI][faceI], 0.0); // weight does not depend on Mr in Dorigon, Hence 2nd parameter is 0.
            for (label bandI = 1; bandI < nBands_; bandI++)
            {
                efields[bandI].boundaryFieldRef()[patchI][faceI] = weightList[bandI];
            }
            efields[0].boundaryFieldRef()[patchI][faceI] = weightList[0];
        }
    }
}

Foam::scalarList Foam::radiation::wsggDorigonAbsorptionEmission::getWeights(
    Foam::scalar T,
    Foam::scalar Xr) const
{
    scalarList weightList {1.0, 0.0, 0.0, 0.0, 0.0};
    for (label bandI = 1; bandI < nBands_; bandI++)
    {
        // equation 11 in Dorigon
        forAll(bi_[bandI - 1], j)
        {
            weightList[bandI] += bi_[bandI - 1][j] * pow(T, j);
        }

        // Correct clear gas weight, sum must be 1
        weightList[0] -= weightList[bandI];
    }
    return weightList;
}
// ************************************************************************* //
