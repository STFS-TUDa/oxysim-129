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

#include "wsggBordbarAbsorptionEmission.H"
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
#include <algorithm>


// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
defineTypeNameAndDebug(wsggBordbarAbsorptionEmission, 0);

addToRunTimeSelectionTable(
    absorptionEmissionModel,
    wsggBordbarAbsorptionEmission,
    dictionary);
}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggBordbarAbsorptionEmission::wsggBordbarAbsorptionEmission(
    const dictionary& dict,
    const fvMesh& mesh)
: wsggAbsorptionEmission(dict, mesh),
  coeffsDict_(dict.subDict(typeName + "Coeffs")),
  d_(coeffsDict_.lookup("dCoeffs")),
  c1_(coeffsDict_.lookup("c1Coeffs")),
  c2_(coeffsDict_.lookup("c2Coeffs")),
  c3_(coeffsDict_.lookup("c3Coeffs")),
  c4_(coeffsDict_.lookup("c4Coeffs")),
  mesh_(mesh),
  // These variables are only required for grey implementation. This is the non-gray implementation, so we don't use these.
  //     pathLength_(coeffsDict_.lookup("pathLength")),
  //     meanBeamPathAutoCalcMode_(coeffsDict_.lookupOrDefault<bool>("meanBeamPathAutoCalcMode", false)),
  //     sector_(coeffsDict_.lookup("sector")),
  nBands_(d_[0].size()),
  efields(nBands_),
  afields(nBands_)

{

    Info << "Number of grey gasses available (Excluding the clear gas): " << nBands_ - 1 << "\n";

    Cijk.append(c1_);
    Cijk.append(c2_);
    Cijk.append(c3_);
    Cijk.append(c4_);
    Info << "Cijk: \n"
         << Cijk << "\n";


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

Foam::radiation::wsggBordbarAbsorptionEmission::~wsggBordbarAbsorptionEmission()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBordbarAbsorptionEmission::aCont(const label bandI) const
{
    return afields[bandI];
}

Foam::tmp<Foam::volScalarField> // This is actually the weight, which is required for the emissitivy. Emissitivy = Weight x Absorptivity --->>>  eCont x aCont
Foam::radiation::wsggBordbarAbsorptionEmission::eCont(const label bandI) const
{
    return efields[bandI];
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggBordbarAbsorptionEmission::ECont(const label bandI) const
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

void Foam::radiation::wsggBordbarAbsorptionEmission::correct // Correct updates the a and e fields. They are kept in memmory and accessed with aCont and eCont
    (
        Foam::volScalarField& aa,
        PtrList<Foam::volScalarField>& aLambda

    ) const
{

    const basicThermo* thermoPtr = this->thermoPtr();
    if (!thermoPtr)
    {
        FatalErrorInFunction
            << "WsggAbsorptionEmission model is not thermo aware and has to be initialised first."
            << exit(FatalError);
    }
    const basicThermo& thermo_ = *thermoPtr;
    const volScalarField& T = thermo_.T();
    const volScalarField& p = thermo_.p();

    if (!isA<basicSpecieMixture>(thermo_))
    // calculate W based on the species
    {
        // not implemented yet
    }
    else
    // access W directly from the table
    {
    }

    // Access mixture properties (T,p,mass fractions)
    const volScalarField& CO2 = mesh_.lookupObject<volScalarField>("CO2");
    const volScalarField& H2O = mesh_.lookupObject<volScalarField>("H2O");
    const volScalarField& W = thermo_.W(); // mean molecular weight

    // Loop all cells
    forAll(afields[0], cellI)
    {
        //--------      Calculate Mr
        scalar XkCO2 = CO2[cellI] * W[cellI] / WCO2_;
        scalar XkH2O = H2O[cellI] * W[cellI] / WH2O_;
        scalar Mr = XkH2O / (XkCO2 + 1e-18);
        //---------     End Mr

        afields[0][cellI] = 0.0;
        scalarList weightList = getWeights(T[cellI], Mr);
        for (label bandI = 1; bandI < nBands_; bandI++)
        {
            afields[bandI][cellI] = 0;
            efields[bandI][cellI] = weightList[bandI];
            // Eq 11 in Bordbar. Absorption coeff
            forAll(d_[bandI - 1], j)
            {
                afields[bandI][cellI] += p[cellI] / 101325.0 * (XkCO2 + XkH2O) * (d_[bandI - 1][j] * pow(Mr, (j)));
            }
        }
        efields[0][cellI] = weightList[0];

    } // Finished looping all cells

    // Loop over all boundary faces, and correct boundary weights, to be used in BC
    forAll(mesh().boundaryMesh(), patchI)
    {

        forAll(efields[0].boundaryFieldRef()[patchI], faceI)
        {
            //--------      Calculate Mr
            scalar XkCO2 = CO2.boundaryField()[patchI][faceI] * W.boundaryField()[patchI][faceI] / WCO2_;
            scalar XkH2O = H2O.boundaryField()[patchI][faceI] * W.boundaryField()[patchI][faceI] / WH2O_;
            scalar Mr = XkH2O / (XkCO2 + 1e-18);
            //---------     End Mr

            scalarList weightList = getWeights(T.boundaryField()[patchI][faceI], Mr);
            for (label bandI = 1; bandI < nBands_; bandI++)
            {
                efields[bandI].boundaryFieldRef()[patchI][faceI] = weightList[bandI];
            }
            efields[0].boundaryFieldRef()[patchI][faceI] = weightList[0];
        }
    }
    // }
}


Foam::scalarList Foam::radiation::wsggBordbarAbsorptionEmission::getWeights(
    Foam::scalar T,
    Foam::scalar Xr) const
{
    scalar dimlessT = min(T / 1200.0, 2400.0 / 1200.0);
    scalarList weightList {1.0, 0.0, 0.0, 0.0, 0.0};
    for (label bandI = 1; bandI < nBands_; bandI++)
    {
        forAll(Cijk[bandI - 1], j)
        {
            double b[Cijk[bandI - 1].size()];
            std::fill_n(b, Cijk[bandI - 1].size(), 0.0);

            forAll(Cijk[bandI - 1][j], k)
            {
                b[j] += Cijk[bandI - 1][j][k] * pow(Xr, (k)); // Fix for CO2=0!!
                // b = b*Mr + Cijk[bandI-1][j][k];  //Horner's method for efficient polynomial calculation
            }
            weightList[bandI] += b[j] * pow(dimlessT, j);
        }
        weightList[0] -= weightList[bandI];
    }
    return weightList;
}


// ************************************************************************* //
