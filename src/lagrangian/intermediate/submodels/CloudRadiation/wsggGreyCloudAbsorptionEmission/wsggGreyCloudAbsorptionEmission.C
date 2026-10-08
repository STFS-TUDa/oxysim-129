/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

#include "wsggGreyCloudAbsorptionEmission.H"
#include "addToRunTimeSelectionTable.H"
#include "radiation.H"
#include "wsggAbsorptionEmission.H"
#include "wsggBinaryAbsorptionEmission.H"
#include "CloudHandler.H"
#include "basicCarbonaceousCloud.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{
defineTypeNameAndDebug(wsggGreyCloudAbsorptionEmission, 0);

addToRunTimeSelectionTable(
    absorptionEmissionModel,
    wsggGreyCloudAbsorptionEmission,
    dictionary);
}
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::wsggGreyCloudAbsorptionEmission::wsggGreyCloudAbsorptionEmission(
    const dictionary& dict,
    const fvMesh& mesh)
: wsggAbsorptionEmission(dict, mesh),
  coeffsDict_(dict.subDict(typeName + "Coeffs")),
  cloudNames_(coeffsDict_.lookup("cloudNames")),
  nBands_(5),
  EDispfields(nBands_),
  eDispfields(nBands_)
{
    // Validate cloudNames contains only carbonaceousCloud
    forAll(cloudNames_, i)
    {
        const word& cloudName = cloudNames_[i];
        CloudHandler& cloudHandler = mesh_.objectRegistry::lookupObjectRef<CloudHandler>("cloudHandler");

        word cloudType(cloudHandler.getCloudType(cloudName));

        if (cloudType != "carbonaceousCloud")
        {
            FatalErrorInFunction
                << "Invalid cloud type: " << cloudType << "for cloud: " << cloudName << nl
                << "Allowed: carbonaceousCloud only"
                << exit(FatalError);
        }
    }

    forAll(EDispfields, i)
    {
        EDispfields.set(
            i,
            new volScalarField(
                IOobject(
                    "EDisp" + name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("EDisp", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)));

        eDispfields.set(
            i,
            new volScalarField(
                IOobject(
                    "eDisp" + name(i),
                    mesh.time().timeName(),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE),
                mesh,
                dimensionedScalar("eDisp", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)));
    }
}
// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::radiation::wsggGreyCloudAbsorptionEmission::~wsggGreyCloudAbsorptionEmission()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggGreyCloudAbsorptionEmission::aDisp(const label) const
{
    tmp<volScalarField> ta(
        new volScalarField(
            IOobject(
                "a",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            mesh_,
            dimensionedScalar(dimless / dimLength, Zero)));

    forAll(cloudNames_, i)
    {
        word cloudName = cloudNames_[i];

        CloudHandler& cloudHandler = mesh_.objectRegistry::lookupObjectRef<CloudHandler>("cloudHandler");

        word cloudType(cloudHandler.getCloudType(cloudName));

        if (cloudType == "carbonaceousCloud")
        {
            const basicCarbonaceousCloud& tc = cloudHandler.getCarbonaceousCloud(cloudName);
            ta.ref() += tc.ap();
        }
    }

    return ta;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggGreyCloudAbsorptionEmission::eDisp(const label bandI) const
{
    tmp<volScalarField> teD(
        new volScalarField(
            IOobject(
                "eDisp",
                mesh_.time().timeName(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE),
            mesh_,
            dimensionedScalar("e", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0)));

    return teD;
}


Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggGreyCloudAbsorptionEmission::EDisp(const label bandI) const
{
    return EDispfields[bandI];
}

Foam::tmp<Foam::volScalarField>
Foam::radiation::wsggGreyCloudAbsorptionEmission::eCont(const label bandI) const
{
    return tmp<volScalarField>::New(
        IOobject(
            "eCont",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh_,
        dimensionedScalar("e", dimensionSet(0, 0, 0, 0, 0, 0, 0), 0.0));
}


template<class T>
void Foam::radiation::wsggGreyCloudAbsorptionEmission::particleSourceHelper(const T& tc) const
{
    // Access gas properties
    const basicThermo* thermoPtr = this->thermoPtr();
    if (!thermoPtr)
    {
        FatalErrorInFunction
            << "WsggAbsorptionEmission model wsggGreyCloudAbsorptionEmission is not thermo aware and has to be initialised first."
            << exit(FatalError);
    }
    const basicThermo& thermo_ = *thermoPtr;
    const volScalarField& W = thermo_.W();
    const volScalarField& CO2 = mesh_.lookupObject<volScalarField>("CO2");
    const volScalarField& H2O = mesh_.lookupObject<volScalarField>("H2O");

    scalar cloudEm_ = tc.constProps().epsilon0();

    // Access gas radiation model for getWeights function
    const radiationModel& radiation_ = mesh_.lookupObject<radiationModel>("radiationProperties");
    // The wsggGreyCloudAbsorptionEmission model can only be used in combination with a wsgg gas absorption emission model
    // Check that it is used together with another model via wsggBinaryAbsorptionEmission
    // In wsggBinaryAbsorptionEmission constructor, it is checked that model 1 is a wsggAbsorptionEmission model
    const wsggBinaryAbsorptionEmission* binaryAbsEmsPtr =
        dynamic_cast<const wsggBinaryAbsorptionEmission*>(&radiation_.absorptionEmission());
    if (!binaryAbsEmsPtr)
    {
        FatalErrorInFunction
            << "The wsggGreyCloudAbsorptionEmission model can only be used in combination with a wsgg gas absorption emission model via"
            << wsggBinaryAbsorptionEmission::typeName << nl
            << "but is of type "
            << radiation_.absorptionEmission().type() << nl
            << "See wsggGreyCloudAbsorptionEmission.H for proper usage." << nl
            << exit(FatalError);
    }
    const wsggBinaryAbsorptionEmission& binaryAbsEms = *binaryAbsEmsPtr;

    const wsggAbsorptionEmission* gasAbsEmsPtr =
        dynamic_cast<const wsggAbsorptionEmission*>(binaryAbsEmsPtr->model1_.get());
    const wsggAbsorptionEmission& gasAbsEms = *gasAbsEmsPtr;

    // Loop over all parcels
    for (const auto& c : tc)
    {
        label cellI = c.cell();
        scalar Tp = c.T();

        scalar XkCO2 = CO2[cellI] * W[cellI] / WCO2_;
        scalar XkH2O = H2O[cellI] * W[cellI] / WH2O_;

        scalar Mr = XkH2O / (XkCO2 + 1e-18);

        scalarList eList = gasAbsEms.getWeights(Tp, Mr);

        for (label bandI = 0; bandI < nBands_; bandI++)
        {
            // add weight*power of individual particles, to the emissivie power, for each band
            EDispfields[bandI][cellI] += 4.0 * eList[bandI] * c.areaP() * pow4(Tp) * c.nParticle() * cloudEm_ * constant::physicoChemical::sigma.value() / (mesh_.V()[cellI]);
        }
    }
}


void Foam::radiation::wsggGreyCloudAbsorptionEmission::correct(
    Foam::volScalarField& aa,
    PtrList<Foam::volScalarField>& aLambda) const
{
    for (label bandI = 0; bandI < nBands_; bandI++)
    {
        EDispfields[bandI] = EDispfields[bandI] * 0; // overwrite all E with 0
    }

    forAll(cloudNames_, i)
    {
        word cloudName = cloudNames_[i];

        CloudHandler& cloudHandler = mesh_.objectRegistry::lookupObjectRef<CloudHandler>("cloudHandler");

        word cloudType(cloudHandler.getCloudType(cloudName));

        if (cloudType == "carbonaceousCloud")
        {
            const basicCarbonaceousCloud& cloud = cloudHandler.getCarbonaceousCloud(cloudName);
            particleSourceHelper<basicCarbonaceousCloud>(cloud);
        }
    }
}

Foam::scalarList Foam::radiation::wsggGreyCloudAbsorptionEmission::getWeights(
    scalar T,
    scalar Yco2h2o) const
{
    FatalErrorInFunction
        << "getWeights is not supported by wsggGreyCloudAbsorptionEmission" << nl
        << "This model uses WSGG weights of the gas absorption emission model."
        << abort(FatalError);

    return scalarList(); // never reached; silences missing-return warnings
}

// ************************************************************************* //
