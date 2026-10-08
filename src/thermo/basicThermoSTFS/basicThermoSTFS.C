/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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

Class
    Foam::basicThermoSTFS

Group
    grpThermo

Description
    Class derived from Thermo that implements an Interface using tabulated chemistry
    with the flameletConfig (hdf5-reader).

Authors
    ...

SourceFiles
    basicThermoSTFS.C

\*---------------------------------------------------------------------------*/

#include "basicThermoSTFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(basicThermoSTFS, 0);
defineRunTimeSelectionTable(basicThermoSTFS, fvMesh);
defineRunTimeSelectionTable(basicThermoSTFS, fvMeshDictPhase);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::basicThermoSTFS::basicThermoSTFS(const fvMesh& mesh, const word& phaseName)
: fluidThermo(mesh, phaseName),
  mesh_(mesh),
  Sct_(1.0),
  Prt_(1.0)
{
    setTurbulenceProperties();
}

Foam::basicThermoSTFS::basicThermoSTFS(const fvMesh& mesh, const dictionary& dict, const word& phaseName)
: fluidThermo(mesh, dict, phaseName),
  mesh_(mesh),
  Sct_(1.0),
  Prt_(1.0)
{
    setTurbulenceProperties();
}

Foam::basicThermoSTFS::basicThermoSTFS(const fvMesh& mesh, const word& phaseName, const word& dictionaryName)
: fluidThermo(mesh, phaseName, dictionaryName),
  mesh_(mesh),
  Sct_(1.0),
  Prt_(1.0)
{
    setTurbulenceProperties();
}

void Foam::basicThermoSTFS::setTurbulenceProperties()
{
    IOdictionary turbulenceProperties(
        IOobject(
            "turbulenceProperties",
            mesh_.time().constant(),
            mesh_,
            IOobject::READ_IF_PRESENT,
            IOobject::NO_WRITE));

    const word simulationType = turbulenceProperties.lookupOrDefault<word>("simulationType", "laminar");
    if (simulationType == "LES")
    {
        Sct_ = turbulenceProperties.subDict("LES").lookupOrDefault<scalar>("Sct", 0.4);
        Prt_ = turbulenceProperties.subDict("LES").lookupOrDefault<scalar>("Prt", Sct_);
        Info << "\nTurbulent Schmidt number Sct = " << Sct_ << endl;
        Info << "Turbulent Prandtl number Prt = " << Prt_ << endl;
    }
    else if (simulationType == "RAS")
    {
        Sct_ = turbulenceProperties.subDict("RAS").lookupOrDefault<scalar>("Sct", 0.7);
        Prt_ = turbulenceProperties.subDict("RAS").lookupOrDefault<scalar>("Prt", Sct_);
        Info << "\nTurbulent Schmidt number Sct = " << Sct_ << endl;
        Info << "Turbulent Prandtl number Prt = " << Prt_ << endl;
    }
    else if (simulationType == "laminar")
    {
        Sct_ = 1.0;
        Prt_ = 1.0;
    }
    else
    {
        FatalErrorIn("Foam::basicThermoSTFS::setTurbulence") << "simulationType in turbulenceProperties is unknown" << endl
                                                             << abort(FatalError);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::basicThermoSTFS> Foam::basicThermoSTFS::New(
    const fvMesh& mesh,
    const word& phaseName)
{
    return basicThermo::New<basicThermoSTFS>(mesh, phaseName);
}

Foam::autoPtr<Foam::basicThermoSTFS> Foam::basicThermoSTFS::New(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName)
{
    return basicThermo::New<basicThermoSTFS>(mesh, phaseName, dictionaryName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::basicThermoSTFS::~basicThermoSTFS()
{
}

//******************** Additional getter functions ********************//

Foam::scalar Foam::basicThermoSTFS::Sct() const
{
    return Sct_;
};

Foam::scalar Foam::basicThermoSTFS::Prt() const
{
    return Prt_;
};
