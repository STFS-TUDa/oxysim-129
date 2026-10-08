/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2016-2020 OpenCFD Ltd.
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

#include "radiationModelSTFS.H"
#include "absorptionEmissionModel.H"
#include "scatterModel.H"
#include "sootModel.H"
#include "fvmSup.H"
#include "basicThermo.H"

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

namespace Foam
{
namespace radiation
{

Foam::word Foam::radiation::radiationModelSTFS::resolveThermoName(
    const dictionary& dict,
    const fvMesh& mesh)
{
    // 1. use entry thermoName in radiationProperties
    word thermoName;
    if (dict.readIfPresent("thermoName", thermoName))
    {
        Info << "radiationModelSTFS: using 'thermoName' entry found in radiationProperties: "
             << thermoName << endl;
        return thermoName;
    }

    // get all basicThermo objects
    const HashTable<const basicThermo*> thermos(mesh.lookupClass<basicThermo>());

    // 2. exactly one basicThermo object exists
    if (thermos.size() == 1)
    {
        thermoName = thermos.begin().key();

        Info << "radiationModelSTFS: no 'thermoName' entry found in "
             << "radiationProperties, using thermo '" << thermoName
             << "' for radiation calculations." << endl;

        return thermoName;
    }

    // 3. no basicThermo object exists
    if (thermos.empty())
    {
        FatalErrorInFunction
            << "radiationModelSTFS constructor requires a basicThermo object."
            << "No basicThermo (or derived) object is registered in the mesh." << nl
            << exit(FatalError);
    }

    // 4. more than one basicThermo object exists
    FatalErrorInFunction
        << "RadiationModelSTFS constructor requires a basicThermo object."
        << "Multiple basicThermo objects are registered in the mesh:" << nl;

    for (auto iter = thermos.cbegin(); iter != thermos.cend(); ++iter)
    {
        FatalError
            << "    " << iter.key() << "  (type: " << iter()->type() << ')' << nl;
    }

    FatalError
        << "Add a 'thermoName' entry to radiationProperties specifying "
        << "which thermo object radiationModelSTFS should use." << nl
        << exit(FatalError);

    return word::null; // unreachable, silences -Wreturn-type
}

}
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::radiation::radiationModelSTFS::radiationModelSTFS(const volScalarField& T)
: radiationModel(T),
  ShGas(
      IOobject(
          "ShGas",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar("ShGas", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)),
  thermoPtr_(&mesh_.lookupObject<basicThermo>(resolveThermoName(*this, T.mesh())))
{
}


Foam::radiation::radiationModelSTFS::radiationModelSTFS(
    const word& type,
    const volScalarField& T)
: radiationModel(type, T),
  ShGas(
      IOobject(
          "ShGas",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar("ShGas", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)),
  thermoPtr_(&mesh_.lookupObject<basicThermo>(resolveThermoName(*this, T.mesh())))
{
}


Foam::radiation::radiationModelSTFS::radiationModelSTFS(
    const word& type,
    const dictionary& dict,
    const volScalarField& T)
: radiationModel(type, dict, T),
  ShGas(
      IOobject(
          "ShGas",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::AUTO_WRITE),
      mesh_,
      dimensionedScalar("ShGas", dimensionSet(1, -1, -3, 0, 0, 0, 0), 0.0)),
  thermoPtr_(&mesh_.lookupObject<basicThermo>(resolveThermoName(*this, T.mesh())))
{
}


// * * * * * * * * * * * * * * * * Destructor    * * * * * * * * * * * * * * //
Foam::radiation::radiationModelSTFS::~radiationModelSTFS()
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::radiation::radiationModelSTFS::correct()
{
    radiationModel::correct();
}

Foam::tmp<Foam::fvScalarMatrix> Foam::radiation::radiationModelSTFS::Sh(
    const basicThermo& thermo,
    const volScalarField& he) const
{
    const volScalarField Cpv(thermo.Cpv());
    const volScalarField T3(pow3(T_));

    tmp<fvScalarMatrix> Sh0 = (Ru()
                               - fvm::Sp(4.0 * Rp() * T3 / Cpv, he)
                               - Rp() * T3 * (T_ - 4.0 * he / Cpv));
    ShGas = (Sh0.ref() & he); // Update tracked source term for gas radiation

    return Sh0;
}


// ************************************************************************* //
