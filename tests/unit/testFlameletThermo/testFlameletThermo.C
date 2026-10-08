/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of an extension to OpenFOAM.

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

Description
    Unit test for the flameletThermo class.

Authors
    Matthias Steinhausen

\*---------------------------------------------------------------------------*/

// Including Catch2
#include <catch2/catch_all.hpp>

// own headers
#include "fvCFD.H"
#include "flameletThermo.H"

// Create a dummy Time object
Foam::Time runTime(
    Foam::Time::controlDictName,
    "",
    "",
    "system",
    "constant");

// Read a single cell mesh from constant/polyMesh
Foam::fvMesh mesh(
    Foam::IOobject(
        fvMesh::defaultRegion,
        runTime.timeName(),
        runTime,
        Foam::IOobject::MUST_READ),
    false);


TEST_CASE("flameletThermo", "[flameletThermo]")
{
    // Suppress warnings and throw exceptions
    // Necessary for REQUIRE_THROWS to work with OpenFOAM errors.
    Foam::Warning.level = 0;
    Foam::FatalError.throwExceptions();
    Foam::FatalIOError.throwExceptions();

    std::cout << "Test flameletThermo" << std::endl
              << std::endl;

    std::system("rm -r 0");
    std::system("ln -s 0.normal 0");
    std::system("cp constant/flameletProperties.org constant/flameletProperties");

    //
    // Test constructor
    //
    Foam::flameletThermo thermo(mesh);

    // Test input-fields are created
    REQUIRE_NOTHROW(thermo.getField("Z"));
    REQUIRE_NOTHROW(thermo.getField("yc"));
    REQUIRE_NOTHROW(thermo.getField("ha"));

    // Test outputFields are created
    // - default fields
    REQUIRE_NOTHROW(thermo.T());
    REQUIRE_NOTHROW(thermo.alpha());
    REQUIRE_NOTHROW(thermo.rho());
    REQUIRE_NOTHROW(thermo.psi());
    // - user-defined fields
    REQUIRE_NOTHROW(thermo.getField("CH4"));
    REQUIRE_NOTHROW(thermo.getField("H2O"));

    //
    // Test thermo update
    //
    REQUIRE_NOTHROW(thermo.update());
    // Check consistency of lookup quantities
    REQUIRE_THAT(thermo.getField("yc").internalField()[0], Catch::Matchers::WithinRel(0.05, 1e-7));
    REQUIRE_THAT(thermo.getField("ha").internalField()[0], Catch::Matchers::WithinRel(-254452.7113, 1e-7));
    // Check values of dependent fields
    REQUIRE_THAT(thermo.rhoRef().internalField()[0], Catch::Matchers::WithinRel(0.235897, 1e-5));
    REQUIRE_THAT(thermo.rhoRef().boundaryField()[0][0], Catch::Matchers::WithinRel(0.235897, 1e-5));
    REQUIRE_THAT(thermo.alpha().internalField()[0], Catch::Matchers::WithinRel(7.45819e-5, 1e-5));
    REQUIRE_THAT(thermo.alpha().boundaryField()[0][0], Catch::Matchers::WithinRel(7.45819e-5, 1e-5));
    REQUIRE_THAT(thermo.psi().internalField()[0], Catch::Matchers::WithinRel(2.32812e-6, 1e-5));
    REQUIRE_THAT(thermo.psi().boundaryField()[0][0], Catch::Matchers::WithinRel(2.32812e-6, 1e-5));
    REQUIRE_THAT(thermo.T().internalField()[0], Catch::Matchers::WithinRel(1394.43, 1e-5));
    REQUIRE_THAT(thermo.T().boundaryField()[0][0], Catch::Matchers::WithinRel(1394.43, 1e-5));

    //
    // Test adding and retrival of fields
    //
    // - additional fields (no lookup quantities)
    REQUIRE_THROWS(thermo.getField("dummyField"));
    REQUIRE_NOTHROW(thermo.appendField("dummyField", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimensionSet(0, -2, 2, 0, 0), 20.0)));
    REQUIRE_NOTHROW(thermo.getField("dummyField"));                                      // Field is found
    REQUIRE(mesh.foundObject<volScalarField>("dummyField") == true);                     // Field is registered to objectRegistry
    REQUIRE(thermo.getField("dummyField").dimensions() == dimensionSet(0, -2, 2, 0, 0)); // Dimensions are correct
    REQUIRE(thermo.getField("dummyField").internalFieldRef()[0] == 20.0);                // Value is correct

    thermo.update();

    // Field is not part of thermo update
    REQUIRE(thermo.getField("dummyField").internalFieldRef()[0] == 20.0);

    // - additional output-fields
    // Check that inputVariables cannot be added as outputFields
    REQUIRE_THROWS(thermo.appendOutputField("yc"));
    // Check standard output-field
    REQUIRE_THROWS(thermo.getField("CO2"));
    REQUIRE_NOTHROW(thermo.appendOutputField("CO2", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 30.0)));
    REQUIRE_NOTHROW(thermo.getField("CO2"));                                     // Field is found
    REQUIRE(mesh.foundObject<volScalarField>("CO2") == true);                    // Field is registered to objectRegistry
    REQUIRE(thermo.getField("CO2").dimensions() == dimensionSet(0, 0, 1, 0, 0)); // Dimensions are correct
    REQUIRE(thermo.getField("CO2").internalFieldRef()[0] == 30.0);               // Value is correct

    // Reload the outputvariable data
    REQUIRE_NOTHROW(thermo.loadOutputVariablesFromFLUT());
    thermo.update();

    // Field is updated by thermo.update()
    // - internalField
    REQUIRE_THAT(thermo.getField("CO2").internalField()[0], Catch::Matchers::WithinRel(thermo.getField("yc").internalField()[0], 1e-5));
    REQUIRE_THAT(thermo.getField("CO2").boundaryField()[0][0], Catch::Matchers::WithinRel(thermo.getField("yc").boundaryField()[0][0], 1e-5));

    // Append outputField that is not in table
    CHECK_THROWS(thermo.appendOutputField("dummyField", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimensionSet(0, 0, 1, 0, 0), 30.0)));

    //
    // Test density correction
    //
    volScalarField deltaRho(
        IOobject(
            "deltaRho",
            runTime.timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE),
        mesh,
        dimensionedScalar("zero", dimDensity, 1e-2));

    // Normal correction
    REQUIRE_NOTHROW(thermo.correctRho(deltaRho));
    REQUIRE_THAT(thermo.rhoRef().internalField()[0], Catch::Matchers::WithinRel(0.235897 + 1e-2, 1e-5));
    // Correction with bounds
    REQUIRE_NOTHROW(thermo.correctRho(deltaRho, dimensionedScalar("min", dimDensity, 0), dimensionedScalar("min", dimDensity, 1)));
    REQUIRE_THAT(thermo.rhoRef().internalField()[0], Catch::Matchers::WithinRel(0.235897 + 2e-2, 1e-5));
    // Correction with lower bound
    REQUIRE_NOTHROW(thermo.correctRho(deltaRho, dimensionedScalar("min", dimDensity, 0), dimensionedScalar("min", dimDensity, 0.2)));
    REQUIRE_THAT(thermo.rhoRef().internalField()[0], Catch::Matchers::WithinRel(0.2, 1e-5));
    // Correction with upper bound
    REQUIRE_NOTHROW(thermo.correctRho(deltaRho, dimensionedScalar("min", dimDensity, 0.3), dimensionedScalar("min", dimDensity, 0.4)));
    REQUIRE_THAT(thermo.rhoRef().internalField()[0], Catch::Matchers::WithinRel(0.3, 1e-5));
}

TEST_CASE("flameletThermo with boundary lookup", "[flameletThermo]")
{
    // Suppress warnings and throw exceptions
    // Necessary for REQUIRE_THROWS to work with OpenFOAM errors.
    Foam::Warning.level = 0;
    Foam::FatalError.throwExceptions();
    Foam::FatalIOError.throwExceptions();

    std::cout << std::endl
              << std::endl
              << "Test flameletThermo with Boundary lookup" << std::endl
              << std::endl;
    // Copy the correct 0 folder for the boundary lookup
    std::system("rm -r 0");
    std::system("ln -s 0.BC 0");
    std::system("cp constant/flameletProperties.org constant/flameletProperties");

    //
    // Test constructor
    //
    Foam::flameletThermo thermo(mesh);

    // Test input-fields are created
    REQUIRE_NOTHROW(thermo.getField("Z"));
    REQUIRE_NOTHROW(thermo.getField("yc"));
    REQUIRE_NOTHROW(thermo.getField("ha"));
    REQUIRE_NOTHROW(thermo.getField("TemperatureBC"));

    // Test outputFields are created
    // - default fields
    REQUIRE_NOTHROW(thermo.T());
    REQUIRE_NOTHROW(thermo.alpha());
    REQUIRE_NOTHROW(thermo.rho());
    REQUIRE_NOTHROW(thermo.psi());
    REQUIRE_NOTHROW(thermo.getField("yc_max"));

    //
    // Test boundary lookup
    //
    // Test thermo update
    REQUIRE_NOTHROW(thermo.update());
    // Enthalpy field is correctly updated
    // - fixedTemperature boundary
    REQUIRE_THAT(thermo.getField("ha").boundaryField()[0][0], Catch::Matchers::WithinRel(-1686266.937, 1e-5));
    REQUIRE_THAT(thermo.T().boundaryField()[0][0], Catch::Matchers::WithinRel(300, 1e-5));
    // - normal / other boundary
    REQUIRE_THAT(thermo.getField("ha").boundaryField()[1][0], Catch::Matchers::WithinRel(-254452.7113, 1e-5));
}

TEST_CASE("flameletThermo with canteraThermo.")
{
    // Suppress warnings and throw exceptions
    // Necessary for REQUIRE_THROWS to work with OpenFOAM errors.
    Foam::Warning.level = 0;
    Foam::FatalError.throwExceptions();
    Foam::FatalIOError.throwExceptions();

    std::cout << std::endl
              << std::endl
              << "Test flameletThermo with canteraThermo." << std::endl
              << std::endl;

    std::system("rm -r 0");
    std::system("ln -s 0.normal 0");
    std::system("cp constant/flameletProperties.org constant/flameletProperties");
    std::system("sed -i 's/useWithCanteraThermo false;/useWithCanteraThermo true;/' constant/flameletProperties");

    //
    // Test constructor
    //
    Foam::flameletThermo thermo(mesh);
    REQUIRE_NOTHROW(thermo.update());
    // Check that default fields are not updated.
    REQUIRE_THAT(thermo.T().internalField()[0], Catch::Matchers::WithinRel(300, 1e-5));
    REQUIRE_THAT(thermo.alpha().internalField()[0], Catch::Matchers::WithinRel(0, 1e-5));
    REQUIRE_THAT(thermo.rhoRef().internalField()[0], Catch::Matchers::WithinRel(0, 1e-5));
}
