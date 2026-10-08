/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
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

#include <exception>

#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/ct_defs.h"
#include "cantera/base/Solution.h"
#include "cantera/thermo/SingleSpeciesTP.h"
#include "cantera/base/Array.h"

#include "canteraSpecieMixture.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(canteraSpecieMixture, 0);
}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const Foam::canteraSpecieMixture::coeffArray&
Foam::canteraSpecieMixture::coeffs(
    const scalar T,
    const label speciei) const
{
    if (T < nasaTMid_[speciei])
    {
        return nasaCoeffsLow_[speciei];
    }
    else
    {
        return nasaCoeffsHigh_[speciei];
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::canteraSpecieMixture::canteraSpecieMixture(
    const dictionary& thermoDict,
    const wordList& specieNames,
    const fvMesh& mesh,
    const word& phaseName)
: basicSpecieMixture(thermoDict, specieNames, mesh, phaseName),
  Rc_(Cantera::GasConstant)
{
    // Create Cantera objects from mechanism file
    const fileName mechanismFile(thermoDict.get<fileName>("mechanismFile"));
    // const std::string phase(thermoDict.get<word>("phaseName"));
    std::shared_ptr<Cantera::Solution> solution(Cantera::newSolution(mechanismFile, phaseName));
    std::shared_ptr<Cantera::ThermoPhase> gas(solution->thermo());
    std::shared_ptr<Cantera::Transport> transport(solution->transport());
    const Cantera::MultiSpeciesThermo& thermo(gas->speciesThermo());

    const label nSpecies = gas->nSpecies();

    W_.reserve(nSpecies);
    Hc_.reserve(nSpecies);
    nasaCoeffsLow_.resize(nSpecies, coeffArray(0.));
    nasaCoeffsHigh_.resize(nSpecies, coeffArray(0.));
    nasaTMid_.resize(nSpecies);
    muCoeffs_.resize(nSpecies, coeffArray5(0.));
    lambdaCoeffsHigh_.resize(nSpecies, coeffArray5(0.));

    // Get molecular weights and enthalpy of formation
    gas->setState_TP(298.15, 101325.0);
    gas->getMolecularWeights(W_.data());
    gas->getEnthalpy_RT_ref(Hc_.data());

    for (label i = 0; i < nSpecies; ++i)
    {
        Hc_[i] *= gas->RT() / W_[i];
    }

    // std::shared_ptr<Cantera::Phase> gasPhase(solution->thermo());
    // std::string specieName = "H2";
    // std::shared_ptr<Cantera::Species> speciesI(gasPhase->species(specieName));
    // Cantera::SingleSpeciesTP thermoI;
    // thermoI.addSpecies(speciesI);

    // Get the viscosity and conductivity polynomial coefficients
    for (label i = 0; i < nSpecies; ++i)
    {
        std::vector<scalar> muPoly(5);
        std::vector<scalar> lambdaPoly(5);

        transport->getViscosityPolynomial(i, muPoly.data());
        transport->getConductivityPolynomial(i, lambdaPoly.data());

        forAll(muCoeffs_[i], j)
        {
            muCoeffs_[i][j] = muPoly[j];
            lambdaCoeffsHigh_[i][j] = lambdaPoly[j];
        }
    }

    // Get the NASA polynomial coefficients
    std::vector<scalar> coeffs(15);
    int type;
    scalar Thigh;
    scalar Tlow;
    scalar pRef;

    for (label i = 0; i < nSpecies; ++i)
    {
        thermo.reportParams(i, type, coeffs.data(), Tlow, Thigh, pRef);
        for (label j = 0; j < 7; ++j)
        {
            nasaCoeffsHigh_[i][j] = coeffs[j + 1] * Rc_;
            nasaCoeffsLow_[i][j] = coeffs[j + 8] * Rc_;
        }
        nasaTMid_[i] = coeffs[0];
    }
}

Foam::scalar Foam::canteraSpecieMixture::W(const label speciei) const
{
    return W_[speciei];
}

Foam::scalar Foam::canteraSpecieMixture::Hc(const label speciei) const
{
    return Hc_[speciei];
}

Foam::scalar Foam::canteraSpecieMixture::Cp(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    const coeffArray& A = coeffs(T, speciei);
    return (A[0] + T * (A[1] + T * (A[2] + T * (A[3] + T * A[4])))) / W_[speciei];
}

Foam::scalar Foam::canteraSpecieMixture::Cv(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    return Cp(speciei, p, T) - Rc_ / W_[speciei];
}

Foam::scalar Foam::canteraSpecieMixture::HE(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    throw std::runtime_error("canteraSpecieMixture::HE not implemented");
}

Foam::scalar Foam::canteraSpecieMixture::Ha(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    const coeffArray& A = coeffs(T, speciei);
    return (T * (A[0] + T * (A[1] / 2. + T * (A[2] / 3. + T * (A[3] / 4. + T * A[4] / 5.)))) + A[5])
         / W_[speciei];
}

Foam::scalar Foam::canteraSpecieMixture::Hs(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    return Ha(speciei, p, T) - Hc(speciei);
}

Foam::scalar Foam::canteraSpecieMixture::S(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    const coeffArray& A = coeffs(T, speciei);
    return (A[0] * log(T) + T * (A[1] + T * (A[2] / 2. + T * (A[3] / 3. + T * A[4] / 4.))) + A[6])
         / W_[speciei];
}

Foam::scalar Foam::canteraSpecieMixture::Es(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    return (Hs(speciei, p, T) - p / rho(speciei, p, T));
}

Foam::scalar Foam::canteraSpecieMixture::G(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    return (Ha(speciei, p, T) - T * S(speciei, p, T));
}

Foam::scalar Foam::canteraSpecieMixture::A(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    throw std::runtime_error("canteraSpecieMixture::A not implemented");
}

Foam::scalar Foam::canteraSpecieMixture::mu(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    const coeffArray5& A = muCoeffs_[speciei];
    scalar lnT = log(T);
    scalar sumAlnT = A[0] + lnT * (A[1] + lnT * (A[2] + lnT * (A[3] + lnT * (A[4]))));

    return sqrt(T) * pow(sumAlnT, 2.);
}

Foam::scalar Foam::canteraSpecieMixture::kappa(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    const coeffArray5& A = lambdaCoeffsHigh_[speciei];
    scalar lnT = log(T);
    scalar sumAlnT = A[0] + lnT * (A[1] + lnT * (A[2] + lnT * (A[3] + lnT * (A[4]))));
    return sqrt(T) * sumAlnT;
}

Foam::scalar Foam::canteraSpecieMixture::alphah(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    return kappa(speciei, p, T) / Cp(speciei, p, T);
}

Foam::scalar Foam::canteraSpecieMixture::rho(
    const label speciei,
    const scalar p,
    const scalar T) const
{
    // Density per species? Not used imo.
    throw std::runtime_error("canteraSpecieMixture::rho not implemented");
}


// ************************************************************************* //
