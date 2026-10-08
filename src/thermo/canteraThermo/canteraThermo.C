/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
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

#include "canteraThermo.H"

#include "gradientEnergyFvPatchScalarField.H"
#include "mixedEnergyFvPatchScalarField.H"
#include "processorFvPatchField.H"

#include "cantera/core.h"
#include "cantera/thermo/ThermoFactory.h"
#include "cantera/transport/TransportFactory.h"
#include "cantera/base/ct_defs.h"
#include "cantera/base/Solution.h"
#include "cantera/base/Array.h"
#include "cantera/base/utilities.h"

#include "cantera/transport/TransportData.h"
#include "cantera/thermo/Species.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(canteraThermo, 0);
defineRunTimeSelectionTable(canteraThermo, fvMesh);
defineRunTimeSelectionTable(canteraThermo, fvMeshDictPhase);
}

Foam::wordList Foam::getSpecieNames(const word& mechanismFile, const word& phaseName)
{
    std::shared_ptr<Cantera::Solution> solution(Cantera::newSolution(mechanismFile, phaseName));

    const auto ctNames(solution->thermo()->speciesNames());
    wordList specieNames(solution->thermo()->nSpecies());
    for (std::size_t i = 0; i < solution->thermo()->nSpecies(); ++i)
    {
        specieNames[i] = ctNames[i];
    }
    return specieNames;
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //
void Foam::canteraThermo::heBoundaryCorrection(volScalarField& he)
{
    volScalarField::Boundary& heBf = he.boundaryFieldRef(false);

    forAll(heBf, patchi)
    {
        if (isA<gradientEnergyFvPatchScalarField>(heBf[patchi]))
        {
            refCast<gradientEnergyFvPatchScalarField>(heBf[patchi]).gradient() = heBf[patchi].fvPatchField::snGrad();
        }
        else if (isA<mixedEnergyFvPatchScalarField>(heBf[patchi]))
        {
            refCast<mixedEnergyFvPatchScalarField>(heBf[patchi]).refGrad() = heBf[patchi].fvPatchField::snGrad();
        }
    }
}

inline std::string Foam::canteraThermo::handleTransportModel(std::string transportModel)
{
    // handling deprecated key-words of transport models for Cantera during initialization
    if (transportModel == "UnityLewis")
    {
        Info << "Replacing deprecated 'UnityLewis' with 'unity-Lewis-number' for Cantera." << endl;
        transportModel = "unity-Lewis-number";
    }
    else if (transportModel == "Mix")
    {
        Info << "Replacing deprecated 'Mix' with 'mixture-averaged' for Cantera." << endl;
        transportModel = "mixture-averaged";
    }
    else if (transportModel == "Multi")
    {
        Info << "Replacing deprecated 'Multi' with 'multicomponent' for Cantera." << endl;
        transportModel = "multicomponent";
    }
    else if (transportModel == "mixture-averaged" || transportModel == "unity-Lewis-number" || transportModel == "multicomponent")
    {
        // do nothing
        ;
    }
    else
    {
        FatalErrorIn("canteraMixture::canteraMixture()")
            << "Transport model " << transportModel
            << "is not supported" << endl
            << "Choose one of: UnityLewis, Mix, Multi, mixture-averaged, unity-Lewis-number, multicomponent! Aborting!"
            << exit(FatalError);
    }
    return transportModel;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //
Foam::canteraThermo::canteraThermo(const fvMesh& mesh, const word& phaseName, const word& dictionaryName, bool initialUpdate)
: basicThermoSTFS(mesh, phaseName, dictionaryName),
  canteraSpecieMixture(mesh.lookupObject<IOdictionary>(dictionaryName), getSpecieNames(get<fileName>("mechanismFile"), get<word>("phaseName")), mesh, phaseName),
  transportModel_(handleTransportModel(get<word>("transportModel"))), // calling 'handleTransportModel' here
  useSoret_(get<bool>("useSoret")),
  sol_(Cantera::newSolution(get<fileName>("mechanismFile"), get<word>("phaseName"), transportModel_)),
  gas_(sol_->thermo()),
  trans_(sol_->transport()),
  nSpecies_(gas_->nSpecies()),
  ctY_(nSpecies_),
  ctD_(transportModel_ == "multicomponent" ? nSpecies_ * nSpecies_ : nSpecies_),
  ctDtherm_(nSpecies_),
  energyType_(get<word>("energyType")),
  he_(
      IOobject(
          energyType_,
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::AUTO_WRITE),
      mesh,
      dimensionedScalar("zero", dimEnergy / dimMass, 0.0),
      this->heBoundaryTypes(),
      this->heBoundaryBaseTypes()),
  hi_(nSpecies_),
  heSensi_(nSpecies_),
  Wi_(nSpecies_),
  heSensRefi_(nSpecies_),
  rho_(
      IOobject(
          (energyType_ == "ea" ? "rho" : "thermo:rho"),
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          energyType_ == "ea" ? IOobject::AUTO_WRITE : IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar("zero", dimMass / dimVolume, 0.0)),
  psi_(
      IOobject(
          "thermo:psi",
          mesh_.time().timeName(),
          mesh_,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar("zero", dimTime * dimTime / dimArea, 0.0)),
  cp_(
      IOobject(
          "cp",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar("zero", dimEnergy / dimMass / dimTemperature, 0.0)),
  cv_(
      IOobject(
          "cv",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar("zero", dimEnergy / dimMass / dimTemperature, 0.0)),
  W_(
      IOobject(
          "W",
          mesh_.time().timeName(),
          mesh_,
          IOobject::NO_READ,
          IOobject::NO_WRITE),
      mesh_,
      dimensionedScalar("zero", dimMass / dimMoles, 0.0)),
  // Transport properties
  D_(transportModel_ == "multicomponent" ? nSpecies_ * nSpecies_ : nSpecies_),
  Dtherm_(nSpecies_),
  nu_(
      IOobject(
          "nu",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("zero", dimArea / dimTime, 0.0)),
  mu_(
      IOobject(
          "mu",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("zero", dimMass / dimLength / dimTime, 0.0)),
  lambda_(
      IOobject(
          "lambda",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("zero", dimMass * dimLength / (dimTime * dimTime * dimTime * dimTemperature), 0.0)),
  speedOfSound_(
      IOobject(
          "speedOfSound",
          mesh.time().timeName(),
          mesh,
          IOobject::READ_IF_PRESENT,
          IOobject::NO_WRITE),
      mesh,
      dimensionedScalar("zero", dimLength / dimTime, 0.0))
{
    Info << "Initializing canteraMixture" << endl;

    // Validate the energyType input
    if (energyType_ == "ha")
    {
        Info << "Use absolute Enthalpy 'ha' as energy type." << endl;
    }
    else if (energyType_ == "ea")
    {
        Info << "Use absolute Energy 'ea' as energy type." << endl;
    }
    else
    {
        FatalErrorIn("canteraMixture::canteraMixture()")
            << "Invalid energy type: " << energyType_
            << " valid options are absolute Enthalpy ('ha') and absolute Energy ('ea')."
            << exit(FatalError);
    }

    Info << "Initializing canteraThermo with "
         << transportModel_ << " transport model\n"
         << endl;

    forAll(hi_, i)
    {
        hi_.set(
            i, new Foam::volScalarField(IOobject("h_" + gas_->speciesName(i), mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedScalar("zero", dimEnergy / dimMass, 0.0)));
    }
    forAll(Wi_, i)
    {
        Wi_.set(
            i, new dimensionedScalar("Wi_" + gas_->speciesName(i), dimMass / dimMoles, gas_->molecularWeight(i)));
    }
    forAll(heSensi_, i)
    {
        heSensi_.set(
            i, new Foam::volScalarField(IOobject("heSens_" + gas_->speciesName(i), mesh_.time().timeName(), mesh_, IOobject::NO_READ, IOobject::NO_WRITE), mesh_, dimensionedScalar("zero", dimEnergy / dimMass, 0.0)));
    }

    // Initialize reference enthalpy value
    const scalar Tstd(constant::standard::Tstd.value());
    gas_->setTemperature(Tstd);

    std::vector<double> tempHSensi(nSpecies_);
    gas_->getEnthalpy_RT_ref(tempHSensi.data());

    forAll(heSensRefi_, i)
    {
        tempHSensi[i] *= gas_->RT() / Wi_[i].value();

        heSensRefi_.set(
            i, new dimensionedScalar("heSensRef_" + gas_->speciesName(i), dimEnergy / dimMass, tempHSensi[i]));
    }

    // Initialze fields
    if (transportModel_ == "multicomponent")
    {
        for (int i = 0; i != nSpecies_; ++i)
        {
            for (int j = 0; j != nSpecies_; ++j)
            {
                D_.set(
                    i + j * nSpecies_, new volScalarField(IOobject("D_" + gas_->speciesName(i) + "_" + gas_->speciesName(j), mesh_.time().timeName(), mesh_, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), mesh_, dimensionedScalar("zero", dimArea / dimTime, 0.0)));
            }
        }
    }
    else if (transportModel_ == "mixture-averaged")
    {
        for (int i = 0; i != nSpecies_; ++i)
        {
            D_.set(
                i, new volScalarField(IOobject("D_" + gas_->speciesName(i), mesh_.time().timeName(), mesh_, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), mesh_, dimensionedScalar("zero", dimArea / dimTime, 0.0)));
        }
    }

    if (useSoret_)
    {
        for (int i = 0; i != nSpecies_; ++i)
        {
            Dtherm_.set(
                i, new volScalarField(IOobject("Dtherm_" + gas_->speciesName(i), mesh_.time().timeName(), mesh_, IOobject::READ_IF_PRESENT, IOobject::NO_WRITE), mesh_, dimensionedScalar("zero", dimMass / dimLength / dimTime, 0.0)));
        }
    }
    // Mixture update
    // This is switch is needed for lagrangian stuff. There we need the cantera
    // thermo for particle coupling, but we don't want it to update properties
    // itself.
    if (initialUpdate)
    {
        updateTPY(
            p_,
            T_,
            psi_,
            cp_,
            cv_,
            W_,
            rho_,
            he_,
            true);
        this->heBoundaryCorrection(he_);
    }
}

// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::canteraThermo> Foam::canteraThermo::New(
    const fvMesh& mesh,
    const word& phaseName)
{
    return basicThermo::New<canteraThermo>(mesh, phaseName);
}

Foam::autoPtr<Foam::canteraThermo> Foam::canteraThermo::New(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName)
{
    return basicThermo::New<canteraThermo>(mesh, phaseName, dictionaryName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::canteraThermo::~canteraThermo()
{
}

// * * * * * * * * * * * * * Public member functions * * * * * * * * * * * //
void Foam::canteraThermo::updateTPY()
{
    updateTPY(
        p_,
        T_,
        psi_,
        cp_,
        cv_,
        W_,
        rho_,
        he_,
        false);
}

void Foam::canteraThermo::updateTPY(
    const volScalarField& p,
    volScalarField& T,
    volScalarField& psi,
    volScalarField& cp,
    volScalarField& cv,
    volScalarField& W,
    volScalarField& rho,
    volScalarField& he,
    const bool doOldTimes)
{
    if (doOldTimes && (p.nOldTimes() || T.nOldTimes()))
    {
        updateTPY(
            p.oldTime(),
            T.oldTime(),
            psi.oldTime(),
            cp.oldTime(),
            cv.oldTime(),
            W.oldTime(),
            rho.oldTime(),
            he.oldTime(),
            true);
    }
    // Get reference to internal fields (-> call to primitiveFieldRef function
    // automatically saves old timestep for time derivative; not needed for p
    // and T as they're solved quantities)
    const scalarField& TCells = T.primitiveField();
    const scalarField& pCells = p.primitiveField();

    scalarField& psiCells = psi.primitiveFieldRef();
    scalarField& cpCells = cp.primitiveFieldRef();
    scalarField& cvCells = cv.primitiveFieldRef();
    scalarField& WCells = W.primitiveFieldRef();
    scalarField& rhoCells = rho.primitiveFieldRef();
    scalarField& heCells = he.primitiveFieldRef();


    forAll(mesh_.C(), cell)
    {
        const scalar ctT = TCells[cell];
        const scalar ctP = pCells[cell];

        forAll(Y_, i)
        {
            ctY_[i] = Y_[i].internalField()[cell];
        }

        // Set state in cantera
        gas_->setState_TPY(ctT, ctP, ctY_.data());

        updateMixtureProperties(cell,
                                psiCells,
                                cpCells,
                                cvCells,
                                WCells,
                                rhoCells);
        updateTransportProperties(ctT, rhoCells[cell], WCells[cell], cell);
        heCells[cell] = getHeFromGasObject();
    }

    // Update boundary fields
    // Get reference to boundary fields (-> call to function automatically saves
    // old timestep for time derivative)
    const volScalarField::Boundary& TBf = T_.boundaryField();
    const volScalarField::Boundary& pBf = p_.boundaryField();

    volScalarField::Boundary& psiBf = psi_.boundaryFieldRef();
    volScalarField::Boundary& cpBf = cp_.boundaryFieldRef();
    volScalarField::Boundary& cvBf = cv_.boundaryFieldRef();
    volScalarField::Boundary& WBf = W_.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();

    forAll(mesh_.boundary(), patch)
    {
        forAll(mesh_.boundary()[patch], face)
        {
            // Set state in cantera
            const scalar ctT = TBf[patch][face];
            const scalar ctP = pBf[patch][face];

            forAll(Y_, i)
            {
                ctY_[i] = Y_[i].boundaryField()[patch][face];
            }
            gas_->setState_TPY(ctT, ctP, ctY_.data());

            updateMixtureProperties(patch,
                                    face,
                                    psiBf,
                                    cpBf,
                                    cvBf,
                                    WBf,
                                    rhoBf);
            updateTransportProperties(ctT, rhoBf[patch][face], WBf[patch][face], patch, face);
            he_.boundaryFieldRef(false)[patch][face] = getHeFromGasObject();
        }
        he_.boundaryFieldRef(false)[patch].useImplicit(T_.boundaryField()[patch].useImplicit());
    }
}

void Foam::canteraThermo::updateHPY()
{
    if (energyType_ == "ea")
    {
        Info << "updateRPY so far incompatible with energy type: absolute energy 'ea'." << endl;
        NotImplemented;
    }

    // Get reference to internal fields (-> call to primitiveFieldRef function
    // automatically saves old timestep for time derivative; not needed for p
    // and he as they're solved quantities)

    const scalarField& heCells = he_.primitiveField();
    const scalarField& pCells = p_.primitiveField();

    scalarField& psiCells = psi_.primitiveFieldRef();
    scalarField& cpCells = cp_.primitiveFieldRef();
    scalarField& cvCells = cv_.primitiveFieldRef();
    scalarField& WCells = W_.primitiveFieldRef();
    scalarField& rhoCells = rho_.primitiveFieldRef();
    scalarField& TCells = T_.primitiveFieldRef();

    forAll(mesh_.C(), cell)
    {
        const scalar ctH = heCells[cell];
        const scalar ctP = pCells[cell];

        forAll(Y_, i)
        {
            ctY_[i] = Y_[i].internalField()[cell];
        }

        // Set state in cantera
        // gas_->setState_PY(ctP, ctY_.data());
        gas_->setMassFractions_NoNorm(ctY_.data());
        gas_->setState_HP(ctH, ctP);

        double Temperature = gas_->RT() / Cantera::GasConstant;
        TCells[cell] = Temperature;

        updateMixtureProperties(cell,
                                psiCells,
                                cpCells,
                                cvCells,
                                WCells,
                                rhoCells);
        updateTransportProperties(Temperature, rhoCells[cell], WCells[cell], cell);
    }

    // Get reference to boundary fields (-> call to function automatically saves
    // old timestep for time derivative)
    volScalarField::Boundary& heBf = he_.boundaryFieldRef();
    const volScalarField::Boundary& pBf = p_.boundaryField();

    volScalarField::Boundary& psiBf = psi_.boundaryFieldRef();
    volScalarField::Boundary& cpBf = cp_.boundaryFieldRef();
    volScalarField::Boundary& cvBf = cv_.boundaryFieldRef();
    volScalarField::Boundary& WBf = W_.boundaryFieldRef();
    volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();

    // Update boundary fields
    forAll(mesh_.boundary(), patch)
    {
        if (TBf[patch].fixesValue())
        {
            forAll(mesh_.boundary()[patch], face)
            {
                // Set state in cantera
                const scalar ctP = pBf[patch][face];
                forAll(Y_, i)
                {
                    ctY_[i] = Y_[i].boundaryField()[patch][face];
                }
                const scalar ctT = TBf[patch][face];
                gas_->setState_TPY(ctT, ctP, ctY_.data());
                heBf[patch][face] = getHeFromGasObject();

                updateMixtureProperties(patch,
                                        face,
                                        psiBf,
                                        cpBf,
                                        cvBf,
                                        WBf,
                                        rhoBf);
                updateTransportProperties(TBf[patch][face], rhoBf[patch][face], WBf[patch][face], patch, face);
            }
        }
        else
        {
            forAll(mesh_.boundary()[patch], face)
            {
                // Set state in cantera
                const scalar ctP = pBf[patch][face];
                forAll(Y_, i)
                {
                    ctY_[i] = Y_[i].boundaryField()[patch][face];
                }
                const scalar ctH = heBf[patch][face];
                gas_->setMassFractions_NoNorm(ctY_.data());
                gas_->setState_HP(ctH, ctP);
                TBf[patch][face] = gas_->RT() / Cantera::GasConstant;

                updateMixtureProperties(patch,
                                        face,
                                        psiBf,
                                        cpBf,
                                        cvBf,
                                        WBf,
                                        rhoBf);
                updateTransportProperties(TBf[patch][face], rhoBf[patch][face], WBf[patch][face], patch, face);
            }
        }
    }
}

void Foam::canteraThermo::updateRPY()
{
    // TODO: Is this necessary here. In general the update should be also valid for absolute enthalpy.
    if (energyType_ == "ha")
    {
        Info << "updateRPY so far incompatible with energy type: absolute enthalpy 'ha'." << endl;
        NotImplemented;
    }

    scalarField& rhoCells = rho_.primitiveFieldRef(); // rho needs to be writable
    const scalarField& pCells = p_.primitiveField();

    scalarField& psiCells = psi_.primitiveFieldRef();
    scalarField& cpCells = cp_.primitiveFieldRef();
    scalarField& cvCells = cv_.primitiveFieldRef();
    scalarField& WCells = W_.primitiveFieldRef();
    scalarField& TCells = T_.primitiveFieldRef();

    forAll(mesh_.C(), cell)
    {
        // Set state in cantera
        const scalar ctRho = rhoCells[cell];
        const scalar ctP = pCells[cell];

        forAll(Y_, i)
        {
            ctY_[i] = Y_[i].internalField()[cell];
        }
        gas_->setMassFractions_NoNorm(ctY_.data());
        gas_->setState_DP(ctRho, ctP);

        double Temperature = gas_->RT() / Cantera::GasConstant;
        TCells[cell] = Temperature;

        updateMixtureProperties(cell,
                                psiCells,
                                cpCells,
                                cvCells,
                                WCells,
                                rhoCells);
        updateTransportProperties(Temperature, rhoCells[cell], WCells[cell], cell);
        // update energy
        he_[cell] = getHeFromGasObject();
    }

    // Update boundary fields
    volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef(); // has to be writable
    volScalarField::Boundary& heBf = he_.boundaryFieldRef();   // has to be writable
    volScalarField::Boundary& pBf = p_.boundaryFieldRef();

    volScalarField::Boundary& psiBf = psi_.boundaryFieldRef();
    volScalarField::Boundary& cpBf = cp_.boundaryFieldRef();
    volScalarField::Boundary& cvBf = cv_.boundaryFieldRef();
    volScalarField::Boundary& WBf = W_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();

    forAll(mesh_.boundary(), patch)
    {
        if (TBf[patch].fixesValue())
        {
            forAll(mesh_.boundary()[patch], face)
            {
                // Set state in cantera
                const scalar ctP = pBf[patch][face];
                const scalar ctT = TBf[patch][face];

                forAll(Y_, i)
                {
                    ctY_[i] = Y_[i].boundaryField()[patch][face];
                }

                gas_->setState_TPY(ctT, ctP, ctY_.data());
                rhoBf[patch][face] = gas_->density();
                heBf[patch][face] = getHeFromGasObject();

                updateMixtureProperties(patch,
                                        face,
                                        psiBf,
                                        cpBf,
                                        cvBf,
                                        WBf,
                                        rhoBf);
                updateTransportProperties(TBf[patch][face], rhoBf[patch][face], WBf[patch][face], patch, face);
            }
        }
        else
        {
            forAll(mesh_.boundary()[patch], face)
            {
                // Set state in cantera
                const scalar ctP = pBf[patch][face];
                const scalar ctRho = rhoBf[patch][face];

                forAll(Y_, i)
                {
                    ctY_[i] = Y_[i].boundaryField()[patch][face];
                }
                gas_->setMassFractions_NoNorm(ctY_.data());
                gas_->setState_DP(ctRho, ctP);
                TBf[patch][face] = gas_->RT() / Cantera::GasConstant;
                heBf[patch][face] = getHeFromGasObject();

                updateMixtureProperties(patch,
                                        face,
                                        psiBf,
                                        cpBf,
                                        cvBf,
                                        WBf,
                                        rhoBf);
                updateTransportProperties(TBf[patch][face], rhoBf[patch][face], WBf[patch][face], patch, face);
            }
        }
    }
}

void Foam::canteraThermo::updateERY()
{
    if (energyType_ == "ha")
    {
        Info << "updateRPY so far incompatible with energy type: absolute enthalpy 'ha'." << endl;
        NotImplemented;
    }

    scalarField& rhoCells = rho_.primitiveFieldRef();
    scalarField& pCells = p_.primitiveFieldRef();
    const scalarField& heCells = he_.primitiveField();

    scalarField& psiCells = psi_.primitiveFieldRef();
    scalarField& cpCells = cp_.primitiveFieldRef();
    scalarField& cvCells = cv_.primitiveFieldRef();
    scalarField& WCells = W_.primitiveFieldRef();
    scalarField& TCells = T_.primitiveFieldRef();

    forAll(mesh_.C(), cell)
    {
        // Set state in cantera
        const scalar e = heCells[cell];
        const scalar v = 1.0 / rhoCells[cell];

        forAll(Y_, i)
        {
            ctY_[i] = Y_[i].internalField()[cell];
        }

        gas_->setMassFractions_NoNorm(ctY_.data());
        gas_->setState_UV(e, v);

        double Temperature = gas_->RT() / Cantera::GasConstant;
        TCells[cell] = Temperature;
        pCells[cell] = gas_->pressure();

        updateMixtureProperties(cell,
                                psiCells,
                                cpCells,
                                cvCells,
                                WCells,
                                rhoCells);
        updateTransportProperties(Temperature, rhoCells[cell], WCells[cell], cell);
    }

    // Update boundary fields
    // Get reference to boundary fields (-> call to function automatically saves
    // old timestep for time derivative)
    volScalarField::Boundary& rhoBf = rho_.boundaryFieldRef();
    volScalarField::Boundary& heBf = he_.boundaryFieldRef();
    volScalarField::Boundary& pBf = p_.boundaryFieldRef();

    volScalarField::Boundary& psiBf = psi_.boundaryFieldRef();
    volScalarField::Boundary& cpBf = cp_.boundaryFieldRef();
    volScalarField::Boundary& cvBf = cv_.boundaryFieldRef();
    volScalarField::Boundary& WBf = W_.boundaryFieldRef();
    volScalarField::Boundary& TBf = T_.boundaryFieldRef();

    forAll(mesh_.boundary(), patch)
    {
        if (TBf[patch].fixesValue())
        {
            forAll(mesh_.boundary()[patch], face)
            {
                // Set state in cantera
                const scalar ctT = TBf[patch][face];
                const scalar v = 1.0 / rhoBf[patch][face];
                const scalar ctP = pBf[patch][face];

                forAll(Y_, i)
                {
                    ctY_[i] = Y_[i].boundaryField()[patch][face];
                }
                gas_->setState_TPY(ctT, ctP, ctY_.data());
                // TODO: Does the pressure need to be updated here?
                // p_[cell] = gas_->pressure();

                heBf[patch][face] = getHeFromGasObject();

                updateMixtureProperties(patch,
                                        face,
                                        psiBf,
                                        cpBf,
                                        cvBf,
                                        WBf,
                                        rhoBf);
                updateTransportProperties(TBf[patch][face], rhoBf[patch][face], WBf[patch][face], patch, face);
            }
        }
        else
        {
            forAll(mesh_.boundary()[patch], face)
            {
                // Set state in cantera
                const scalar e = heBf[patch][face];
                const scalar v = 1.0 / rhoBf[patch][face];

                forAll(Y_, i)
                {
                    ctY_[i] = Y_[i].boundaryField()[patch][face];
                }

                gas_->setMassFractions_NoNorm(ctY_.data());
                gas_->setState_UV(e, v);

                // TODO: Does the pressure need to be updated here?
                pBf[patch][face] = gas_->pressure();

                TBf[patch][face] = gas_->RT() / Cantera::GasConstant;
                updateMixtureProperties(patch,
                                        face,
                                        psiBf,
                                        cpBf,
                                        cvBf,
                                        WBf,
                                        rhoBf);
                updateTransportProperties(TBf[patch][face], rhoBf[patch][face], WBf[patch][face], patch, face);
            }
        }
    }
}

//
// Calculation of mixture properties
//
double Foam::canteraThermo::getHeFromGasObject() const
{
    if (energyType_ == "ha")
    {
        return gas_->enthalpy_mass();
    }
    else if (energyType_ == "ea")
    {
        return gas_->intEnergy_mass();
    }
    else
    {
        // ToDo: Make this pretty
        NotImplemented;
        return 0;
    }
}

void Foam::canteraThermo::updateMixtureProperties(
    const label cell,
    scalarField& psiCells,
    scalarField& cpCells,
    scalarField& cvCells,
    scalarField& WCells,
    scalarField& rhoCells)
{
    const scalar RT = gas_->RT();
    const scalar W = gas_->meanMolecularWeight();
    psiCells[cell] = W / RT;
    cpCells[cell] = gas_->cp_mass();
    cvCells[cell] = gas_->cv_mass();
    WCells[cell] = W;
    rhoCells[cell] = gas_->density();
    speedOfSound_[cell] = gas_->soundSpeed();
    // Specific anbsolute enthalpy of species at current temperature
    std::vector<scalar> tempHSensi(nSpecies_);
    gas_->getEnthalpy_RT_ref(tempHSensi.data());
    forAll(hi_, i)
    {
        heSensi_[i][cell] = tempHSensi[i];
        hi_[i][cell] = RT / Wi_[i].value()
                     * heSensi_[i].internalField()[cell];
    }
}


void Foam::canteraThermo::updateMixtureProperties(
    const label patch,
    const label face,
    volScalarField::Boundary& psiBf,
    volScalarField::Boundary& cpBf,
    volScalarField::Boundary& cvBf,
    volScalarField::Boundary& WBf,
    volScalarField::Boundary& rhoBf)
{
    const scalar RT = gas_->RT();
    const scalar W = gas_->meanMolecularWeight();
    psiBf[patch][face] = W / RT;
    cpBf[patch][face] = gas_->cp_mass();
    cvBf[patch][face] = gas_->cv_mass();
    WBf[patch][face] = W;
    rhoBf[patch][face] = gas_->density();
    speedOfSound_.boundaryFieldRef(false)[patch][face] = gas_->soundSpeed();
    std::vector<scalar> tempHSensi(nSpecies_);
    gas_->getEnthalpy_RT_ref(tempHSensi.data());
    forAll(hi_, i)
    {
        heSensi_[i].boundaryFieldRef(false)[patch][face] = tempHSensi[i];
        hi_[i].boundaryFieldRef(false)[patch][face] = RT / Wi_[i].value()
                                                    * heSensi_[i].boundaryField()[patch][face];
    }
}


//
// Calculation of transport properties
//
void Foam::canteraThermo::updateTransportProperties(
    const scalar T,
    const scalar rho,
    const scalar W,
    const label cell)
{
    lambda_[cell] = trans_->thermalConductivity();
    alpha_[cell] = trans_->thermalConductivity() / gas_->cp_mass();
    mu_[cell] = trans_->viscosity();
    nu_[cell] = trans_->viscosity() / gas_->density();

    // Calculate the diffusion coefficients.
    if (transportModel_ == "unity-Lewis-number")
    {
    }
    else if (transportModel_ == "mixture-averaged")
    {
        calculateDMix(cell);
        if (useSoret_)
        {
            calculateDthermMix(T, rho, W, cell);
        }
    }
    else if (transportModel_ == "multicomponent")
    {
        calculateDMulti(cell);
        if (useSoret_)
        {
            calculateDthermMulti(cell);
        }
    }
}

void Foam::canteraThermo::updateTransportProperties(
    const scalar T,
    const scalar rho,
    const scalar W,
    const label patch,
    const label face)
{
    lambda_.boundaryFieldRef(false)[patch][face] = trans_->thermalConductivity();
    alpha_.boundaryFieldRef(false)[patch][face] = trans_->thermalConductivity() / gas_->cp_mass();
    mu_.boundaryFieldRef(false)[patch][face] = trans_->viscosity();
    nu_.boundaryFieldRef(false)[patch][face] = trans_->viscosity() / gas_->density();

    // Calculate the diffusion coefficients.
    if (transportModel_ == "unity-Lewis-number")
    {
    }
    else if (transportModel_ == "mixture-averaged")
    {
        calculateDMix(patch, face);
        if (useSoret_)
        {
            calculateDthermMix(T, rho, W, patch, face);
        }
    }
    else if (transportModel_ == "multicomponent")
    {
        calculateDMulti(patch, face);
        if (useSoret_)
        {
            calculateDthermMulti(patch, face);
        }
    }
}

void Foam::canteraThermo::calculateDMix(
    const label cell)
{
    trans_->getMixDiffCoeffs(ctD_.data());
    for (label i = 0; i < nSpecies_; ++i)
    {
        D_[i][cell] = ctD_[i];
    }
}

void Foam::canteraThermo::calculateDMix(
    const label patch,
    const label face)
{
    trans_->getMixDiffCoeffs(ctD_.data());
    for (label i = 0; i < nSpecies_; ++i)
    {
        D_[i].boundaryFieldRef(false)[patch][face] = ctD_[i];
    }
}

void Foam::canteraThermo::calculateDMulti(
    const label cell)
{
    // Note:
    // Cantera returns the array in column-major fashion (Fortran-like)
    trans_->getMultiDiffCoeffs(nSpecies_, ctD_.data());

    for (label i = 0; i < nSpecies_; ++i)
    {
        for (label j = 0; j < nSpecies_; ++j)
        {
            // Todo: Check if conversion from column to row major is correct
            D_[nSpecies_ * i + j][cell] = ctD_[nSpecies_ * j + i];
        }
    }
}

void Foam::canteraThermo::calculateDMulti(
    const label patch,
    const label face)
{
    // Note:
    // Cantera returns the array in column-major fashion (Fortran-like)
    trans_->getMultiDiffCoeffs(nSpecies_, ctD_.data());
    for (label i = 0; i < nSpecies_; ++i)
    {
        for (label j = 0; j < nSpecies_; ++j)
        {
            // Todo: Check if conversion from column to row major is correct
            D_[nSpecies_ * i + j].boundaryFieldRef(false)[patch][face] = ctD_[nSpecies_ * j + i];
        }
    }
}

void Foam::canteraThermo::calculateDthermMix(
    const scalar T,
    const scalar rho,
    const scalar W,
    const label cell)
{
    std::vector<scalar> k_Tivec(nSpecies_);
    calculateKTi(T, rho, W, k_Tivec);
    for (label i = 0; i < nSpecies_; ++i)
    {
        Dtherm_[i][cell] = rho * (Wi_[i].value() / W) * D_[i][cell] * k_Tivec[i];
    }
}

void Foam::canteraThermo::calculateDthermMix(
    const scalar T,
    const scalar rho,
    const scalar W,
    const label patch,
    const label face)
{
    std::vector<scalar> k_Tivec(nSpecies_);
    calculateKTi(T, rho, W, k_Tivec);
    for (label i = 0; i < nSpecies_; ++i)
    {
        Dtherm_[i].boundaryFieldRef(false)[patch][face] = rho * (Wi_[i].value() / W) * D_[i].boundaryField()[patch][face] * k_Tivec[i];
    }
}

void Foam::canteraThermo::calculateKTi(const scalar T, const scalar rho, const scalar W, std::vector<scalar>& k_Tivec) const
{
    if (trans_->CKMode())
    {
        // In cantera there are special considerations for these cases, which are not implemented here
        NotImplemented;
    }

    size_t viscosityPolyDegree = 5;
    // COLL_INT_POLY_DEGREE is defined as 8 in cantera
    size_t collIntPolyDegree = 8;

    const scalar logT = log(T);

    std::vector<scalar> logTvec(viscosityPolyDegree);

    logTvec[0] = 1.0;
    logTvec[1] = logT;
    logTvec[2] = logT * logT;
    logTvec[3] = logT * logT * logT;
    logTvec[4] = logT * logT * logT * logT;

    std::vector<scalar> X(nSpecies_);
    gas_->getMoleFractions(X.data());
    std::vector<scalar> DBinary_ij(nSpecies_ * nSpecies_);
    trans_->getBinaryDiffCoeffs(nSpecies_, DBinary_ij.data());

    std::vector<scalar> nu(nSpecies_);
    std::vector<scalar> nu_i_coeffs(viscosityPolyDegree);
    std::vector<scalar> Phi_ij(nSpecies_ * nSpecies_);
    std::vector<scalar> a(nSpecies_);

    std::vector<scalar> aStar_coeffs(collIntPolyDegree);
    std::vector<scalar> bStar_coeffs(collIntPolyDegree);
    std::vector<scalar> cStar_coeffs(collIntPolyDegree);

    // needed for polar correction
    std::vector<double> sigma(nSpecies_);
    std::vector<double> alpha(nSpecies_);
    std::vector<double> epsilon(nSpecies_);
    std::vector<bool> polar(nSpecies_);
    std::vector<std::vector<double>> dipole(nSpecies_, std::vector<double>(nSpecies_));

    for (label i = 0; i < nSpecies_; ++i)
    {
        trans_->getViscosityPolynomial(i, nu_i_coeffs.data());
        // for degree == 5: the polynomial fit is done for sqrt(visc/sqrt(T))
        nu[i] = sqrt(sqrt(T)) * Cantera::dot5(logTvec, nu_i_coeffs) * sqrt(sqrt(T)) * Cantera::dot5(logTvec, nu_i_coeffs);

        std::shared_ptr<Cantera::Species> si = gas_->species(i);
        const Cantera::GasTransportData* sptrani = dynamic_cast<Cantera::GasTransportData*>(si->transport.get());
        epsilon[i] = sptrani->well_depth;
        sigma[i] = sptrani->diameter;
        alpha[i] = sptrani->polarizability;
        dipole[i][i] = sptrani->dipole;
        polar[i] = (sptrani->dipole > 0);
    }

    for (label i = 0; i < nSpecies_; ++i)
    {
        // there are different definitions for Phi_ij in the literature: sometimes there is a 8 in the denominator, sometimes not
        scalar X_j_Phi_ij = 0.0;
        scalar k_Ti = 0;
        for (label j = 0; j < nSpecies_; ++j)
        {
            if (j != i)
            {
                const scalar factor1 = (1 + sqrt(nu[i] / nu[j]) * sqrt(sqrt((Wi_[j].value() / Wi_[i].value()))));
                // it needs to be considered, if sqrt(8.0) is needed in the denominator or not!
                Phi_ij[nSpecies_ * j + i] = (factor1 * factor1) / (sqrt(1 + (Wi_[i].value() / Wi_[j].value())));
                X_j_Phi_ij += X[j] * Phi_ij[nSpecies_ * j + i];

                // This is at the moment not needed here, since only diagnonal elements are used in polar correction
                // dipole[i][j] = sqrt(dipole[i][i]*dipole[j][j]);
            }
        }
        const scalar lambda_iMon = (15 / 4) * (Cantera::GasConstant * nu[i] / Wi_[i].value());

        if (X[i] > 1e-18)
        {
            a[i] = lambda_iMon * (1 / (1 + (1.065 / (2 * sqrt(2.0) * X[i])) * X_j_Phi_ij));
        }
        else
        {
            a[i] = 0;
        }

        for (label j = 0; j < nSpecies_; ++j)
        {
            trans_->getCollisionIntegralPolynomial(i, j, aStar_coeffs.data(), bStar_coeffs.data(), cStar_coeffs.data());

            scalar eps_ij = sqrt(epsilon[i] * epsilon[j]);
            // Corrections of well depth
            double f_eps;
            // no correction if both are nonpolar, or both are polar
            if (polar[i] == polar[j])
            {
                f_eps = 1.0;
            }
            else
            {
                // corrections to the effective well depth
                // if one is polar and one is non-polar
                label kp = (polar[i] ? i : j); // the polar one
                label knp = (i == kp ? j : i); // the nonpolar one
                double d3np, d3p, alpha_star, mu_p_star, xi;
                d3np = std::pow(sigma[knp], 3);
                d3p = std::pow(sigma[kp], 3);
                alpha_star = alpha[knp] / d3np;
                mu_p_star = dipole[kp][kp] / sqrt(4 * Cantera::Pi * Cantera::epsilon_0 * d3p * epsilon[kp]);
                xi = 1.0 + 0.25 * alpha_star * mu_p_star * mu_p_star * sqrt(epsilon[kp] / epsilon[knp]);
                f_eps = xi * xi;
            }

            eps_ij *= f_eps;

            scalar logTstar = log(T * Cantera::Boltzmann / eps_ij);
            const scalar C_ij = Cantera::poly8(logTstar, cStar_coeffs.data());

            k_Ti += ((1.2 * C_ij - 1) / (DBinary_ij[nSpecies_ * j + i])) * ((ctY_[i] * a[j] - ctY_[j] * a[i]) / (Wi_[i].value() + Wi_[j].value()));
        }
        k_Tivec[i] = k_Ti * W * W / (Cantera::GasConstant * rho);
    }
}

void Foam::canteraThermo::calculateDthermMulti(
    const label cell)
{
    trans_->getThermalDiffCoeffs(ctDtherm_.data());
    for (label i = 0; i < nSpecies_; ++i)
    {
        Dtherm_[i][cell] = ctDtherm_[i];
    }
}

void Foam::canteraThermo::calculateDthermMulti(
    const label patch,
    const label face)
{
    trans_->getThermalDiffCoeffs(ctDtherm_.data());
    for (label i = 0; i < nSpecies_; ++i)
    {
        Dtherm_[i].boundaryFieldRef(false)[patch][face] = ctDtherm_[i];
    }
}

//
// Getter and setter
//
Foam::label Foam::canteraThermo::nSpecies() const
{
    return nSpecies_;
}

const std::vector<std::string> Foam::canteraThermo::speciesNames() const
{
    return gas_->speciesNames();
}

const std::vector<std::string> Foam::canteraThermo::elementNames() const
{
    return gas_->elementNames();
}

const size_t Foam::canteraThermo::speciesIndex(const std::string& name) const
{
    return gas_->speciesIndex(name);
}

const size_t Foam::canteraThermo::elementIndex(const std::string& name) const
{
    return gas_->elementIndex(name);
}

double Foam::canteraThermo::nAtoms(size_t speciesIndex, size_t elementIndex) const
{
    return gas_->nAtoms(speciesIndex, elementIndex);
}

Foam::word Foam::canteraThermo::transportModel() const
{
    return transportModel_;
}

bool Foam::canteraThermo::useSoret() const
{
    return useSoret_;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::rho() const
{
    return rho_;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::rho(const label patchi) const
{
    return rho_.boundaryField()[patchi];
}

void Foam::canteraThermo::correctRho(
    const Foam::volScalarField& deltaRho,
    const dimensionedScalar& rhoMin,
    const dimensionedScalar& rhoMax)
{
    rho_ += deltaRho;
    rho_ = max(rho_, rhoMin);
    rho_ = min(rho_, rhoMax);
}

void Foam::canteraThermo::correctRho(const Foam::volScalarField& deltaRho)
{
    rho_ += deltaRho;
}

Foam::volScalarField& Foam::canteraThermo::rhoRef()
{
    return rho_;
}

Foam::volScalarField& Foam::canteraThermo::he()
{
    return he_;
}

const Foam::volScalarField& Foam::canteraThermo::he() const
{
    return he_;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::he(
    const scalarField& p,
    const scalarField& T,
    const labelList& cells) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    std::vector<scalar> ctY(nSpecies_);

    forAll(T, celli)
    {
        forAll(Y_, i)
        {
            ctY[i] = Y_[i].internalField()[cells[celli]];
        }

        gas_->setState_TPY(T[celli], p[celli], ctY.data());
        he[celli] = getHeFromGasObject();
    }
    return the;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::he(
    const scalarField& p,
    const scalarField& T,
    const label patchi) const
{
    tmp<scalarField> the(new scalarField(T.size()));
    scalarField& he = the.ref();

    std::vector<scalar> ctY(nSpecies_);

    forAll(T, facei)
    {
        forAll(Y_, i)
        {
            ctY[i] = Y_[i].boundaryField()[patchi][facei];
        }

        gas_->setState_TPY(T[facei], p[facei], ctY.data());
        he[facei] = getHeFromGasObject();
    }
    return the;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::Cp() const
{
    return cp_;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::Cv() const
{
    return cv_;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::gamma() const
{
    return cp_ / cv_;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::Cpv(
    const scalarField& p,
    const scalarField& T,
    const label patchi) const
{
    tmp<scalarField> tCpv(new scalarField(T.size()));
    scalarField& Cpv = tCpv.ref();

    std::vector<scalar> ctY(nSpecies_);

    forAll(T, facei)
    {
        forAll(Y_, i)
        {
            ctY[i] = Y_[i].boundaryField()[patchi][facei];
        }
        gas_->setState_TPY(T[facei], p[facei], ctY.data());
        if (energyType_ == "ha")
        {
            Cpv[facei] = gas_->cp_mass();
        }
        else if (energyType_ == "ea")
        {
            Cpv[facei] = gas_->cv_mass();
        }
    }
    return tCpv;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::Cpv() const
{
    if (energyType_ == "ha")
    {
        return cp_;
    }
    else if (energyType_ == "ea")
    {
        return cv_;
    }
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::CpByCpv() const
{
    return cp_ / Cpv();
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::CpByCpv(
    const scalarField& p,
    const scalarField& T,
    const label patchi) const
{
    return cp_.boundaryField()[patchi] / Cpv(p, T, patchi);
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::W() const
{
    return W_;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::kappa() const
{
    Foam::tmp<Foam::volScalarField> kappa(this->Cp() * this->alpha_);
    kappa.ref().rename("kappa");
    return kappa;
};

Foam::tmp<Foam::scalarField> Foam::canteraThermo::kappa(
    const Foam::label patchi) const
{
    return this->cp_.boundaryField()[patchi] * this->alpha_.boundaryField()[patchi];
};

const Foam::volScalarField& Foam::canteraThermo::psi() const
{
    return psi_;
}

Foam::PtrList<Foam::volScalarField>& Foam::canteraThermo::Y()
{
    return Y_;
}

const Foam::PtrList<Foam::volScalarField>& Foam::canteraThermo::Y() const
{
    return Y_;
}

Foam::volScalarField& Foam::canteraThermo::hi(const label i)
{
    return hi_[i];
}

Foam::dimensionedScalar& Foam::canteraThermo::Wi(const label i)
{
    return Wi_[i];
}

Foam::scalar Foam::canteraThermo::Wi(const label i) const
{
    return Wi_[i].value();
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::mu() const
{
    return mu_;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::mu(const label patchi) const
{
    return mu_.boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::nu() const
{
    return nu_;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::nu(const label patchi) const
{
    return nu_.boundaryField()[patchi];
}

const Foam::volScalarField& Foam::canteraThermo::D(const label i) const
{
    return D_[i];
}

const Foam::volScalarField& Foam::canteraThermo::Dtherm(const label i) const
{
    return Dtherm_[i];
}

const Foam::volScalarField& Foam::canteraThermo::lambda() const
{
    return lambda_;
}
const Foam::volScalarField& Foam::canteraThermo::speedOfSound() const
{
    return speedOfSound_;
}

Foam::tmp<Foam::volScalarField> Foam::canteraThermo::alphahe() const
{
    tmp<Foam::volScalarField> alphaEff(this->CpByCpv() * this->alpha_);
    alphaEff.ref().rename("alphahe");
    return alphaEff;
}

Foam::tmp<Foam::scalarField> Foam::canteraThermo::alphahe(const label patchi) const
{
    return this->CpByCpv(
               this->p_.boundaryField()[patchi],
               this->T_.boundaryField()[patchi],
               patchi)
         * this->alpha_.boundaryField()[patchi];
}

Foam::volScalarField& Foam::canteraThermo::alpha()
{
    return alpha_;
}
