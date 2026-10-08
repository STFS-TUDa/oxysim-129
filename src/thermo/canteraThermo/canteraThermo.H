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

Class
    Foam::canteraThermo

Group
    grpThermo

Description
    Class implementing the interface of fluidThermo (and basicSpecieMixture) as
    well as providing access to additional thermo and transport properties
    through Cantera.

Authors

SourceFiles
    canteraThermo.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <string>

#include "basicThermoSTFS.H"
#include "canteraSpecieMixture.H"

#include "constants.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Cantera
{
class Solution;
class ThermoPhase;
class Transport;
}
namespace Foam
{

//- Global auxiliary function to get the list of species names from Cantera
wordList getSpecieNames(const word& mechanismFile, const word& phaseName);

/*---------------------------------------------------------------------------*\
                         Class canteraThermo Declaration
\*---------------------------------------------------------------------------*/
/** @addtogroup canteraThermo
 *  @{
 */
class canteraThermo
: public basicThermoSTFS,
  public canteraSpecieMixture
{
private:
    static inline std::string handleTransportModel(std::string);

protected:
    // Protected data

    //
    // Cantera related properties and helper fields
    //
    //- Transport model
    const std::string transportModel_;
    //- Switch for soret diffusion
    const bool useSoret_;
    // //- Switch for Thermo type
    // const bool psiThermo_;
    //- Cantera solution object
    std::shared_ptr<Cantera::Solution> sol_;
    //- Cantera gas object
    std::shared_ptr<Cantera::ThermoPhase> gas_;
    //- Cantera transport object
    std::shared_ptr<Cantera::Transport> trans_;
    //- Number of species
    const label nSpecies_;
    // Auxillary arrays
    std::vector<scalar> ctY_;
    std::vector<scalar> ctD_;
    std::vector<scalar> ctDtherm_;

    //
    // Mixture properties
    //
    //- Energy type to use inside the thermo
    // (current options are ha: absolute enthalpy; ea: absolute energy)
    const word energyType_;
    //- Energy field [J/kg]
    volScalarField he_;
    //- Species specific enthalpies / internal energies [J/kg]
    PtrList<volScalarField> hi_;
    //- Sensible enthalpy / internal energy of species [J/kg]
    PtrList<volScalarField> heSensi_;
    //- Species molar weights [kg/kmol]
    PtrList<dimensionedScalar> Wi_;
    //- Sensible enthalpy / internal energy of species at reference temperature [J/kg]
    PtrList<dimensionedScalar> heSensRefi_;
    //- Thermodynamic density [kg/m^3]
    volScalarField rho_;
    //- Compressibility [s^2/m^2]
    volScalarField psi_;
    //- Mean specific heat capacity at constant pressure [J/kg/K]
    volScalarField cp_;
    //- Mean specific heat capacity at constant volume [J/kg/K]
    volScalarField cv_;
    //- Mean molecular weight [kg/kmol]
    volScalarField W_;

    //
    // Transport properties
    //
    //- Diffusion coefficients [m^2/s]
    PtrList<volScalarField> D_;
    //- Thermal diffusion coefficients [kg/m/s]
    PtrList<volScalarField> Dtherm_;
    //- Kinematic visosity [m^2/s]
    volScalarField nu_;
    //- Dynamic visosity [kg/m/s]
    volScalarField mu_;
    //- Thermal conductivity [kg*m/s^3/K]
    volScalarField lambda_;
    //- SpeedOfSound [m/s]
    volScalarField speedOfSound_;

    // Protected member functions
    //- Correct the absolute enthalpy/internal energy field boundaries
    void heBoundaryCorrection(volScalarField& h);

    // Protected member functions

    //- get the energy corresponding to energyType_ from gas_
    double getHeFromGasObject() const;

    //- Mixture properties update (psi, cp, cv, W, heSensi, hi)
    // for index cell in internal field
    void updateMixtureProperties(const label cell,
                                 scalarField& psiCells,
                                 scalarField& cpCells,
                                 scalarField& cvCells,
                                 scalarField& WBf,
                                 scalarField& rhoCells);
    //- Mixture properties update (psi, cp, cv, W, heSensi, hi)
    // for patch and face of a boundary field
    void updateMixtureProperties(
        const label patch,
        const label face,
        volScalarField::Boundary& psiBf,
        volScalarField::Boundary& cpBf,
        volScalarField::Boundary& cvBf,
        volScalarField::Boundary& WBf,
        volScalarField::Boundary& rhoBf);

    //- Transport properties update (lambda, alpha, mu, nu, diffusion coefficients)
    // for index cell in internal field
    inline void updateTransportProperties(
        const scalar T,
        const scalar rho,
        const scalar W,
        const label cell);
    //- Transport properties update (lambda, alpha, mu, nu, diffusion coefficients)
    // for patch and face of a boundary field
    inline void updateTransportProperties(
        const scalar T,
        const scalar rho,
        const scalar W,
        const label patch,
        const label face);

    //- Calculate species diffusion coefficients (mixture-average)
    // for index cell in internal field
    inline void calculateDMix(
        const label cell);
    //- Calculate species diffusion coefficients (mixture-average)
    // for patch and face of a boundary field
    inline void calculateDMix(
        const label patch,
        const label face);

    //- Calculate species diffusion coefficients (multi-component)
    // for index cell in internal field
    inline void calculateDMulti(
        const label cell);
    //- Calculate species diffusion coefficients (multi-component)
    // for patch and face of a boundary field
    inline void calculateDMulti(
        const label patch,
        const label face);

    //- Calculate thermal diffusion coefficients (mixutre-average)
    // for index cell in internal field
    inline void calculateDthermMix(
        const scalar T,
        const scalar rho,
        const scalar W,
        const label cell);
    //- Calculate thermal diffusion coefficients (mixutre-average)
    // for patch and face of a boundary field
    inline void calculateDthermMix(
        const scalar T,
        const scalar rho,
        const scalar W,
        const label patch,
        const label face);

    //- Helper function for the calculation of the thermal diffusion coefficients (mixture-average)
    inline void calculateKTi(
        const scalar T,
        const scalar rho,
        const scalar W,
        std::vector<scalar>& k_Tivec) const;

    //- Calculate thermal diffusion coefficients (multi-component)
    // for index cell in internal field
    inline void calculateDthermMulti(
        const label cell);
    //- Calculate thermal diffusion coefficients (multi-component)
    // for patch and face of a boundary field
    inline void calculateDthermMulti(
        const label patch,
        const label face);

public:
    //- Runtime type information
    TypeName("canteraThermo");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        canteraThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName));

    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable(
        autoPtr,
        canteraThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName));

    // Constructors
    //- This constructor is not implemented.
    canteraThermo(
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName);

    //- Construct from mesh, phase and dictionary name
    canteraThermo(
        const fvMesh&,
        const word& phaseName = "",
        const word& dictionaryName = "thermophysicalProperties",
        bool initialUpdate = true);

    //- Construct as copy
    canteraThermo(const canteraThermo&) = delete;

    //- Selector
    static autoPtr<canteraThermo> New(
        const fvMesh& mesh,
        const word& phaseName = word::null);

    //- Selector
    static autoPtr<canteraThermo> New(
        const fvMesh& mesh,
        const word& phaseName = word::null,
        const word& dictionaryName = word::null);

    ~canteraThermo();

    // Member Functions
    //- Update mixture and transport properties
    //- from T, p and Y
    void updateTPY();

    void updateTPY(
        const volScalarField& p,
        volScalarField& T,
        volScalarField& psi,
        volScalarField& cp,
        volScalarField& cv,
        volScalarField& W,
        volScalarField& rho,
        volScalarField& he,
        const bool doOldTimes);
    // normal gas phase solver (low-mach)
    //- from h, p and Y
    void updateHPY();
    // compressible
    //- from rho, p and Y
    void updateRPY();
    //- from e, p and Y
    void updateERY();

    //- Return information of thermo and transport models
    label nSpecies() const;
    const std::vector<std::string> speciesNames() const;
    const std::vector<std::string> elementNames() const;
    const size_t speciesIndex(const std::string& name) const;
    const size_t elementIndex(const std::string& name) const;
    double nAtoms(size_t speciesIndex, size_t elementIndex) const;

    word transportModel() const;
    bool useSoret() const;

    // Return access to fields
    // Thermodynamic properties inherited from basicThermo
    //- Density [kg/m^3]
    tmp<volScalarField> rho() const;
    //- Density for patch [kg/m^3]
    tmp<scalarField> rho(const label patchi) const;

    //- Add the given density correction to the density field.
    //  Used to update the density field following pressure solution
    //  Limit thermo rho between rhoMin and rhoMax
    // ToDo: Add to unit test
    virtual void correctRho(
        const volScalarField& deltaRho,
        const dimensionedScalar& rhoMin,
        const dimensionedScalar& rhoMax);

    //- Add the given density correction to the density field.
    //  Used to update the density field following pressure solution
    // ToDo: Add to unit test
    virtual void correctRho(
        const volScalarField& deltaRho);

    //- Enthalpy/Internal energy [J/kg]
    tmp<scalarField> he(
        const scalarField& p,
        const scalarField& T,
        const labelList& cells) const override;
    tmp<scalarField> he(
        const scalarField& p,
        const scalarField& T,
        const label patch) const override;
    //  Non-const access allowed for transport equations
    volScalarField& he() override;
    //- Enthalpy/Internal energy [J/kg]
    const volScalarField& he() const override;
    //- Heat capacity at constant pressure [J/kg/K]
    tmp<volScalarField> Cp() const override;
    //- Heat capacity at constant volume [J/kg/K]
    tmp<volScalarField> Cv() const override;
    //- Gamma = Cp/Cv []
    tmp<volScalarField> gamma() const override;
    //- Heat capacity at constant pressure/volume for patch [J/kg/K]
    tmp<scalarField> Cpv(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const override;
    //- Heat capacity at constant pressure/volume [J/kg/K]
    tmp<volScalarField> Cpv() const override;
    //- Heat capacity ratio []
    tmp<volScalarField> CpByCpv() const override;
    tmp<scalarField> CpByCpv(const scalarField& p, const scalarField& T, const label patchi) const override;
    //- Molecular weight [kg/kmol]
    tmp<volScalarField> W() const override;
    //- Thermal diffusivity for temperature of mixture [J/m/s/K]
    tmp<volScalarField> kappa() const override;
    //- Speed pf Sound [m/s]
    [[nodiscard]] const volScalarField& speedOfSound() const;
    //- Thermal diffusivity for temperature
    //  of mixture for patch [J/m/s/K]
    tmp<scalarField> kappa(
        const label patchi) const override;
    //- Compressibility [s^2/m^2]
    const volScalarField& psi() const override;

    // Thermodynamic properties introduced here
    //- Species mass fractions [-]
    volScalarField& rhoRef();
    [[nodiscard]] PtrList<volScalarField>& Y();
    [[nodiscard]] const PtrList<volScalarField>& Y() const;
    //- Specific internal energy / specific enthalpy of species i [J/kg]
    [[nodiscard]] volScalarField& hi(const label i);
    //- Molecular weight of species i [kg/mol]
    [[nodiscard]] dimensionedScalar& Wi(const label i);
    scalar Wi(const label i) const;

    // Transport properties inherited from basicThermo
    //- Dynamic viscosity of mixture [kg/m/s]
    tmp<volScalarField> mu() const override;
    tmp<scalarField> mu(const label patchi) const override;

    //- Kinematic viscosity of mixture [m^2/s]
    tmp<volScalarField> nu() const override;
    tmp<scalarField> nu(const label patchi) const override;

    // Transport properties introduced here
    //- Diffusion coefficients [m^2/s]
    [[nodiscard]] const volScalarField& D(const label i) const;
    //- Thermal diffusion coefficients [m^2/s]
    [[nodiscard]] const volScalarField& Dtherm(const label i) const;
    //- Thermal conductivity [W/m/K]
    [[nodiscard]] const volScalarField& lambda() const;

    tmp<volScalarField> alphahe() const override;
    tmp<scalarField> alphahe(const label patchi) const override;

    // Non-const access to alpha
    volScalarField& alpha();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**@} */
} // End namespace Foam

// ************************************************************************* //
