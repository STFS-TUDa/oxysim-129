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

#pragma once

#include "fluidThermo.H"
#include "runTimeSelectionTables.H"
#include "dimensionedScalar.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicThermoSTFS Declaration
\*---------------------------------------------------------------------------*/
/** @defgroup basicThermo
 *   @ingroup basic
 *  @{
 */

class basicThermoSTFS : public fluidThermo
{
protected:
    // Protected data

    const fvMesh& mesh_;

    // Turbulence related properties
    scalar Sct_;
    scalar Prt_;

public:
    // Constructors & Destructors
    //- Runtime type information
    TypeName("basicThermoSTFS");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        basicThermoSTFS,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName));

    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable(
        autoPtr,
        basicThermoSTFS,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName));

    // Constructors

    //- Construct from mesh and phase name
    basicThermoSTFS(
        const fvMesh&,
        const word& phaseName);


    basicThermoSTFS(
        const fvMesh& mesh,
        const dictionary& dict,
        const word& phaseName);

    //- Construct from mesh,dictionary,phase name with a single temperature
    basicThermoSTFS(

        const fvMesh&,
        const word& phaseName,
        const word& dictionaryName);

    //- Construct as copy
    basicThermoSTFS(const basicThermoSTFS&) = delete;

    //- Selector
    static autoPtr<basicThermoSTFS> New(
        const fvMesh& mesh,
        const word& phaseName = word::null);

    //- Selector
    static autoPtr<basicThermoSTFS> New(
        const fvMesh& mesh,
        const word& phaseName = word::null,
        const word& dictionaryName = word::null);

    //- Destructor
    virtual ~basicThermoSTFS();

    // Member Functions

    //******************** Additional getter functions ********************//

    void setTurbulenceProperties();
    scalar Sct() const;
    scalar Prt() const;

    //******************** Basic Thermo Functions ********************//

    //- Update properties
    virtual void correct() { NotImplemented; };

    //- Return the name of the thermo physics
    virtual word thermoName() const
    {
        NotImplemented;
        return "dummy";
    };

    //- Return true if the equation of state is incompressible
    //  i.e. rho != f(p)
    virtual bool incompressible() const
    {
        NotImplemented;
        return true;
    };

    //- Return true if the equation of state is isochoric
    //  i.e. rho = const
    virtual bool isochoric() const
    {
        NotImplemented;
        return false;
    };

    // Access to thermodynamic state variables

    //- Density [kg/m^3]
    virtual tmp<volScalarField> rho() const
    {
        NotImplemented;
        return T_;
    };

    //- Density for patch [kg/m^3]
    virtual tmp<scalarField> rho(const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Enthalpy/Internal energy [J/kg]
    //  Non-const access allowed for transport equations
    virtual volScalarField& he()
    {
        NotImplemented;
        return T_;
    };

    //- Enthalpy/Internal energy [J/kg]
    virtual const volScalarField& he() const
    {
        NotImplemented;
        return T_;
    };

    //- Enthalpy/Internal energy
    //  for given pressure and temperature [J/kg]
    virtual tmp<volScalarField> he(
        const volScalarField& p,
        const volScalarField& T) const
    {
        NotImplemented;
        return T_;
    };

    //- Enthalpy/Internal energy for cell-set [J/kg]
    virtual tmp<scalarField> he(
        const scalarField& p,
        const scalarField& T,
        const labelList& cells) const
    {
        NotImplemented;
        return T_;
    };

    //- Enthalpy/Internal energy for patch [J/kg]
    virtual tmp<scalarField> he(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Chemical enthalpy [J/kg]
    virtual tmp<volScalarField> hc() const
    {
        NotImplemented;
        return T_;
    };

    //- Temperature from enthalpy/internal energy for cell-set
    virtual tmp<scalarField> THE(
        const scalarField& h,
        const scalarField& p,
        const scalarField& T0, // starting temperature
        const labelList& cells) const
    {
        NotImplemented;
        return T_;
    };

    //- Temperature from enthalpy/internal energy for patch
    virtual tmp<scalarField> THE(
        const scalarField& h,
        const scalarField& p,
        const scalarField& T0, // starting temperature
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };


    // Fields derived from thermodynamic state variables

    //- Heat capacity at constant pressure [J/kg/K]
    virtual tmp<volScalarField> Cp() const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity at constant pressure for patch [J/kg/K]
    virtual tmp<scalarField> Cp(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity using pressure and temperature [J/kg/K]
    virtual tmp<scalarField> Cp(
        const scalarField& p,
        const scalarField& T,
        const labelList& cells) const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity at constant volume [J/kg/K]
    virtual tmp<volScalarField> Cv() const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity at constant volume for patch [J/kg/K]
    virtual tmp<scalarField> Cv(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Density from pressure and temperature from EoS
    virtual tmp<scalarField> rhoEoS(
        const scalarField& p,
        const scalarField& T,
        const labelList& cells) const
    {
        NotImplemented;
        return T_;
    };

    //- Gamma = Cp/Cv []
    virtual tmp<volScalarField> gamma() const
    {
        NotImplemented;
        return T_;
    };

    //- Gamma = Cp/Cv for patch []
    virtual tmp<scalarField> gamma(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity at constant pressure/volume [J/kg/K]
    virtual tmp<volScalarField> Cpv() const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity at constant pressure/volume for patch [J/kg/K]
    virtual tmp<scalarField> Cpv(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity ratio []
    virtual tmp<volScalarField> CpByCpv() const
    {
        NotImplemented;
        return T_;
    };

    //- Heat capacity ratio for patch []
    virtual tmp<scalarField> CpByCpv(
        const scalarField& p,
        const scalarField& T,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Molecular weight [kg/kmol]
    virtual tmp<volScalarField> W() const
    {
        NotImplemented;
        return T_;
    };


    // Fields derived from transport state variables

    //- Thermal diffusivity for temperature of mixture [J/m/s/K]
    virtual tmp<volScalarField> kappa() const
    {
        NotImplemented;
        return T_;
    };

    //- Thermal diffusivity for temperature
    //  of mixture for patch [J/m/s/K]
    virtual tmp<scalarField> kappa(
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Thermal diffusivity for energy of mixture [kg/m/s]
    virtual tmp<volScalarField> alphahe() const
    {
        NotImplemented;
        return T_;
    };

    //- Thermal diffusivity for energy of mixture for patch [kg/m/s]
    virtual tmp<scalarField> alphahe(const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Effective thermal diffusivity for temperature
    //  of mixture [J/m/s/K]
    virtual tmp<volScalarField> kappaEff(
        const volScalarField&) const
    {
        NotImplemented;
        return T_;
    };

    //- Effective thermal diffusivity for temperature
    //  of mixture for patch [J/m/s/K]
    virtual tmp<scalarField> kappaEff(
        const scalarField& alphat,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Effective thermal diffusivity of mixture [kg/m/s]
    virtual tmp<volScalarField> alphaEff(
        const volScalarField& alphat) const
    {
        NotImplemented;
        return T_;
    };

    //- Effective thermal diffusivity of mixture for patch [kg/m/s]
    virtual tmp<scalarField> alphaEff(
        const scalarField& alphat,
        const label patchi) const
    {
        NotImplemented;
        return T_;
    };


    //- Read thermophysical properties dictionary
    virtual bool read()
    {
        NotImplemented;
        return false;
    };

    //******************** Fluid Thermo Functions ********************//

    // Access to thermodynamic state variables

    //- Add the given density correction to the density field.
    //  Used to update the density field following pressure solution
    //  Limit thermo rho between rhoMin and rhoMax
    virtual void correctRho(
        const volScalarField& deltaRho,
        const dimensionedScalar& rhoMin,
        const dimensionedScalar& rhoMax) { NotImplemented; };


    //- Add the given density correction to the density field.
    //  Used to update the density field following pressure solution
    virtual void correctRho(
        const volScalarField& deltaRho) { NotImplemented; };

    //- Compressibility [s^2/m^2]
    virtual const volScalarField& psi() const
    {
        NotImplemented;
        return T_;
    };


    // Access to transport state variables

    //- Dynamic viscosity of mixture [kg/m/s]
    virtual tmp<volScalarField> mu() const
    {
        NotImplemented;
        return T_;
    };

    //- Dynamic viscosity of mixture for patch [kg/m/s]
    virtual tmp<scalarField> mu(const label patchi) const
    {
        NotImplemented;
        return T_;
    };

    //- Kinematic viscosity of mixture [m^2/s]
    virtual tmp<volScalarField> nu() const
    {
        NotImplemented;
        return T_;
    };

    //- Kinematic viscosity of mixture for patch [m^2/s]
    virtual tmp<scalarField> nu(const label patchi) const
    {
        NotImplemented;
        return T_;
    };
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**@}*/
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
