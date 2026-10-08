/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2014-2017 OpenFOAM Foundation
    Copyright (C) 2020-2021 OpenCFD Ltd.
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
    Foam::canteraSpecieMixture

Group
    grpMixtures

Description
    Cantera specific specialization of basicSpecieMixture for a mixture
    consisting of a number for molecular species.

SourceFiles
    canteraSpecieMixture.C

\*---------------------------------------------------------------------------*/

#ifndef canteraSpecieMixture_H
#define canteraSpecieMixture_H

#include "basicSpecieMixture.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                 Class canteraSpecieMixture Declaration
\*---------------------------------------------------------------------------*/
/** @addtogroup canteraThermo
 *  @{
 */
class canteraSpecieMixture
: public basicSpecieMixture
{
    //- Private members

    //- Cantera gas constant
    const scalar Rc_;

    //- Polynomial array
    static constexpr int nCoeffs_ = 7;
    typedef FixedList<scalar, nCoeffs_> coeffArray;

    typedef FixedList<scalar, 5> coeffArray5;

    //- Nasa Enthalpy Low coeffs
    List<coeffArray> nasaCoeffsLow_;

    //- Nasa Enthalpy High coeffs
    List<coeffArray> nasaCoeffsHigh_;

    //- Nasa medium Temperature
    List<scalar> nasaTMid_;

    //- Species Molar Weight
    std::vector<scalar> W_;

    //- Species Chemical Enthalpy
    std::vector<scalar> Hc_;

    //- Species Viscosity coeffs
    List<coeffArray5> muCoeffs_;

    //- Species Thermal conductivity coeffs
    List<coeffArray5> lambdaCoeffsHigh_;

    //- Return the coefficients corresponding to the given temperature
    const coeffArray& coeffs(const scalar T, const label speciei) const;

public:
    //- Run time type information
    TypeName("canteraSpecieMixture");

    //- The base class of the mixture
    typedef canteraSpecieMixture basicMixtureType;

    // Constructors

    //- Construct from dictionary, species names, mesh and phase name
    canteraSpecieMixture(
        const dictionary& thermoDict,
        const wordList& specieNames,
        const fvMesh& mesh,
        const word& phaseName);


    //- Destructor
    virtual ~canteraSpecieMixture() = default;


    // Member Functions

    // Per specie properties

    //- Molecular weight of the given specie [kg/kmol]
    scalar W(const label speciei) const override;

    //- Chemical enthalpy [J/kg]
    scalar Hc(const label speciei) const override;


    // Per specie thermo properties

    //- Heat capacity at constant pressure [J/(kg K)]
    scalar Cp(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Heat capacity at constant volume [J/(kg K)]
    scalar Cv(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Enthalpy/Internal energy [J/kg]
    scalar HE(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Absolute enthalpy [J/kg]
    scalar Ha(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Sensible enthalpy [J/kg]
    scalar Hs(
        const label speciei,
        const scalar p,
        const scalar T) const override;


    //- Entropy [J/(kg K)]
    scalar S(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Sensible internal energy [J/kg]
    scalar Es(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Gibbs free energy [J/kg]
    scalar G(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Helmholtz free energy [J/kg]
    scalar A(
        const label speciei,
        const scalar p,
        const scalar T) const override;


    // Per specie transport properties

    //- Dynamic viscosity [kg/m/s]
    scalar mu(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Thermal conductivity [W/m/K]
    scalar kappa(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Thermal diffusivity of enthalpy [kg/m/s]
    scalar alphah(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Density [kg/m3]
    scalar rho(
        const label speciei,
        const scalar p,
        const scalar T) const override;

    //- Species composition
    autoPtr<speciesCompositionTable> specieComposition() const override
    {
        return autoPtr<speciesCompositionTable>::New();
    }
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**@} */
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //