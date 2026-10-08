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
    Foam::COxidationKineticDiffusionLimitedRateSTFS

Group
    grpLagrangianIntermediateSurfaceReactionSubModels

Description
    Kinetic/diffusion limited rate surface reaction model for carbonaceous fuel parcels.
    Limited to:

        C(s) + Sb*O2 -> CO

    where Sb is the stoichiometry of the reaction

    The char is oxidized to CO.

    Usage
        COxidationKineticDiffusionLimitedRateSTFSCoeffs
        {
            Sb              1;
            C1              7.43e-13; //Cdiff
            C2              1.886e2; // A
            E               1.286e8; // Ea
        dummySpecies    AR;
        }
        The dummy species gets the mass released to CO. In tabulated chemistry cases it can make sense to use a

\*---------------------------------------------------------------------------*/

#pragma once

#include "SurfaceReactionModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward class declarations
template<class CloudType>
class COxidationKineticDiffusionLimitedRateSTFS;

/*---------------------------------------------------------------------------*\
            Class COxidationKineticDiffusionLimitedRateSTFS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class COxidationKineticDiffusionLimitedRateSTFS
: public SurfaceReactionModel<CloudType>
{
    // Private data

    // Model constants

    //- Stoichiometry of reaction
    const scalar Sb_;

    //- Mass diffusion limited rate constant, C1
    const scalar C1_;

    //- Kinetics limited rate pre-exponential constant, C2
    const scalar C2_;

    //- Kinetics limited rate activation energy
    const scalar E_;


    // Addressing

    //- Cs positions in global/local lists
    label CsLocalId_;

    //- O2 position in global list
    label O2GlobalId_;

    //- CO2 positions in global list
    label CO2GlobalId_;

    label COGlobalId_;
    label COcharGlobalId_;
    // Local copies of thermo properties

    //- Molecular weight of C [kg/kmol]
    scalar WC_;

    //- Molecular weight of O2 [kg/kmol]
    scalar WO2_;

    //- Formation enthalpy for CO2 [J/kg]
    scalar HcCO2_;

    //- Formation enthalpy for CO [J/kg]
    scalar HcCO_;


public:
    //- Runtime type information
    TypeName("COxidationKineticDiffusionLimitedRateSTFS");


    // Constructors

    //- Construct from dictionary
    COxidationKineticDiffusionLimitedRateSTFS(
        const dictionary& dict,
        CloudType& owner);

    //- Construct copy
    COxidationKineticDiffusionLimitedRateSTFS(
        const COxidationKineticDiffusionLimitedRateSTFS<CloudType>& srm);

    //- Construct and return a clone
    virtual autoPtr<SurfaceReactionModel<CloudType>> clone() const
    {
        return autoPtr<SurfaceReactionModel<CloudType>>(
            new COxidationKineticDiffusionLimitedRateSTFS<CloudType>(*this));
    }


    //- Destructor
    virtual ~COxidationKineticDiffusionLimitedRateSTFS() = default;


    // Member Functions

    //- Update surface reactions
    virtual scalar calculate(
        const scalar dt,
        const scalar Re,
        const scalar nu,
        const label celli,
        const scalar d,
        const scalar T,
        const scalar Tc,
        const scalar pc,
        const scalar rhoc,
        const scalar mass,
        const scalarField& YGas,
        const scalarField& YLiquid,
        const scalarField& YSolid,
        const scalarField& YMixture,
        const scalar N,
        scalarField& dMassGas,
        scalarField& dMassLiquid,
        scalarField& dMassSolid,
        scalarField& dMassSRCarrier) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "COxidationKineticDiffusionLimitedRateSTFS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
