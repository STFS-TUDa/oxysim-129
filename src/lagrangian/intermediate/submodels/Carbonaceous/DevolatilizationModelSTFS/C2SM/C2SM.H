/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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

Class
    Foam::C2SM

Group
    grpLagrangianIntermediateDevolatilisationSubModels

Description
    Competing Two Step Model
    - acts on a per-specie basis
    - Rates given by Arrhenius eqns

        kappa1 = A1.exp(- E1/R.T)
        kappa2 = A2.exp(- E2/R.T)

        dMassDV = dt.massVolatile.(kappa1.a1 + kappa2.a2)
        dMassC = dt.massVolatile.(kappa1.(1-a1) + kappa2.(1-a2))

    Where:
        kappa = rate constant
        Ai    = activation energy (user input)
        Ei    = pre-exponential factor (user input)
        ai    = weighing between reactions
        R     = universal gas constant
        T     = temperature

    Usage:

        C2SMCoeffs
        {
            volatileData
            (
            (H2 0.484 1.65e5 4.8988e7 0.752 4.312e8 1.324e8)      // (name alpha1 A1 E1 alpha2 A2 E2)
            (H2O 0.484 1.65e5 4.8988e7 0.752 4.312e8 1.324e8)     // (name alpha1 A1 E1 alpha2 A2 E2)
            (CO 0.484 1.65e5 4.8988e7 0.752 4.312e8 1.324e8)      // (name alpha1 A1 E1 alpha2 A2 E2)
            (CO2 0.484 1.65e5 4.8988e7 0.752 4.312e8 1.324e8)     // (name alpha1 A1 E1 alpha2 A2 E2)
            (CH4 0.484 1.65e5 4.8988e7 0.752 4.312e8 1.324e8)     // (name alpha1 A1 E1 alpha2 A2 E2)
            (A1XC6H6 0.484 1.65e5 4.8988e7 0.752 4.312e8 1.324e8) // (name alpha1 A1 E1 alpha2 A2 E2)
            );

            residualCoeff 1e-6;
        }

\*---------------------------------------------------------------------------*/

#pragma once

#include "DevolatilisationModelSTFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class C2SM Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class C2SM
: public DevolatilisationModelSTFS<CloudType>
{
    //- Helper class to store specie-local volatile data
    class volatileData
    {
        // Private Data

        //- Specie name
        word name_;

        scalar a1_;
        //- Activation energy
        scalar A11_;

        //- Pre-exponential factor
        scalar E1_;

        scalar a2_;
        //- Activation energy
        scalar A12_;

        //- Pre-exponential factor
        scalar E2_;

    public:
        // Constructors

        //- Default construct
        volatileData()
        : name_(),
          a1_(0),
          A11_(0),
          E1_(0),
          a2_(0),
          A12_(0),
          E2_(0)
        {
        }

        //- Construct from Istream
        explicit volatileData(Istream& is)
        : name_(is),
          a1_(readScalar(is)),
          A11_(readScalar(is)),
          E1_(readScalar(is)),
          a2_(readScalar(is)),
          A12_(readScalar(is)),
          E2_(readScalar(is))
        {
        }


        // Member Functions

        // Access

        //- Return const access to the name
        const word& name() const
        {
            return name_;
        }

        scalar a1() const
        {
            return a1_;
        }

        scalar A11() const
        {
            return A11_;
        }

        scalar E1() const
        {
            return E1_;
        }

        scalar a2() const
        {
            return a2_;
        }

        scalar A12() const
        {
            return A12_;
        }

        scalar E2() const
        {
            return E2_;
        }

        // IOstream Operators

        //- Read from Istream
        friend Istream& operator>>(Istream& is, volatileData& vd)
        {
            is.readBegin("volatileData");
            is >> vd.name_ >> vd.a1_ >> vd.A11_ >> vd.E1_ >> vd.a2_ >> vd.A12_ >> vd.E2_;
            is.readEnd("volatileData");

            return is;
        }

        //- Write to Ostream
        friend Ostream& operator<<(Ostream& os, const volatileData& vd)
        {
            os << token::BEGIN_LIST
               << vd.name_ << token::SPACE
               << vd.a1_ << token::SPACE
               << vd.A11_ << token::SPACE
               << vd.E1_ << token::SPACE
               << vd.a2_ << token::SPACE
               << vd.A12_ << token::SPACE
               << vd.E2_
               << token::END_LIST;

            return os;
        }
    };


    // Private Data

    // Model Constants

    //- List of volatile data - (name A1 E)
    List<volatileData> volatileData_;

    //- List of initial volatile mass fractions
    List<scalar> YVolatile0_;

    //- Mapping between local and cloud gaseous species
    List<label> volatileToGasMap_;

    //- Volatile residual coefficient (0-1)
    //  When the fraction of volatiles are depleted below this
    //  threshold, combustion can occur
    const scalar residualCoeff_;


public:
    //- Declare type-name, virtual type (with debug switch)
    TypeName("C2SM");


    // Constructors

    //- Construct from dictionary
    C2SM(
        const dictionary& dict,
        CloudType& owner);

    //- Construct and return a clone
    virtual autoPtr<DevolatilisationModelSTFS<CloudType>> clone() const
    {
        return autoPtr<DevolatilisationModelSTFS<CloudType>>(
            new C2SM<CloudType>(*this));
    }


    //- Destructor
    virtual ~C2SM() = default;


    // Member Functions

    //- Update model
    virtual void calculate(
        const scalar dt,
        const scalar age,
        const scalar mass0,
        const scalar mass,
        const scalar T,
        const scalarField& YGasEff,
        const scalarField& YLiquidEff,
        const scalarField& YSolidEff,
        label& canCombust,
        scalarField& dMassDV,
        scalarField& dMassC) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "C2SM.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
