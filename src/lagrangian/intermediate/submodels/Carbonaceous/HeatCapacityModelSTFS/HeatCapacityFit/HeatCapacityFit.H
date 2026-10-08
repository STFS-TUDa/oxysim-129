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
    Foam::HeatCapacityFit

Group
    grpLagrangianIntermediateHeatCapacitySubModels

Description
    Heat capacity is described with a polynomial and is a function of temperature.
    Cp = a + b * T + c * T^2 + d * T^3 + e * T^4

    Usage:
    heatCapacityFitCoeffs
    {
        coefficients
        (-1.198e3 10.465 -0.0124 5.888e-6 -9.60e-10); // a b c d e
    }

\*---------------------------------------------------------------------------*/

#pragma once

#include "HeatCapacityModelSTFS.H"
#include "Polynomial.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
              Class HeatCapacityFit Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HeatCapacityFit
: public HeatCapacityModelSTFS<CloudType>
{

public:
    //- Declare type-name, virtual type (with debug switch)
    TypeName("heatCapacityFit");


    // Constructors

    //- Construct from dictionary
    HeatCapacityFit(const dictionary& dict, CloudType& owner);

    //- Construct copy
    HeatCapacityFit(const HeatCapacityFit<CloudType>& dm);

    //- Construct and return a clone
    virtual autoPtr<HeatCapacityModelSTFS<CloudType>> clone() const
    {
        return autoPtr<HeatCapacityModelSTFS<CloudType>>(
            new HeatCapacityFit<CloudType>(*this));
    }


    //- Destructor
    virtual ~HeatCapacityFit() = default;


    // Member Functions

    //- Update model
    virtual void calculate(
        const scalar Cp0,
        const scalar T,
        const scalarField& YGas,
        scalar& Cp) const;

private:
    List<scalar> coeffs_;
    Polynomial<5> poly;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "HeatCapacityFit.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
