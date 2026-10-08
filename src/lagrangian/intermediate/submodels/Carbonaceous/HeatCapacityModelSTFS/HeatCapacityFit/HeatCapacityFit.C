/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "HeatCapacityFit.H"
// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType>
Foam::HeatCapacityFit<CloudType>::HeatCapacityFit(
    const dictionary& dict,
    CloudType& owner)
: HeatCapacityModelSTFS<CloudType>(dict, owner, typeName),
  coeffs_(this->coeffDict().lookup("coefficients"))
{
    poly = Polynomial<5>(coeffs_);
}

template<class CloudType>
Foam::HeatCapacityFit<CloudType>::HeatCapacityFit(
    const HeatCapacityFit<CloudType>& dm)
: HeatCapacityModelSTFS<CloudType>(dm.owner_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::HeatCapacityFit<CloudType>::calculate(
    const scalar Cp0,
    const scalar T,
    const scalarField& YGas,
    scalar& Cp) const
{
    Cp = poly.value(T);
}


// ************************************************************************* //
