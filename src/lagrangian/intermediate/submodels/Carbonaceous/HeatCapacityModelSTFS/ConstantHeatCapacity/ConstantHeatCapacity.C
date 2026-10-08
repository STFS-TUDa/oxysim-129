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

#include "ConstantHeatCapacity.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //


template<class CloudType>
Foam::ConstantHeatCapacity<CloudType>::ConstantHeatCapacity(
    const dictionary&,
    CloudType& owner)
: HeatCapacityModelSTFS<CloudType>(owner)
{
}

template<class CloudType>
Foam::ConstantHeatCapacity<CloudType>::ConstantHeatCapacity(
    const ConstantHeatCapacity<CloudType>& dm)
: HeatCapacityModelSTFS<CloudType>(dm.owner_)
{
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ConstantHeatCapacity<CloudType>::calculate(
    const scalar Cp0,
    const scalar T,
    const scalarField& YGas,
    scalar& Cp) const
{
    Cp = Cp0;
}


// ************************************************************************* //
