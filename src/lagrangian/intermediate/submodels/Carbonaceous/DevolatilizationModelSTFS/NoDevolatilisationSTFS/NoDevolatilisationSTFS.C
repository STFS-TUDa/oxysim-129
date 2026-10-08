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

#include "NoDevolatilisationSTFS.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoDevolatilisationSTFS<CloudType>::NoDevolatilisationSTFS(
    const dictionary&,
    CloudType& owner)
: DevolatilisationModelSTFS<CloudType>(owner)
{
}


template<class CloudType>
Foam::NoDevolatilisationSTFS<CloudType>::NoDevolatilisationSTFS(
    const NoDevolatilisationSTFS<CloudType>& dm)
: DevolatilisationModelSTFS<CloudType>(dm.owner_)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::NoDevolatilisationSTFS<CloudType>::~NoDevolatilisationSTFS()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
bool Foam::NoDevolatilisationSTFS<CloudType>::active() const
{
    return false;
}


template<class CloudType>
void Foam::NoDevolatilisationSTFS<CloudType>::calculate(
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalar,
    const scalarField&,
    const scalarField&,
    const scalarField&,
    label& canCombust,
    scalarField&,
    scalarField&) const
{
    // Model does not stop combustion taking place
    canCombust = true;
}


// ************************************************************************* //
