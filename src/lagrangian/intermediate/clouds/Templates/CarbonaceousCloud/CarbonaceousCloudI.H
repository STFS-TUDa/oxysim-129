/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2013 OpenFOAM Foundation
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

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
inline const Foam::CarbonaceousCloud<CloudType>&
Foam::CarbonaceousCloud<CloudType>::cloudCopy() const
{
    return *cloudCopyPtr_;
}

template<class CloudType>
inline const typename CloudType::particleType::constantProperties&
Foam::CarbonaceousCloud<CloudType>::constProps() const
{
    return constProps_;
}


template<class CloudType>
inline typename CloudType::particleType::constantProperties&
Foam::CarbonaceousCloud<CloudType>::constProps()
{
    return constProps_;
}


template<class CloudType>
inline const Foam::DevolatilisationModelSTFS<
    Foam::CarbonaceousCloud<CloudType>>&
Foam::CarbonaceousCloud<CloudType>::devolatilisationSTFS() const
{
    return devolatilisationModelSTFS_;
}


template<class CloudType>
inline Foam::DevolatilisationModelSTFS<
    Foam::CarbonaceousCloud<CloudType>>&
Foam::CarbonaceousCloud<CloudType>::devolatilisationSTFS()
{
    return *devolatilisationModelSTFS_;
}

template<class CloudType>
inline const Foam::HeatCapacityModelSTFS<
    Foam::CarbonaceousCloud<CloudType>>&
Foam::CarbonaceousCloud<CloudType>::heatCapacitySTFS() const
{
    return heatCapacityModelSTFS_;
}


template<class CloudType>
inline Foam::HeatCapacityModelSTFS<
    Foam::CarbonaceousCloud<CloudType>>&
Foam::CarbonaceousCloud<CloudType>::heatCapacitySTFS()
{
    return *heatCapacityModelSTFS_;
}
// ************************************************************************* //
