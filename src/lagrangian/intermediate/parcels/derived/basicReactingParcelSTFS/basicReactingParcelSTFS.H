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
    Foam::basicReactingParcelSTFS

Description
    Definition of basic reacting parcel

SourceFiles
    basicReactingParcelSTFS.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "contiguous.H"
#include "particle.H"
#include "KinematicParcel.H"
#include "ThermoParcelSTFS.H"
#include "ReactingParcelSTFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
typedef ReactingParcelSTFS<ThermoParcelSTFS<KinematicParcel<particle>>>
    basicReactingParcelSTFS;

//- Non-contiguous data for basicReactingParcelSTFS
template<>
struct is_contiguous<basicReactingParcelSTFS> : std::false_type
{
};
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
