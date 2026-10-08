/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::NoHeatCapacity

Group
    grpLagrangianIntermediateHeatCapacitySubModels

Description
    Dummy devolatilisation model for 'none'

\*---------------------------------------------------------------------------*/

#pragma once

#include "HeatCapacityModelSTFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class NoHeatCapacity Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoHeatCapacity
: public HeatCapacityModelSTFS<CloudType>
{

public:
    //- Runtime type information
    TypeName("none");


    // Constructors

    //- Construct from dictionary
    NoHeatCapacity(const dictionary&, CloudType& owner);

    //- Construct copy
    NoHeatCapacity(const NoHeatCapacity<CloudType>& dm);

    //- Construct and return a clone
    virtual autoPtr<HeatCapacityModelSTFS<CloudType>> clone() const
    {
        return autoPtr<HeatCapacityModelSTFS<CloudType>>(
            new NoHeatCapacity<CloudType>(*this));
    }


    //- Destructor
    virtual ~NoHeatCapacity();


    // Member Functions

    //- Flag to indicate whether model activates devolatilisation model
    virtual bool active() const;

    //- Update model
    virtual void calculate(
        const scalar,
        const scalar,
        const scalarField&,
        scalar&) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "NoHeatCapacity.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
