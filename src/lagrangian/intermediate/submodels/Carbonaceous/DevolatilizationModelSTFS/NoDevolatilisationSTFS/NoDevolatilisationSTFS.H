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
    Foam::NoDevolatilisationSTFS

Group
    grpLagrangianIntermediateDevolatilisationSubModels

Description
    Dummy devolatilisation model for 'none'

\*---------------------------------------------------------------------------*/

#pragma once

#include "DevolatilisationModelSTFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class NoDevolatilisationSTFS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class NoDevolatilisationSTFS
: public DevolatilisationModelSTFS<CloudType>
{

public:
    //- Runtime type information
    TypeName("none");


    // Constructors

    //- Construct from dictionary
    NoDevolatilisationSTFS(const dictionary&, CloudType& owner);

    //- Construct copy
    NoDevolatilisationSTFS(const NoDevolatilisationSTFS<CloudType>& dm);

    //- Construct and return a clone
    virtual autoPtr<DevolatilisationModelSTFS<CloudType>> clone() const
    {
        return autoPtr<DevolatilisationModelSTFS<CloudType>>(
            new NoDevolatilisationSTFS<CloudType>(*this));
    }


    //- Destructor
    virtual ~NoDevolatilisationSTFS();


    // Member Functions

    //- Flag to indicate whether model activates devolatilisation model
    virtual bool active() const;

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
#include "NoDevolatilisationSTFS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
