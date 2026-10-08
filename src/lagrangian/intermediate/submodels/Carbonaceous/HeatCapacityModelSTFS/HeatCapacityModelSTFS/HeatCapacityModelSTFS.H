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
    Foam::HeatCapacityModelSTFS

Group
    grpLagrangianIntermediateHeatCapacitySubModels

Description
    Templated devolatilisation model class

SourceFiles
    HeatCapacityModelSTFS.C
    HeatCapacityModelNewSTFS.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "IOdictionary.H"
#include "autoPtr.H"
#include "runTimeSelectionTables.H"
#include "CloudSubModelBase.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                    Class HeatCapacityModelSTFS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class HeatCapacityModelSTFS
: public CloudSubModelBase<CloudType>
{
protected:
    // Protected data

    //- Mass of lagrangian phase converted


public:
    //- Runtime type information
    TypeName("heatCapacityModelSTFS");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        HeatCapacityModelSTFS,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner),
        (dict, owner));


    // Constructors

    //- Construct null from owner
    HeatCapacityModelSTFS(CloudType& owner);

    //- Construct from dictionary
    HeatCapacityModelSTFS(
        const dictionary& dict,
        CloudType& owner,
        const word& type);

    //- Construct copy
    // HeatCapacityModelSTFS(const HeatCapacityModelSTFS<CloudType>& dm); //??

    //- Construct and return a clone
    virtual autoPtr<HeatCapacityModelSTFS<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~HeatCapacityModelSTFS();


    //- Selector
    static autoPtr<HeatCapacityModelSTFS<CloudType>> New(
        const dictionary& dict,
        CloudType& owner);


    // Member Functions

    //- Update model
    virtual void calculate(
        const scalar Cp0,
        const scalar T,
        const scalarField& YGas,
        scalar& Cp) const = 0;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeHeatCapacityModelSTFS(CloudType)                \
                                                            \
    typedef Foam::CloudType::carbonaceousCloudType          \
        carbonaceousCloudType;                              \
    defineNamedTemplateTypeNameAndDebug(                    \
        Foam::HeatCapacityModelSTFS<carbonaceousCloudType>, \
        0);                                                 \
    namespace Foam                                          \
    {                                                       \
    defineTemplateRunTimeSelectionTable(                    \
        HeatCapacityModelSTFS<carbonaceousCloudType>,       \
        dictionary);                                        \
    }


#define makeHeatCapacityModelSTFSType(SS, CloudType)                         \
                                                                             \
    typedef Foam::CloudType::carbonaceousCloudType                           \
        carbonaceousCloudType;                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<carbonaceousCloudType>, 0); \
                                                                             \
    Foam::HeatCapacityModelSTFS<carbonaceousCloudType>::                     \
        adddictionaryConstructorToTable<Foam::SS<carbonaceousCloudType>>     \
            add##SS##CloudType##carbonaceousCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "HeatCapacityModelSTFS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
