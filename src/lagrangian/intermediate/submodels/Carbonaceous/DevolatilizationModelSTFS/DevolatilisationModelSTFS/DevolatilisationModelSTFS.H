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
    Foam::DevolatilisationModelSTFS

Group
    grpLagrangianIntermediateDevolatilisationSubModels

Description
    Templated devolatilisation model class

SourceFiles
    DevolatilisationModelSTFS.C
    DevolatilisationModelNewSTFS.C

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
                    Class DevolatilisationModelSTFS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class DevolatilisationModelSTFS
: public CloudSubModelBase<CloudType>
{
protected:
    // Protected data

    //- Mass of lagrangian phase converted
    scalar dMass_;


public:
    //- Runtime type information
    TypeName("devolatilisationModelSTFS");

    //- Declare runtime constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        DevolatilisationModelSTFS,
        dictionary,
        (
            const dictionary& dict,
            CloudType& owner),
        (dict, owner));


    // Constructors

    //- Construct null from owner
    DevolatilisationModelSTFS(CloudType& owner);

    //- Construct from dictionary
    DevolatilisationModelSTFS(
        const dictionary& dict,
        CloudType& owner,
        const word& type);

    //- Construct copy
    DevolatilisationModelSTFS(const DevolatilisationModelSTFS<CloudType>& dm);

    //- Construct and return a clone
    virtual autoPtr<DevolatilisationModelSTFS<CloudType>> clone() const = 0;


    //- Destructor
    virtual ~DevolatilisationModelSTFS();


    //- Selector
    static autoPtr<DevolatilisationModelSTFS<CloudType>> New(
        const dictionary& dict,
        CloudType& owner);


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
        scalarField& dMassC) const = 0;

    //- Add to devolatilisation mass
    void addToDevolatilisationMass(const scalar dMass);


    // I-O

    //- Write injection info to stream
    virtual void info(Ostream& os);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeDevolatilisationModelSTFS(CloudType)                \
                                                                \
    typedef Foam::CloudType::carbonaceousCloudType              \
        carbonaceousCloudType;                                  \
    defineNamedTemplateTypeNameAndDebug(                        \
        Foam::DevolatilisationModelSTFS<carbonaceousCloudType>, \
        0);                                                     \
    namespace Foam                                              \
    {                                                           \
    defineTemplateRunTimeSelectionTable(                        \
        DevolatilisationModelSTFS<carbonaceousCloudType>,       \
        dictionary);                                            \
    }


#define makeDevolatilisationModelSTFSType(SS, CloudType)                     \
                                                                             \
    typedef Foam::CloudType::carbonaceousCloudType                           \
        carbonaceousCloudType;                                               \
    defineNamedTemplateTypeNameAndDebug(Foam::SS<carbonaceousCloudType>, 0); \
                                                                             \
    Foam::DevolatilisationModelSTFS<carbonaceousCloudType>::                 \
        adddictionaryConstructorToTable<Foam::SS<carbonaceousCloudType>>     \
            add##SS##CloudType##carbonaceousCloudType##ConstructorToTable_;


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "DevolatilisationModelSTFS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
