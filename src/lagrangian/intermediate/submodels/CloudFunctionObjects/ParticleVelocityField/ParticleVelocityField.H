/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2017 OpenFOAM Foundation
    Copyright (C) 2020 OpenCFD Ltd.
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
    Foam::ParticleVelocityField

Group
    grpLagrangianIntermediateFunctionObjects

Description
    Writes particle counts and velocities to a field at a given interval (every
    n times steps). They can also be averaged and particle count can be
    integrated over time.

    cloudFunctions
    {
      particleVelocityField
      {
        type particleVelocityField;
    interval 10;
    average true;
    dMin 10e-6; // default 0.0
    dMax 20e-6; // default 1.0
    suffix _10_20; //suffix of field
      }
    }

SourceFiles
    ParticleVelocityField.C

\*---------------------------------------------------------------------------*/

#ifndef ParticleVelocityField_H
#define ParticleVelocityField_H

#include "CloudFunctionObject.H"
#include "volFields.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                       Class ParticleVelocityField Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ParticleVelocityField
: public CloudFunctionObject<CloudType>
{
    // Private Data

    // Typedefs

    //- Convenience typedef for parcel type
    typedef typename CloudType::parcelType parcelType;


    //- Void fraction field
    autoPtr<volVectorField> UParticlePtr_;

    autoPtr<volScalarField> nParticlePtr_;
    autoPtr<volScalarField> nParcelPtr_;

    autoPtr<volVectorField> UParticleIntPtr_;

    autoPtr<volScalarField> nParticleIntPtr_;
    autoPtr<volScalarField> nParcelIntPtr_;

    autoPtr<volVectorField> UParticleMeanPtr_;

    scalar dmax_;

    scalar dmin_;

    scalar interval_;

    word suffix_;

    label counter_;

    bool average_;


protected:
    // Protected Member Functions

    //- Write post-processing info
    virtual void write();


public:
    //- Runtime type information
    TypeName("particleVelocityField");


    // Constructors

    //- Construct from dictionary
    ParticleVelocityField(
        const dictionary& dict,
        CloudType& owner,
        const word& modelName);

    //- Construct copy
    ParticleVelocityField(const ParticleVelocityField<CloudType>& vf);

    //- Construct and return a clone
    virtual autoPtr<CloudFunctionObject<CloudType>> clone() const
    {
        return autoPtr<CloudFunctionObject<CloudType>>(
            new ParticleVelocityField<CloudType>(*this));
    }


    //- Destructor
    virtual ~ParticleVelocityField() = default;


    // Member Functions

    // Evaluation

    //- Pre-evolve hook
    virtual void preEvolve(
        const typename parcelType::trackingData& td);

    //- Post-evolve hook
    virtual void postEvolve(
        const typename parcelType::trackingData& td);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "ParticleVelocityField.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
