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

\*---------------------------------------------------------------------------*/

#include "ParticleVelocityField.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleVelocityField<CloudType>::write()
{
    if (UParticlePtr_)
    {
        UParticlePtr_->write();
    }
    else
    {
        FatalErrorInFunction
            << "thetaPtr not valid" << abort(FatalError);
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ParticleVelocityField<CloudType>::ParticleVelocityField(
    const dictionary& dict,
    CloudType& owner,
    const word& modelName)
: CloudFunctionObject<CloudType>(dict, owner, modelName, typeName),
  UParticlePtr_(nullptr),
  nParticlePtr_(nullptr),
  nParcelPtr_(nullptr),
  UParticleIntPtr_(nullptr),
  nParticleIntPtr_(nullptr),
  nParcelIntPtr_(nullptr),
  UParticleMeanPtr_(nullptr),
  dmax_(this->coeffDict().getOrDefault("dMax", 1.0)),
  dmin_(this->coeffDict().getOrDefault("dMin", 0.0)),
  interval_(this->coeffDict().getOrDefault("interval", 10)),
  suffix_(this->coeffDict().template getOrDefault<word>("suffix", "")),
  average_(this->coeffDict().template getOrDefault<bool>("average", false)),
  counter_(0)
{
    Info << "Initializing particleVelocityField with dmin=" << dmin_ << " m, dmax=" << dmax_ << " m, interval=" << interval_ << ", and suffix=" << suffix_ << "." << endl;
}


template<class CloudType>
Foam::ParticleVelocityField<CloudType>::ParticleVelocityField(
    const ParticleVelocityField<CloudType>& vf)
: CloudFunctionObject<CloudType>(vf),
  UParticlePtr_(nullptr),
  nParticlePtr_(nullptr),
  nParcelPtr_(nullptr),
  UParticleIntPtr_(nullptr),
  nParticleIntPtr_(nullptr),
  nParcelIntPtr_(nullptr),
  UParticleMeanPtr_(nullptr),
  dmax_(this->coeffDict().getOrDefault("dMax", 1.0)),
  dmin_(this->coeffDict().getOrDefault("dMin", 0.0)),
  interval_(this->coeffDict().getOrDefault("interval", 10)),
  suffix_(this->coeffDict().template getOrDefault<Foam::word>("suffix", "")),
  average_(this->coeffDict().template getOrDefault<bool>("average", false)),
  counter_(0)
{
    Info << "Initializing particleVelocityField with Dmin = " << dmin_ << " m, dMax = " << dmax_ << " m, interval = " << interval_ << ", and suffix=" << suffix_ << "." << endl;
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::ParticleVelocityField<CloudType>::preEvolve(
    const typename parcelType::trackingData& td)
{
    if (counter_ == interval_)
    {
        counter_ = 0;
    }
    if (counter_ == 0)
    {
        if (UParticlePtr_)
        {
            UParticlePtr_->primitiveFieldRef() = vector(0.0, 0.0, 0.0);
            nParticlePtr_->primitiveFieldRef() = 0.0;
            nParcelPtr_->primitiveFieldRef() = 0.0;
        }
        else
        {
            const fvMesh& mesh = this->owner().mesh();

            UParticlePtr_.reset(
                new volVectorField(
                    IOobject(
                        this->owner().name() + "UParticle" + suffix_,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimVelocity, vector::zero)));
            nParticlePtr_.reset(
                new volScalarField(
                    IOobject(
                        this->owner().name() + "nParticle" + suffix_,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedScalar(dimless, Zero)));
            nParcelPtr_.reset(
                new volScalarField(
                    IOobject(
                        this->owner().name() + "nParcel" + suffix_,
                        mesh.time().timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedScalar(dimless, Zero)));

            if (average_)
            {
                UParticleIntPtr_.reset(
                    new volVectorField(
                        IOobject(
                            this->owner().name() + "UParticleInt" + suffix_,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE),
                        mesh,
                        dimensionedVector(dimVelocity, vector::zero)));
                nParticleIntPtr_.reset(
                    new volScalarField(
                        IOobject(
                            this->owner().name() + "nParticleInt" + suffix_,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE),
                        mesh,
                        dimensionedScalar(dimless, Zero)));
                nParcelIntPtr_.reset(
                    new volScalarField(
                        IOobject(
                            this->owner().name() + "nParcelInt" + suffix_,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::READ_IF_PRESENT,
                            IOobject::AUTO_WRITE),
                        mesh,
                        dimensionedScalar(dimless, Zero)));
                UParticleMeanPtr_.reset(
                    new volVectorField(
                        IOobject(
                            this->owner().name() + "UParticleMean" + suffix_,
                            mesh.time().timeName(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE),
                        mesh,
                        dimensionedVector(dimVelocity, vector::zero)));
            }
        }
    }
}


template<class CloudType>
void Foam::ParticleVelocityField<CloudType>::postEvolve(
    const typename parcelType::trackingData& td)
{
    if (counter_ == 0)
    {
        Info << "Updating UParticle field" << endl;
        volVectorField& UParticle = UParticlePtr_();
        volScalarField& nParticle = nParticlePtr_();
        volScalarField& nParcel = nParcelPtr_();


        const auto& c = this->owner();

        label parceli = 0;
        forAllConstIters(c, parcelIter)
        {
            const parcelType& p = parcelIter();
            if ((dmin_ < p.d()) && (p.d() < dmax_))
            {
                UParticle[p.cell()] += p.nParticle() * p.U();
                nParticle[p.cell()] += p.nParticle();
                nParcel[p.cell()]++;
            }
        }

        UParticle.primitiveFieldRef(false) /= nParticle.primitiveFieldRef(false) + 1e-10;

        if (average_)
        {
            volVectorField& UParticleInt = UParticleIntPtr_();
            volScalarField& nParticleInt = nParticleIntPtr_();
            volScalarField& nParcelInt = nParcelIntPtr_();

            volVectorField& UParticleMean = UParticleMeanPtr_();
            UParticleInt.primitiveFieldRef(false) += UParticle.primitiveFieldRef(false) * nParticle.primitiveFieldRef(false);
            nParticleInt.primitiveFieldRef(false) += nParticle.primitiveFieldRef(false);
            nParcelInt.primitiveFieldRef(false) += nParcel.primitiveFieldRef(false);

            UParticleMean.primitiveFieldRef(false) = UParticleInt.primitiveFieldRef(false) / (nParticleInt.primitiveFieldRef(false) + 1e-10);
        }
    }
    counter_++;
}

// ************************************************************************* //
