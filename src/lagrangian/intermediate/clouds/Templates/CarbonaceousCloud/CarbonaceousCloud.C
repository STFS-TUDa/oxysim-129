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

#include "CarbonaceousCloud.H"

// * * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * //

template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::setModels()
{
    //    devolatilisationModel_.reset(
    //        DevolatilisationModel<CarbonaceousCloud<CloudType>>::New(
    //            this->subModelProperties(),
    //            *this)
    //            .ptr());

    devolatilisationModelSTFS_.reset(
        DevolatilisationModelSTFS<CarbonaceousCloud<CloudType>>::New(
            this->subModelProperties(),
            *this)
            .ptr());

    //    surfaceReactionModel_.reset(
    //        SurfaceReactionModel<CarbonaceousCloud<CloudType>>::New(
    //            this->subModelProperties(),
    //            *this)
    //            .ptr());
    heatCapacityModelSTFS_.reset(
        HeatCapacityModelSTFS<CarbonaceousCloud<CloudType>>::New(
            this->subModelProperties(),
            *this)
            .ptr());
}

template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::cloudReset(
    CarbonaceousCloud<CloudType>& c)
{
    CloudType::cloudReset(c);

    //    devolatilisationModel_.reset(c.devolatilisationModel_.ptr());
    devolatilisationModelSTFS_.reset(c.devolatilisationModelSTFS_.ptr());
    //    surfaceReactionModel_.reset(c.surfaceReactionModel_.ptr());
    heatCapacityModelSTFS_.reset(c.heatCapacityModelSTFS_.ptr());

    //    dMassDevolatilisation_ = c.dMassDevolatilisation_;
    //    dMassSurfaceReaction_ = c.dMassSurfaceReaction_;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CarbonaceousCloud<CloudType>::CarbonaceousCloud(
    const word& cloudName,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& thermo,
    bool readFields)
: CloudType(cloudName, rho, U, g, thermo, false),
  carbonaceousCloud(),
  cloudCopyPtr_(nullptr),
  constProps_(this->particleProperties()),
  // devolatilisationModel_(nullptr),
  devolatilisationModelSTFS_(nullptr),
  //  surfaceReactionModel_(nullptr),
  heatCapacityModelSTFS_(nullptr)
//  dMassDevolatilisation_(0.0),
//  dMassSurfaceReaction_(0.0)
{
    if (this->solution().active())
    {
        setModels();

        if (readFields)
        {
            parcelType::readFields(*this, this->composition());
            this->deleteLostParticles();
        }
    }

    if (this->solution().resetSourcesOnStartup())
    {
        resetSourceTerms();
    }
}


template<class CloudType>
Foam::CarbonaceousCloud<CloudType>::CarbonaceousCloud(
    CarbonaceousCloud<CloudType>& c,
    const word& name)
: CloudType(c, name),
  carbonaceousCloud(),
  cloudCopyPtr_(nullptr),
  constProps_(c.constProps_),
  //  devolatilisationModel_(c.devolatilisationModel_->clone()),
  devolatilisationModelSTFS_(c.devolatilisationModelSTFS_->clone()),
  //  surfaceReactionModel_(c.surfaceReactionModel_->clone()),
  heatCapacityModelSTFS_(c.heatCapacityModelSTFS_->clone())
//  dMassDevolatilisation_(c.dMassDevolatilisation_),
//  dMassSurfaceReaction_(c.dMassSurfaceReaction_)
{
}


template<class CloudType>
Foam::CarbonaceousCloud<CloudType>::CarbonaceousCloud(
    const fvMesh& mesh,
    const word& name,
    const CarbonaceousCloud<CloudType>& c)
: CloudType(mesh, name, c),
  carbonaceousCloud(),
  cloudCopyPtr_(nullptr),
  constProps_(),
  //  devolatilisationModel_(nullptr),
  devolatilisationModelSTFS_(nullptr),
  //  surfaceReactionModel_(nullptr),
  heatCapacityModelSTFS_(nullptr)
//  dMassDevolatilisation_(0.0),
//  dMassSurfaceReaction_(0.0)
{
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CarbonaceousCloud<CloudType>::~CarbonaceousCloud()
{
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::setParcelThermoProperties(
    parcelType& parcel,
    const scalar lagrangianDt)
{
    CloudType::setParcelThermoProperties(parcel, lagrangianDt);
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::checkParcelProperties(
    parcelType& parcel,
    const scalar lagrangianDt,
    const bool fullyDescribed)
{
    CloudType::checkParcelProperties(parcel, lagrangianDt, fullyDescribed);
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::storeState()
{
    cloudCopyPtr_.reset(
        static_cast<CarbonaceousCloud<CloudType>*>(
            clone(this->name() + "Copy").ptr()));
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::restoreState()
{
    cloudReset(cloudCopyPtr_());
    cloudCopyPtr_.clear();
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::resetSourceTerms()
{
    CloudType::resetSourceTerms();
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::evolve()
{
    if (this->solution().canEvolve())
    {
        typename parcelType::trackingData td(*this);

        this->solve(*this, td);
    }
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::autoMap(
    const mapPolyMesh& mapper)
{
    Cloud<parcelType>::autoMap(mapper);

    this->updateMesh();
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::info()
{
    CloudType::info();
}


template<class CloudType>
void Foam::CarbonaceousCloud<CloudType>::writeFields() const
{
    if (this->compositionModel_)
    {
        CloudType::particleType::writeFields(*this, this->composition());
    }
}


// ************************************************************************* //
