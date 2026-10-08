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
    Foam::CarbonaceousCloud

Group
    grpLagrangianIntermediateClouds

Description
    Templated base class for multiphase reacting cloud

    - Adds to reacting cloud
      - multiphase composition
      - devolatilisatsion
      - surface reactions

SourceFiles
    CarbonaceousCloudI.H
    CarbonaceousCloud.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "Cloud.H"
#include "SLGThermo.H"
#include "carbonaceousCloud.H"
#include "DevolatilisationModel.H"
#include "DevolatilisationModelSTFS.H"
#include "SurfaceReactionModel.H"
#include "HeatCapacityModelSTFS.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

// template<class CloudType>
// class DevolatilisationModel;

template<class CloudType>
class DevolatilisationModelSTFS;

// template<class CloudType>
// class SurfaceReactionModel;

template<class CloudType>
class HeatCapacityModelSTFS;

/*---------------------------------------------------------------------------*\
                   Class CarbonaceousCloud Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CarbonaceousCloud
: public CloudType,
  public carbonaceousCloud
{
public:
    // Public typedefs

    //- Type of cloud this cloud was instantiated for
    typedef CloudType cloudType;

    //- Type of parcel the cloud was instantiated for
    typedef typename CloudType::particleType parcelType;

    //- Convenience typedef for this cloud type
    typedef CarbonaceousCloud<CloudType> carbonaceousCloudType;


private:
    // Private data

    //- Cloud copy pointer
    autoPtr<CarbonaceousCloud<CloudType>> cloudCopyPtr_;


    // Private member functions

    //- No copy construct
    CarbonaceousCloud(const CarbonaceousCloud&) = delete;

    //- No copy assignment
    void operator=(const CarbonaceousCloud&) = delete;


protected:
    // Protected data

    //- Parcel constant properties
    typename parcelType::constantProperties constProps_;


    // References to the cloud sub-models

    //- Devolatilisation model
    //    autoPtr<
    //        DevolatilisationModel<CarbonaceousCloud<CloudType>>>
    //        devolatilisationModel_;

    autoPtr<
        DevolatilisationModelSTFS<CarbonaceousCloud<CloudType>>>
        devolatilisationModelSTFS_;

    //- Surface reaction model
    //    autoPtr<
    //        SurfaceReactionModel<CarbonaceousCloud<CloudType>>>
    //        surfaceReactionModel_;

    //- Heat capacity model
    autoPtr<
        HeatCapacityModelSTFS<CarbonaceousCloud<CloudType>>>
        heatCapacityModelSTFS_;

    // Check

    //- Total mass transferred to continuous phase via devolatilisation
    scalar dMassDevolatilisation_;

    //- Total mass transferred to continuous phase via surface
    //  reactions
    scalar dMassSurfaceReaction_;


    // Protected Member Functions

    // Initialisation

    //- Set cloud sub-models
    void setModels();


    // Cloud evolution functions

    //- Reset state of cloud
    void cloudReset(CarbonaceousCloud<CloudType>& c);


public:
    // Constructors

    //- Construct given carrier gas fields
    CarbonaceousCloud(
        const word& cloudName,
        const volScalarField& rho,
        const volVectorField& U,
        const dimensionedVector& g,
        const SLGThermo& thermo,
        bool readFields = true);


    //- Copy constructor with new name
    CarbonaceousCloud(
        CarbonaceousCloud<CloudType>& c,
        const word& name);

    //- Copy constructor with new name - creates bare cloud
    CarbonaceousCloud(
        const fvMesh& mesh,
        const word& name,
        const CarbonaceousCloud<CloudType>& c);

    //- Construct and return clone based on (this) with new name
    virtual autoPtr<Cloud<parcelType>> clone(const word& name)
    {
        return autoPtr<Cloud<parcelType>>(
            new CarbonaceousCloud(*this, name));
    }

    //- Construct and return bare clone based on (this) with new name
    virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
    {
        return autoPtr<Cloud<parcelType>>(
            new CarbonaceousCloud(this->mesh(), name, *this));
    }


    //- Destructor
    virtual ~CarbonaceousCloud();


    // Member Functions

    // Access

    //- Return a reference to the cloud copy
    inline const CarbonaceousCloud& cloudCopy() const;

    //- Return the constant properties
    inline const typename parcelType::constantProperties&
    constProps() const;

    //- Return access to the constant properties
    inline typename parcelType::constantProperties& constProps();


    // Sub-models

    //- Return const access to devolatilisation model
    //    inline const DevolatilisationModel<
    //        CarbonaceousCloud<CloudType>>&
    //    devolatilisation() const;

    //- Return reference to devolatilisation model
    //    inline DevolatilisationModel<
    //        CarbonaceousCloud<CloudType>>&
    //    devolatilisation();


    //- Return const access to devolatilisation model
    inline const DevolatilisationModelSTFS<
        CarbonaceousCloud<CloudType>>&
    devolatilisationSTFS() const;

    //- Return reference to devolatilisation model
    inline DevolatilisationModelSTFS<
        CarbonaceousCloud<CloudType>>&
    devolatilisationSTFS();

    //- Return const access to reacting surface reaction model
    //    inline const SurfaceReactionModel<
    //        CarbonaceousCloud<CloudType>>&
    //    surfaceReaction() const;

    //- Return reference to reacting surface reaction model
    //    inline SurfaceReactionModel<
    //        CarbonaceousCloud<CloudType>>&
    //    surfaceReaction();

    inline const HeatCapacityModelSTFS<
        CarbonaceousCloud<CloudType>>&
    heatCapacitySTFS() const;

    //- Return reference to heat capacity model
    inline HeatCapacityModelSTFS<
        CarbonaceousCloud<CloudType>>&
    heatCapacitySTFS();


    // Cloud evolution functions

    //- Set parcel thermo properties
    void setParcelThermoProperties(
        parcelType& parcel,
        const scalar lagrangianDt);

    //- Check parcel properties
    void checkParcelProperties(
        parcelType& parcel,
        const scalar lagrangianDt,
        const bool fullyDescribed);

    //- Store the current cloud state
    void storeState();

    //- Reset the current cloud to the previously stored state
    void restoreState();

    //- Reset the cloud source terms
    void resetSourceTerms();

    //- Evolve the cloud
    void evolve();


    // Mapping

    //- Remap the cells of particles corresponding to the
    //  mesh topology change with a default tracking data object
    virtual void autoMap(const mapPolyMesh&);


    // I-O

    //- Print cloud information
    void info();

    //- Write the field data for the cloud
    virtual void writeFields() const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "CarbonaceousCloudI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "CarbonaceousCloud.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
