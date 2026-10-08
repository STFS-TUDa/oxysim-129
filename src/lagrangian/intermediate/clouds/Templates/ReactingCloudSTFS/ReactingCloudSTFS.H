/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2016 OpenCFD Ltd.
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
    Foam::ReactingCloudSTFS

Group
    grpLagrangianIntermediateClouds

Description
    Templated base class for reacting cloud

    - Adds to thermodynamic cloud
      - Variable composition (single phase)
      - Phase change

SourceFiles
    ReactingCloudSTFSI.H
    ReactingCloudSTFS.C

\*---------------------------------------------------------------------------*/

#pragma once

#include "reactingCloud.H"

#include "ReactingCloudSTFS.H"
#include "CompositionModel.H"
#include "PhaseChangeModel.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of classes

template<class CloudType>
class CompositionModel;

template<class CloudType>
class PhaseChangeModel;

/*---------------------------------------------------------------------------*\
                      Class ReactingCloudSTFS Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class ReactingCloudSTFS
: public CloudType,
  public reactingCloud
{
public:
    // Public typedefs

    //- Type of cloud this cloud was instantiated for
    typedef CloudType cloudType;

    //- Type of parcel the cloud was instantiated for
    typedef typename CloudType::particleType parcelType;

    //- Convenience typedef for this cloud type
    typedef ReactingCloudSTFS<CloudType> reactingCloudType;


private:
    // Private data

    //- Cloud copy pointer
    autoPtr<ReactingCloudSTFS<CloudType>> cloudCopyPtr_;


    // Private member functions

    //- No copy construct
    ReactingCloudSTFS(const ReactingCloudSTFS&) = delete;

    //- No copy assignment
    void operator=(const ReactingCloudSTFS&) = delete;


protected:
    // Protected data

    //- Parcel constant properties
    typename parcelType::constantProperties constProps_;

    // References to the cloud sub-models

    //- Reacting composition model
    autoPtr<CompositionModel<ReactingCloudSTFS<CloudType>>>
        compositionModel_;

    //- Reacting phase change model
    autoPtr<PhaseChangeModel<ReactingCloudSTFS<CloudType>>>
        phaseChangeModel_;


    // Sources

    //- Mass transfer fields - one per carrier phase specie
    PtrList<volScalarField::Internal> rhoTrans_;

    //- Absolute enthalpy transfer [J/kg]
    autoPtr<volScalarField::Internal> haTrans_;

    //- Coefficient for carrier phase ha equation [W/K]
    autoPtr<volScalarField::Internal> haCoeff_;

    //- Kinetic Energy transfer [J/kg]
    autoPtr<volScalarField::Internal> KTrans_;


    // Protected Member Functions

    // New parcel helper functions

    //- Check that size of a composition field is valid
    void checkSuppliedComposition(
        const scalarField& YSupplied,
        const scalarField& Y,
        const word& YName);


    // Initialisation

    //- Set cloud sub-models
    void setModels();


    // Cloud evolution functions

    //- Reset state of cloud
    void cloudReset(ReactingCloudSTFS<CloudType>& c);


public:
    // Constructors

    //- Construct given carrier gas fields
    ReactingCloudSTFS(
        const word& cloudName,
        const volScalarField& rho,
        const volVectorField& U,
        const dimensionedVector& g,
        const SLGThermo& thermo,
        bool readFields = true);

    //- Copy constructor with new name
    ReactingCloudSTFS(ReactingCloudSTFS<CloudType>& c, const word& name);

    //- Copy constructor with new name - creates bare cloud
    ReactingCloudSTFS(
        const fvMesh& mesh,
        const word& name,
        const ReactingCloudSTFS<CloudType>& c);

    //- Construct and return clone based on (this) with new name
    virtual autoPtr<Cloud<parcelType>> clone(const word& name)
    {
        return autoPtr<Cloud<parcelType>>(
            new ReactingCloudSTFS(*this, name));
    }

    //- Construct and return bare clone based on (this) with new name
    virtual autoPtr<Cloud<parcelType>> cloneBare(const word& name) const
    {
        return autoPtr<Cloud<parcelType>>(
            new ReactingCloudSTFS(this->mesh(), name, *this));
    }


    //- Destructor
    virtual ~ReactingCloudSTFS() = default;


    // Member Functions

    // Access

    //- Return a reference to the cloud copy
    inline const ReactingCloudSTFS& cloudCopy() const;

    //- Return the constant properties
    inline const typename parcelType::constantProperties&
    constProps() const;

    //- Return access to the constant properties
    inline typename parcelType::constantProperties& constProps();


    // Sub-models

    //- Return const access to reacting composition model
    inline const CompositionModel<ReactingCloudSTFS<CloudType>>&
    composition() const;

    //- Return const access to reacting phase change model
    inline const PhaseChangeModel<ReactingCloudSTFS<CloudType>>&
    phaseChange() const;

    //- Return reference to reacting phase change model
    inline PhaseChangeModel<ReactingCloudSTFS<CloudType>>&
    phaseChange();


    // Sources

    //- Mass

    //- Return reference to mass source for field i
    inline volScalarField::Internal&
    rhoTrans(const label i);

    //- Return const access to mass source fields
    inline const PtrList<volScalarField::Internal>&
    rhoTrans() const;

    //- Return reference to mass source fields
    inline PtrList<volScalarField::Internal>&
    rhoTrans();

    //- Return mass source term for specie i - specie eqn
    inline tmp<fvScalarMatrix> SYi(
        const label i,
        volScalarField& Yi) const;

    //- Return tmp mass source for field i - fully explicit
    inline tmp<volScalarField::Internal>
    Srho(const label i) const;

    //- Return tmp total mass source for carrier phase
    //  - fully explicit
    inline tmp<volScalarField::Internal> Srho() const;

    //- Return total mass source term [kg/m3/s]
    inline tmp<fvScalarMatrix> Srho(volScalarField& rho) const;

    //- Absolute Enthalpy

    //- Creation enthalpy transfer [J/kg]
    inline volScalarField::Internal& haTrans();

    //- Creation enthalpy transfer [J/kg]
    inline const volScalarField::Internal&
    haTrans() const;

    //- Return coefficient for carrier phase ha equation
    inline volScalarField::Internal& haCoeff();

    //- Return const coefficient for carrier phase ha equation
    inline const volScalarField::Internal&
    haCoeff() const;

    //- Return absolute enthalpy source term [J/kg/m3/s]
    inline tmp<fvScalarMatrix> Sha(volScalarField& ha) const;

    //- Kinetic Energy

    //- Kinetic Energy transfer [J/kg]
    inline volScalarField::Internal& KTrans();

    //- Kinetic Energy transfer [J/kg]
    inline const volScalarField::Internal&
    KTrans() const;

    //- Return creation enthalpy source term [J/kg/m3/s]
    inline tmp<fvScalarMatrix> SK(volScalarField& K) const;

    // Check

    //- Total absolute enthalpy of system
    inline scalar absoluteEnthalpyOfSystem() const;

    //- Max diameter
    inline scalar Dmin() const;

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

    //- Apply relaxation to (steady state) cloud sources
    void relaxSources(const ReactingCloudSTFS<CloudType>& cloudOldTime);

    //- Apply scaling to (transient) cloud sources
    void scaleSources();

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

    //- Write particle fields as objects into the obr registry
    virtual void writeObjects(objectRegistry& obr) const;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ReactingCloudSTFSI.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
#include "ReactingCloudSTFS.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// ************************************************************************* //
