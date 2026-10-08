/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
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
    Foam::FlameletThermo

Group
    grpThermo

Description
    Class derived from Thermo that implements an Interface using tabulated chemistry
    with the flameletConfig (hdf5-reader).

Authors
    ...

SourceFiles
    FlameletThermo.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <string>
#include <vector>

#include "fvMesh.H"
#include "IOdictionary.H"
#include "volFields.H"
#include "dimensionedScalar.H"
#include "hashedWordList.H"
#include "HashPtrTable.H"

#include "runTimeSelectionTables.H"
#include "basicThermoSTFS.H"

// Include the newFlameletConfig
#include "lookup_table.hpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                Structs to interact with the flameletConfig.
\*---------------------------------------------------------------------------*/
/** @defgroup flameletConfig
 *   @ingroup flamelet
 *  @{
 */
class Flut_object
{

protected:
    std::string name_ = "none";
    Foam::volScalarField& field_;

public:
    Flut_object(
        std::string name,
        Foam::volScalarField& f)
    : name_(name), field_(f) {};

    // Function to check if the object has a specific name
    bool operator==(const std::string& targetName) const
    {
        return name_ == targetName;
    }

    const std::string getName() const
    {
        return name_;
    }

    const Foam::volScalarField& getField() const { return field_; }
    Foam::volScalarField& getField() { return field_; }
};

template<typename T>
bool contains(const std::vector<T>& vec, const std::string& value)
{
    return std::find(vec.begin(), vec.end(), value) != vec.end();
}

template<typename T>
bool contains(const Foam::List<T>& vec, const std::string& value)
{
    return std::find(vec.begin(), vec.end(), value) != vec.end();
}
/**@} */
/*---------------------------------------------------------------------------*\
                           Class flameletThermo Declaration
\*---------------------------------------------------------------------------*/
/** @addtogroup flamelet
 *  @{
 */
class flameletThermo : public basicThermoSTFS
{

protected:
    // Protected data
    const fvMesh& mesh_;

    // Custom flameletThermo properties
    const IOdictionary thermoDict_;

    // Custom functionalities of flameletThermo
    dimensionedScalar pRef_;

    std::shared_ptr<FLUT::LookupTable> Flut_ {nullptr};
    std::shared_ptr<FLUT::LookupTable> FlutBoundary_ {nullptr};

    // Collection of fields that are either updated by thermo or required by the equationHandler.
    // Includes all fields owned by the thermo-class
    // Includes additional field that are created outside the thermo-class
    //      e.g. created by canteraThermo or inside the solver.
    std::map<Foam::word, std::reference_wrapper<volScalarField>> fields_;
    std::vector<Foam::word> ownedFieldNames_; // Names of the fields owned by the flameletThermo in fields_

    std::vector<Flut_object> inputVariables_;       // Control variables (used for table access)
    std::vector<double> inputVariable_values_;      // Serves as input for the table lookup.
    std::vector<Flut_object> outputVariables_;      // output variables (looked up quantities)
    std::vector<std::string> outputVariable_names_; // names of output variables

    std::vector<Flut_object> inputVariablesBoundary_;  // Control variables (used for boundary table access)
    std::vector<double> inputVariableBoundary_values_; // Serves as input for the table lookup.
    std::vector<Flut_object> outputVariablesBoundary_; // output variables (look up quantities)

    // Treatment for variablesToSolve
    hashedWordList variablesToSolve_;

    mutable std::string boundaryLookupType_ = "none";
    mutable std::vector<label> boundaryLookupPatchIDs_;

    const bool useWithCanteraThermo_;

    // Member functions (private)
    void init();

    // Standard initialization functions.
    void init_fc();
    void init_inputVariables();
    void init_variablesToSolve();
    void init_outputVariables();

    // Initialization functions for nonAdiabatic boundary lookup for enthalpy tables.
    void init_fcBoundary();
    void init_inputVariablesBoundary();
    void init_outputVariablesBoundary();

    // Update functions.
    void saveOldTimes();
    void updateFields();
    void updateBoundaries();
    void updateBoundary(label patchID);


    /**
     * @brief Register a field as outputVariable.
     *
     * This function registers a field as outputVariable, i.e.
     * the field is interpolated based on the FLUT during thermo.correct()
     * Variables that cannot be registered as outputVariables are: [variablesToSolve, transportedSpecies, inputVariables]
     *
     * @param field_name The name of the field.
     */
    void registerOutputField(word field_name);

    // Functions for boundaryLookup initialization
    void setBoundaryLookupType();
    void setBoundary(std::string type, label patchID) const;
    void boundaryLookUp();

    // Protected Member Functions

public:
    //- Runtime type information
    TypeName("flameletThermo");

    //- Declare run-time constructor selection table
    declareRunTimeSelectionTable(
        autoPtr,
        flameletThermo,
        fvMesh,
        (const fvMesh& mesh, const word& phaseName),
        (mesh, phaseName));

    //- Declare run-time constructor selection table for dictionary based
    declareRunTimeSelectionTable(
        autoPtr,
        flameletThermo,
        fvMeshDictPhase,
        (const fvMesh& mesh, const word& phaseName, const word& dictName),
        (mesh, phaseName, dictName));

    // Constructors

    //- Construct from mesh,dictionary,phase name with a single temperature
    flameletThermo(

        const fvMesh&,
        const word& phaseName = "",
        const word& dictionaryName = "flameletProperties");

    //- Construct as copy
    flameletThermo(const flameletThermo&) = delete;

    //- Selector
    static autoPtr<flameletThermo> New(
        const fvMesh& mesh,
        const word& phaseName = word::null);

    //- Selector
    static autoPtr<flameletThermo> New(
        const fvMesh& mesh,
        const word& phaseName = word::null,
        const word& dictionaryName = word::null);

    ~flameletThermo();

    // Member Functions

    /**
     * @brief Append a output field to the thermo-class.
     *
     * This function appends a field to the thermo-class (see appendField).
     * Then the field is registered as outputVariable (see registerField).
     * outputVariables are looked-up from the FLUT.
     *
     * @param field_name The name of the field.
     * @param readOpt    The read option for IOobject (default: NO_READ).
     * @param writeOpt   The write option for IOobject (default: AUTO_WRITE).
     * @param dim        The dimensionedScalar representing the field's dimensions (default: dimensionedScalar("zero", dimless, 0.0)).
     */
    void appendOutputField(word field_name, IOobject::readOption readOpt = IOobject::NO_READ, IOobject::writeOption writeOpt = IOobject::AUTO_WRITE, dimensionedScalar dim = dimensionedScalar("zero", dimless, 0.0));

    /**
     * @brief Append a field to the thermo-class.
     *
     * This function appends a field to the thermo-class. If the field is not found in the
     * object registry, it constructs the field and adds it to the thermo-class's fields.
     * If the field is found in the object registry but not in the thermo-class's fields,
     * it registers the field in the thermo-class's fields.
     *
     * @param field_name The name of the field.
     * @param readOpt    The read option for IOobject (default: NO_READ).
     * @param writeOpt   The write option for IOobject (default: AUTO_WRITE).
     * @param dim        The dimensionedScalar representing the field's dimensions (default: dimensionedScalar("zero", dimless, 0.0)).
     */
    void appendField(word field_name, IOobject::readOption readOpt = IOobject::NO_READ, IOobject::writeOption writeOpt = IOobject::AUTO_WRITE, dimensionedScalar dim = dimensionedScalar("zero", dimless, 0.0));

    // Function to load the outputVariable data into the FLUTs.
    // This function should only be called onced in the solver after
    // the initialization is completed.
    void loadOutputVariablesFromFLUT();

    // Update function
    void update();

    // Return access
    // Custom functionalities
    // ToDo: Add to unit test
    [[nodiscard]] const fvMesh& getMesh() const;

    // Return access
    // ToDo: Add to unit test
    const dictionary& dict() const;
    // dictionary& dict();

    // Custom getter for aribary fields
    // Get an outputField of the thermoClass.
    const volScalarField& getField(word name) const;
    volScalarField& getField(word name);

    // OpenFOAM defaults

    //- Laminar thermal diffusivity [kg/m/s]
    //  Non-const access allowed for overwritten of field by canteraThermo
    volScalarField& alpha();

    // Necessary for fluidThermo interface
    //- Density [kg/m^3]
    virtual tmp<volScalarField> rho() const;
    //- Density for patch [kg/m^3]
    virtual tmp<scalarField> rho(const label patchi) const;

    //- Return reference to the density field [kg/m^3]
    virtual volScalarField& rhoRef();

    //- Add the given density correction to the density field.
    //  Used to update the density field following pressure solution
    //  Limit thermo rho between rhoMin and rhoMax
    // ToDo: Add to unit test
    virtual void correctRho(
        const volScalarField& deltaRho,
        const dimensionedScalar& rhoMin,
        const dimensionedScalar& rhoMax);

    //- Add the given density correction to the density field.
    //  Used to update the density field following pressure solution
    // ToDo: Add to unit test
    virtual void correctRho(
        const volScalarField& deltaRho);

    //- Compressibility [s^2/m^2]
    //  pure virtual function from fluidThermo
    virtual const volScalarField& psi() const;

    // Transport properties inherited from basicThermo
    //- Dynamic viscosity of mixture [kg/m/s]
    tmp<volScalarField> mu() const override;
    tmp<scalarField> mu(const label patchi) const override;

    //- Heat capacity at constant pressure/volume [J/kg/K] - Calculated on the fly
    tmp<volScalarField> Cpv() const override;

    //- Molecular weight [kg/kmol]. Requires "W" to be listed in
    //  outputVariablesFlameletTable.
    virtual tmp<volScalarField> W() const override;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**@} */
} // End namespace Foam

// ************************************************************************* //
