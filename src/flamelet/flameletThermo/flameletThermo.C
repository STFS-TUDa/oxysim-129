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

\*---------------------------------------------------------------------------*/

#include <iostream>
#include <unordered_set>

#include "flameletThermo.H"
#include "boundaryLookUpFvPatchScalarField.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
defineTypeNameAndDebug(flameletThermo, 0);
defineRunTimeSelectionTable(flameletThermo, fvMesh);
defineRunTimeSelectionTable(flameletThermo, fvMeshDictPhase);
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::flameletThermo::flameletThermo(const fvMesh& mesh, const word& phaseName, const word& dictionaryName)
: basicThermoSTFS(mesh, phaseName, dictionaryName),
  mesh_(mesh),
  thermoDict_(
      IOobject(
          dictionaryName,
          mesh.time().constant(),
          mesh,
          IOobject::MUST_READ,
          IOobject::NO_WRITE)),
  // TODO: Check the default dict for pRef_ definition and use this one.
  pRef_(dimPressure, mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<scalar>("pRef", 101325)),
  useWithCanteraThermo_(lookupOrDefault("useWithCanteraThermo", false))
{
    init();
}


// * * * * * * * * * * * * * * * * Selectors * * * * * * * * * * * * * * * * //

Foam::autoPtr<Foam::flameletThermo> Foam::flameletThermo::New(
    const fvMesh& mesh,
    const word& phaseName)
{
    return basicThermo::New<flameletThermo>(mesh, phaseName);
}

Foam::autoPtr<Foam::flameletThermo> Foam::flameletThermo::New(
    const fvMesh& mesh,
    const word& phaseName,
    const word& dictionaryName)
{
    return basicThermo::New<flameletThermo>(mesh, phaseName, dictionaryName);
}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::flameletThermo::~flameletThermo()
{
    // checkout fields owned by flameletThermo
    for (const auto& fieldName : ownedFieldNames_)
    {
        mesh_.thisDb().checkOut(fieldName);
    }
}

// * * * * * * * * * * * * * * * Member Functions for initialization  * * * * * * * * * * * * * //

void Foam::flameletThermo::init()
{
    Info << "flameletThermo::flameletThermo() called." << endl;

    // * * * Initialize flamelet-config * * * //
    init_fc();
    // * * * Read in variables to solve (these cannot be used as outputVariables) * * * //
    variablesToSolve_ = lookup("variablesToSolve");

    // * * * Initialize inputVariable Fields * * * //
    init_inputVariables();

    // * * * Initialize outputVariable Fields * * * //
    init_outputVariables();

    setBoundaryLookupType();
    if (boundaryLookupType_ != "none")
    {
        init_fcBoundary();
        init_inputVariablesBoundary();
        init_outputVariablesBoundary();
    }

    Info << "variablesToSolve_: " << endl;
    Info << variablesToSolve_ << endl;
}

void Foam::flameletThermo::init_fc()
{
    Info << endl
         << "Create FLUT:" << endl;

    std::string flutFile = lookupOrDefault<fileName>("flutFile", "FLUT.h5");

    // Initialize the FLUT object
    Flut_ = std::make_unique<FLUT::LookupTable>(flutFile);
    // Print information of the FLUT
    if (Pstream::master())
    {
        Flut_->printInfo();
    }
}

void Foam::flameletThermo::init_inputVariables()
{
    std::vector<std::string> inputVariable_names = Flut_->getInputVariables();

    // Check if inputVariables from flameletProperties match the inputVariables from the Flut
    hashedWordList inputVariable_names_flameletProperties = get<hashedWordList>("inputVariablesFlameletTable");


    wordList wl(inputVariable_names.size());
    for (label i = 0; i < label(inputVariable_names.size()); ++i)
    {
        wl[i] = word(inputVariable_names[i]);
    }
    hashedWordList inputVariable_names_Flut(std::move(wl));

    if (inputVariable_names_flameletProperties != inputVariable_names_Flut)
    {
        FatalErrorInFunction
            << "inputVariables from flameletProperties don't match inputVariables from FLUT!\n"
            << "  flameletProperties: " << inputVariable_names_flameletProperties << "\n"
            << "  FLUT: " << inputVariable_names_Flut << nl
            << abort(FatalError);
    }

    inputVariable_values_.resize(inputVariable_names.size());

    Info << "Create fields for inputVariables: " << endl;
    for (auto field_name : inputVariable_names)
    {
        Info << "Register " << field_name << " as inputVariable for the flut." << endl;

        appendField(field_name, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
        // Add the field to inputVariable_fields
        inputVariables_.push_back(Flut_object(field_name, getField(field_name)));
    }
}

void Foam::flameletThermo::init_outputVariables()
{
    // Creation of default outputVariables
    // Initialize field that are part of the fluidThermo class.
    Info << "Create default outputVariables: " << endl;

    if (useWithCanteraThermo_)
    {
        appendField("thermo:rho", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimDensity, 0.0));
    }
    else
    {
        // these fields (T_, alpha_) are already part of the basicThermo --> no appendField needed!
        outputVariables_.push_back(Flut_object("T", T_));
        outputVariable_names_.push_back("T");
        outputVariables_.push_back(Flut_object("alpha", alpha_));
        outputVariable_names_.push_back("alpha");
        // Register additional required fields.
        appendOutputField("thermo:rho", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimDensity, 0.0));
        appendOutputField("mu", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimMass / dimLength / dimTime, 0.0));
    }
    appendField("thermo:psi", IOobject::NO_READ, IOobject::NO_WRITE, dimensionedScalar("zero", dimensionSet(0, -2, 2, 0, 0), 0.0));

    // Creation of additional variables requested by the user.
    Info << "Create fields for user-defined outputVariablesFlameletTable: " << endl;
    hashedWordList user_fields(lookup("outputVariablesFlameletTable"), true);
    for (const auto& field_name : user_fields)
    {
        // treatment of inputVariables and variablesToSolve in registerOutputField
        appendOutputField(field_name, IOobject::NO_READ, IOobject::AUTO_WRITE);
    }
}

void Foam::flameletThermo::loadOutputVariablesFromFLUT()
{
    Flut_->readOutputVariableData(outputVariable_names_);
}


void Foam::flameletThermo::init_fcBoundary()
{
    Info << endl
         << "Create Boundary FLUT:" << endl;

    std::string flutFile = lookupOrDefault<fileName>("flutFileBC", "FLUT_BC.h5");

    // Initialize the FLUT object
    FlutBoundary_ = std::make_unique<FLUT::LookupTable>(flutFile);
    // Print information of the FLUT
    if (Pstream::master())
    {
        FlutBoundary_->printInfo();
    }
};

void Foam::flameletThermo::init_inputVariablesBoundary()
{
    std::vector<std::string> inputVariable_names = FlutBoundary_->getInputVariables();
    inputVariableBoundary_values_.resize(inputVariable_names.size());

    Info << "Create fields for inputVariables: " << endl;
    for (auto field_name : inputVariable_names)
    {
        Info << "Register " << field_name << " as inputVariable for the boundary flut." << endl;

        if (field_name == "ccBoundary")
        {
            appendField(field_name, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE);
        }
        else
        {
            appendField(field_name, IOobject::MUST_READ, IOobject::AUTO_WRITE);
        }
        // Add the field to inputVariable_fields
        inputVariablesBoundary_.push_back(Flut_object(field_name, getField(field_name)));
    }
}

void Foam::flameletThermo::init_outputVariablesBoundary()
{
    // output-field from the Flut_ (not FlutBoundary_).
    if (boundaryLookupType_ == "specialBoundaryLookup")
    {
        appendOutputField("yc_max", IOobject::NO_READ, IOobject::NO_WRITE);
    }

    outputVariablesBoundary_.push_back(Flut_object("ha", getField("ha")));
    std::vector<std::string> outputVariableBoundary_names;
    outputVariableBoundary_names.push_back("ha");
    FlutBoundary_->readOutputVariableData(outputVariableBoundary_names);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::flameletThermo::saveOldTimes()
{
    // Write outputVariables to field.
    for (auto& outputVariable : outputVariables_)
    {
        outputVariable.getField().storeOldTimes();
    }
}

void Foam::flameletThermo::update()
{
    // Check if the outputData has already been initialized in the solver
    // If not, initialize the output variables.
    if (!Flut_->hasInitializedOutputData())
    {
        loadOutputVariablesFromFLUT();
    }

    saveOldTimes();
    updateFields();
    updateBoundaries();
    if (boundaryLookupType_ != "none")
    {
        boundaryLookUp();
    }

    // calculate psi
    if (Foam::contains(inputVariables_, "p"))
    {
        // If pressure is an inputVariable use the current pressure
        getField("thermo:psi") = getField("thermo:rho") / p_;
    }
    else
    {
        // If pressure is not an inputVariable use pRef
        // i.e. pressure at which the table was created
        getField("thermo:psi") = getField("thermo:rho") / pRef_;
    }
}

void Foam::flameletThermo::updateFields()
{
    forAll(mesh_.cells(), cellI)
    {
        // Set input variable values
        for (size_t i = 0; i < inputVariables_.size(); ++i)
        {
            inputVariable_values_[i] = inputVariables_[i].getField()[cellI];
        }

        const std::vector<double>& outputVariable_values = Flut_->lookup(inputVariable_values_);

        // Write outputVariables to field.
        for (size_t i = 0; i < outputVariables_.size(); ++i)
        {
            outputVariables_[i].getField()[cellI] = outputVariable_values[i];
        }
    }

    // Update event number
    for (auto& outputVariable : outputVariables_)
    {
        outputVariable.getField().setUpToDate();
    }
}

void Foam::flameletThermo::updateBoundaries()
{
    forAll(mesh_.boundary(), patchID)
    {
        updateBoundary(patchID);
    }
}

void Foam::flameletThermo::updateBoundary(label patchID)
{
    forAll(mesh_.boundary()[patchID], facei)
    {
        // Set input Variables to field.
        for (size_t i = 0; i < inputVariables_.size(); ++i)
        {
            inputVariable_values_[i] = inputVariables_[i].getField().boundaryFieldRef(false)[patchID][facei];
        }

        const std::vector<double>& outputVariable_values = Flut_->lookup(inputVariable_values_);

        // Write outputVariables to field.
        for (size_t i = 0; i < outputVariables_.size(); ++i)
        {
            outputVariables_[i].getField().boundaryFieldRef(false)[patchID][facei] = outputVariable_values[i];
        }
    }
}

//
// Non-adiabatic lookup for enthalpy tables.
//
void Foam::flameletThermo::boundaryLookUp()
{
    if (Foam::contains(inputVariablesBoundary_, "ccBoundary") && (boundaryLookupType_ == "specialBoundaryLookup"))
    {
        // Calculate ccBoundary for the species boundary lookup.
        // This is ignored otherwise. Then, the normalization need to be performed
        // directly in the FlutBoundary_
        // ccBoundary is the normalized quantity of the boundary FLUT and needs to be
        // named accordingly for the lookup to work.
        volScalarField& ccBoundary = getField("ccBoundary");
        // Calculate ccBoundary for the table lookup
        ccBoundary = getField("yc") / getField("yc_max");
    }

    for (const auto& patchID : boundaryLookupPatchIDs_)
    {
        forAll(mesh_.boundary()[patchID], facei)
        {
            // Set input variable values
            for (size_t i = 0; i < inputVariablesBoundary_.size(); ++i)
            {
                inputVariableBoundary_values_[i] = inputVariablesBoundary_[i].getField().boundaryFieldRef(false)[patchID][facei];
            }

            std::vector<double> outputVariable_values = FlutBoundary_->lookup(inputVariableBoundary_values_);

            // Write outputVariables to field.
            for (size_t i = 0; i < outputVariablesBoundary_.size(); ++i)
            {
                outputVariablesBoundary_[i].getField().boundaryFieldRef(false)[patchID][facei] = outputVariable_values[i];
            }
        }
        updateBoundary(patchID);
    }
}

void Foam::flameletThermo::registerOutputField(word field_name)
{
    if (variablesToSolve_.found(field_name))
    {
        FatalErrorInFunction << "Field " + field_name + " cannot be used as outputVariable, because it is a variableToSolve." << exit(FatalError);
    }
    else if (Foam::contains(inputVariables_, field_name))
    {
        FatalErrorInFunction << "Field " + field_name + " cannot be used as outputVariable, because it is an inputVariable." << exit(FatalError);
    }
    else if (Foam::contains(outputVariables_, field_name))
    {
        Info << "Field " << field_name << " already registered as outputVariable. Nothing is done." << endl;
    }
    else
    {
        // Exception for fields with naming add-on in OpenFOAM.
        // These are special cases due to double field definition in the solvers.
        std::vector<std::string> THERMO_FIELDS = {"thermo:rho", "thermo:psi"};
        std::string flut_name = field_name;
        std::string prefix = "thermo:";
        if (Foam::contains(THERMO_FIELDS, field_name))
        {
            flut_name = field_name.substr(prefix.length());
        }

        Info << "Register " << field_name << " as outputVariable (field); Flut-reader name: " << flut_name << endl;
        Flut_->containsOutputVariable(flut_name);
        // List of field names with flut-reader names
        outputVariables_.push_back(Flut_object(field_name, getField(field_name)));
        outputVariable_names_.push_back(flut_name);
    }
}

void Foam::flameletThermo::appendOutputField(word field_name, IOobject::readOption readOpt, IOobject::writeOption writeOpt, dimensionedScalar dim)
{
    appendField(field_name, readOpt, writeOpt, dim);
    registerOutputField(field_name);
}

void Foam::flameletThermo::appendField(word field_name, IOobject::readOption readOpt, IOobject::writeOption writeOpt, dimensionedScalar dim)
{

    // Attempt to get the pointer to the field from the object registry
    auto* ptr = mesh_.objectRegistry::getObjectPtr<volScalarField>(field_name);

    // Field not found in objectRegistry -> field is constructed and added to fields_
    if (!ptr)
    {
        Info << "Field " << field_name << " is created and registered in fields_." << endl;

        // Create the field based on read option
        if (readOpt == IOobject::MUST_READ)
        {
            ptr = new volScalarField(
                IOobject(
                    field_name,
                    mesh_.time().timeName(),
                    mesh_,
                    readOpt,
                    writeOpt),
                mesh_);
        }
        else
        {
            ptr = new volScalarField(
                IOobject(
                    field_name,
                    mesh_.time().timeName(),
                    mesh_,
                    readOpt,
                    writeOpt),
                mesh_,
                dim);
        }

        // Transfer ownership of this object to the objectRegistry
        ptr->store();

        // Add name of field to the list of ownedFields of the thermoClass.
        ownedFieldNames_.push_back(field_name);

        // Register the field in the thermo-class's fields_
        fields_.emplace(field_name, *ptr);
    }
    // Field is part of object registry
    else
    {
        // Field is already registerd in fields_
        if (fields_.contains(field_name))
        {
            Info << "Field " << field_name << " already exists in fields_. Nothing is done." << endl;
        }
        // Field is in objectRegistry, but not in fields_ -> field gets registered.
        else
        {
            Info << "Existing Field " << field_name << " from object registry is registered in fields_." << endl;
            fields_.emplace(field_name, *ptr);
        }
    }
}

void Foam::flameletThermo::setBoundaryLookupType()
{
    // Setup of the boundary lookup (special treatment).
    if (Foam::contains(inputVariables_, "ha"))
    {
        // Check if a boundary lookup needs to be performed
        volScalarField& ha = getField("ha");
        forAll(ha.boundaryField(), patchI)
        {
            const fvPatchScalarField& patch = ha.boundaryField()[patchI];

            if (isA<boundaryLookUpFvPatchScalarField>(patch))
            {
                const boundaryLookUpFvPatchScalarField& bounday_patch = dynamic_cast<const boundaryLookUpFvPatchScalarField&>(patch);
                const word lookupType = bounday_patch.getLookupType();
                setBoundary(lookupType, patchI);
            }
        }
    }
}

void Foam::flameletThermo::setBoundary(std::string type, label patchID) const
{
    if (boundaryLookupType_ == "none")
    {
        if (type == "specialBoundaryLookup")
        {
            Info << "Use special Lookup for boundary look-up similar to Ketelheun2013. (Default behaviour before 04/2022)" << endl;
            boundaryLookupType_ = type;
        }
        else if (type == "standardBoundaryLookup")
        {
            Info << "Use flut reader handled boundary lookup. This was used for the HOQ-EGR Fluts in Luo2023." << endl;
            boundaryLookupType_ = type;
        }
        else
        {
            FatalErrorInFunction << "The given boundaryType: " + type + " is unknown use: specialBoundaryLookup or boundaryLookup." << exit(FatalError);
        }
    }
    else if (boundaryLookupType_ != type)
    {
        FatalErrorInFunction << "The given boundaryType: " + type + " does not match the already set boundaryType: " + boundaryLookupType_ + ". Check your boundary conditions for ha consistency." << exit(FatalError);
    }

    // Add patchID to patches to the boundary patches to update.
    auto it = std::find(boundaryLookupPatchIDs_.begin(), boundaryLookupPatchIDs_.end(), patchID);
    if (it == boundaryLookupPatchIDs_.end())
    {
        Info << "Add patch " << patchID << " to boundaryFaces." << endl;
        boundaryLookupPatchIDs_.push_back(patchID);
    }
}

//
// Return access (custom functions)
//
const Foam::dictionary& Foam::flameletThermo::dict() const
{
    return thermoDict_;
};

const Foam::volScalarField& Foam::flameletThermo::getField(word name) const
{
    // Attempt to retrieve the field from the map
    try
    {
        return fields_.at(name);
    }
    catch (const std::out_of_range& e)
    {
        // Handle the case where the field does not exist
        FatalErrorInFunction << "Field " << name << " not found in fields_ in flameletThermo class. The field can be registered using thermo.appendField(FIELD_NAME)." << exit(FatalError);
    }
}

Foam::volScalarField& Foam::flameletThermo::getField(word name)
{
    // Attempt to retrieve the field from the map
    try
    {
        return fields_.at(name);
    }
    catch (const std::out_of_range& e)
    {
        // Handle the case where the field does not exist
        FatalErrorInFunction << "Field " << name << " not found in fields_ in flameletThermo class. The field can be registered using thermo.appendField(FIELD_NAME)." << exit(FatalError);
    }
}

Foam::volScalarField& Foam::flameletThermo::alpha()
{
    return alpha_;
}

const Foam::fvMesh& Foam::flameletThermo::getMesh() const
{
    return mesh_;
}

//
// Return access (OF-specific)
//
Foam::tmp<Foam::volScalarField> Foam::flameletThermo::rho() const
{
    return getField("thermo:rho");
}

Foam::tmp<Foam::scalarField> Foam::flameletThermo::rho(const label patchi) const
{
    return getField("thermo:rho").boundaryField()[patchi];
}

Foam::volScalarField& Foam::flameletThermo::rhoRef()
{
    return getField("thermo:rho");
}

void Foam::flameletThermo::correctRho(
    const Foam::volScalarField& deltaRho,
    const dimensionedScalar& rhoMin,
    const dimensionedScalar& rhoMax)
{

    getField("thermo:rho") += deltaRho;
    getField("thermo:rho") = max(getField("thermo:rho"), rhoMin);
    getField("thermo:rho") = min(getField("thermo:rho"), rhoMax);
}

void Foam::flameletThermo::correctRho(const Foam::volScalarField& deltaRho)
{
    getField("thermo:rho") += deltaRho;
}

const Foam::volScalarField& Foam::flameletThermo::psi() const
{
    return getField("thermo:psi");
}

Foam::tmp<Foam::volScalarField> Foam::flameletThermo::mu() const
{
    return getField("mu");
}

Foam::tmp<Foam::scalarField> Foam::flameletThermo::mu(const label patchi) const
{
    return getField("mu").boundaryField()[patchi];
}

Foam::tmp<Foam::volScalarField> Foam::flameletThermo::Cpv() const
{
    // convert cp to cv, based on gas law
    //    Info << "Cpv is: " << getField("cp") - p_ / ( getField("thermo:rho") * T_ ) << endl;
    return getField("cp") - p_ / (getField("thermo:rho") * T_);
}
Foam::tmp<Foam::volScalarField> Foam::flameletThermo::W() const
{
    return getField("W");
}
// ************************************************************************* //
