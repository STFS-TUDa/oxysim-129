#include <iostream>
#include "EquationHandler.H"

// Constructor
InputVariableEquations::EquationHandler::EquationHandler()
{
    // Create functionMap
    std::map<std::string, std::pair<InputVariableEquations::InitFunctionType, InputVariableEquations::FunctionType>> functionMap;
    // Add functions to functionMap
    appendFunctionMap();
}

InputVariableEquations::EquationHandler::EquationHandler(flameletThermo& thermo)
{
    // Create functionMap
    std::map<std::string, std::pair<InputVariableEquations::InitFunctionType, InputVariableEquations::FunctionType>> functionMap;
    // Add functions to functionMap
    appendFunctionMap();

    initializeEquations(thermo);
}

// Function map to hold your equations.
void InputVariableEquations::EquationHandler::appendFunctionMap()
{
    // Here we append functions to the functionMap
    // 1. Function for the initalization of the necessary fields.
    // 2. Equation to solve in OpenFOAM.
    // These are given in the pairs: {INIT_FUNCTION, FUNCTION}
    functionMap["noUpdate"] = {init_empty, noUpdate};
    functionMap["unityLewis"] = {init_sourceTerm, unityLewis};
    functionMap["unityLewis_noSource"] = {init_name, unityLewis_noSource};
    functionMap["varianceAlgebraic"] = {init_name, varianceAlgebraic};
    // Lagrangian (particle–gas coupling)
    functionMap["unityLewis_particleSource"] = {init_particleSource, unityLewis_particleSource};
    functionMap["unityLewis_sourceTerm_particleSource"] = {init_sourceTerm_particleSource, unityLewis_sourceTerm_particleSource};
    functionMap["particleSource"] = {init_empty, particleSource};
    functionMap["sum"] = {init_empty, sum};
    functionMap["fraction"] = {init_name, fraction};
    functionMap["unityLewis_particleEnthalpySource"] = {init_empty, unityLewis_particleEnthalpySource};
    functionMap["scaleSourceTerm"] = {init_empty, scaleSourceTerm};
}

// Initializes functions (e.g., adds necessary fields) and registers run-time selected update equations into function pointer.
void InputVariableEquations::EquationHandler::initializeAndRegisterEquations(
    flameletThermo& thermo,
    compressible::turbulenceModel& turbulence,
    volScalarField& rho,
    surfaceScalarField& phi,
    fv::options& fvOptions)
{
    initializeEquations(thermo);
    registerEquations(thermo, turbulence, rho, phi, fvOptions);
}

void InputVariableEquations::EquationHandler::initializeEquations(
    flameletThermo& thermo)
{
    const Foam::dictionary& flameletProperties = thermo.dict();

    // Get function and variable names (and additional Settings) from the flameletProperties
    std::vector<InputVariableEquations::VariableEquationSettings> funcVarnameSettings = getFunctionAndVariableNames(flameletProperties);


    for (const auto& funcVarnameSetting : funcVarnameSettings)
    {
        std::string functionName = funcVarnameSetting.updateType;
        std::string controlVariableName = funcVarnameSetting.variable;
        dictionary settings = funcVarnameSetting.settings;

        // Treatment for controlVariables
        auto it = functionMap.find(functionName);
        if (it != functionMap.end())
        {
            // Initialization
            Info << "Initialize field for function: " << functionName << endl;
            it->second.first(controlVariableName, settings, thermo);
        }
        else
        {
            throw std::invalid_argument("Initialization function " + functionName + " for variable: " + controlVariableName + " not found.");
        }
    }
}


void InputVariableEquations::EquationHandler::registerEquations(
    flameletThermo& thermo,
    compressible::turbulenceModel& turbulence,
    volScalarField& rho,
    surfaceScalarField& phi,
    fv::options& fvOptions)
{

    const Foam::dictionary& flameletProperties = thermo.dict();

    // Get function and variable names (and additional Settings) from the flameletProperties
    std::vector<InputVariableEquations::VariableEquationSettings> funcVarnameSettings = getFunctionAndVariableNames(flameletProperties);


    for (const auto& funcVarnameSetting : funcVarnameSettings)
    {
        std::string functionName = funcVarnameSetting.updateType;
        std::string controlVariableName = funcVarnameSetting.variable;
        dictionary settings = funcVarnameSetting.settings;

        // Treatment for controlVariables
        auto it = functionMap.find(functionName);
        if (it != functionMap.end())
        {
            // Update Functions (added to the function vectors)
            Info << "Add update function for variable " << controlVariableName << " to function pointer: " << functionName << endl;
            variableEquations.push_back([&, it, controlVariableName, settings]()
                                        { it->second.second(controlVariableName, settings, thermo, turbulence, rho, phi, fvOptions); });
        }
        else
        {
            throw std::invalid_argument("Update function " + functionName + " for variable: " + controlVariableName + " not found.");
        }
    }
}


void InputVariableEquations::EquationHandler::solve()
{
    for (const auto& function : variableEquations)
    {
        function(); // Invoke the function through the std::function
    }
}


// Helper functions
std::vector<InputVariableEquations::VariableEquationSettings> InputVariableEquations::getFunctionAndVariableNames(const dictionary& flameletProperties)
{
    // Create struct to store updateTyp, variable names and additional settings
    std::vector<InputVariableEquations::VariableEquationSettings> funcVarnameSettings;

    // Get equations to register from flameletProperties
    hashedWordList variablesToSolve(flameletProperties.lookup("variablesToSolve"), true);

    // Iterate over the variablesToSolve from flameletProperties
    for (auto variable : variablesToSolve)
    {
        {
            initializeAndAppendFunction(variable, variable, flameletProperties, funcVarnameSettings);
        }
    }

    return funcVarnameSettings;
}

void InputVariableEquations::initializeAndAppendFunction(std::string variable,
                                                         std::string variable_name_in_dict,
                                                         const dictionary& flameletProperties,
                                                         std::vector<InputVariableEquations::VariableEquationSettings>& funcVarnameSettings)
{
    const Foam::dictionary& varEq = flameletProperties.subDict("variableEquations").subDict(variable_name_in_dict);

    InputVariableEquations::VariableEquationSettings s;
    s.updateType = varEq.get<word>("updateType");
    s.variable = variable;

    List<keyType> varEqKeyList = varEq.keys();
    int n_keys = varEqKeyList.size();

    if (n_keys > 1)
    {
        Info << "Settings found for variable " << variable << " with updateType " << s.updateType << ":" << endl;
        s.settings = varEq;
        s.settings.remove("updateType");
        Info << s.settings << endl;
    }
    else
    {
        s.settings.clear();
    }

    funcVarnameSettings.push_back(std::move(s));
}
