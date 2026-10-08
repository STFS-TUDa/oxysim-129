
#pragma once

#include <map>
#include <string> // Include this for std::string

#include "flameletThermo.H"
#include "fvCFD.H"
#include "fvOptions.H"
#include "hashedWordList.H"
#include "tmp.H"
#include "flameletTypeDefs.H"

#include "inputVariableEquations.H"

namespace InputVariableEquations
{
struct VariableEquationSettings
{
    word updateType;     // e.g. "mix", "calcFromTransportedSpecies"
    word variable;       // e.g. "Z", "yc"
    dictionary settings; // additional Settings from the dict
};

// ------------------------------------------------------------------ //
// Class definition: EquationHandler
// ------------------------------------------------------------------ //
/** @addtogroup equations
 *  @{
 */
class EquationHandler
{
    std::vector<InputVariableEquations::FlameletFunction> variableEquations;

protected:
    std::map<std::string, std::pair<InitFunctionType, FunctionType>> functionMap;

public:
    // Constructor
    EquationHandler();
    EquationHandler(flameletThermo& thermo);
    virtual ~EquationHandler() = default;

    virtual void appendFunctionMap();
    void solve();

    void initializeAndRegisterEquations(
        flameletThermo& thermo,
        compressible::turbulenceModel& turbulence,
        volScalarField& rho,
        surfaceScalarField& phi,
        fv::options& fvOptions);
    void initializeEquations(
        flameletThermo& thermo);
    void registerEquations(
        flameletThermo& thermo,
        compressible::turbulenceModel& turbulence,
        volScalarField& rho,
        surfaceScalarField& phi,
        fv::options& fvOptions);
};

// Helper equations
std::vector<VariableEquationSettings> getFunctionAndVariableNames(const dictionary& flameletProperties);
void initializeAndAppendFunction(std::string variable, std::string variable_name_in_dict, const dictionary& flameletProperties, std::vector<VariableEquationSettings>& funcVarnameSettings);
/**@} */
} // end namespace
