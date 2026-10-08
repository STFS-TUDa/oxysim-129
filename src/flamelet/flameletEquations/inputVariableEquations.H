// ------------------------------------------------------------------ //
// Scalar transport equations for lookup variables
// ------------------------------------------------------------------ //

// Info: When you add an equation, you also need to add the respective option to the function
// 'InputVariableEquations::EquationHandler::appendFunctionMap()'

#include "flameletTypeDefs.H"
#include "flameletThermo.H"

namespace InputVariableEquations
{
/** @addtogroup equations
 *  @{
 */
// Initialization functions to create necessary fields
void init_sourceTerm(std::string name, dictionary settings, flameletThermo& thermo);
void init_empty(std::string name, dictionary settings, flameletThermo& thermo);
void init_name(std::string name, dictionary settings, flameletThermo& thermo);

// Update functions
void noUpdate(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
void unityLewis(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
void unityLewis_noSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);

void varianceAlgebraic(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);

// Lagrangian (particle–gas coupling) equations
void init_particleSource(std::string name, dictionary settings, flameletThermo& thermo);
void init_sourceTerm_particleSource(std::string name, dictionary settings, flameletThermo& thermo);

void unityLewis_particleSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
void unityLewis_sourceTerm_particleSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);

template<class T>
void particleSourceHelper(const T& cloud, List<Tuple2<word, scalar>> speciesWeights, volScalarField& particleSourceTerm);

void particleSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
void sum(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
void fraction(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);

template<class T>
Foam::fvScalarMatrix unityLewis_particleEnthalpySourceHelper(const T& cloud, volScalarField& controlVariable);
void unityLewis_particleEnthalpySource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
void scaleSourceTerm(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions);
/**@} */
} // end namespace
