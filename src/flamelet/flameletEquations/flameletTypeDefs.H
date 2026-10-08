
#include <functional>
#include "fvCFD.H"
#include "fvOptions.H"
#include "flameletThermo.H"
#include "turbulentFluidThermoModel.H"

namespace InputVariableEquations
{
/** @addtogroup equations
 *  @{
 */
typedef std::function<void(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)> FunctionType;
typedef std::function<void(std::string name, dictionary settings, flameletThermo& thermo)> InitFunctionType;
typedef std::function<void()> FlameletFunction;
/**@} */
} // end namespace
