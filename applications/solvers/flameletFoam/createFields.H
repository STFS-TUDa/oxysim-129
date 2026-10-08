// Create velocity field
Info << "Reading field U\n"
     << endl;
volVectorField U(
    IOobject(
        "U",
        runTime.timeName(),
        mesh,
        IOobject::MUST_READ,
        IOobject::AUTO_WRITE),
    mesh);

// -
// Dummy T field for flamelet thermo class
volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimTemperature, 0.0));

// - Create flameletThermo class.
flameletThermo thermo(mesh);

// - Initialization of the equationHandler.
// - needs to be done after thermo initalisation and before the initial thermo update
auto equationHandler = InputVariableEquations::EquationHandler(thermo);

// - Initial thermo update
Info << "Load output-variables from FLUT." << endl;
thermo.loadOutputVariablesFromFLUT();

Info << "Performing initial thermo update: \n"
     << endl;
thermo.update();

// - Create Fields
volScalarField& p(thermo.p());
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    thermo.rho());

const volScalarField& psi(thermo.psi());

// determine if low mach formulation of pressure or thermodynamik pressure is taken
const bool lowMachPressure(pimple.dict().lookupOrDefault<bool>("lowMach", true));
std::unique_ptr<volScalarField> pd;
if (lowMachPressure)
{
    // creating new pd object when low mach formulation is used
    Info << "Taking the low mach formulation of the pressure." << endl;
    pd = std::make_unique<volScalarField>(
        IOobject(
            "pd",
            runTime.timeName(),
            mesh,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE),
        mesh);
}
else
{
    // resetting the pd pointer to the field p when "normal formulation is used"
    Info << "Taking the compressible pressure equation from reactingFoam." << endl;
    pd.reset(&p);
}

pressureControl pressureControl(*pd, rho, pimple.dict(), false);
mesh.setFluxRequired(pd->name());


#include "compressibleCreatePhi.H"
#include "createDpdt.H"
#include "createK.H"
#include "createMRF.H"
#include "createFvOptions.H"

Info << "Create turbulence model.\n"
     << nl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        thermo));

benchmarkUtility timer(runTime);