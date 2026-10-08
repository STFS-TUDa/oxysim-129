// - Register database.
#include "readGravitationalAcceleration.H"
const objectRegistry& db = mesh.thisDb();

volScalarField T(
    IOobject(
        "T",
        runTime.timeName(),
        mesh,
        IOobject::READ_IF_PRESENT,
        IOobject::AUTO_WRITE),
    mesh,
    dimensionedScalar("zero", dimTemperature, 0.0));


Foam::canteraThermo thermo(mesh, "", "thermophysicalProperties", false);
SLGThermo slgThermo(mesh, thermo);

// volScalarField ha_init(
//     IOobject(
//         "ha",
//         runTime.timeName(),
//         mesh,
//         IOobject::MUST_READ,
//         IOobject::NO_WRITE),
//     mesh);
//
// Info << thermo.he() << endl;
//
// volScalarField& he = thermo.he();
// he = ha_init;
// forAll(he.boundaryField(), patchi)
//{
//     fvPatchScalarField& phe = he.boundaryFieldRef()[patchi];
//     fvPatchScalarField& pha_init = ha_init.boundaryFieldRef()[patchi];
//     forAll(phe, facei)
//     {
//         phe[facei] = pha_init[facei];
//     }
// }
//
// Info << he << endl;
//
Foam::flameletThermo flameletThermo(mesh);
//
// Info << flameletThermo.getField("ha") << endl;
// Info << thermo.he() << endl;

auto equationHandler = InputVariableEquations::EquationHandler();
equationHandler.initializeEquations(flameletThermo);

wordList variablesToSolve_ = flameletThermo.lookup("variablesToSolve");
for (word specieName : thermo.species())
{
    if (!variablesToSolve_.found(specieName))
    {
        flameletThermo.appendOutputField(specieName);
    }
    else
    {
        flameletThermo.appendField(specieName);
    }
}
flameletThermo.appendOutputField("cp");

//// Flamelet thermo needs to be created after canteraThermo.
//// This way the fields inside the flameletThermo are created as references to the canteraThermo class
//// Exception: thermo.alpha(). This fields needs to be manually overwritten.
//// This cannot be avoided, becaue both thermo clases inherit from fluidThermo, where the alpha field is created.
//// - Create flameletThermo class.
// Foam::flameletThermo flameletThermo(mesh);

Info << "Initial flamelet update" << endl;
flameletThermo.update();
thermo.alpha() = flameletThermo.alpha();

// - Create Fields
volScalarField rho(
    IOobject(
        "rho",
        runTime.timeName(),
        mesh,
        IOobject::NO_READ,
        IOobject::AUTO_WRITE),
    flameletThermo.rho());

const volScalarField& psi = thermo.psi();
// Create canteraThermo
// This thermo is named thermo, because it is needed for some of the #include functionalities.
// They want to have a thermo called "thermo"
// The canteraThermo is the better one to use, because it implements all thermo interfaces from OpenFOAM
// What the canteraThermo does:
// - createFields
// - h = f(T,p,Yi)
// Problem: Yi might be incorrect -> enthalpy is incorrect

// volScalarField& he = thermo.he();
// he = ha_init;

// - Lookup if species are transported
// - if yes, reconstruct species from the given variables in the input_dict
// - needs to be done after thermo initalisation and before thermo update
IOdictionary flameletProperties(
    IOobject(
        "flameletProperties",
        mesh.time().constant(),
        mesh,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

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

volScalarField& p(flameletThermo.p());

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

scalar pRef_dimless(
    mesh.solutionDict().subDict("PIMPLE").lookupOrDefault<scalar>("pRef", 101325));

dimensionedScalar pRef("pRef", dimensionSet(1, -1, -2, 0, 0, 0, 0), pRef_dimless);

#include "compressibleCreatePhi.H"

// pressureControl pressureControl(p, rho, pimple.dict(), false);
pressureControl pressureControl(*pd, rho, pimple.dict(), false);
mesh.setFluxRequired(pd->name());

Info << "Creating turbulence model.\n"
     << nl;
autoPtr<compressible::turbulenceModel> turbulence(
    compressible::turbulenceModel::New(
        rho,
        U,
        phi,
        flameletThermo));

#include "createDpdt.H"

#include "createK.H"

#include "createMRF.H"

#include "createClouds.H"

#include "createRadiationModel.H"

const dimensionedScalar rhoMax("rhoMax", dimDensity, GREAT, pimple.dict());
const dimensionedScalar rhoMin("rhoMin", dimDensity, Zero, pimple.dict());

#include "createFvOptions.H"
