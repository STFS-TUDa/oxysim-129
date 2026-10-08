#include "inputVariableEquations.H"
#include "CloudHandler.H"
#include "basicCarbonaceousCloud.H"
#include "SLGThermo.H"
#include <numeric>
#include "cellSet.H"
#include "radiationModel.H"

// Initialization functions
void InputVariableEquations::init_sourceTerm(std::string name, dictionary settings, flameletThermo& thermo)
{
    thermo.appendField(name, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE, dimensionedScalar(name, dimensionSet(0, 0, 0, 0, 0), 0.0));
    thermo.appendOutputField("omega_" + name, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE, dimensionedScalar("omega_", dimDensity / dimTime, 0.0));
}

void InputVariableEquations::init_empty(std::string name, dictionary settings, flameletThermo& thermo) {}

void InputVariableEquations::init_name(std::string name, dictionary settings, flameletThermo& thermo)
{
    thermo.appendField(name, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE, dimensionedScalar(name, dimensionSet(0, 0, 0, 0, 0), 0.0));
}


// Update functions
void InputVariableEquations::noUpdate(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
}


void InputVariableEquations::unityLewis(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    // Get input fields
    volScalarField& controlVariable = thermo.getField(name);
    // Get dependent fields from FLUT.
    volScalarField& sourceTerm = thermo.getField("omega_" + name);

    fvScalarMatrix Eqn(
        fvm::ddt(rho, controlVariable)
            + fvm::div(phi, controlVariable)
            - fvm::laplacian(thermo.alpha() + turbulence.mut() / thermo.Sct(), controlVariable)
        == sourceTerm
               + fvOptions(rho, controlVariable));

    Eqn.relax();
    fvOptions.constrain(Eqn);
    Eqn.solve();
    fvOptions.correct(controlVariable);

    (controlVariable).correctBoundaryConditions();
}

void InputVariableEquations::unityLewis_noSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    // Get input fields
    volScalarField& controlVariable = thermo.getField(name);

    fvScalarMatrix Eqn(
        fvm::ddt(rho, controlVariable)
            + fvm::div(phi, controlVariable)
            - fvm::laplacian(thermo.alpha() + turbulence.mut() / thermo.Sct(), controlVariable)
        == fvOptions(rho, controlVariable));

    Eqn.relax();
    fvOptions.constrain(Eqn);
    Eqn.solve();
    fvOptions.correct(controlVariable);

    (controlVariable).correctBoundaryConditions();
}

void InputVariableEquations::varianceAlgebraic(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{

    const volScalarField& delta = thermo.getMesh().objectRegistry::lookupObject<volScalarField>("delta");

    scalar Cvar = settings.get<scalar>("Cvar");

    word variableName = settings.get<word>("variableName");

    const volScalarField& variable = thermo.getField(variableName);

    volScalarField& varianceVariable = thermo.getField(name);

    varianceVariable = Cvar * magSqr(delta) * magSqr((fvc::grad(variable)));
}

// -------------------- Lagrangian (particle–gas coupling) -------------------- //

void InputVariableEquations::init_particleSource(std::string name, dictionary settings, flameletThermo& thermo)
{
    thermo.appendField(name, IOobject::READ_IF_PRESENT, IOobject::AUTO_WRITE, dimensionedScalar(name, dimensionSet(0, 0, 0, 0, 0), 0.0));
    thermo.appendField("S" + name, IOobject::NO_READ, IOobject::AUTO_WRITE, dimensionedScalar("S", dimDensity / dimTime, 0.0));
}

void InputVariableEquations::init_sourceTerm_particleSource(std::string name, dictionary settings, flameletThermo& thermo)
{
    thermo.appendField("S" + name, IOobject::NO_READ, IOobject::AUTO_WRITE, dimensionedScalar("S", dimDensity / dimTime, 0.0));
    thermo.appendOutputField("omega_" + name, IOobject::NO_READ, IOobject::AUTO_WRITE, dimensionedScalar("omega", dimDensity / dimTime, 0.0));
}

void InputVariableEquations::unityLewis_particleSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    volScalarField& controlVariable = thermo.getField(name);
    volScalarField& particleSourceTerm = thermo.getField("S" + name);

    fvScalarMatrix Eqn(
        fvm::ddt(rho, controlVariable)
            + fvm::div(phi, controlVariable)
            - fvm::laplacian(thermo.alpha() + turbulence.mut() / thermo.Sct(), controlVariable)
        == particleSourceTerm
               + fvOptions(rho, controlVariable));

    Eqn.relax();
    fvOptions.constrain(Eqn);
    Eqn.solve();
    fvOptions.correct(controlVariable);

    (controlVariable).correctBoundaryConditions();
}

void InputVariableEquations::unityLewis_sourceTerm_particleSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    volScalarField& controlVariable = thermo.getField(name);
    volScalarField& sourceTerm = thermo.getField("omega_" + name);
    volScalarField& particleSourceTerm = thermo.getField("S" + name);

    fvScalarMatrix Eqn(
        fvm::ddt(rho, controlVariable)
            + fvm::div(phi, controlVariable)
            - fvm::laplacian(thermo.alpha() + turbulence.mut() / thermo.Sct(), controlVariable)
        == sourceTerm + particleSourceTerm
               + fvOptions(rho, controlVariable));

    Eqn.relax();
    fvOptions.constrain(Eqn);
    Eqn.solve();
    fvOptions.correct(controlVariable);

    (controlVariable).correctBoundaryConditions();
}

template<class T>
void InputVariableEquations::particleSourceHelper(const T& cloud, List<Tuple2<word, scalar>> speciesWeights, volScalarField& particleSourceTerm)
{
    if (cloud.solution().active())
    {
        forAll(speciesWeights, i)
        {
            const word& speciesName = speciesWeights[i].first();
            const scalar& speciesWeight = speciesWeights[i].second();
            label gasIdxi = cloud.thermo().carrierId(speciesName);
            particleSourceTerm.internalFieldRef() += speciesWeight * cloud.Srho(gasIdxi);
        }
    }
}

void InputVariableEquations::particleSource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    Info << "Updating " << name << endl;
    volScalarField& particleSourceTerm = thermo.getField(name);
    particleSourceTerm *= 0.0;

    List<Tuple2<word, scalar>> speciesWeights = settings.lookup("speciesWeights");
    List<word> cloudNames(settings.get<List<word>>("cloudNames"));

    for (word cloudName : cloudNames)
    {
        const fvMesh& mesh = thermo.getMesh();
        CloudHandler& cloudHandler = mesh.objectRegistry::lookupObjectRef<CloudHandler>("cloudHandler");
        word cloudType(cloudHandler.getCloudType(cloudName));
        if (cloudType == "carbonaceousCloud")
        {
            const basicCarbonaceousCloud& cloud = cloudHandler.getCarbonaceousCloud(cloudName);
            particleSourceHelper<basicCarbonaceousCloud>(cloud, speciesWeights, particleSourceTerm);
        }
    }
    Info << name << " min/max: " << min(particleSourceTerm).value() << ", " << max(particleSourceTerm).value() << endl;
}

void InputVariableEquations::sum(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    volScalarField& controlVariable = thermo.getField(name);
    List<word> summandList = settings.lookup("summands");
    controlVariable = std::accumulate(summandList.begin(), summandList.end(), controlVariable * 0, [&](volScalarField a, word b)
                                      { return a + thermo.getField(b); });
}

void InputVariableEquations::fraction(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    volScalarField& controlVariable = thermo.getField(name);
    List<word> numeratorList = settings.lookup("numerators");
    List<word> denominatorList = settings.lookup("denominators");
    controlVariable = std::accumulate(numeratorList.begin(), numeratorList.end(), controlVariable * 0, [&](volScalarField a, word b)
                                      { return a + thermo.getField(b); })
                    / (std::accumulate(denominatorList.begin(), denominatorList.end(), controlVariable * 0, [&](volScalarField a, word b)
                                       { return a + thermo.getField(b); })
                       + 1e-15);
}

template<class T>
Foam::fvScalarMatrix InputVariableEquations::unityLewis_particleEnthalpySourceHelper(const T& cloud, volScalarField& controlVariable)
{
    return cloud.Sha(controlVariable);
}

void InputVariableEquations::unityLewis_particleEnthalpySource(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    volScalarField& controlVariable = thermo.getField(name);
    Foam::tmp<Foam::fvScalarMatrix> haMatrix(new fvScalarMatrix(controlVariable, dimEnergy / dimTime));
    Foam::fvScalarMatrix& haMatrixRef = haMatrix.ref();

    const fvMesh& mesh = thermo.getMesh();

    List<word> cloudNames(settings.get<List<word>>("cloudNames"));
    for (word cloudName : cloudNames)
    {
        CloudHandler& cloudHandler = mesh.objectRegistry::lookupObjectRef<CloudHandler>("cloudHandler");
        word cloudType(cloudHandler.getCloudType(cloudName));
        if (cloudType == "carbonaceousCloud")
        {
            const basicCarbonaceousCloud& cloud = cloudHandler.getCarbonaceousCloud(cloudName);
            haMatrixRef += InputVariableEquations::unityLewis_particleEnthalpySourceHelper(cloud, controlVariable);
        }
    }

    Foam::radiation::radiationModel& radiation =
        mesh.objectRegistry::lookupObjectRef<Foam::radiation::radiationModel>("radiationProperties");

    fvScalarMatrix Eqn(
        fvm::ddt(rho, controlVariable)
            + fvm::div(phi, controlVariable)
            - fvm::laplacian(thermo.alpha() + turbulence.mut() / thermo.Sct(), controlVariable)
        == haMatrix
               + radiation.Sh(thermo, controlVariable)
               + fvOptions(rho, controlVariable));

    Eqn.relax();
    fvOptions.constrain(Eqn);
    Eqn.solve();
    fvOptions.correct(controlVariable);

    (controlVariable).correctBoundaryConditions();
}

void InputVariableEquations::scaleSourceTerm(std::string name, dictionary settings, flameletThermo& thermo, compressible::turbulenceModel& turbulence, volScalarField& rho, surfaceScalarField& phi, fv::options& fvOptions)
{
    word variableName(settings.getWord("variableName"));
    volScalarField& controlVariable = thermo.getField(variableName);
    word setName(settings.getWord("cellSetName"));
    scalar factor(settings.getScalar("factor"));
    scalar summand(settings.getScalar("summand"));

    cellSet cells(thermo.getMesh(), setName);
    const labelList& cells_ = cells.toc();

    forAll(cells_, i)
    {
        label cellID = cells_[i];
        controlVariable[cellID] *= factor;
        controlVariable[cellID] += summand;
    }
}
