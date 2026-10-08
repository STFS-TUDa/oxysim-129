rho = thermo.rho();

// Thermodynamic density needs to be updated by psi*d(p) after the
// pressure solution
// const volScalarField psip0(psi*pd);

volScalarField rAU(1.0 / UEqn.A());
surfaceScalarField rhorAUf("rhorAUf", fvc::interpolate(rho* rAU));
volVectorField HbyA(constrainHbyA(rAU* UEqn.H(), U, *pd));

if (pimple.nCorrPISO() <= 1)
{
    tUEqn.clear();
}

surfaceScalarField phiHbyA(
    "phiHbyA",
    fvc::interpolate(rho) * fvc::flux(HbyA) // generally it is a good idea to interpolate individually, however we should check what is the difference
        + rhorAUf * fvc::ddtCorr(rho, U, phi));

// Update the pressure BCs to ensure flux consistency
constrainPressure(*pd, rho, U, phiHbyA, rhorAUf, MRF);

while (pimple.correctNonOrthogonal())
{
    fvScalarMatrix pdEqn(
        fvc::ddt(rho)
        // Unsure if needed: should be zero as ddt(p) is 0
        // + psi*correction(fvm::ddt(p))
        + fvc::div(phiHbyA)
        - fvm::laplacian(rhorAUf, *pd)
        //==
        //   fvOptions(psi, *pd, rho.name()) // not sure what psi is for us -> so commented for now
    );

    pdEqn.solve(mesh.solver(pd->select(pimple.finalInnerIter())));

    if (pimple.finalNonOrthogonalIter())
    {
        // Calculate the conservative fluxes
        phi = phiHbyA + pdEqn.flux();

        // Explicitly relax pressure for momentum corrector
        pd->relax();

        // Correct the momentum source with the pressure gradient flux
        // calculated from the relaxed pressure
        U = HbyA + rAU * fvc::reconstruct(pdEqn.flux() / rhorAUf);
        U.correctBoundaryConditions();
        fvOptions.correct(U);
    }
}

#include "lowMachContinuityErrs.H"
