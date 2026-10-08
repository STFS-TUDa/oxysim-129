// Solve the Momentum equation
timer.startBenchmark("UEqn");

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn(
    fvm::ddt(rho, U)
        + fvm::div(phi, U)
        + MRF.DDt(rho, U)
        //- fvm::laplacian(mu, U)
        //- fvc::div(mu * dev2(Foam::T(fvc::grad(U)())))
        + turbulence->divDevRhoReff(U)
    == fvOptions(rho, U) + rho * g);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    // depending on lowMach choice pd is either pd from the low mach formulation (separate field)
    // or this is reference to p when the "compressible" formulations are considered
    solve(UEqn == -fvc::grad(*pd));

    fvOptions.correct(U);
    K = 0.5 * magSqr(U);
}

timer.stopBenchmark("UEqn");
