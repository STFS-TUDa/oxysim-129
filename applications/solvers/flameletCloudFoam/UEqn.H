// Solve the Momentum equation

MRF.correctBoundaryVelocity(U);

tmp<fvVectorMatrix> tUEqn(
    fvm::ddt(rho, U) + fvm::div(phi, U)
        + MRF.DDt(rho, U)
        + turbulence->divDevRhoReff(U)
    == fvOptions(rho, U) + cloudHandler.SUClouds(U, false) + g * rho);
fvVectorMatrix& UEqn = tUEqn.ref();

UEqn.relax();

fvOptions.constrain(UEqn);

if (pimple.momentumPredictor())
{
    solve(UEqn == -fvc::grad(*pd));

    fvOptions.correct(U);
    K = 0.5 * magSqr(U);
}
