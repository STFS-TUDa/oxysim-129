{
    dimensionedScalar totalMass = fvc::domainIntegrate(rho);
    volScalarField continuityErr(fvc::ddt(rho) + fvc::div(phi));

    scalar sumLocalContErr =
        runTime.deltaTValue() * (fvc::domainIntegrate(mag(continuityErr)) / totalMass).value();

    scalar globalContErr =
        runTime.deltaTValue() * (fvc::domainIntegrate(continuityErr) / totalMass).value();

    cumulativeContErr += globalContErr;

    Info << "time step continuity errors : sum local = " << sumLocalContErr
         << ", global = " << globalContErr
         << ", cumulative = " << cumulativeContErr
         << ", total Mass in system = " << totalMass.value()
         << endl;
}

// ************************************************************************* //
