equationHandler.registerEquations(flameletThermo, turbulence.ref(), rho, phi, fvOptions);

// get inital fields from the thermo class (e.g. Cp, alpha, ...) after init of the equations
// this thermo.update() is repeated here, because the one in createFields is done before (especially important for the transportSpecies option)
// Info << "Performing initial lookup:\n"
//      << endl;
// thermo.update();
// Info << "Finished initial lookup\n"
//      << endl;
