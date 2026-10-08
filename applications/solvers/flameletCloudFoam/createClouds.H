IOdictionary controlDict(
    IOobject(
        "controlDict",
        runTime.system(),
        runTime,
        IOobject::MUST_READ,
        IOobject::NO_WRITE));

dictionary cloudsDict = controlDict.subDict("clouds");


CloudHandler cloudHandler(mesh, cloudsDict, rho, U, g, slgThermo);

// word cloudName("coalCloud");
// const basicCarbonaceousCloud& carbonaceousCloud = cloudHandler.getCarbonaceousCloud(cloudName);
//
// const CloudHandler& cloudHandler_const = mesh.objectRegistry::lookupObject<CloudHandler>("cloudHandler");
// CloudHandler& cloudHandler_ref = mesh.objectRegistry::lookupObjectRef<CloudHandler>("cloudHandler");
