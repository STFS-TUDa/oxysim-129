#include "CloudHandler.H"

defineTypeNameAndDebug(CloudHandler, 0);

CloudHandler::CloudHandler(
    fvMesh& i_mesh,
    dictionary cloudDict,
    const volScalarField& rho,
    const volVectorField& U,
    const dimensionedVector& g,
    const SLGThermo& slgThermo) : regIOobject(IOobject("cloudHandler",
                                                       i_mesh.time().name(),
                                                       i_mesh,
                                                       IOobject::NO_READ,
                                                       IOobject::NO_WRITE,
                                                       true)),
                                  mesh(i_mesh)
{
    int idx_carb(0);
    int idx(0);

    List<keyType> cloudKeyList = cloudDict.keys();
    for (keyType key : cloudKeyList)
    {
        word cloudName(key);
        word cloudType(cloudDict.lookup(cloudName));

        if (cloudType == "carbonaceousCloud")
        {
            Info << "\nConstructing carbonaceous cloud: " << cloudName << endl;

            carbonaceousCloudList.push_back(
                std::make_unique<basicCarbonaceousCloud>(
                    cloudName,
                    rho,
                    U,
                    g,
                    slgThermo));

            Info << "Construction done" << endl;
            idx = idx_carb;
            idx_carb++;
        }
        else
        {
            FatalErrorInFunction
                << "Unknown cloud type: " << cloudType
                << ". Valid cloud type is carbonaceousCloud."
                << abort(FatalError);
        }

        cloudMap[cloudName] = {cloudType, idx};
    }
}

void CloudHandler::evolveClouds()
{
    evolve<basicCarbonaceousCloud>(carbonaceousCloudList);
}

template<class T>
void CloudHandler::evolve(std::vector<std::unique_ptr<T>>& cloudTypeList)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        cloud->evolve();
    }
}

void CloudHandler::writeClouds()
{
    write<basicCarbonaceousCloud>(carbonaceousCloudList);
}

template<class T>
void CloudHandler::write(std::vector<std::unique_ptr<T>>& cloudTypeList)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        cloud->write();
    }
}

tmp<DimensionedField<scalar, volMesh>> CloudHandler::SrhoClouds()
{
    tmp<volScalarField::Internal> trhoTrans(
        new volScalarField::Internal(
            IOobject(
                "rhoTransAll",
                mesh.thisDb().time().timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            mesh,
            dimensionedScalar(
                dimMass / dimTime / dimVolume, Zero)));

    scalarField& sourceField = trhoTrans.ref();

    Srho<basicCarbonaceousCloud>(carbonaceousCloudList, sourceField);

    sourceField /= mesh.thisDb().time().deltaTValue() * mesh.V();

    return trhoTrans;
}

template<class T>
void CloudHandler::Srho(std::vector<std::unique_ptr<T>>& cloudTypeList, scalarField& sourceField)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        if (cloud->solution().coupled())
        {
            PtrList<DimensionedField<scalar, volMesh>>& rhoTrans = cloud->rhoTrans();
            for (int i = 0; i < rhoTrans.size(); i++)
            {
                sourceField += rhoTrans[i];
            }
        }
    }
}

tmp<fvScalarMatrix> CloudHandler::SrhoClouds(volScalarField& rho)
{
    tmp<volScalarField::Internal> trhoTrans(
        new volScalarField::Internal(
            IOobject(
                "rhoTransAll",
                mesh.thisDb().time().timeName(),
                mesh.thisDb(),
                IOobject::NO_READ,
                IOobject::NO_WRITE,
                false),
            mesh,
            dimensionedScalar(
                dimMass / dimTime / dimVolume, Zero)));

    tmp<fvScalarMatrix> tmatrix(new fvScalarMatrix(rho, dimMass / dimTime));
    fvScalarMatrix& matrix = tmatrix.ref();

    Srho<basicCarbonaceousCloud>(carbonaceousCloudList, trhoTrans, matrix, rho);

    return tmatrix;
}

template<class T>
void CloudHandler::Srho(std::vector<std::unique_ptr<T>>& cloudTypeList, tmp<volScalarField::Internal>& trhoTrans, fvScalarMatrix& matrix, volScalarField& rho)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        if (cloud->solution().coupled())
        {
            scalarField& sourceField = trhoTrans.ref();
            sourceField *= 0;
            PtrList<DimensionedField<scalar, volMesh>>& rhoTrans = cloud->rhoTrans();

            for (int i = 0; i < rhoTrans.size(); i++)
            {
                sourceField += rhoTrans[i];
            }

            if (cloud->solution().semiImplicit("rho"))
            {
                sourceField /= mesh.thisDb().time().deltaTValue() * mesh.V();
                matrix += fvm::SuSp(trhoTrans() / rho, rho);
            }
            else
            {
                matrix.source() += -trhoTrans() / mesh.thisDb().time().deltaT();
            }
        }
    }
}

tmp<fvScalarMatrix> CloudHandler::SYiClouds(const label i, volScalarField& Yi)
{
    tmp<fvScalarMatrix> tmatrix(new fvScalarMatrix(Yi, dimMass / dimTime));
    fvScalarMatrix& matrix = tmatrix.ref();

    SYi<basicCarbonaceousCloud>(carbonaceousCloudList, matrix, Yi, i);

    return tmatrix;
}

template<class T>
void CloudHandler::SYi(std::vector<std::unique_ptr<T>>& cloudTypeList, fvScalarMatrix& matrix, volScalarField& Yi, const label i)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        matrix += cloud->SYi(i, Yi);
    }
}

tmp<fvScalarMatrix> CloudHandler::ShaClouds(volScalarField& ha)
{
    tmp<fvScalarMatrix> tmatrix(new fvScalarMatrix(ha, dimEnergy / dimTime));
    fvScalarMatrix& matrix = tmatrix.ref();

    Sha<basicCarbonaceousCloud>(carbonaceousCloudList, matrix, ha);

    return tmatrix;
}

template<class T>
void CloudHandler::Sha(std::vector<std::unique_ptr<T>>& cloudTypeList, fvScalarMatrix& matrix, volScalarField& ha)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        matrix += cloud->Sha(ha);
    }
}

tmp<Foam::fvVectorMatrix> CloudHandler::SUClouds(volVectorField& U, bool incompressible)
{
    dimensionSet dim(dimForce);
    if (incompressible)
    {
        dim.reset(dimForce / dimDensity);
    }

    tmp<fvVectorMatrix> tmatrix(new fvVectorMatrix(U, dim));
    fvVectorMatrix& matrix = tmatrix.ref();

    SU<basicCarbonaceousCloud>(carbonaceousCloudList, matrix, U, incompressible, dim);

    return tmatrix;
}

template<class T>
void CloudHandler::SU(std::vector<std::unique_ptr<T>>& cloudTypeList, fvVectorMatrix& matrix, volVectorField& U, bool incompressible, dimensionSet dim)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        if (cloud->solution().coupled())
        {
            if (cloud->solution().semiImplicit("U"))
            {
                volScalarField::Internal Vdt(mesh.V() * mesh.thisDb().time().deltaT());

                if (incompressible)
                {
                    Vdt.dimensions() *= dimDensity;
                }

                matrix += cloud->UTrans() / Vdt - fvm::Sp(cloud->UCoeff() / Vdt, U) + cloud->UCoeff() / Vdt * U;
            }
            else
            {
                matrix.source() += -cloud->UTrans() / (mesh.thisDb().time().deltaT());
            }
        }
    }
}

const basicCarbonaceousCloud& CloudHandler::getCarbonaceousCloud(word cloudName)
{
    if (!cloudMap.contains(cloudName))
    {
        Info << "Cloud not found. Solver is going to crash." << endl;
    }
    else
    {
        std::pair<std::string, int>& pair = cloudMap[cloudName];
        int idx = pair.second;
        return *carbonaceousCloudList[idx];
    }
}

word CloudHandler::getCloudType(word cloudName)
{
    return cloudMap[cloudName].first;
}

bool CloudHandler::writeData(Ostream& os) const
{
    return os.good();
}

void CloudHandler::deleteParticles(dictionary deleteDict)
{
    List<keyType> cloudKeyList = deleteDict.keys();
    for (keyType key : cloudKeyList)
    {
        word groupName(key);
        Info << groupName << ": ";
        dictionary actions(deleteDict.subDict(groupName));
        word cloudName(actions.getWord("cloudName"));
        actions.remove("cloudName");

        std::pair<std::string, int>& pair = cloudMap[cloudName];
        std::string cloudType = pair.first;
        int idx = pair.second;

        if (cloudType == "carbonaceousCloud")
        {
            deleteParticlesCloudType<basicCarbonaceousCloud>(carbonaceousCloudList, idx, actions, cloudName);
        }
        else
        {
            FatalErrorInFunction
                << "Unknown cloud type for particle deletion: " << cloudType
                << abort(FatalError);
        }
    }
}

template<class T>
void CloudHandler::deleteParticlesCloudType(std::vector<std::unique_ptr<T>>& cloudTypeList, int idx, dictionary actions, word cloudName)
{
    List<keyType> actionsKeyList = actions.keys();
    int n_actions = actionsKeyList.size();
    std::vector<word> actionNames(n_actions);
    std::vector<word> field(n_actions);
    std::vector<word> op(n_actions);
    std::vector<scalar> value(n_actions);

    for (int i = 0; i < n_actions; i++)
    {
        dictionary& action(actions.subDict(actionsKeyList[i]));
        actionNames[i] = action.getWord("action");
        field[i] = action.getWord("field");
        op[i] = action.getWord("operator");
        value[i] = action.getScalar("value");
    }

    T& cloud = *cloudTypeList[idx];
    using ParticleType = decltype(*std::begin(cloud));

    std::unordered_map<std::string, std::function<double(ParticleType&)>> getParticleValue;

    getParticleValue["x"] = [](ParticleType& p) { return p.position()[0]; };
    getParticleValue["y"] = [](ParticleType& p) { return p.position()[1]; };
    getParticleValue["z"] = [](ParticleType& p) { return p.position()[2]; };
    getParticleValue["rz"] = [](ParticleType& p)
    { return std::sqrt(std::pow(p.position()[0], 2) + std::pow(p.position()[1], 2)); };
    getParticleValue["Uz"] = [](ParticleType& p) { return p.U()[2]; };
    getParticleValue["Umag"] = [](ParticleType& p) { return mag(p.U()); };
    getParticleValue["T"] = [](ParticleType& p) { return p.T(); };

    int j = 0;
    int oldTot(returnReduce(cloud.nParcels(), sumOp<label>()));

    std::unordered_map<std::string, bool> addMap;
    addMap["add"] = true;
    addMap["subtract"] = false;

    std::unordered_map<std::string, std::function<void(bool&, double, double&, bool&)>> operatorMap;
    operatorMap["greater"] = [](bool& selection, double particleValue, double& value, bool& r)
    { selection = (particleValue > value) ? r : selection; };
    operatorMap["less"] = [](bool& selection, double particleValue, double& value, bool& r)
    { selection = (particleValue < value) ? r : selection; };

    for (auto& p : cloud)
    {
        bool selection = false;

        for (label i = 0; i < n_actions; i++)
        {
            bool r = addMap[actionNames[i]];
            operatorMap[op[i]](selection, getParticleValue[field[i]](p), value[i], r);
        }

        if (selection)
        {
            cloud.deleteParticle(p);
            ++j;
        }
    }
    Info << "Reduced number of particles in " << cloudName << "  by "
         << returnReduce(j, sumOp<int>()) << ", from " << oldTot
         << " to " << returnReduce(cloud.nParcels(), sumOp<label>()) << "." << endl;
}

void CloudHandler::particleProcessorDistributionInfo()
{
    Info << "--- Cloud processor distribution info ---" << endl;
    particleProcessorDistributionInfoCloudType<basicCarbonaceousCloud>(carbonaceousCloudList);
    Info << "------" << endl;
}

template<class T>
void CloudHandler::particleProcessorDistributionInfoCloudType(std::vector<std::unique_ptr<T>>& cloudTypeList)
{
    for (std::unique_ptr<T>& cloud : cloudTypeList)
    {
        int total(returnReduce(cloud->nParcels(), sumOp<label>()));
        Info << "Total number of particles: " << total << endl;
        labelList parcelsProc(UPstream::listGatherValues(cloud->nParcels()));
        Info << "Homogeneous distribution: "
             << total / (parcelsProc.size() + 1e-16)
             << setprecision(IOstream::defaultPrecision()) << ", "
             << 100 / (parcelsProc.size() + 1e-16) << setprecision(3) << "%" << endl;
        Info << "Min number of particles: " << min(parcelsProc) << " on processor " << findMin(parcelsProc) << endl;
        Info << "Max number of particles: " << max(parcelsProc) << " on processor " << findMax(parcelsProc) << endl;
        for (int i = 0; i < parcelsProc.size(); i++)
        {
            Info << "Processor " << i << ": "
                 << setprecision(IOstream::defaultPrecision()) << parcelsProc[i] << ", "
                 << setprecision(3) << parcelsProc[i] / (total + 1e-16) * 100 << "%" << endl;
        }
    }
    Info << setprecision(IOstream::defaultPrecision());
}
