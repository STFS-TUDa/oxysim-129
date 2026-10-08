/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2011-2016 OpenFOAM Foundation
    Copyright (C) 2019-2020 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Class
    CloudHandler

Group
    src

Description

The CloudHandler manages carbonaceousCloud instances. It evolves them each
time step and accumulates exchange source terms for the solver equations.

Multiple clouds of the same type can be combined:

clouds
{
        biomassCloud    carbonaceousCloud;
        coalCloud       carbonaceousCloud;
}

The CloudHandler can also delete particles based on position, temperature,
or velocity criteria (see deleteParticlesDict).

\*---------------------------------------------------------------------------*/

#pragma once
#include "fvCFD.H"
#include "fvOptions.H"
#include <string>
#include "basicCarbonaceousCloud.H"
#include "basicReactingCloudSTFS.H"
#include "SLGThermo.H"

class CloudHandler
: public regIOobject
{

    std::vector<std::unique_ptr<basicCarbonaceousCloud>> carbonaceousCloudList;
    fvMesh& mesh;
    std::map<std::string, std::pair<std::string, int>> cloudMap;

    template<class T>
    void evolve(std::vector<std::unique_ptr<T>>& cloudTypeList);

    template<class T>
    void write(std::vector<std::unique_ptr<T>>& cloudTypeList);

    template<class T>
    void Srho(std::vector<std::unique_ptr<T>>& cloudTypeList, scalarField& sourceField);

    template<class T>
    void Srho(std::vector<std::unique_ptr<T>>& cloudTypeList, tmp<volScalarField::Internal>& trhoTrans, fvScalarMatrix& matrix, volScalarField& rho);

    template<class T>
    void SYi(std::vector<std::unique_ptr<T>>& cloudTypeList, fvScalarMatrix& matrix, volScalarField& Yi, const label i);

    template<class T>
    void Sha(std::vector<std::unique_ptr<T>>& cloudTypeList, fvScalarMatrix& matrix, volScalarField& ha);

    template<class T>
    void SU(std::vector<std::unique_ptr<T>>& cloudTypeList, fvVectorMatrix& matrix, volVectorField& U, bool incompressible, dimensionSet dim);

    template<class T>
    void deleteParticlesCloudType(std::vector<std::unique_ptr<T>>& cloudList, int idx, dictionary actions, word cloudName);

    template<class T>
    void particleProcessorDistributionInfoCloudType(std::vector<std::unique_ptr<T>>& cloudTypeList);

public:
    TypeName("CloudHandler");
    // Constructor
    // CloudHandler();
    CloudHandler(
        fvMesh& i_mesh,
        dictionary cloudDict,
        const volScalarField& rho,
        const volVectorField& U,
        const dimensionedVector& g,
        const SLGThermo& slgThermo);

    void evolveClouds();

    void writeClouds();

    tmp<DimensionedField<scalar, volMesh>> SrhoClouds();

    tmp<fvScalarMatrix> SrhoClouds(volScalarField& rho);

    tmp<fvScalarMatrix> SYiClouds(const label i, volScalarField& Yi);

    tmp<fvScalarMatrix> ShaClouds(volScalarField& ha);

    tmp<Foam::fvVectorMatrix> SUClouds(volVectorField& U, bool incompressible);

    const basicCarbonaceousCloud& getCarbonaceousCloud(word cloudName);

    word getCloudType(word cloudName);

    bool writeData(Ostream& os) const;

    void deleteParticles(dictionary dict);

    void particleProcessorDistributionInfo();
};
