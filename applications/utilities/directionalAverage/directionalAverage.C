/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | foam-extend: Open Source CFD
   \\    /   O peration     | Version:     4.0
    \\  /    A nd           | Web:         http://www.foam-extend.org
     \\/     M anipulation  | For copyright notice see file Copyright
-------------------------------------------------------------------------------
License
    This file is part of foam-extend.

    foam-extend is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    foam-extend is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with foam-extend.  If not, see <http://www.gnu.org/licenses/>.

Application
    Average a axial symmetric domain in circumferential direction

Description
    For each time:
    The velocity is displayed: (U_r,U_theat,U_axis)

\*---------------------------------------------------------------------------*/


#include "mathematicalConstants.H"
#include "fvCFD.H"
#include <vector>
#include "HashSet.H"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

int main(int argc, char* argv[])
{ // main paranthesis
    // adds the option of the -constant and the -noZero flag
    timeSelector::addOptions();
    // adds the fields flag which needs to be defined
    argList::validOptions.insert("fields", "fields");
    argList::validOptions.insert("cellsize", "Specify a cell size as a double value");

#include "setRootCase.H"
#include "createTime.H"

    instantList timeDirs = timeSelector::select0(runTime, args);

#include "createMesh.H"

    HashSet<wordRe> selectedFields;
    scalar cellSize = 0.0; // Default value

    // this gets the names of the current fields
    if (args.found("fields"))
    {
        selectedFields = args.getList<wordRe>("fields");
    }
    // this checks if fields were specified
    if (!selectedFields.size())
    {
        Info << "Specify the fields to average in x and z direction with the option -fields " << endl;
    }

    if (args.found("cellsize"))
    {
        cellSize = readScalar(args["cellsize"]);
    }
    else
    {
        Info << "Specify the cellsize of the equidistant grid with the option -cellsize " << endl;
        FatalErrorInFunction << "Invalid cell size! Must be > 0." << exit(FatalError);
    }
    Info << "Using cell size: " << cellSize << endl;


    // define directions in which to average in -> dont avergage in y (1) direction, but in x and z
    int averagingIndices[2] = {0, 2};
    int independentCord = 1;

    // int nCells= mesh.nCells();
    Foam::Vector<int> nCells = mesh.solutionD(); // Get mesh dimensions

    vector minCoords = gMin(mesh.C()); // Get minimum (x_min, y_min, z_min)
    vector maxCoords = gMax(mesh.C()); // Get maximum (x_max, y_max, z_max)


    const double dy = cellSize;                                     //(maxCoords.y()-minCoords.y())/ny;
    const int ny = (int)((maxCoords.y() - minCoords.y()) / dy) + 1; // Number of cells in y-direction
    Info << "Assuming structured grid. This grid has " << ny << " cells in y-direction " << endl;

    // Iterating over all time vals
    forAll(timeDirs, timeI)
    { // time loop
        runTime.setTime(timeDirs[timeI], timeI);

        Info << "Time = " << runTime.timeName() << endl;

        // iterate over all fields
        for (HashSet<wordRe>::iterator iter = selectedFields.begin(); iter != selectedFields.end(); iter++)
        { // fields loop paranthesis

            // check if its not a vector field
            if (selectedFields.found(iter.key()) && iter.key() != "U" && iter.key() != "UMean" && iter.key() != "vorticity" && iter.key() != "UPrime")
            {

                Info << "Averaging scalar field  " << iter.key() << endl;
                // sample Field is called Z here
                volScalarField ZMean(
                    IOobject(
                        iter.key(),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE),
                    mesh);

                volScalarField ZMeanAve(
                    IOobject(
                        iter.key() + "AveXZ",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    ZMean);

                typedef std::vector<double> single_vector;
                single_vector weight(ny);
                single_vector averagedVal(ny);

                // loop for averaging
                forAll(ZMean, cellI)
                {                                             // loop over all cells
                    scalar currentyPos = mesh.C()[cellI].y(); // Get the y-coordinate of the cell center
                    // get the vector position on where to put this in
                    // assumption: grid starts at 0
                    int yVectorPos = (int)((currentyPos / dy) - 0.5); // current y position in vector, e.g. lowest cell center has position 0
                    weight[yVectorPos] += 1;

                    averagedVal[yVectorPos] = 1. / weight[yVectorPos] * (ZMean[cellI] + (weight[yVectorPos] - 1) * averagedVal[yVectorPos]);
                }
                // loop for writing average
                forAll(ZMean, cellI)
                {

                    scalar currentyPos = mesh.C()[cellI].y(); // Get the y-coordinate of the cell center
                    // get the vector position on where to put this in
                    // assumption: grid starts at 0
                    int yVectorPos = (int)((currentyPos / dy) - 0.5); // current y position in vector, e.g. lowest cell center has positio

                    ZMeanAve[cellI] = max(averagedVal[yVectorPos], 1.0E-20);
                }

                ZMeanAve.write();
                Info << "Writing field  " << iter.key() + "AveXZ field" << endl;
            } // paranthesis vector field
            // Averaging vector fields
            else if (selectedFields.found(iter.key()) && (iter.key() == "U" || iter.key() == "UMean" || iter.key() == "vorticity" || iter.key() == "UPrime"))
            {
                Info << "Averaging vector field  " << iter.key() << endl;
                volVectorField UMean(
                    IOobject(
                        iter.key(),
                        runTime.timeName(),
                        mesh,
                        IOobject::MUST_READ,
                        IOobject::NO_WRITE),
                    mesh);

                volVectorField UMeanAve(
                    IOobject(
                        iter.key() + "AveXZ",
                        runTime.timeName(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    UMean);

                typedef std::vector<double> single_vector;
                typedef std::vector<single_vector> double_vector;
                single_vector weight(ny);
                double_vector U(ny, single_vector(3));

                // calc new Values
                forAll(UMean, cellI)
                {

                    scalar currentyPos = mesh.C()[cellI].y(); // Get the y-coordinate of the cell center
                    // get the vector position on where to put this in
                    // assumption: grid starts at 0
                    int yVectorPos = (int)((currentyPos / dy) - 0.5); // current y position in vector, e.g. lowest cell center has positio
                    weight[yVectorPos] += 1;

                    U[yVectorPos][0] = 1. / weight[yVectorPos] * (UMean[cellI].x() + (weight[yVectorPos] - 1) * U[yVectorPos][0]);
                    U[yVectorPos][1] = 1. / weight[yVectorPos] * (UMean[cellI].y() + (weight[yVectorPos] - 1) * U[yVectorPos][1]);
                    U[yVectorPos][2] = 1. / weight[yVectorPos] * (UMean[cellI].z() + (weight[yVectorPos] - 1) * U[yVectorPos][2]);
                }
                // write new vals
                forAll(UMean, cellI)
                {
                    scalar currentyPos = mesh.C()[cellI].y(); // Get the y-coordinate of the cell center
                    // get the vector position on where to put this in
                    // assumption: grid starts at 0
                    int yVectorPos = (int)((currentyPos / dy) - 0.5); // current y position in vector, e.g. lowest cell center has positio
                    UMeanAve[cellI].z() = U[yVectorPos][2];
                    UMeanAve[cellI].x() = U[yVectorPos][0];
                    UMeanAve[cellI].y() = U[yVectorPos][1];
                }

                UMeanAve.write();
                Info << "Writing field  " << iter.key() + "AveXZ field" << endl;
            } // velocity fields
        } // loop over all fields
    } // loop over all time vals
} // main paranthesis
