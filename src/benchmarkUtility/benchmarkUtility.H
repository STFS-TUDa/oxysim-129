/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | www.openfoam.com
     \\/     M anipulation  |
-------------------------------------------------------------------------------
    Copyright (C) 2013-2016 OpenFOAM Foundation
    Copyright (C) 2019 OpenCFD Ltd.
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
    Foam::benchmarkUtility

Group
    basic

Description
    Class used for benchmarking of code parts.

Authors
     Vinzenz Schuh

SourceFiles
    benchmarkUtility.C

\*---------------------------------------------------------------------------*/

#pragma once

#include <map>
#include <string>
#include <chrono>

#include "fvCFD.H"

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class basicThermoSTFS Declaration
\*---------------------------------------------------------------------------*/
/** @defgroup benchmarkUtility
 *   @ingroup basic
 *  @{
 */

class benchmarkUtility
{
private:
    // Private data

    const Time& runTime_;

    std::map<std::string, std::chrono::time_point<std::chrono::high_resolution_clock>> start_data;
    std::map<std::string, float> sum_data;

    const label print_interval;
    const bool print_every_time_step;

public:
    // Constructors

    //- Construct from mesh and phase name
    benchmarkUtility(const Time&);

    // Member Functions
    void startBenchmark(const std::string&);

    void stopBenchmark(const std::string&);

    void printAndResetBenchmark();
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/**@}*/
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
