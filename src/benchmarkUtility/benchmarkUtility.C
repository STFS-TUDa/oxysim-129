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

#include "benchmarkUtility.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

// * * * * * * * * * * * * * * * *  Constructors  * * * * * * * * * * * * * * //

Foam::benchmarkUtility::benchmarkUtility(const Time& Time)
: runTime_(Time),
  start_data(),
  sum_data(),
  print_interval(runTime_.controlDict().lookupOrDefault<label>("benchmarkInterval", 0)),
  print_every_time_step(runTime_.controlDict().lookupOrDefault<bool>("benchmarkPrintEveryTimeStep", true))
{
    if (print_interval > 0)
    {
        Info << "Benchmarking with an averaging interval of " << print_interval << endl;
    }
}

// * * * * * * * * * * * * * * * *  Member Functions  * * * * * * * * * * * * * * //

void Foam::benchmarkUtility::startBenchmark(const std::string& Key)
{
    if (print_interval == 0)
    {
        return;
    }
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time = std::chrono::high_resolution_clock::now();
    // store the start time in a map
    static_cast<void>(start_data.insert_or_assign(Key, start_time));
}

void Foam::benchmarkUtility::stopBenchmark(const std::string& Key)
{
    if (print_interval == 0)
    {
        return;
    }
    // check if timer has been started for this quantity
    if (!start_data.contains(Key))
    {
        FatalErrorIn("Foam::benchmarkUtility::stopBenchmark(const std::string& Key)") << "Benchmarking for Key '" << Key << "' not started. Please start the benchmark for the key first before trying to access the end time.";
    }
    // get stored start time
    std::chrono::time_point<std::chrono::high_resolution_clock> start = start_data[Key];
    // get end time and calculate duration
    std::chrono::time_point<std::chrono::high_resolution_clock> end = std::chrono::high_resolution_clock::now();
    // getting duration, casting it to milliseconds dividing and then to float -> duration is float in milliseconds
    float duration = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();

    // get max of all processors
    reduce(duration, maxOp<float>());

    if (sum_data.contains(Key))
    {
        // sum and store duration in map (for some reason += does not work)
        sum_data[Key] = sum_data[Key] + duration;
    }
    else
    {
        // no data for summing available -> just inserting
        static_cast<void>(sum_data.insert_or_assign(Key, duration));
    }
    // print every time step if wanted
    if (print_every_time_step)
    {
        Info << "TIME: Duration of " << Key << ": " << duration << " [ms]" << endl;
    }
}

void Foam::benchmarkUtility::printAndResetBenchmark()
{
    if (print_interval == 0)
    {
        return;
    }
    if (sum_data.empty())
    {
        Info << "\nNo averaged benchmarking data found." << endl;
        return;
    }

    if (runTime_.timeIndex() % print_interval == 0)
    {
        Info << "TIME: These are the benchmarking data averaged over " << print_interval << " timesteps.\n";
        for (const auto& [key_, val_] : sum_data)
        {
            Info << "TIME: Averaged duration of " << key_ << ": " << (val_ / static_cast<float>(print_interval)) << " [ms]" << endl; //
            // clean map
            sum_data[key_] = 0.0;
        }
        Info << "TIME: Resetting average benchmarking data.\n"
             << endl;
    }
}