/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */

#include <cxxtest/TestSuite.h>
#include "CheckpointArchiveTypes.hpp"

#include <iostream>
#include <iomanip>

#include <boost/lexical_cast.hpp>

#include "SimulationTime.hpp"
#include "RandomNumberGenerator.hpp"


#include "Debug.hpp"
#include "ExecutableSupport.hpp"

/*
 * These headers handle passing parameters to the executable.
 *
 * We need to do the following to the hostconfig file:
 *  + boost_libs = ['boost_serialization', 'boost_filesystem', 'boost_boost::program_options']
 *  - boost_libs = ['boost_serialization', 'boost_filesystem']
 */
#include <boost/program_options/options_description.hpp>
#include <boost/program_options/variables_map.hpp>
#include <boost/program_options/parsers.hpp>

void SetupSingletons();
void DestroySingletons();
void SetupAndRunSimulation(unsigned simulation_id);

int main(int argc, char *argv[])
{
    // This sets up PETSc and prints out copyright information, etc.
    ExecutableSupport::StandardStartup(&argc, &argv);

    // Define command line options
    boost::program_options::options_description general_options("This is a  sample chaste executable .\n\n");
    general_options.add_options()
                    ("help", "produce help message")
                    ("ID", boost::program_options::value<unsigned>()->default_value(0),"ID of the simulation (for output)");

    // define parse command line into variables_map
    boost::program_options::variables_map variables_map;
    boost::program_options::store(parse_command_line(argc, argv, general_options), variables_map);

    // print help message if wanted
    if (variables_map.count("help"))
    {
//        std::cout << setprecision(3) << general_options << "\n";
        std::cout << general_options << "\n";
        return 1;
    }

    // get id and name from command line
    unsigned simulation_id=variables_map["ID"].as<unsigned>();


    SetupSingletons();
    SetupAndRunSimulation(simulation_id);
    DestroySingletons();
}

void SetupSingletons()
{
    // Set up what the test suite would do
    SimulationTime::Instance()->SetStartTime(0.0);

    // we want every realisation of the simulation to be different!
    RandomNumberGenerator::Instance()->Reseed(time(NULL));
//    CellPropertyRegistry::Instance()->Clear();
//    CellId::ResetMaxCellId();
}

void DestroySingletons()
{
    // this is from the tearDown method of the test suite
    SimulationTime::Destroy();
    RandomNumberGenerator::Destroy();
//    CellPropertyRegistry::Instance()->Clear();
}

void SetupAndRunSimulation(unsigned simulation_id)
{
    PRINT_VARIABLE(simulation_id);
}
