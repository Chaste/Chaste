/*

Copyright (c) 2005-2020, University of Oxford.
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

#ifndef _TESTMETADATASUBMODULE_HPP_
#define _TESTMETADATASUBMODULE_HPP_

#include <cxxtest/TestSuite.h>

#include <array>
#include <cstdio>
#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>

#include <boost/algorithm/string.hpp>

#include "FileFinder.hpp"

//This test is always run sequentially (never in parallel)
#include "FakePetscSetup.hpp"

/**
 * There is a git submodule in python/pycml/ontologies
 * which is shared with other projects, and is now imported into
 * the Chaste source via a git submodule instead of living here.
 * 
 * This test suite makes sure that the submodule has been initialised and is up to date.
 * 
 * The second test should be updated by the person who decides that the
 * Chaste copy of the metadata should be updated to match the remote ontology
 * at  https://github.com/ModellingWebLab/ontologies/
 * 
 * This update is done manually, make sure you are on the develop branch and do:
 * 
 * cd $CHASTE_SRC/python/pycml/ontologies
 * git pull origin master
 * cd $CHASTE_SRC
 * git add python/pycml/ontologies
 * git commit -m "Update CellML Metadata ontology to latest remote version."
 * git push
 * 
 * Then type 
 * git submodule
 * and copy the commit hash into the member variable in the below test. If anyone runs the
 * latest version of the code, then it will fail to remind them to do a `git submodule update`.
 * 
 */
class TestMetadataSubmodule : public CxxTest::TestSuite
{
    std::string GetOutputOfCommand(const char* cmd)
    {
        std::array<char, 256> buffer;
        std::string result;
        std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
        if (!pipe)
        {
            throw std::runtime_error("popen() failed!");
        }
        while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr)
        {
            result += buffer.data();
        }
        return result;
    }

public:
    void TestForPresenceOfTheOntologyFiles()
    {
        std::string base_filename = "oxford-metadata";
        FileFinder ontology_ttl("python/pycml/ontologies/" + base_filename + ".ttl", RelativeTo::ChasteSourceRoot);
        FileFinder ontology_rdf("python/pycml/ontologies/" + base_filename + ".rdf", RelativeTo::ChasteSourceRoot);

        bool ttl_present = ontology_ttl.IsFile();
        bool rdf_present = ontology_rdf.IsFile();

        TS_ASSERT_EQUALS(ttl_present, true);
        TS_ASSERT_EQUALS(rdf_present, true);

        if (!(ttl_present && rdf_present))
        {
            std::cout << "This test has failed because the ontology files in the Chaste source folder python/pycml/ontologies are not present."
                         "These are provided by a git submodule.\n"
                         "Please ensure that you initialise the submodule by typing :\n"
                         "'git submodule update --init'\n"
                         "into the terminal in the Chaste source directory."
                      << std::endl;
        }
    }

    void TestForPendingUpdatesToSubmodule()
    {
        std::string latest_commit_hash = "e5376c254ae3ea1642fa5cc3838f150ebd2b7b29";

        FileFinder chaste_source("", RelativeTo::ChasteSourceRoot);
        std::stringstream command;
        command << "git --git-dir=" << chaste_source.GetAbsolutePath() << ".git submodule";
        std::string command_string = command.str();
        //std::cout << command_string << "\n" << output << std::endl;

        std::string output = GetOutputOfCommand(command_string.c_str());

        boost::trim_left(output);

        unsigned first_space_idx = output.find(" ");
        std::string current_commit_hash = output.substr(0, first_space_idx);
        std::string rest_of_output = output.substr(first_space_idx, output.length());
        boost::trim_left(rest_of_output);
        std::string submodule_foldername = rest_of_output.substr(0, rest_of_output.find(" "));

        TS_ASSERT_EQUALS(submodule_foldername, "python/pycml/ontologies");
        TS_ASSERT_EQUALS(current_commit_hash, latest_commit_hash);

        if (current_commit_hash != latest_commit_hash)
        {
            std::cout << "This test has failed because the submodule containing"
                         " the CellML metadata ontology is out of date.\n"
                         "Please run the command:\n"
                         "'git submodule update'\n"
                         "in the terminal in the Chaste source directory."
                      << std::endl;
        }
    }
};

#endif //_TESTMETADATASUBMODULE_HPP_
