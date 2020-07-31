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

//#include <array>
//#include <cstdio>
//#include <iostream>
//#include <memory>
//#include <stdexcept>
//#include <string>

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
class TestCodegenPresent : public CxxTest::TestSuite
{
public:
    void TestChaste_codegenInstalled()
    {
        FileFinder codegen("codegen_python3_venv/bin/chaste_codegen", RelativeTo::ChasteSourceRoot);
        FileFinder python("codegen_python3_venv/bin/python", RelativeTo::ChasteSourceRoot);

        TS_ASSERT_EQUALS(codegen.IsFile(), true);
        TS_ASSERT_EQUALS(python.IsFile(), true);
    }

};

#endif //_TESTMETADATASUBMODULE_HPP_
