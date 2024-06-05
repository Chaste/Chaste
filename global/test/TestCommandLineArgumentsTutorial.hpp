/*

Copyright (c) 2005-2024, University of Oxford.
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

#ifndef TESTCOMMANDLINEARGUMENTSTUTORIAL_HPP_
#define TESTCOMMANDLINEARGUMENTSTUTORIAL_HPP_

#include <cxxtest/TestSuite.h>
/* Most Chaste code uses PETSc to solve linear algebra problems.  This involves starting PETSc at the beginning of a test-suite
 * and closing it at the end.  (If you never run code in parallel then it is safe to replace PetscSetupAndFinalize.hpp with FakePetscSetup.hpp)
 */
#include "PetscSetupAndFinalize.hpp"
#include "CommandLineArguments.hpp"

/**
 * @file
 *
 * This is a tutorial for providing command line arguements to a test from a bash script.
 * The example bash script also includes a for loop as an example of how such input could
 * be used to run multiple tests with varying parameter inputs.
 */

 // NOTE: This test will not work if directlly executed from terminal due to requiring command line arguements.
 // Investiagte the accompanying runcommandlinetutorial.sh bash script to see how this tutorial should be executed.
class TestCommandLineArgumentsTutorial : public CxxTest::TestSuite
{
public:
    void TestCommandLineTutorial()
    {
        // Here we utilise CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption() to take in our command line arguements.
        int outp1 = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt1");
        int outp2 = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt2");
        int outp3 = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-opt3");

        int sum = outp1 + outp2 + outp3;
        
        std::cout << "When we add "<< outp1 << " ,"<<  outp2 << " and "<< outp3 <<" we get " << sum <<"\n";
    }
};

#endif /*TESTHELLO_ECADTURNOVERMODEL_HPP_*/
