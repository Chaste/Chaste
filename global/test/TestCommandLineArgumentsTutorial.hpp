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
