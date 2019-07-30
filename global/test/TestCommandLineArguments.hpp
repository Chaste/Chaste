/*

Copyright (c) 2005-2019, University of Oxford.
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

#ifndef TESTCOMMANDLINEARGUMENTS_HPP_
#define TESTCOMMANDLINEARGUMENTS_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <cassert>
#include "CommandLineArguments.hpp"
#include "CommandLineArgumentsMocker.hpp"
#include "PetscSetupAndFinalize.hpp"

/* HOW_TO_TAG General
 * Read and use parameters from the command line
 *
 * If you want to use parameters that are supplied in the command line, then
 *  (i) add lines such as "double x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-myparam");" below
 *  (ii) compile but do not run the test (see ChasteGuides/RunningBinariesFromCommandLine)
 *  (iii) run the compiled executable from the command line (see ChasteGuides/RunningBinariesFromCommandLine), with your parameter.
 *        If, at this step, you "undefined symbol:" errors then set your LD_LIBRARY_PATH (see ChasteGuides/RunningBinariesFromCommandLine)
 *
 *
 *
 *  For example:
 *  scons co=1 projects/you/TestBlah.hpp
 *  ./projects/you/build/debug/TestBlahRunner -myparam 10.4
 *
 *  Alternatively, you can add the arguments to SCons so that you can compile and run in one go:
 *  scons run_time_flags="--verbose true" global/test/TestCommandLineArguments.hpp
 *  This should produce "You have successfully set --verbose to take the value 1." for this test suite.
 *
 * Note: error messages such as
 *   WARNING! There are options you set that were not used!
 *   WARNING! could be spelling mistake, etc!
 * are due to PETSc thinking the parameter must have been for it.
 *
 */
class TestCommandLineArguments : public CxxTest::TestSuite
{
public:

    void TestCommandLineArgumentsSingleton()
    {
        // Test that argc and argv are populated
        int argc = *(CommandLineArguments::Instance()->p_argc);
        TS_ASSERT_LESS_THAN(0, argc); // argc should always be 1 or greater

        /* argv[0] will be equal to global/build/debug/TestCommandLineArguments
         * or global/build/optimised/TestCommandLineArguments, etc
         * Variations on Windows (CMake) include .../TestCommandLineArguments.exe
         * Test executables built with SCons (deprecated) have the word Runner
         * * .../TestCommandLineArgumentsRunner
         * * .../TestCommandLineArgumentsRunner.exe
         */
        char** argv = *(CommandLineArguments::Instance()->p_argv);
        assert(argv != NULL);
        std::string arg_as_string(argv[0]);
        size_t pos = arg_as_string.find("TestCommandLineArguments");
        // If TestCommandLineArguments is not a substring of the commandline then pos==std::string::npos
        TS_ASSERT_DIFFERS(pos, std::string::npos);

        // Now test OptionExists() and GetValueCorrespondingToOption()
        //
        // The following tests would require the following arguments to be passed
        // in:
        // ./global/build/debug/TestCommandLineArguments -myoption -myintval 24 -mydoubleval 3.14 -3.14 -m2intval -42 -mystrings Baboons Monkeys Gibbons -mystring more_baboons
        //
        // To test the methods we overwrite the arg_c and arg_v contained in the
        // singleton with the arguments that were needed.
        int new_argc = 25;
        char** new_argv = new char*[25];

        char new_argv0[] = "..";
        char new_argv1[] = "-myoption";
        char new_argv2[] = "-myintval";
        char new_argv3[] = "24";
        char new_argv4[] = "-mydoubleval";
        char new_argv5[] = "3.14";
        char new_argv6[] = "-3.14";
        char new_argv7[] = "-m2intval";
        char new_argv8[] = "-42";
        char new_argv9[] = "-mystrings";
        char new_argv10[] = "Baboons";
        char new_argv11[] = "Monkeys";
        char new_argv12[] = "Gibbons";
        char new_argv13[] = "-mystring";
        char new_argv14[] = "more_baboons";
        char new_argv15[] = "-mybool";
        char new_argv16[] = "0";
        char new_argv17[] = "-mybool2";
        char new_argv18[] = "fAlSe";
        char new_argv19[] = "-mybool3";
        char new_argv20[] = "1";
        char new_argv21[] = "-mybool4";
        char new_argv22[] = "TrUe";
        char new_argv23[] = "-mybool5";
        char new_argv24[] = "doDgy";

        new_argv[0] = new_argv0;
        new_argv[1] = new_argv1;
        new_argv[2] = new_argv2;
        new_argv[3] = new_argv3;
        new_argv[4] = new_argv4;
        new_argv[5] = new_argv5;
        new_argv[6] = new_argv6;
        new_argv[7] = new_argv7;
        new_argv[8] = new_argv8;
        new_argv[9] = new_argv9;
        new_argv[10] = new_argv10;
        new_argv[11] = new_argv11;
        new_argv[12] = new_argv12;
        new_argv[13] = new_argv13;
        new_argv[14] = new_argv14;
        new_argv[15] = new_argv15;
        new_argv[16] = new_argv16;
        new_argv[17] = new_argv17;
        new_argv[18] = new_argv18;
        new_argv[19] = new_argv19;
        new_argv[20] = new_argv20;
        new_argv[21] = new_argv21;
        new_argv[22] = new_argv22;
        new_argv[23] = new_argv23;
        new_argv[24] = new_argv24;

        // Save the real args to be restored at the end
        int* p_real_argc = CommandLineArguments::Instance()->p_argc;
        char*** p_real_argv = CommandLineArguments::Instance()->p_argv;

        // Overwrite the args
        CommandLineArguments::Instance()->p_argc = &new_argc;
        CommandLineArguments::Instance()->p_argv = &new_argv;

        // Test OptionExists()
        TS_ASSERT_EQUALS(CommandLineArguments::Instance()->OptionExists("-myoption"), true);

        TS_ASSERT_EQUALS(CommandLineArguments::Instance()->OptionExists("-asddsgijdfgokgfgurgher"), false);

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->OptionExists("-42"),
                "A command line option must begin with '-' followed by a non-numeric character.");

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetStringsCorrespondingToOption("-myoption"),
                "No value(s) given after command line option '-myoption'");

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetStringsCorrespondingToOption("-mynonsense"),
                "Command line option '-mynonsense' does not exist");

        // Test GetNumberOfArgumentsForOption()
        TS_ASSERT_EQUALS(CommandLineArguments::Instance()->GetNumberOfArgumentsForOption("-myoption", false), 0);

        // Test GetValueCorrespondingToOption()
        char* val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("-myintval");
        unsigned i = atol(val);
        TS_ASSERT_EQUALS(i, 24u);

        val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("-m2intval");
        int j = atol(val);
        TS_ASSERT_EQUALS(j, -42);

        j = CommandLineArguments::Instance()->GetIntCorrespondingToOption("-m2intval");
        TS_ASSERT_EQUALS(j, -42);

        TS_ASSERT_THROWS_THIS(i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-m2intval"),
                              "Option is a negative number and cannot be converted to unsigned.");

        i = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-myintval");
        TS_ASSERT_EQUALS(i, 24u);

        val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mydoubleval");
        double x = atof(val);
        TS_ASSERT_EQUALS(x, 3.14);

        x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-mydoubleval");
        TS_ASSERT_EQUALS(x, 3.14);

        // Test exceptions in GetValueCorrespondingToOption()
        TS_ASSERT_THROWS_CONTAINS(CommandLineArguments::Instance()->GetValueCorrespondingToOption("-rwesdb"), "does not exist");


        val = CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mystring");
        std::string argument(val);
        TS_ASSERT_EQUALS(argument, "more_baboons");

        std::string string_argument = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-mystring");
        TS_ASSERT_EQUALS(string_argument, "more_baboons");

        std::vector<std::string> string_arguments = CommandLineArguments::Instance()->GetStringsCorrespondingToOption("-mystrings");
        TS_ASSERT_EQUALS(string_arguments.size(), 3u);
        TS_ASSERT_EQUALS(string_arguments[0], "Baboons");
        TS_ASSERT_EQUALS(string_arguments[1], "Monkeys");
        TS_ASSERT_EQUALS(string_arguments[2], "Gibbons");
        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mystrings",4),
                              "Index=4 requested for '-mystrings', but only 3 given.");

        std::vector<double> double_arguments = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption("-mydoubleval");
        TS_ASSERT_EQUALS(double_arguments.size(), 2u);
        TS_ASSERT_DELTA(double_arguments[0], 3.14, 1e-6);
        TS_ASSERT_DELTA(double_arguments[1], -3.14, 1e-6);

        std::vector<unsigned> unsigned_args = CommandLineArguments::Instance()->GetUnsignedsCorrespondingToOption("-myintval");
        TS_ASSERT_EQUALS(unsigned_args.size(), 1u);
        TS_ASSERT_EQUALS(unsigned_args[0],24u);

        std::vector<int> int_args = CommandLineArguments::Instance()->GetIntsCorrespondingToOption("-m2intval");
        TS_ASSERT_EQUALS(int_args.size(), 1u);
        TS_ASSERT_EQUALS(int_args[0],-42);

        // Check the bool methods
        // (NB alternatively just check for the presence of an option,
        // but this is nicer for checking all arguments have been specified).
        bool result = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-mybool");
        TS_ASSERT_EQUALS(result, false);
        result = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-mybool2");
        TS_ASSERT_EQUALS(result, false);
        result = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-mybool3");
        TS_ASSERT_EQUALS(result, true);
        result = CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-mybool4");
        TS_ASSERT_EQUALS(result, true);

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetBoolCorrespondingToOption("-mybool5"),
                "The option '-mybool5' argument 'doDgy' cannot be converted to a bool.");

        // Fool the arguments into thinking that the options end at 5 entries to check an exception:
        new_argc = 5;
        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mydoubleval"),
                "No value(s) given after command line option '-mydoubleval'");

        delete[] new_argv;

        // Restore the real args
        CommandLineArguments::Instance()->p_argc = p_real_argc;
        CommandLineArguments::Instance()->p_argv = p_real_argv;
    }

    void TestCommandLineArgumentsMocker()
    {
        {
            /* HOW_TO_TAG General
             * Use mock (pretend) command line arguments
             */
            CommandLineArgumentsMocker wrapper("--option1 choice1 --option2 2 --option3 1.0 2.0 3.0");

            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("--option1"));
            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("--option2"));
            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("--option3"));

            std::vector<double> some_doubles = CommandLineArguments::Instance()->GetDoublesCorrespondingToOption("--option3");
            TS_ASSERT_EQUALS(some_doubles.size(), 3u);
            TS_ASSERT_DELTA(some_doubles[0], 1.0, 1e-12);
            TS_ASSERT_DELTA(some_doubles[1], 2.0, 1e-12);
            TS_ASSERT_DELTA(some_doubles[2], 3.0, 1e-12);

            unsigned a_number = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("--option2");
            TS_ASSERT_EQUALS(a_number, 2u);

            std::string a_string = CommandLineArguments::Instance()->GetStringCorrespondingToOption("--option1");
            TS_ASSERT_EQUALS(a_string,"choice1");
        }

        // Check the original arguments are still working...
        char** argv = *(CommandLineArguments::Instance()->p_argv);
        assert(argv != NULL);
        std::string arg_as_string(argv[0]);
        size_t pos = arg_as_string.find("TestCommandLineArguments");
        // If TestCommandLineArguments is not a substring of the commandline then pos==std::string::npos
        TS_ASSERT_DIFFERS(pos, std::string::npos);
    }

    /* A test which a user can run in order to check that they are passing command line arguments correctly*/
    void TestCommandLineArgumentsParrotting()
    {
        std::string verb = "--verbose";
        if (CommandLineArguments::Instance()->OptionExists(verb))
        {
            bool is_verbose = CommandLineArguments::Instance()->GetBoolCorrespondingToOption(verb);
            std::cout << "You have successfully set "<< verb << " to take the value "<< is_verbose<<".\n";
        }
    }
};

#endif /*TESTCOMMANDLINEARGUMENTS_HPP_*/
