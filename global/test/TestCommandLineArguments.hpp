/*

Copyright (C) University of Oxford, 2005-2012

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#ifndef TESTCOMMANDLINEARGUMENTS_HPP_
#define TESTCOMMANDLINEARGUMENTS_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>
#include <cassert>
#include "CommandLineArguments.hpp"


class TestCommandLineArguments : public CxxTest::TestSuite
{
public:

    void TestCommandLineArgumentsSingleton() throw(Exception)
    {
        // Test that argc and argv are populated
        int argc = *(CommandLineArguments::Instance()->p_argc);
        TS_ASSERT_LESS_THAN(0, argc); // argc should always be 1 or greater

        // argv[0] will be equal to global/build/debug/TestCommandLineArgumentsRunner
        // or global/build/optimised/TestCommandLineArgumentsRunner, etc
        char** argv = *(CommandLineArguments::Instance()->p_argv);
        assert(argv != NULL);
        std::string arg_as_string(argv[0]);
        std::string final_part_of_string = arg_as_string.substr(arg_as_string.length()-30,arg_as_string.length());
        TS_ASSERT_EQUALS("TestCommandLineArgumentsRunner",final_part_of_string);

        // Now test OptionExists() and GetValueCorrespondingToOption()
        //
        // The following tests would require the following arguments to be passed
        // in:
        // ./global/build/debug/TestCommandLineArgumentsRunner -myoption -myintval 24 -mydoubleval 3.14 -3.14 -m2intval -42 -mystrings Baboons Monkeys Gibbons -mystring more_baboons
        //
        // To test the methods we overwrite the arg_c and arg_v contained in the
        // singleton with the arguments that were needed.
        int new_argc = 15;
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

        char** new_argv = new char*[15];
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

        // Save the real args to be restored at the end
        int* p_real_argc = CommandLineArguments::Instance()->p_argc;
        char*** p_real_argv = CommandLineArguments::Instance()->p_argv;

        // Overwrite the args
        CommandLineArguments::Instance()->p_argc = &new_argc;
        CommandLineArguments::Instance()->p_argv = &new_argv;

        // Test OptionExists()
        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-myoption"));

        TS_ASSERT( ! CommandLineArguments::Instance()->OptionExists("-asddsgijdfgokgfgurgher"));

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->OptionExists("-42"),
                "A command line option must begin with '-' followed by a non-numeric character.");

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetStringsCorrespondingToOption("-myoption"),
                "No value(s) given after command line option '-myoption'");

        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetStringsCorrespondingToOption("-mynonsense"),
                "Command line option '-mynonsense' does not exist");


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

        // Fool the arguments into thinking that the options end at 5 entries to check an exception:
        new_argc = 5;
        TS_ASSERT_THROWS_THIS(CommandLineArguments::Instance()->GetValueCorrespondingToOption("-mydoubleval"),
                "No value(s) given after command line option '-mydoubleval'");

        delete new_argv;

        // Restore the real args
        CommandLineArguments::Instance()->p_argc = p_real_argc;
        CommandLineArguments::Instance()->p_argv = p_real_argv;
    }
};

#endif /*TESTCOMMANDLINEARGUMENTS_HPP_*/
