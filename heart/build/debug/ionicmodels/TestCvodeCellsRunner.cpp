/* Generated file, do not edit */

#ifndef CXXTEST_RUNNING
#define CXXTEST_RUNNING
#endif

#define _CXXTEST_HAVE_STD
#define _CXXTEST_HAVE_EH
#include <cxxtest/TestListener.h>
#include <cxxtest/TestTracker.h>
#include <cxxtest/TestRunner.h>
#include <cxxtest/RealDescriptions.h>
#include <cxxtest/ErrorPrinter.h>

#include "CommandLineArguments.hpp"
int main( int argc, char *argv[] ) {
 CommandLineArguments::Instance()->p_argc = &argc;
 CommandLineArguments::Instance()->p_argv = &argv;
 return CxxTest::ErrorPrinter().run();
}
#include "heart/test/ionicmodels/TestCvodeCells.hpp"

static TestCvodeCells suite_TestCvodeCells;

static CxxTest::List Tests_TestCvodeCells = { 0, 0 };
CxxTest::StaticSuiteDescription suiteDescription_TestCvodeCells( "heart/test/ionicmodels/TestCvodeCells.hpp", 125, "TestCvodeCells", suite_TestCvodeCells, Tests_TestCvodeCells );

static class TestDescription_TestCvodeCells_TestLuoRudyCvodeCell : public CxxTest::RealTestDescription {
public:
 TestDescription_TestCvodeCells_TestLuoRudyCvodeCell() : CxxTest::RealTestDescription( Tests_TestCvodeCells, suiteDescription_TestCvodeCells, 129, "TestLuoRudyCvodeCell" ) {}
 void runTest() { suite_TestCvodeCells.TestLuoRudyCvodeCell(); }
} testDescription_TestCvodeCells_TestLuoRudyCvodeCell;

static class TestDescription_TestCvodeCells_TestShannon2004 : public CxxTest::RealTestDescription {
public:
 TestDescription_TestCvodeCells_TestShannon2004() : CxxTest::RealTestDescription( Tests_TestCvodeCells, suiteDescription_TestCvodeCells, 318, "TestShannon2004" ) {}
 void runTest() { suite_TestCvodeCells.TestShannon2004(); }
} testDescription_TestCvodeCells_TestShannon2004;

static class TestDescription_TestCvodeCells_TestArchivingCvodeCells : public CxxTest::RealTestDescription {
public:
 TestDescription_TestCvodeCells_TestArchivingCvodeCells() : CxxTest::RealTestDescription( Tests_TestCvodeCells, suiteDescription_TestCvodeCells, 487, "TestArchivingCvodeCells" ) {}
 void runTest() { suite_TestCvodeCells.TestArchivingCvodeCells(); }
} testDescription_TestCvodeCells_TestArchivingCvodeCells;

#include <cxxtest/Root.cpp>
