#!/usr/bin/env python

"""Copyright (c) 2005-2019, University of Oxford.
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
"""

"""
Tester for Chaste Python code.  Runs unittest tests and makes the output
follow the same format as the cxxtest shipped with Chaste.

The first argument should be the path to a Python source file.  If this
file defines a MakeTestSuite method, this will be called to obtain a test
suite to run.  Otherwise, we scan the file for TestCase subclasses, and
run them all.

The optional second argument can restrict which test cases to run, as for
unittest.main.  For example, it can be 'MyTestCase' or
'MyTestCase.TestMethod'.

If no arguments are supplied it will run the example tests defined in this
module, class TestTest.

The module-scope variable CHASTE_TEST_OUTPUT will be set in test source
files to the path to the output folder, which is guaranteed to exist.

Also, the module-scope variable CHASTE_NUM_PROCS will be set to the maximum
number of processes to use for any test that can be run in parallel.
"""

import imp
import os
import sys
import time

try:
    import unittest2 as unittest
except ImportError:
    import unittest

if sys.version_info[:2] < (2,7):
    # Add subprocess.check_output to earlier Python versions
    import subprocess
    def check_output(*popenargs, **kwargs):
        """Run command with arguments and return its output as a string. Copied from Python 2.7"""
        if 'stdout' in kwargs:
            raise ValueError('stdout argument not allowed, it will be overridden.')
        process = subprocess.Popen(stdout=subprocess.PIPE, *popenargs, **kwargs)
        output, unused_err = process.communicate()
        retcode = process.poll()
        if retcode:
            cmd = kwargs.get("args")
            if cmd is None:
                cmd = popenargs[0]
            print >>sys.stderr, "Called process failed; output was:"
            print >>sys.stderr, output
            raise subprocess.CalledProcessError(retcode, cmd)
        return output
    subprocess.check_output = check_output


class TestTest(unittest.TestCase):
    """A simple test case for testing the framework."""
    def TestOk(self):
        self.failUnless(True)

    def TestFail(self):
        self.failIf(True)

    def TestError(self):
        self.assertEqual(1, 1/0)

#def MakeTestSuite():
#    """An example of how to create a suite of tests."""
#    return unittest.TestSuite(map(TestTest, ['TestOk', 'TestFail', 'TestError']))


class ChasteTestResult(unittest.TestResult):
    """A test result class that can print cxxtest-formatted text results to a stream.

    Used by _ChasteTestRunner.
    """
    separator1 = '=' * 70
    separator2 = '-' * 70

    def __init__(self, stream, descriptions=True):
        unittest.TestResult.__init__(self)
        self.stream = stream
        self.descriptions = descriptions
        self.debug = bool(int(os.environ.get('CHASTE_DEBUG', 0)))
        if self.debug:
            import pdb
            self.pdb = pdb
    
    def _Debug(self, err):
        type, value, traceback = err
        if (hasattr(sys, 'ps1')
            or not sys.stderr.isatty() or not sys.stdin.isatty() or not sys.stdout.isatty()
            or type == SyntaxError):
            # Debugging is not possible
            return
        self.pdb.post_mortem(traceback)

    def getDescription(self, test):
        if self.descriptions:
            return test.shortDescription() or str(test)
        else:
            return str(test)

    def startTest(self, test):
        unittest.TestResult.startTest(self, test)
        self.stream.write("Entering " + self.getDescription(test) + "\n")

    def addSuccess(self, test):
        unittest.TestResult.addSuccess(self, test)
        self.stream.write("Passed\n")

    def addError(self, test, err):
        unittest.TestResult.addError(self, test, err)
        self.stream.write(self.errors[-1][1])
        self.stream.write("Failed\n")
        if self.debug:
            self._Debug(err)

    def addFailure(self, test, err):
        unittest.TestResult.addFailure(self, test, err)
        self.stream.write(self.failures[-1][1])
        self.stream.write("Failed\n")
        if self.debug:
            self._Debug(err)

class ChasteTestRunner:
    """A test runner class that displays results in Chaste's cxxtest format."""
    def __init__(self, stream=sys.stdout, descriptions=True, profile=False, lineProfile=False):
        self.stream = stream
        self.descriptions = descriptions
        self.profiler = None
        if profile:
            import cProfile
            self.profiler = cProfile.Profile()
            self.profiler_type = 'func'
        if lineProfile:
            import line_profiler
            self.profiler = line_profiler.LineProfiler()
            self.profiler_type = 'line'
            import __builtin__
            __builtin__.__dict__['line_profile'] = self.profiler

    def _makeResult(self):
        return ChasteTestResult(self.stream, descriptions=self.descriptions)

    def run(self, test):
        "Run the given test case or test suite."
        result = self._makeResult()
        if self.profiler:
            self.profiler.enable()
        start_time = time.time()
        test(result)
        stop_time = time.time()
        if self.profiler:
            self.profiler.disable()
        time_taken = stop_time - start_time
        num_run = result.testsRun
        self.stream.write("\nRan %d test%s in %.3fs\n\n" %
                            (num_run, num_run != 1 and "s" or "", time_taken))
        if not result.wasSuccessful():
            failed, errored = map(len, (result.failures, result.errors))
            num_bad = failed + errored
            self.stream.write("Failed %d of %d tests\n" % (num_bad, num_run))
        else:
            self.stream.write("OK!\n")
        if self.profiler:
            # Append a profile report to the output
            if self.profiler_type == 'line':
                self.profiler.print_stats(self.stream)
            else:
                import pstats
                stats = pstats.Stats(self.profiler, stream=self.stream)
                self.stream.write('\n\nProfile report:\n\n')
                stats.sort_stats('time')
                stats.print_stats(50)
                self.stream.write('\n\nMethod callees:\n\n')
                stats.print_callees(.2, 30)
                self.stream.write('\n\nMethod callers:\n\n')
                stats.print_callers(.2, 30)
        return result

class ChasteTestLoader(unittest.TestLoader):
    """Allow test methods to start with either Test or test."""
    def getTestCaseNames(self, testCaseClass):
        """Return a sorted sequence of method names found within testCaseClass."""
        def isTestMethod(attrname, testCaseClass=testCaseClass):
            return (callable(getattr(testCaseClass, attrname)) and
                    (attrname.startswith('Test') or attrname.startswith('test')))
        test_names = filter(isTestMethod, dir(testCaseClass))
        for base_class in testCaseClass.__bases__:
            for test_name in self.getTestCaseNames(base_class):
                if test_name not in test_names:  # handle overridden methods
                    test_names.append(test_name)
        if self.sortTestMethodsUsing:
            test_names.sort(self.sortTestMethodsUsing)
        return test_names

def SetTestOutput(module):
    """Set the CHASTE_TEST_OUTPUT attribute in the given test module, and ensure the folder exists."""
    module.CHASTE_TEST_OUTPUT = os.getenv('CHASTE_TEST_OUTPUT', 'testoutput')
    try:
        os.makedirs(module.CHASTE_TEST_OUTPUT)
    except os.error:
        pass

def main(filepath, profile=False, lineProfile=False, numProcs=1):
    """Run tests defined in the given Python file.

    :param profile: whether to enable profiling of the test execution using cProfile.
    :param lineProfile: whether to enable profiling of the test execution using line_profiler.
    :param numProcs: maximum number of processes to use, if a test supports running in parallel.
    """
    if not os.path.isfile(filepath):
        raise ValueError(filepath + ' is not a file')
    base, ext = os.path.splitext(os.path.basename(filepath))
    if ext != '.py':
        raise ValueError(filepath + ' is not a Python source file')
    # Load Python file
    dirpath = os.path.dirname(filepath)
    (file, pathname, desc) = imp.find_module(base, [dirpath])
    try:
        module = imp.load_module(base, file, pathname, desc)
    except ImportError:
        print "Python module search path:", sys.path, os.environ.get('PYTHONPATH')
        raise
    finally:
        file.close()
    # Extract and run its tests
    SetTestOutput(module)
    module.CHASTE_NUM_PROCS = numProcs
    runner = ChasteTestRunner(profile=profile, lineProfile=lineProfile)
    if hasattr(module, 'MakeTestSuite') and callable(module.MakeTestSuite):
        suite = module.MakeTestSuite()
        result = runner.run(suite)
        sys.exit(not result.wasSuccessful())
    else:
        unittest.main(module=module, argv=[sys.argv[0]], testRunner=runner, testLoader=ChasteTestLoader())

if __name__ == '__main__':
    if len(sys.argv) > 1:
        profile = '--profile' in sys.argv
        if profile:
            sys.argv.remove('--profile')
        line_profile = '--line-profile' in sys.argv
        if line_profile:
            sys.argv.remove('--line-profile')
        try:
            i = sys.argv.index('--num-procs')
            num_procs = int(sys.argv[i+1])
            sys.argv[i:i+2] = []
        except ValueError:
            num_procs = 1
        original_argv = sys.argv[:]
        sys.argv[0:1] = [] # Remove this file from list
        main(sys.argv[0], profile=profile, lineProfile=line_profile, numProcs=num_procs)
        sys.argv = original_argv
    else:
        # Default test of this file
        unittest.main(testRunner=ChasteTestRunner(), testLoader=ChasteTestLoader())
