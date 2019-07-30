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

#ifndef _TESTIMPEDANCEPOSTPROCESSOR_HPP_
#define _TESTIMPEDANCEPOSTPROCESSOR_HPP_

#include <cxxtest/TestSuite.h>

#include "ImpedancePostProcessor.hpp"


class TestImpedancePostProcessor : public CxxTest::TestSuite
{
private:

public:
    void TestExceptions()
    {
        std::vector<double> freqs;
        freqs.push_back(4);
        freqs.push_back(7);
        freqs.push_back(20);

        std::vector<std::complex<double> > imps;
        imps.push_back(std::complex<double>(5, 2));
        imps.push_back(std::complex<double>(4, 4));
        imps.push_back(std::complex<double>(3.5, 4.5));

        TS_ASSERT_THROWS_CONTAINS(ImpedancePostProcessor(freqs, imps), "Impedance post processor requires data points at 5 Hz & 20 Hz");
    }

    void TestCalculatePropertiesWithNiceData()
    {
        std::vector<double> freqs;
        freqs.push_back(4);
        freqs.push_back(5);
        freqs.push_back(7);
        freqs.push_back(10);
        freqs.push_back(20);
        freqs.push_back(30);

        std::vector<std::complex<double> > imps;
        imps.push_back(std::complex<double>(8,-10));
        imps.push_back(std::complex<double>(7,-5));
        imps.push_back(std::complex<double>(6, 0));
        imps.push_back(std::complex<double>(5, 2));
        imps.push_back(std::complex<double>(4, 4));
        imps.push_back(std::complex<double>(3.5, 4.5));

        ImpedancePostProcessor processor(freqs, imps);

        TS_ASSERT_DELTA(processor.GetR5MinusR20(), 3.0, 1e-4);

        TS_ASSERT_DELTA(processor.GetResonantFrequency(), 7.0, 1e-4);
        TS_ASSERT_DELTA(processor.GetAx(), -5.0, 1e-4);

        TS_ASSERT_DELTA(processor.GetRrs(), 4.48, 1e-2);
    }

    void TestCalculatePropertiesWithHorridData()
    {
        //Corner cases
        std::vector<double> freqs;
        freqs.push_back(5);
        freqs.push_back(6);
        freqs.push_back(7);
        freqs.push_back(10);
        freqs.push_back(20);
        freqs.push_back(30);

        std::vector<std::complex<double> > imps;
        imps.push_back(std::complex<double>(6,-2));
        imps.push_back(std::complex<double>(6,-1));
        imps.push_back(std::complex<double>(6, 2));
        imps.push_back(std::complex<double>(6, 2));
        imps.push_back(std::complex<double>(6, 4));
        imps.push_back(std::complex<double>(6, 4.5));

        ImpedancePostProcessor processor(freqs, imps);

        TS_ASSERT_DELTA(processor.GetR5MinusR20(), 0.0, 1e-4);

        TS_ASSERT_DELTA(processor.GetResonantFrequency(), 6.3333, 1e-4);
        TS_ASSERT_DELTA(processor.GetAx(), -1.66666, 1e-4);

        TS_ASSERT_DELTA(processor.GetRrs(), 6.0, 1e-2);
    }
};

#endif /*_TESTIMPEDANCEPOSTPROCESSOR_HPP_*/
