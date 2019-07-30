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

#ifndef IMPEDANCEPOSTPROCESSOR_HPP_
#define IMPEDANCEPOSTPROCESSOR_HPP_

#include <complex>
#include <vector>

/**
 * A class to facilitate easy calculation of clinically relevant impedance metrics from impedance simulations.
 *
 * Note that this class will attempt to interpolate values if insufficient frequency data is available, but
 * for increased accuracy it is advisable that the frequency list contains 5Hz, 20Hz and a few entries close
 * to the actual resonant frequency (~10Hz in a healthy human, less in an asthmatic)
 */
class ImpedancePostProcessor
{

public:
    /**
     * Constructor
     *
     * @param rFrequencies A vector of test frequencies. Must be monotonically increasing.
     * @param rImpedances A corresponding vector of impedances.
     */
    ImpedancePostProcessor(std::vector<double>& rFrequencies, std::vector<std::complex<double> >& rImpedances);

    /**
     * Returns the resistance at 5Hz minus the resistance at 20Hz
     *
     * @return R5 - R20
     */
    double GetR5MinusR20();

    /**
     * Gets the resonant frequency of the lung (FRes)
     *
     * @return The resonant frequency of the lung
     */
    double GetResonantFrequency();

    /**
     * Gets the integral of the reactance between 5Hz and FRes
     *
     * @return The integral of the reactance between 5Hz and FRes
     */
    double GetAx();

    /**
     * Gets the average of the resistance from 5Hz onwards
     *
     * @return The average of the resistance from 5Hz onwards
     */
    double GetRrs();


    /**
     * Gets the std deviation of the resistance from 5Hz onwards
     *
     * @return The std deviation of the resistance from 5Hz onwards
     */
//    double GetStdDevRrs();

private:
    /**
     * The set of test frequencies
     */
    std::vector<double>& mrFrequencies;

    /**
     * The set of impedances
     */
    std::vector<std::complex<double> >& mrImpedances;
};

#endif /* IMPEDANCEPOSTPROCESSOR_HPP_ */
