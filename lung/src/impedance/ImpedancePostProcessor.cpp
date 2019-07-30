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

#include "ImpedancePostProcessor.hpp"
#include "Exception.hpp"
#include <assert.h>
#include <iostream>
#include <complex>
#include <algorithm>

ImpedancePostProcessor::ImpedancePostProcessor(std::vector<double>& rFrequencies,
                                               std::vector<std::complex<double> >& rImpedances) : mrFrequencies(rFrequencies),
                                                                                                  mrImpedances(rImpedances)
{
    assert(mrFrequencies.size() == mrImpedances.size());

    if ((std::find(mrFrequencies.begin(), mrFrequencies.end(), 5.0 ) == mrFrequencies.end()) ||
        (std::find(mrFrequencies.begin(), mrFrequencies.end(), 20.0) == mrFrequencies.end()))
    {
        EXCEPTION("Impedance post processor requires data points at 5 Hz & 20 Hz");
    }
}


double ImpedancePostProcessor::GetR5MinusR20()
{
    double R5 = 0.0;
    double R20 = 0.0;

    for (unsigned freq_index = 0; freq_index < mrFrequencies.size(); ++freq_index)
    {
        if (mrFrequencies[freq_index] == 5.0)
        {
            R5 = real(mrImpedances[freq_index]);
        }

        if (mrFrequencies[freq_index] == 20.0)
        {
            R20 = real(mrImpedances[freq_index]);
        }
    }

    return R5 - R20;
}


double ImpedancePostProcessor::GetResonantFrequency()
{
    unsigned lower_bound_index = 0u;
    unsigned upper_bound_index = 0u;

    for (unsigned freq_index = 0; freq_index < mrFrequencies.size(); ++freq_index)
    {
        if (imag(mrImpedances[freq_index]) <= 0.0)
        {
            lower_bound_index = freq_index;
        }
        else if (imag(mrImpedances[freq_index]) > 0.0)
        {
            upper_bound_index = freq_index;
            break;
        }
    }

    double slope = (imag(mrImpedances[upper_bound_index]) - imag(mrImpedances[lower_bound_index]))/
            (mrFrequencies[upper_bound_index] - mrFrequencies[lower_bound_index]);

    double b = imag(mrImpedances[lower_bound_index]) - slope*mrFrequencies[lower_bound_index];

    return -b/slope;
}

double ImpedancePostProcessor::GetAx()
{
    double ax_integral = 0.0;
    double fres = GetResonantFrequency();

    unsigned freq_index = 0u;

    while (mrFrequencies[freq_index] <= fres)
    {
        if (mrFrequencies[freq_index] >= 5)
        {
            if (mrFrequencies[freq_index+1] <= fres)
            {
                ax_integral += (mrFrequencies[freq_index+1] - mrFrequencies[freq_index])*(imag(mrImpedances[freq_index+1]) + imag(mrImpedances[freq_index]))/2.0;
            }
            else
            {
                ax_integral += (fres - mrFrequencies[freq_index])*(0.0 + imag(mrImpedances[freq_index]))/2.0;
            }
        }

        freq_index++;
    }

    return ax_integral;
}

double ImpedancePostProcessor::GetRrs()
{
    double rrs_integral = 0.0;

    for (unsigned freq_index = 0; freq_index < mrFrequencies.size(); ++freq_index)
    {
        if (mrFrequencies[freq_index] >= 5 && freq_index+1 < mrFrequencies.size())
        {
            rrs_integral += (mrFrequencies[freq_index+1] - mrFrequencies[freq_index])*(real(mrImpedances[freq_index+1]) + real(mrImpedances[freq_index]))/2.0;
        }
    }

    return rrs_integral / (mrFrequencies[mrFrequencies.size() - 1] - 5);
}
