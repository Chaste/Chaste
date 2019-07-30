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


#ifndef CHECKMONOLR91VARS_HPP_
#define CHECKMONOLR91VARS_HPP_

#include "MonodomainProblem.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void CheckMonoLr91Vars(MonodomainProblem<ELEMENT_DIM, SPACE_DIM>& problem)
{

    DistributedVector voltage = problem.rGetMesh().GetDistributedVectorFactory()->CreateDistributedVector(problem.GetSolution());
    for (DistributedVector::Iterator index = voltage.Begin();
         index != voltage.End();
         ++index)
    {
        // assuming LR model has Ena = 54.4 and Ek = -77
        double Ena   =  54.4;
        double Ek    = -77.0;

        TS_ASSERT_LESS_THAN_EQUALS( voltage[index] , Ena +  30);
        TS_ASSERT_LESS_THAN_EQUALS(-voltage[index] + (Ek-30), 0);

        std::vector<double> ode_vars = problem.GetMonodomainTissue()->GetCardiacCell(index.Global)->GetStdVecStateVariables();
        for (int j=0; j<8; j++)
        {
            // if not voltage or calcium ion conc, test whether between 0 and 1
            if ((j!=0) && (j!=7))
            {
                TS_ASSERT_LESS_THAN_EQUALS(  ode_vars[j], 1.0);
                TS_ASSERT_LESS_THAN_EQUALS( -ode_vars[j], 0.0);
            }
        }
    }
}

#endif /*CHECKMONOLR91VARS_HPP_*/
