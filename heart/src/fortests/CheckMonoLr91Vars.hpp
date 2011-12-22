/*

Copyright (C) University of Oxford, 2005-2011

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

        std::vector<double> odeVars = problem.GetMonodomainTissue()->GetCardiacCell(index.Global)->rGetStateVariables();
        for (int j=0; j<8; j++)
        {
            // if not voltage or calcium ion conc, test whether between 0 and 1
            if ((j!=0) && (j!=7))
            {
                TS_ASSERT_LESS_THAN_EQUALS(  odeVars[j], 1.0);
                TS_ASSERT_LESS_THAN_EQUALS( -odeVars[j], 0.0);
            }
        }
    }
}

#endif /*CHECKMONOLR91VARS_HPP_*/
