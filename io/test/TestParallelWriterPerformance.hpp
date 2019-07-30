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

#ifndef TESTPARALLELWRITERPERFORMANCE_HPP_
#define TESTPARALLELWRITERPERFORMANCE_HPP_

#include <cxxtest/TestSuite.h>

#include "ParallelColumnDataWriter.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include <petsc.h>
#include "PetscSetupAndFinalize.hpp"

class TestParallelWriterPerformance : public CxxTest::TestSuite
{
public:

    void Test1()
    {
        const unsigned SIZE = 1000;
        const unsigned REPETITIONS = 10;

        // Create a distibuted vector
        DistributedVectorFactory factory(SIZE);

        Vec petsc_vec = factory.CreateVec();
        DistributedVector distributed_vector = factory.CreateDistributedVector(petsc_vec);
        for (DistributedVector::Iterator index = distributed_vector.Begin();
             index!= distributed_vector.End();
             ++index)
        {
            distributed_vector[index] = -(double)(index.Local*index.Global);
        }
        distributed_vector.Restore();

        // Set up a parallel writer
        ParallelColumnDataWriter parallel_writer("TestParallelWriterPerformance","ParallelColumnWriter");
        unsigned time_var_id = parallel_writer.DefineUnlimitedDimension("Time","msecs");
        parallel_writer.DefineFixedDimension("Node","dimensionless", SIZE);
        unsigned var1_id = parallel_writer.DefineVariable("Var1","LightYears");
        parallel_writer.EndDefineMode();

        // Write multiple times
        for (unsigned i=0; i<REPETITIONS; i++)
        {
            double time = (double)i;
            parallel_writer.PutVariable(time_var_id, time);
            parallel_writer.PutVector(var1_id, petsc_vec);
            parallel_writer.AdvanceAlongUnlimitedDimension();
        }

        PetscTools::Destroy(petsc_vec);
    }
};

#endif /*TESTPARALLELWRITERPERFORMANCE_HPP_*/
