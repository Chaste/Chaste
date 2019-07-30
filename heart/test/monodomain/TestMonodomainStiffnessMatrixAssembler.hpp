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

#ifndef TESTMONODOMAINSTIFFNESSMATRIX_HPP_
#define TESTMONODOMAINSTIFFNESSMATRIX_HPP_

#include "TetrahedralMesh.hpp"
#include "MonodomainStiffnessMatrixAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "LuoRudy1991.hpp"
#include "PlaneStimulusCellFactory.hpp"
#include "MonodomainTissue.hpp"
#include "PetscMatTools.hpp"

class TestMonodomainStiffnessMatrixAssembler : public CxxTest::TestSuite
{
public:

    void TestMonodomainStiffnessMatrixAssembler1d()
    {
        TetrahedralMesh<1,1> mesh;
        double h = 0.1;
        mesh.ConstructRegularSlabMesh(h, 0.5);

        Mat mat;
        PetscTools::SetupMat(mat, mesh.GetNumNodes(), mesh.GetNumNodes(), 3);

        PlaneStimulusCellFactory<CellLuoRudy1991FromCellML, 1> cell_factory;
        cell_factory.SetMesh(&mesh);
        MonodomainTissue<1> monodomain_tissue( &cell_factory );
        MonodomainStiffnessMatrixAssembler<1,1> assembler(&mesh, &monodomain_tissue);

        assembler.SetMatrixToAssemble(mat);
        assembler.Assemble();

        PetscMatTools::Finalise(mat);

        double sigma = 1.75;

        int lo, hi;
        MatGetOwnershipRange(mat, &lo, &hi);

        for (unsigned i=lo; i<(unsigned)hi; i++)
        {
            for (unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                if (i>0 && i<mesh.GetNumNodes()-1)
                {
                    // All rows except first and last should look like
                    // [0, .., 0, -sigma/h, 2*sigma/h, -sigma/h, 0, .., 0]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, 2*sigma/h, 1e-5);
                    }
                    else if (j==i+1 || j+1==i)
                    {
                        TS_ASSERT_DELTA(value, -sigma/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if (i==0)
                {
                    // top row: [sigma/h, -sigma/h, 0, .., 0]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, sigma/h, 1e-5);
                    }
                    else if (j==i+1)
                    {
                        TS_ASSERT_DELTA(value, -sigma/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if (i+1==mesh.GetNumNodes())
                {
                    // bottom row: [0, .., 0, -sigma/h, sigma/h]
                    if (j==i)
                    {
                        TS_ASSERT_DELTA(value, sigma/h, 1e-5);
                    }
                    else if (j+1==i)
                    {
                        TS_ASSERT_DELTA(value, -sigma/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }

            }
        }

        PetscTools::Destroy(mat);
    }
};

#endif /* TESTMONODOMAINSTIFFNESSMATRIX_HPP_ */
