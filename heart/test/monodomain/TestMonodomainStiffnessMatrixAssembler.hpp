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

    void TestMonodomainStiffnessMatrixAssembler1d() throw(Exception)
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

        for(unsigned i=lo; i<(unsigned)hi; i++)
        {
            for(unsigned j=0; j<mesh.GetNumNodes(); j++)
            {
                double value = PetscMatTools::GetElement(mat,i,j);
                if(i>0 && i<mesh.GetNumNodes()-1)
                {
                    // All rows except first and last should look like
                    // [0, .., 0, -sigma/h, 2*sigma/h, -sigma/h, 0, .., 0]
                    if(j==i)
                    {
                        TS_ASSERT_DELTA(value, 2*sigma/h, 1e-5);
                    }
                    else if(j==i+1 || j+1==i)
                    {
                        TS_ASSERT_DELTA(value, -sigma/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if(i==0)
                {
                    // top row: [sigma/h, -sigma/h, 0, .., 0]
                    if(j==i)
                    {
                        TS_ASSERT_DELTA(value, sigma/h, 1e-5);
                    }
                    else if(j==i+1)
                    {
                        TS_ASSERT_DELTA(value, -sigma/h, 1e-5);
                    }
                    else
                    {
                        TS_ASSERT_DELTA(value, 0.0, 1e-5);
                    }
                }
                if(i+1==mesh.GetNumNodes())
                {
                    // bottom row: [0, .., 0, -sigma/h, sigma/h]
                    if(j==i)
                    {
                        TS_ASSERT_DELTA(value, sigma/h, 1e-5);
                    }
                    else if(j+1==i)
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

        MatDestroy(mat);
    }

};

#endif /* TESTMONODOMAINSTIFFNESSMATRIX_HPP_ */
