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

#ifndef TESTABSTRACTFUNCTIONALCALCULATOR_HPP_
#define TESTABSTRACTFUNCTIONALCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>

#include "AbstractFunctionalCalculator.hpp"
#include "TetrahedralMesh.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "DistributedVector.hpp"

/*
 * Returns 1.0 everywhere so that the total integral over the mesh of
 * this integrand is just the volume of the mesh. For testing.
 */
template<unsigned DIM>
class VolumeCalculator : public AbstractFunctionalCalculator<DIM,DIM,1>
{
    double GetIntegrand(ChastePoint<DIM>& rX,
                        c_vector<double,1>& rU,
                        c_matrix<double,1,DIM>& rGradU)
    {
        return 1.0;
    }
};

// Check x and u are interpolated correctly
class ExampleFunctionalOne : public AbstractFunctionalCalculator<2,2,2>
{
    double GetIntegrand(ChastePoint<2>& rX,
                        c_vector<double,2>& rU,
                        c_matrix<double,2,2>& rGradU)
    {
        return rX[0]*rU[0] + rX[1]*rU[1];
    }
};

// Check grad_u is interpolated correctly
class ExampleFunctionalTwo : public AbstractFunctionalCalculator<2,2,2>
{
    double GetIntegrand(ChastePoint<2>& rX,
                        c_vector<double,2>& rU,
                        c_matrix<double,2,2>& rGradU)
    {
        return rX[0]*rU[0] + rX[1]*rU[1] + 0.5*(rGradU(0,0)+rGradU(0,1)+rGradU(1,0)+rGradU(1,1));
    }
};

class TestAbstractFunctionalCalculator : public CxxTest::TestSuite
{
public:

    void TestWithVolumeCalculator()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);
        DistributedTetrahedralMesh<2,2> distributed_mesh;
        distributed_mesh.ConstructFromMeshReader(reader);

        VolumeCalculator<2> volume_calculator;

        Vec vec = PetscTools::CreateAndSetVec(mesh.GetNumNodes(), 0.0);

        double result = volume_calculator.Calculate(mesh,vec);
        TS_ASSERT_DELTA(result, mesh.GetVolume(), 1e-12);
        double distributed_result = volume_calculator.Calculate(distributed_mesh,vec);
        TS_ASSERT_DELTA(result, distributed_result, 1e-12);

        Vec bad_vec = PetscTools::CreateAndSetVec(mesh.GetNumNodes()+1, 0.0);
        TS_ASSERT_THROWS_THIS(volume_calculator.Calculate(mesh,bad_vec),"The solution size does not match the mesh");

        VecDestroy(vec);
        VecDestroy(bad_vec);
    }

    void TestWithExampleFunctionals()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        DistributedTetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Test interpolation of x and u
        // Integrate x^2 + 2y over the unit square
        // = 4/3
        ExampleFunctionalOne calculator;

        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        Vec petsc_vec = p_factory->CreateVec(2);
        DistributedVector vec1 = p_factory->CreateDistributedVector(petsc_vec);
        DistributedVector::Stripe u1(vec1, 0);
        DistributedVector::Stripe v1(vec1, 1);
        for (DistributedVector::Iterator index = vec1.Begin();
             index != vec1.End();
             ++index)
        {
            Node<2>* p_node = mesh.GetNode(index.Global);
            u1[index] = p_node->rGetLocation()[0];
            v1[index] = 2.0;
        }
        vec1.Restore();

        double result = calculator.Calculate(mesh, petsc_vec);
        TS_ASSERT_DELTA(result, 4.0/3.0, 1e-6);

        // Test interpolation of grad_u
        // Integrate x^2 + y^2 + 1 over the unit square
        // = 5/3
        ExampleFunctionalTwo other_calculator;

        DistributedVector vec2 = p_factory->CreateDistributedVector(petsc_vec);
        DistributedVector::Stripe u2(vec2, 0);
        DistributedVector::Stripe v2(vec2, 1);
        for (DistributedVector::Iterator index = vec2.Begin();
             index != vec2.End();
             ++index)
        {
            Node<2>* p_node = mesh.GetNode(index.Global);
            u2[index] = p_node->rGetLocation()[0];
            v2[index] = p_node->rGetLocation()[1];
        }
        vec2.Restore();

        result = other_calculator.Calculate(mesh, petsc_vec);
        TS_ASSERT_DELTA(result, 1 + 2.0/3.0, 1e-6);

        VecDestroy(petsc_vec);
    }

    void TestWithExampleFunctionalsNonDistributed()
    {
        TrianglesMeshReader<2,2> reader("mesh/test/data/square_128_elements");
        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Test interpolation of x and u
        // Integrate x^2 + 2y over the unit square
        // = 4/3
        ExampleFunctionalOne calculator;

        DistributedVectorFactory* p_factory = mesh.GetDistributedVectorFactory();
        Vec petsc_vec = p_factory->CreateVec(2);
        DistributedVector vec1 = p_factory->CreateDistributedVector(petsc_vec);
        DistributedVector::Stripe u1(vec1, 0);
        DistributedVector::Stripe v1(vec1, 1);
        for (DistributedVector::Iterator index = vec1.Begin();
             index != vec1.End();
             ++index)
        {
            Node<2>* p_node = mesh.GetNode(index.Global);
            u1[index] = p_node->rGetLocation()[0];
            v1[index] = 2.0;
        }
        vec1.Restore();

        double result = calculator.Calculate(mesh, petsc_vec);
        TS_ASSERT_DELTA(result, 4.0/3.0, 1e-6);

        // Test interpolation of grad_u
        // Integrate x^2 + y^2 + 1 over the unit square
        // = 5/3
        ExampleFunctionalTwo other_calculator;

        DistributedVector vec2 = p_factory->CreateDistributedVector(petsc_vec);
        DistributedVector::Stripe u2(vec2, 0);
        DistributedVector::Stripe v2(vec2, 1);
        for (DistributedVector::Iterator index = vec2.Begin();
             index != vec2.End();
             ++index)
        {
            Node<2>* p_node = mesh.GetNode(index.Global);
            u2[index] = p_node->rGetLocation()[0];
            v2[index] = p_node->rGetLocation()[1];
        }
        vec2.Restore();

        result = other_calculator.Calculate(mesh, petsc_vec);
        TS_ASSERT_DELTA(result, 1 + 2.0/3.0, 1e-6);

        VecDestroy(petsc_vec);
    }
};

#endif /*TESTABSTRACTFUNCTIONALCALCULATOR_HPP_*/
