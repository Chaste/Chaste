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


#ifndef TESTPSEUDOECGCALCULATOR_HPP_
#define TESTPSEUDOECGCALCULATOR_HPP_

#include <cxxtest/TestSuite.h>
#include <iostream>

#include "TetrahedralMesh.hpp" //must be first, it gets UblasIncludes from the mesh classes (ChastePoint.hpp)
#include "DistributedTetrahedralMesh.hpp"
#include "PseudoEcgCalculator.hpp"
#include "ReplicatableVector.hpp"
#include "Hdf5DataReader.hpp"
#include "Hdf5DataWriter.hpp"
#include "OutputFileHandler.hpp"
#include "TrianglesMeshReader.hpp"
#include "TrianglesMeshReader.hpp"
#include "TetrahedralMesh.hpp"
#include "HeartConfig.hpp"
#include "PetscSetupAndFinalize.hpp"

class TestPseudoEcgCalculator : public CxxTest::TestSuite
{

public:

    void TestCalculator1DLinearGradient() throw (Exception)
    {

        //read in the 1D mesh, from 0 to 1
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        ////////////////////////////////////////////////////
        //First we write an hdf5 file with a gradient of V
        // i.e. V(x) = x; We write 4 time steps.
        ///////////////////////////////////////////////////
        unsigned number_nodes = mesh.GetNumNodes();

        DistributedVectorFactory factory(number_nodes);
        Hdf5DataWriter writer(factory, "hdf5", "gradient_V", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        //we call the variable like the default for the pseudoecg calculator
        int V_id = writer.DefineVariable("V", "mV");
        writer.DefineUnlimitedDimension("Time", "msec");
        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        //4 time steps
        unsigned number_of_time_steps = 4;
        for (unsigned time_step=0; time_step<number_of_time_steps; time_step++)
        {
            // Write some values, the value of V is the same as the x coordinate, i.e. a gradient
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
                distributed_vector_2[index] = index.Global/100.0;
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();

            // Write the vector
            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(V_id, petsc_data_2);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);

        ///////////////////////////////////////////////////
        // Now we compute the pseudo ECG. We set an electrode at x=15.
        ///////////////////////////////////////////////////
        ChastePoint<1> probe_electrode(15.0);

        PseudoEcgCalculator<1,1,1> calculator (mesh, probe_electrode, "hdf5", "gradient_V");
        double pseudo_ecg;

        // The expected result is the integral of: - d(gradV)*dgrad(1/r) in dx
        // Because in this simple case d(gradV)=1, the result is simply -1/(15-x) evaluated
        // between 0 and 1
        double expected_result = -(1/14.0-1/15.0);
        for (unsigned k = 0; k< number_of_time_steps; k++)
        {
            pseudo_ecg = calculator.ComputePseudoEcgAtOneTimeStep(k);
            TS_ASSERT_DELTA(pseudo_ecg, expected_result,1e-6);
        }

        //now we test the writer method
        calculator.WritePseudoEcg();

        std::string output_dir = "ChasteResults/output";//default value
        std::string command;
        command = "diff -a -I \"Created by Chaste\" " + OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/PseudoEcgFromElectrodeAt_15_0_0.dat "
                  + "heart/test/data/ValidPseudoEcg1D.dat";
        TS_ASSERT_EQUALS(system(command.c_str()), 0);

        ChastePoint<1> bad_probe_electrode(0.0021132486540519);
        PseudoEcgCalculator<1,1,1> bad_calculator (mesh, bad_probe_electrode, "hdf5", "gradient_V");
        TS_ASSERT_THROWS_THIS(bad_calculator.ComputePseudoEcgAtOneTimeStep(0), "Probe is on a mesh Gauss point.");

    }

    void TestCalculator1DParabolic() throw (Exception)
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        ////////////////////////////////////////////////////
        //First we write an hdf5 file with V as a parabolic function of x
        // i.e. V(x) = x^2; We write 4 time steps.
        ///////////////////////////////////////////////////
        unsigned number_nodes = mesh.GetNumNodes();

        DistributedVectorFactory factory(number_nodes);

        Hdf5DataWriter writer(factory, "hdf5", "parabolic_V", false);
        writer.DefineFixedDimension(number_nodes);

        int node_id = writer.DefineVariable("Node", "dimensionless");
        //we call the variable like the default for the pseudoecg calculator
        int V_id = writer.DefineVariable("V", "mV");
        writer.DefineUnlimitedDimension("Time", "msec");
        writer.EndDefineMode();

        Vec petsc_data_1 = factory.CreateVec();
        DistributedVector distributed_vector_1 = factory.CreateDistributedVector(petsc_data_1);

        Vec petsc_data_2 = factory.CreateVec();
        DistributedVector distributed_vector_2 = factory.CreateDistributedVector(petsc_data_2);

        //4 time steps
        unsigned number_of_time_steps = 4;
        for (unsigned time_step=0; time_step<number_of_time_steps; time_step++)
        {
            // Write some values, the value of V is x^2, i.e. a parabola
            for (DistributedVector::Iterator index = distributed_vector_1.Begin();
                 index!= distributed_vector_1.End();
                 ++index)
            {
                distributed_vector_1[index] = index.Global;
                distributed_vector_2[index] = (index.Global/100.0)*(index.Global/100.0);
            }
            distributed_vector_1.Restore();
            distributed_vector_2.Restore();

            // Write the vector
            writer.PutVector(node_id, petsc_data_1);
            writer.PutVector(V_id, petsc_data_2);
            writer.PutUnlimitedVariable(time_step);
            writer.AdvanceAlongUnlimitedDimension();
        }

        writer.Close();

        VecDestroy(petsc_data_1);
        VecDestroy(petsc_data_2);

        ///////////////////////////////////////////////////
        // Now we compute the pseudo ECG. We set an electrode at x=15.
        ///////////////////////////////////////////////////

        ChastePoint<1> probe_electrode(15.0);

        PseudoEcgCalculator<1,1,1> calculator (mesh, probe_electrode, "hdf5", "parabolic_V", "V", true);

        double pseudo_ecg; //stores the results

        calculator.SetDiffusionCoefficient(1.0);

        // The expected result is the integral of: - d(gradV)*dgrad(1/r) in dx
        // In this case d(gradV)/dx=2x, hence the result is the integral of: - (2x * d(1/(15-x))/dx)
        // Integrating by parts with f = 2x and g = 1/(15-x) and evaluating between 0 and 1

        double expected_result = -( (2/14.0) - 2.0*log(1/14.0) + 2.0*log(1/15.0));
        for (unsigned k = 0; k< number_of_time_steps; k++)
        {
            pseudo_ecg = calculator.ComputePseudoEcgAtOneTimeStep(k);
            TS_ASSERT_DELTA(pseudo_ecg, expected_result,1e-6);
        }

        // Now try a different diffusion coefficient
        double diff_coeff = 2.0;
        calculator.SetDiffusionCoefficient(diff_coeff);

        //since we are assuming D to be constant, the expected result is just mulitplied by diff_coeff
        for (unsigned k = 0; k< number_of_time_steps; k++)
        {
            pseudo_ecg = calculator.ComputePseudoEcgAtOneTimeStep(k);
            TS_ASSERT_DELTA(pseudo_ecg, diff_coeff*expected_result,1e-6);
        }

    }

 };


#endif /*TESTPSEUDOECGCALCULATOR_HPP_*/
