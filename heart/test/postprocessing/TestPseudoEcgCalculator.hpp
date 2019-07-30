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
#include "FileComparison.hpp"
#include "SimpleBathProblemSetup.hpp"
#include "BidomainWithBathProblem.hpp"

/* HOW_TO_TAG Cardiac/Post-processing
 * Compute pseudo-ECGs
 */
class TestPseudoEcgCalculator : public CxxTest::TestSuite
{

public:

    void TestCalculator1DLinearGradient()
    {

        //read in the 1D mesh, from 0 to 1
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        for (AbstractTetrahedralMesh<1,1>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            // A correct calculation depends on this not breaking!
            TS_ASSERT(!HeartRegionCode::IsRegionBath((*it).GetUnsignedAttribute()));
        }

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

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);

        ///////////////////////////////////////////////////
        // Now we compute the pseudo ECG. We set an electrode at x=15.
        ///////////////////////////////////////////////////
        ChastePoint<1> probe_electrode(15.0);

        PseudoEcgCalculator<1,1,1> calculator (mesh, probe_electrode,
                                               FileFinder("hdf5",RelativeTo::ChasteTestOutput),
                                               "gradient_V");
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
        {
            FileComparison comparer(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/PseudoEcgFromElectrodeAt_15_0_0.dat",
                                    "heart/test/data/ValidPseudoEcg1D.dat");
            TS_ASSERT(comparer.CompareFiles());
        }

        {
            // Time-step striding
            PseudoEcgCalculator<1,1,1> course_calculator (mesh, probe_electrode,
                                                   FileFinder("hdf5",RelativeTo::ChasteTestOutput),
                                                   "gradient_V", "V" /*Voltage name*/, 2 /*Time stride*/);
            course_calculator.WritePseudoEcg();
            FileComparison comparer(OutputFileHandler::GetChasteTestOutputDirectory() + output_dir + "/PseudoEcgFromElectrodeAt_15_0_0.dat",
                    "heart/test/data/CoarsePseudoEcg1D.dat");
            TS_ASSERT(comparer.CompareFiles());
        }

        ChastePoint<1> bad_probe_electrode(0.0021132486540519);
        PseudoEcgCalculator<1,1,1> bad_calculator (mesh,
                                                   bad_probe_electrode,
                                                   FileFinder("hdf5",RelativeTo::ChasteTestOutput),
                                                   "gradient_V");
        TS_ASSERT_THROWS_THIS(bad_calculator.ComputePseudoEcgAtOneTimeStep(0), "Probe is on a mesh Gauss point.");

    }

    void TestCalculator1DParabolic()
    {
        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        DistributedTetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        for (AbstractTetrahedralMesh<1,1>::ElementIterator it = mesh.GetElementIteratorBegin();
             it != mesh.GetElementIteratorEnd();
             ++it)
        {
            // A correct calculation depends on this not breaking!
            TS_ASSERT(!HeartRegionCode::IsRegionBath((*it).GetUnsignedAttribute()));
        }

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

        PetscTools::Destroy(petsc_data_1);
        PetscTools::Destroy(petsc_data_2);

        ///////////////////////////////////////////////////
        // Now we compute the pseudo ECG. We set an electrode at x=15.
        ///////////////////////////////////////////////////

        ChastePoint<1> probe_electrode(15.0);

        PseudoEcgCalculator<1,1,1> calculator (mesh, probe_electrode,
                FileFinder("hdf5",RelativeTo::ChasteTestOutput), "parabolic_V", "V");

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

    /**
     * Test for BidomainWithBath problems.
     */
    void TestBathEcgCalculations()
    {
        HeartConfig::Instance()->SetSimulationDuration(2.0);  //ms - set to 500 to see whole AP trace.
        HeartConfig::Instance()->SetOutputDirectory("BidomainBath1d_PseudoEcg");
        HeartConfig::Instance()->SetOutputFilenamePrefix("bidomain_bath_1d");

        c_vector<double,1> centre;
        centre(0) = 0.5;
        BathCellFactory<1> cell_factory(-1e6, centre); // stimulates x=0.5 node

        BidomainWithBathProblem<1> bidomain_problem( &cell_factory );

        TrianglesMeshReader<1,1> reader("mesh/test/data/1D_0_to_1_100_elements");
        TetrahedralMesh<1,1> mesh;
        mesh.ConstructFromMeshReader(reader);

        // Set the x<0.25 and x>0.75 regions as the bath region
        for (unsigned i=0; i<mesh.GetNumElements(); i++)
        {
            double x = mesh.GetElement(i)->CalculateCentroid()[0];
            if ((x<0.25) || (x>0.75))
            {
                mesh.GetElement(i)->SetAttribute(HeartRegionCode::GetValidBathId());
            }
        }

        bidomain_problem.SetMesh(&mesh);
        bidomain_problem.Initialise();

        bidomain_problem.Solve();

        Vec sol = bidomain_problem.GetSolution();
        ReplicatableVector sol_repl(sol);

        // Test some Pseudo ECG calculation methods #2107
        ChastePoint<1> point1(0.0); // Point at which to calculate a pseudo-ECG - this is allowed.
        PseudoEcgCalculator<1,1,1>  ecg_calculator1(mesh, point1, FileFinder("BidomainBath1d_PseudoEcg",RelativeTo::ChasteTestOutput), "bidomain_bath_1d");

        // At the beginning we should have zero as voltage is uniform (in the tissue, this failed before #2107)
        double pseudo_ecg1_at_t_0 = ecg_calculator1.ComputePseudoEcgAtOneTimeStep(0);
        TS_ASSERT_DELTA(pseudo_ecg1_at_t_0, 0.0, 1e-9);

        ChastePoint<1> point2(0.3); // Point at which to calculate a pseudo-ECG - this should look decent
        PseudoEcgCalculator<1,1,1> ecg_calculator2(mesh, point2, FileFinder("BidomainBath1d_PseudoEcg",RelativeTo::ChasteTestOutput), "bidomain_bath_1d");

        // At the beginning we should have zero as voltage is uniform (in the tissue, this failed before #2107)
        double pseudo_ecg2_at_t_0 = ecg_calculator2.ComputePseudoEcgAtOneTimeStep(0);
        TS_ASSERT_DELTA(pseudo_ecg2_at_t_0, 0.0, 1e-9);

        // For completeness
        ecg_calculator1.WritePseudoEcg();
        ecg_calculator2.WritePseudoEcg();
    }

 };


#endif /*TESTPSEUDOECGCALCULATOR_HPP_*/
