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

#ifndef _TESTIMPEDANCEPROBLEM_HPP_
#define _TESTIMPEDANCEPROBLEM_HPP_

#include <cxxtest/TestSuite.h>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "SimpleImpedanceProblem.hpp"


class TestSimpleImpedanceProblem : public CxxTest::TestSuite
{
private:

public:
    void TestAcinarImpedance()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/y_branch_3d_mesh");
        mesh.ConstructFromMeshReader(mesh_reader);

        TS_ASSERT_THROWS_CONTAINS(SimpleImpedanceProblem(mesh, 1u), "Outlet node is not a boundary node");

        SimpleImpedanceProblem problem(mesh, 0u);
        problem.rGetMesh(); //for coverage

        unsigned node_index = 3; //Arbitrary terminal node

        problem.SetElastance(2*M_PI/2.0);

        TS_ASSERT_DELTA(real(problem.CalculateAcinusImpedance(mesh.GetNode(node_index), 0.0)), 0.0, 1e-6);
        TS_ASSERT_DELTA(real(problem.CalculateAcinusImpedance(mesh.GetNode(node_index), 1.0)), 0.0, 1e-6);
        TS_ASSERT_DELTA(imag(problem.CalculateAcinusImpedance(mesh.GetNode(node_index), 1.0)), -1.0, 1e-6);
    }

    void TestElementReactanceAndInertance()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/y_branch_3d_mesh");
        mesh.ConstructFromMeshReader(mesh_reader);

        SimpleImpedanceProblem problem(mesh, 0u);

        problem.SetRho(M_PI);
        problem.SetMu(M_PI);

        double l = 2.0;
        double r = 2.0;

        TS_ASSERT_DELTA(problem.CalculateElementResistance(r, l), 8*l/(r*r*r*r), 1e-6);
        TS_ASSERT_DELTA(problem.CalculateElementInertance(r, l), l/(r*r), 1e-6);
    }

    void TestSolveTreeImpedance()
    {
        TetrahedralMesh<1,3> mesh;
        TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/y_branch_3d_mesh");
        mesh.ConstructFromMeshReader(mesh_reader);

        SimpleImpedanceProblem problem(mesh, 0u);

        //Pure resistance calculation
        problem.SetRho(0.0);
        problem.SetMu(M_PI);
        problem.SetElastance(0.0);
        problem.SetFrequency(0.0);

        problem.Solve();

        double small_r = (mesh.GetNode(2u)->rGetNodeAttributes()[0] + mesh.GetNode(1u)->rGetNodeAttributes()[0])/2.0;
        double small_l = sqrt(2.0);
        double big_r = (mesh.GetNode(1u)->rGetNodeAttributes()[0] + mesh.GetNode(0u)->rGetNodeAttributes()[0])/2.0;
        double big_l = 1.0;

        double small_Raw = problem.CalculateElementResistance(small_r, small_l);
        double big_Raw = problem.CalculateElementResistance(big_r, big_l);

        TS_ASSERT_DELTA(real(problem.GetImpedance()), big_Raw + 1.0/(2.0/small_Raw), 1e-6);
        TS_ASSERT_DELTA(imag(problem.GetImpedance()), 0.0, 1e-6);

        //Pure airway inertance calculation (no acinus elastance)
        problem.SetRho(1);
        problem.SetMu(0.0);
        problem.SetElastance(0.0);
        problem.SetFrequency(1.0);

        problem.Solve();

        double small_Ini = problem.CalculateElementInertance(small_r, small_l);
        double big_Ini = problem.CalculateElementInertance(big_r, big_l);

        TS_ASSERT_DELTA(real(problem.GetImpedance()), 0.0, 1e-6);
        TS_ASSERT_DELTA(imag(problem.GetImpedance()), 2*M_PI*big_Ini + 1.0/(2.0/(2*M_PI*small_Ini)), 1e-6);

        //Pure acinar elastance calculation
        problem.SetRho(0.0);
        problem.SetMu(0.0);
        problem.SetElastance(4*M_PI);
        problem.SetFrequency(1.0);

        problem.Solve();

        std::complex<double> i(0, 1);
        std::complex<double> acinus_elastance = -i;
        std::complex<double> total_elastance = 2.0/(1.0/acinus_elastance);

        TS_ASSERT_DELTA(real(problem.GetImpedance()), 0.0, 1e-6);
        TS_ASSERT_DELTA(imag(problem.GetImpedance()), imag(total_elastance), 1e-6);

        //Pure zero flow calculation - check that divide by zero errors are delt with.
        problem.SetRho(0.0);
        problem.SetMu(0.0);
        problem.SetElastance(0.0);
        problem.SetFrequency(0.0);

        problem.Solve();

        TS_ASSERT_DELTA(real(problem.GetImpedance()), 0.0, 1e-6);
        TS_ASSERT_DELTA(imag(problem.GetImpedance()), 0.0, 1e-6);
    }

    void TestMultipleFrequencies()
    {
        TetrahedralMesh<1,3> mesh;
        //TrianglesMeshReader<1,3> mesh_reader("mesh/test/data/y_branch_3d_mesh");
        TrianglesMeshReader<1,3> mesh_reader("lung/test/data/TestSubject002");
        mesh.ConstructFromMeshReader(mesh_reader);

        //Scale all radii by 0.7 to give an FRC equivalent lung
        for (TetrahedralMesh<1,3>::NodeIterator node_iter = mesh.GetNodeIteratorBegin();
             node_iter != mesh.GetNodeIteratorEnd();
             ++node_iter)
        {
            node_iter->rGetNodeAttributes()[0] *= 0.7;
        }

        std::vector<double> test_frequencies;
        test_frequencies.push_back(1.0);
        test_frequencies.push_back(2.0);
        test_frequencies.push_back(3.0);
        test_frequencies.push_back(5.0);
        test_frequencies.push_back(10.0);
        test_frequencies.push_back(20.0);
        test_frequencies.push_back(30.0);

        SimpleImpedanceProblem problem(mesh, 0u);
        problem.SetMeshInMilliMetres();
        problem.SetFrequencies(test_frequencies);               //Set & get frequencies for coverage
        std::vector<double>& freqs = problem.rGetFrequencies();

        TS_ASSERT_EQUALS(freqs.size(), 7u);

        problem.Solve();

        std::vector<std::complex<double> > impedances = problem.rGetImpedances();

        TS_ASSERT_EQUALS(impedances.size(), 7u);

        //These are hard coded from previous runs, but are as expected for
        //a patient with moderate to severe asthma
        TS_ASSERT_DELTA(real(impedances[0])*1e-3/98, 8.45, 1e-2);
        TS_ASSERT_DELTA(imag(impedances[0])*1e-3/98, -3.65, 1e-2);
        TS_ASSERT_DELTA(real(impedances[6])*1e-3/98, 5.77, 1e-2);
        TS_ASSERT_DELTA(imag(impedances[6])*1e-3/98, 4.12, 1e-2);
    }
};

#endif /*_TESTIMPEDANCEPROBLEM_HPP_*/
