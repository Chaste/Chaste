/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef _TESTVENTILATIONPROBLEM_HPP_
#define _TESTVENTILATIONPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include <queue>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "ReplicatableVector.hpp"
#include "VentilationProblem.hpp"

// Pressures read from file.  Pressures labelled P10, P20, P21, P30, P31, P32, P33 map to the mesh
// nodes 1, 2, 3, 4, 5, 6, 7 (respectively).
std::vector<double> pressureAt1;
std::vector<double> pressureAt2;
std::vector<double> pressureAt3;
std::vector<double> pressureAt4;
std::vector<double> pressureAt5;
std::vector<double> pressureAt6;
std::vector<double> pressureAt7;

void LinearTimeBCs(VentilationProblem* pProblem, double time)
{
    pProblem->SetConstantInflowPressures(15*time);
}

void SineBCs(VentilationProblem* pProblem, double time)
{
    pProblem->SetConstantInflowPressures(15.0*sin(time*5.0/(2.0*M_PI)));
}
void FileBCs(VentilationProblem* pProblem, double time)
{
    unsigned timestep= (unsigned) floor(time*100.0+0.5);
    pProblem->SetConstantInflowPressures(pressureAt7[timestep]);
}
class TestVentilationProblem : public CxxTest::TestSuite
{
private:

    void ReadDataFromFile(const std::string& filePath)
    {
        pressureAt1.clear();
        pressureAt2.clear();
        pressureAt3.clear();
        pressureAt4.clear();
        pressureAt5.clear();
        pressureAt6.clear();
        pressureAt7.clear();
        std::ifstream filestream(filePath.c_str());
        TS_ASSERT(filestream.is_open());
        //Skip header line
        std::string line;
        getline(filestream, line);
        while (filestream.good())
        {
            double data;
            filestream >> data; //Time
            if (filestream.eof())
            {
                break;
            }
            filestream >> data; //Pressure at root
            TS_ASSERT_DELTA(data, 0.0, 1e-15);
            filestream >> data;
            pressureAt1.push_back(data);
            filestream >> data;
            pressureAt2.push_back(data);
            filestream >> data;
            pressureAt3.push_back(data);
            filestream >> data;
            pressureAt4.push_back(data);
            filestream >> data;
            pressureAt5.push_back(data);
            filestream >> data;
            pressureAt6.push_back(data);
            filestream >> data;
            pressureAt7.push_back(data);
        }
        filestream.close();
    }
public:
    void TestVentilationProblemOnBranch() throw(Exception)
    {
        VentilationProblem problem("mesh/test/data/y_branch_3d_mesh", 0u);
        TS_ASSERT_DELTA(problem.GetViscosity(), 1.92e-8, 1e-12);
        problem.SetViscosity(1.0);
        problem.SetOutflowPressure(100);
        problem.SetConstantInflowPressures(0.0);
        problem.Solve();

        ReplicatableVector solution_vector_repl( problem.GetSolution());
        TS_ASSERT_DELTA(solution_vector_repl[0], 92813610, 1);
        TS_ASSERT_DELTA(solution_vector_repl[1], 46406805, 1);
        TS_ASSERT_DELTA(solution_vector_repl[2], 46406805, 1);
        TS_ASSERT_DELTA(solution_vector_repl[3], 100.0, 1e-6); //BC
        TS_ASSERT_DELTA(solution_vector_repl[4], 84.5107, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[5], 0.0, 1e-6); //BC
        TS_ASSERT_DELTA(solution_vector_repl[6], 0.0, 1e-6); //BC

#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "small_conical");
#endif
    }
    void TestVentilationProblemCylindrical() throw(Exception)
    {
        VentilationProblem problem("mesh/test/data/y_branch_3d_mesh");
        //These values are equivalent to Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
        TS_ASSERT_DELTA(problem.GetViscosity(), 1.92e-8, 1e-12);
        problem.SetViscosity(1.0);
        TS_ASSERT_DELTA(problem.GetDensity(), 1.511e-9, 1e-12);
        problem.SetDensity(1.0);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(100);
        problem.SetConstantInflowPressures(0.0);
        problem.Solve();

        ReplicatableVector solution_vector_repl( problem.GetSolution());
        TS_ASSERT_DELTA(solution_vector_repl[0], 19932007, 1);
        TS_ASSERT_DELTA(solution_vector_repl[1], 9966003, 1);
        TS_ASSERT_DELTA(solution_vector_repl[2], 9966003, 1);
        TS_ASSERT_DELTA(solution_vector_repl[3], 100.0, 1e-6); //BC
        TS_ASSERT_DELTA(solution_vector_repl[4], 91.8789, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[5], 0.0, 1e-6); //BC
        TS_ASSERT_DELTA(solution_vector_repl[6], 0.0, 1e-6); //BC

#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "small_cylindrical");
#endif
    }
private:
public:
    void TestThreeBifurcationsWithRadiusOnEdgeFile() throw (Exception)
    {
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(0.00148608);
        problem.Solve();
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations");
#endif
        const unsigned num_edge = problem.rGetMesh().GetNumElements();
        ReplicatableVector solution_vector_repl( problem.GetSolution());
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+1], 0.0006,   1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+2], 0.001274, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+3], 0.001274, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+4], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+5], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+6], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+7], 0.00148608, 1e-8); //BC
    }

    void TestThreeBifurcationsWithRadiusOnNodeFile() throw (Exception)
    {
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(15);
        problem.Solve();
        const unsigned num_edge = problem.rGetMesh().GetNumElements();
        ReplicatableVector solution_vector_repl( problem.GetSolution());
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+1], 6.6666,   1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+2], 12.2222, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+3], 12.2222, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+4], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+5], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+6], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+7], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[0], -284.0705, 1e-4); // (Outflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[3],  -71.0176, 1e-4); // (Inflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[4],  -71.0176, 1e-4); // (Inflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[5],  -71.0176, 1e-4); // (Inflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[6],  -71.0176, 1e-4); // (Inflow flux)
    }

    void TestThreeBifurcationsWithRadiusOnNodeFileFluxBoundaries() throw (Exception)
    {
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowFluxes(-71.0176);
        problem.Solve();
        const unsigned num_edge = problem.rGetMesh().GetNumElements();
        ReplicatableVector solution_vector_repl( problem.GetSolution());
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+1], 6.6666,   1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+2], 12.2222, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+3], 12.2222, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+4], 15, 1e-4); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+5], 15, 1e-4); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+6], 15, 1e-4); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+7], 15, 1e-4); //BC
        TS_ASSERT_DELTA(solution_vector_repl[0], -284.0705, 1e-4); // (Outflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[3],  -71.0176, 1e-8); // BC (Inflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[4],  -71.0176, 1e-8); // BC (Inflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[5],  -71.0176, 1e-8); // BC (Inflow flux)
        TS_ASSERT_DELTA(solution_vector_repl[6],  -71.0176, 1e-8); // BC (Inflow flux)
    }

    void TestThreeBifurcationsWithDynamicResistance() throw (Exception)
    {
        /*
         * HOW_TO_TAG Continuum mechanics/Ventilation
         * Solve a simple ventilation problem with no time variation.
         *
         * Pressure is prescribed at the root (trachea) and at the leaves (acini).
         * Dynamic (Pedley) resistance is used at higher Reynolds numbers.
         * There is no coupled acinus compliance model in this version.
         */
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(15);
        problem.SetDynamicResistance();
        problem.Solve();
        const unsigned num_edge = problem.rGetMesh().GetNumElements();
        ReplicatableVector solution_vector_repl( problem.GetSolution());
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+1], 6.6878,   1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+2], 12.2292, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+3], 12.2292, 1e-4);
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+4], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+5], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+6], 15, 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+7], 15, 1e-8); //BC
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_pedley");
#endif
    }

    void TestTimeVaryingThreeBifurcations() throw (Exception)
    {
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 1.0, 0.1);
        problem.Solve(stepper, &LinearTimeBCs, "TestVentilation", "three_bifurcations_time");
    }

    void TestSineThreeBifurcations() throw (Exception)
    {
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 25.0, 0.1);
        problem.Solve(stepper, &SineBCs, "TestVentilation", "three_bifurcations_sine");
    }
    /*
     * HOW_TO_TAG Continuum mechanics/Ventilation
     * Solve a simple ventilation problem defined in a file.
     *
     * The file prescribes pressure at the root (trachea) and at the leaves (acini).
     * The file sets density and viscosity
     *
     * Output file names are given to the function (but might instead be given in the file)
     *
     */
    void TestSolveProblemDefinedInFile() throw (Exception)
    {
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations");
        problem.SolveProblemFromFile("continuum_mechanics/test/data/ChasteVentilationInput.txt", "VentilationOutput", "3_bifurcations");
    }
    void TestReadFile()  throw (Exception)
    {
        ReadDataFromFile("continuum_mechanics/test/data/Pdata-Normal.txt");
        TS_ASSERT_EQUALS(pressureAt1.size(), 2501u);
        TS_ASSERT_DELTA(pressureAt1[2500], 0.703592, 1e-6);
        TS_ASSERT_DELTA(pressureAt2[2500], 1.05539, 1e-5);
        TS_ASSERT_DELTA(pressureAt4[2500], 1.23129, 1e-5);
    }
    void TestNormalAgainstFile() throw (Exception)
    {
        ReadDataFromFile("continuum_mechanics/test/data/Pdata-Normal.txt");
        TS_ASSERT_EQUALS(pressureAt1.size(), 2501u);
        TS_ASSERT_DELTA(pressureAt1[2500], 0.703592, 1e-6);
        TS_ASSERT_DELTA(pressureAt2[2500], 1.05539, 1e-5);
        TS_ASSERT_DELTA(pressureAt4[2500], 1.23129, 1e-5);
        TS_ASSERT_DELTA(pressureAt7[2500], 1.23129, 1e-5);
        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 25.0, 0.01);
        problem.Solve(stepper, &FileBCs, "TestVentilation", "three_bifurcations_normal");
        const unsigned num_edge = problem.rGetMesh().GetNumElements();
        ReplicatableVector solution_vector_repl( problem.GetSolution());//check pressure at time @ 25
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+0], 0.0, 1e-8); //BC
        //Note: the two implementations are using different resistance models, so are not yet the same
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+1], pressureAt1[2500], 0.2); ///\todo #2300
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+2], pressureAt2[2500], 0.1); ///\todo #2300
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+3], pressureAt3[2500], 0.1); ///\todo #2300
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+4], pressureAt4[2500], 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+5], pressureAt5[2500], 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+6], pressureAt6[2500], 1e-8); //BC
        TS_ASSERT_DELTA(solution_vector_repl[num_edge+7], pressureAt7[2500], 1e-8); //BC
    }
//    void TestAnnotateGenerations() throw (Exception)
//    {
//        EXIT_IF_PARALLEL;
//        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
//        Node<3>* p_root_node = problem.rGetMesh().GetNode(0);
//        TS_ASSERT_EQUALS(p_root_node->rGetNodeAttributes().size(), 1u);
//       // p_root_node->Set
//        std::map<int, int> max_branch_at_generation;
//        std::queue<Node<3>*> dfs_queue;
//        dfs_queue.push(p_root_node);
//
//        while (!dfs_queue.empty())
//        {
//            std::cout<<dfs_queue.front()->GetIndex()<<"\n";
//            dfs_queue.pop();
//        }
//    }
    void TestExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(VentilationProblem bad_problem("mesh/test/data/y_branch_3d_mesh", 1u),
                "Outlet node is not a boundary node");

        VentilationProblem problem("continuum_mechanics/test/data/three_bifurcations");
        TS_ASSERT_THROWS_THIS(problem.SolveProblemFromFile("DoesNotExist.txt", "out", "out"), "Could not open file DoesNotExist.txt");

    }

};

#endif /*_TESTVENTILATIONPROBLEM_HPP_*/
