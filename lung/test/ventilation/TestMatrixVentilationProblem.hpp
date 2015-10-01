/*

Copyright (c) 2005-2015, University of Oxford.
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

#ifndef _TESTMATRIXVENTILATIONPROBLEM_HPP_
#define _TESTMATRIXVENTILATIONPROBLEM_HPP_

#include <cxxtest/TestSuite.h>
#include <queue>

#include "TetrahedralMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "LinearSystem.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "MatrixVentilationProblem.hpp"

// Pressures read from file.  Pressures labelled P10, P20, P21, P30, P31, P32, P33 map to the mesh
// nodes 1, 2, 3, 4, 5, 6, 7 (respectively).
std::vector<double> pressureAt1;
std::vector<double> pressureAt2;
std::vector<double> pressureAt3;
std::vector<double> pressureAt4;
std::vector<double> pressureAt5;
std::vector<double> pressureAt6;
std::vector<double> pressureAt7;

void LinearTimeBCs(MatrixVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double time = rTimeStepper.GetTime();
    pProblem->SetPressureAtBoundaryNode(rNode, 15*time);
}

void SineBCs(MatrixVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double time = rTimeStepper.GetTime();
    pProblem->SetPressureAtBoundaryNode(rNode, 15.0*sin(time*5.0/(2.0*M_PI)));
}

void FileBCs(MatrixVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double time = rTimeStepper.GetTime();

    unsigned timestep= (unsigned) floor(time*100.0+0.5);
    pProblem->SetPressureAtBoundaryNode(rNode, pressureAt7[timestep]);
}

void GravitationalBCs(MatrixVentilationProblem* pProblem, TimeStepper& rTimeStepper, const Node<3>& rNode)
{
    double x_max =  6.0;
    double x_min = -6.0;
    double delta_p = 0.5;

    double x = rNode.rGetLocation()[0];

    double pressure = delta_p * (x - x_min)/(x_max - x_min);

    pProblem->SetPressureAtBoundaryNode(rNode, pressure);
}

class TestMatrixVentilationProblem : public CxxTest::TestSuite
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
    void TestMatrixVentilationProblemOnBranch() throw(Exception)
    {
        MatrixVentilationProblem problem("mesh/test/data/y_branch_3d_mesh", 0u);
        problem.SetMeshInMilliMetres();
        TS_ASSERT_DELTA(problem.GetViscosity(), 1.92e-5, 1e-12);
        problem.SetViscosity(1.0);
        problem.SetOutflowPressure(100);
        problem.SetConstantInflowPressures(0.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(flux[0], 0.092813610, 1e-8);
        TS_ASSERT_DELTA(flux[1], 0.046406805, 1e-8);
        TS_ASSERT_DELTA(flux[2], 0.046406805, 1e-8);
        TS_ASSERT_DELTA(pressure[0], 100.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[1], 84.5107, 1e-4);
        TS_ASSERT_DELTA(pressure[2], 0.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[3], 0.0, 1e-6); //BC
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "small_conical");
#endif
    }
    void TestMatrixVentilationProblemCylindrical() throw(Exception)
    {
        MatrixVentilationProblem problem("mesh/test/data/y_branch_3d_mesh");
        problem.SetMeshInMilliMetres();
        //These values are equivalent to Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
        TS_ASSERT_DELTA(problem.GetViscosity(), 1.92e-5, 1e-12);
        problem.SetViscosity(1.0);
        TS_ASSERT_DELTA(problem.GetDensity(), 1.15, 1e-10);
        problem.SetDensity(1.0);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(100);
        problem.SetConstantInflowPressures(0.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(flux[0], 0.019932007, 1e-8);
        TS_ASSERT_DELTA(flux[1], 0.009966003, 1e-8);
        TS_ASSERT_DELTA(flux[2], 0.009966003, 1e-8);
        TS_ASSERT_DELTA(pressure[0], 100.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[1], 91.8789, 1e-4);
        TS_ASSERT_DELTA(pressure[2], 0.0, 1e-6); //BC
        TS_ASSERT_DELTA(pressure[3], 0.0, 1e-6); //BC

#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "small_cylindrical");
#endif
    }

    void TestThreeBifurcationsWithRadiusOnEdgeFile() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(0.00148608);
        problem.Solve();
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations");
#endif
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 0.0006,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 0.001274, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 0.001274, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 0.00148608, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 0.00148608, 1e-8); //BC

        // Total flux
        // TS_ASSERT_DELTA(flux[0], -2.8143e-5, 1e-8); // Mesh in meters
        TS_ASSERT_DELTA(flux[0], -2.8143e-14, 1e-15); //Mesh in millimeters
    }

    void TestThreeBifurcationsWithRadiusOnNodeFile() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(15.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.66666,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 12.22223, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 12.22222, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 15.0, 1e-8); //BC
        TS_ASSERT_DELTA(flux[0], -2.8407e-10 , 1e-13); // (Outflow flux)
        TS_ASSERT_DELTA(flux[3], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[4], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[5], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[6], -7.102e-11, 1e-13); // (Inflow flux)
    }

    void TestThreeBifurcationsExtraLinkWithRadiusOnNodeFile() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations_extra_links", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0 + 1.0);
        problem.SetConstantInflowPressures(15.0 + 1.0);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0 + 1.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.66666 + 1.0,   1e-4);
        TS_ASSERT_DELTA(pressure[2], 12.22223 + 1.0, 1e-4);
        TS_ASSERT_DELTA(pressure[3], 12.22223 + 1.0, 1e-4);
        TS_ASSERT_DELTA(pressure[4], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(pressure[5], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(pressure[6], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(pressure[7], 15 + 1.0, 1e-7); //BC
        TS_ASSERT_DELTA(flux[0], 2.8407e-10, 1e-4); // (Outflow flux)
        TS_ASSERT_DELTA(flux[10], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[11], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[12], -7.102e-11, 1e-13); // (Inflow flux)
        TS_ASSERT_DELTA(flux[13], -7.102e-11, 1e-13); // (Inflow flux)
        //This is the extra node at the Trachea
        TS_ASSERT_DELTA(pressure[8], 2.22225 + 1.0, 1e-4); //Between root and first bifurcation

#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_extra_links");
#endif
    }

    void TestThreeBifurcationsWithRadiusOnNodeFileFluxBoundaries() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowFluxes(-7.10176e-11);
        problem.Solve();

        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 6.6666,   2e-3);
        TS_ASSERT_DELTA(pressure[2], 12.2222, 2e-4);
        TS_ASSERT_DELTA(pressure[3], 12.2222, 2e-4);
        TS_ASSERT_DELTA(pressure[4], 15, 2e-4);
        TS_ASSERT_DELTA(pressure[5], 15, 2e-4);
        TS_ASSERT_DELTA(pressure[6], 15, 2e-4);
        TS_ASSERT_DELTA(pressure[7], 15, 2e-4);
        TS_ASSERT_DELTA(flux[0], -2.840704e-10, 1e-13); // (Outflow flux)
        TS_ASSERT_DELTA(flux[3],  -7.10176e-11, 1e-16); // BC (Inflow flux)
        TS_ASSERT_DELTA(flux[4],  -7.10176e-11, 1e-16); // BC (Inflow flux)
        TS_ASSERT_DELTA(flux[5],  -7.10176e-11, 1e-16); // BC (Inflow flux)
        TS_ASSERT_DELTA(flux[6],  -7.10176e-11, 1e-16); // BC (Inflow flux)
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
        EXIT_IF_PARALLEL; ///\todo #2300 There is a problem with the Windows parallel implementation
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetMeshInMilliMetres();
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(150000); //Needed to increase the resistance in these artificial airways
        problem.SetDynamicResistance();
        problem.Solve();
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[1], 91108.7409,   1e-1);
        TS_ASSERT_DELTA(pressure[2], 132694.0014, 1e-2);
        TS_ASSERT_DELTA(pressure[3], 132694.0014, 1e-2);
        TS_ASSERT_DELTA(pressure[4], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], 1.5e5, 1e-8); //BC
        TS_ASSERT_DELTA(flux[6], -4.424511e-7, 1e-11);
#ifdef CHASTE_VTK
        problem.WriteVtk("TestVentilation", "three_bifurcations_pedley");
#endif
    }

    void TestTimeVaryingThreeBifurcations() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 1.0, 0.1);
        problem.Solve(stepper, &LinearTimeBCs, "TestVentilation", "three_bifurcations_time");
    }

    void TestSineThreeBifurcations() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 25.0, 0.1);
        problem.Solve(stepper, &SineBCs, "TestVentilation", "three_bifurcations_sine");
    }

    void TestGravitationalVaryingThreeBifurcations() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 1.0, 0.1);
        problem.Solve(stepper, &GravitationalBCs, "TestVentilation", "three_bifurcations_gravity");
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
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations");
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
        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations", 0u);
        problem.SetRadiusOnEdge();
        problem.SetOutflowPressure(0.0);
        TimeStepper stepper(0.0, 25.0, 0.01);
        problem.Solve(stepper, &FileBCs, "TestVentilation", "three_bifurcations_normal");
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure); //check pressure at time @ 25
        TS_ASSERT_DELTA(pressure[0], 0.0, 1e-8); //BC
        //Note: the two implementations are using different resistance models, so are not yet the same
        TS_ASSERT_DELTA(pressure[1], pressureAt1[2500], 0.2); ///\todo #2300
        TS_ASSERT_DELTA(pressure[2], pressureAt2[2500], 0.1); ///\todo #2300
        TS_ASSERT_DELTA(pressure[3], pressureAt3[2500], 0.1); ///\todo #2300
        TS_ASSERT_DELTA(pressure[4], pressureAt4[2500], 1e-8); //BC
        TS_ASSERT_DELTA(pressure[5], pressureAt5[2500], 1e-8); //BC
        TS_ASSERT_DELTA(pressure[6], pressureAt6[2500], 1e-8); //BC
        TS_ASSERT_DELTA(pressure[7], pressureAt7[2500], 1e-8); //BC
    }

    void TestTopOfAirwaysPatientData() throw (Exception)
    {
        MatrixVentilationProblem problem("continuum_mechanics/test/data/top_of_tree", 0u);
        PetscOptionsSetValue("-ksp_monitor", "");
        problem.SetOutflowPressure(0.0);
        problem.SetConstantInflowPressures(50.0);
        //problem.SetConstantInflowFluxes(100.0);
        TetrahedralMesh<1, 3>& r_mesh=problem.rGetMesh();
        TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 31u);
        TS_ASSERT_EQUALS(r_mesh.GetNumElements(), 30u);
        problem.Solve();
        problem.Solve();
        // For debugging...
        std::vector<double> flux, pressure;
        problem.GetSolutionAsFluxesAndPressures(flux, pressure);
    }

    void OnlyWorksWithUMFPACKTestPatientData() throw (Exception)
    {
        MatrixVentilationProblem problem("notforrelease_lung/test/data/Novartis002", 0u);
        PetscOptionsSetValue("-ksp_monitor", "");

        problem.SetOutflowPressure(0.);
        problem.SetConstantInflowPressures(50.0);
//        problem.SetConstantInflowFluxes(100.0);
        //TetrahedralMesh<1, 3>& r_mesh=problem.rGetMesh();
       // TS_ASSERT_EQUALS(r_mesh.GetNumNodes(), 31u);
     //   TS_ASSERT_EQUALS(r_mesh.GetNumElements(), 30u);
        //problem.SetRadiusOnEdge();
        problem.SetMeshInMilliMetres();
        //problem.SetDynamicResistance();
        problem.Solve();

//        std::vector<double> flux, pressure;
//        problem.GetSolutionAsFluxesAndPressures(flux, pressure); //check pressure at time @ 25
//
//        VtkMeshWriter<1, 3> vtk_writer("TestMatrixVentilation", "flow_solution", false);
//        vtk_writer.AddCellData("Flux", flux);
//        vtk_writer.AddPointData("Pressure", pressure);
//
//        std::vector<double> radii(problem.rGetMesh().GetNumElements());
//        for(TetrahedralMesh<1,3>::ElementIterator iter = problem.rGetMesh().GetElementIteratorBegin();
//            iter != problem.rGetMesh().GetElementIteratorEnd();
//            ++iter)
//        {
//            radii[(*iter).GetIndex()] = (*iter).GetAttribute();
//        }
//
//        vtk_writer.AddCellData("Radii", radii);
//        vtk_writer.WriteFilesUsingMesh(problem.rGetMesh());
    }
    void TestExceptions() throw(Exception)
    {
        TS_ASSERT_THROWS_THIS(MatrixVentilationProblem bad_problem("mesh/test/data/y_branch_3d_mesh", 1u),
                "Outlet node is not a boundary node");

        MatrixVentilationProblem problem("continuum_mechanics/test/data/three_bifurcations");
        TS_ASSERT_THROWS_THIS(problem.SolveProblemFromFile("DoesNotExist.txt", "out", "out"), "Could not open file DoesNotExist.txt");

        TS_ASSERT_THROWS_THIS(problem.SetPressureAtBoundaryNode(3u, 0.0), "Boundary conditions cannot be set at internal nodes");
        TS_ASSERT_THROWS_THIS(problem.SetFluxAtBoundaryNode(3u, 0.0), "Boundary conditions cannot be set at internal nodes");

    }
};

#endif /*_TESTMATRIXVENTILATIONPROBLEM_HPP_*/
