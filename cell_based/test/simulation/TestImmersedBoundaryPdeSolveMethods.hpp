/*

Copyright (c) 2005-2023, University of Oxford.
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


#ifndef TESTIMMERSEDBOUNDARYPDESOLVEMETHODS_HPP_
#define TESTIMMERSEDBOUNDARYPDESOLVEMETHODS_HPP_

// Needed for test framework
#include <cxxtest/TestSuite.h>
#include "AbstractCellBasedTestSuite.hpp"

#define BOOST_DISABLE_ASSERTS
#include "boost/multi_array.hpp"

// Includes from trunk
#include "CellsGenerator.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "UniformCellCycleModel.hpp"

// Includes from projects/ImmersedBoundary
#include "ImmersedBoundaryCellPopulation.hpp"
#include "ImmersedBoundaryMesh.hpp"
#include "ImmersedBoundaryMeshWriter.hpp"
#include "ImmersedBoundaryMeshReader.hpp"
#include "ImmersedBoundarySimulationModifier.hpp"
#include "ImmersedBoundaryPalisadeMeshGenerator.hpp"
#include "SuperellipseGenerator.hpp"
#include "ImmersedBoundaryLinearMembraneForce.hpp"

// This test is never run in parallel
#include "FakePetscSetup.hpp"

///\todo why are all these tests commented out?
class TestImmersedBoundaryPdeSolveMethods : public AbstractCellBasedTestSuite
{
public:

    void TestForcePropagation()
    {
//        // Create a vector of nodes forming a rectangle in (0,1)x(0,1)
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, true, 0.45, 0.25));
//        nodes.push_back(new Node<2>(1, true, 0.65, 0.25));
//        nodes.push_back(new Node<2>(2, true, 0.65, 0.75));
//        nodes.push_back(new Node<2>(3, true, 0.65, 0.85));
//        nodes.push_back(new Node<2>(4, true, 0.45, 0.85));
//
//        // Create a vector of immersed boundary elements and create an element with the nodes above
//        std::vector<ImmersedBoundaryElement<2,2>*> ib_element;
//        ib_element.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));
//
//        // Create a mesh with the nodes and elements vectors
//        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_element, 10, 10);
//
//        // Set element parameters
//        ImmersedBoundaryElement<2,2>* p_elem = p_mesh->GetElement(0u);
//        p_elem->SetMembraneSpringConstant(1.0);
//        p_elem->SetMembraneRestLength(0.4);
//
//        p_mesh->SetNumGridPtsX(10);
//        p_mesh->SetNumGridPtsY(10);
//
//        std::vector<CellPtr> cells;
//        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
//        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
//
//        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
//
//        OffLatticeSimulation<2> simulator(cell_population);
//
//        // Add main immersed boundary simulation modifier
//        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_mod);
//        simulator.AddSimulationModifier(p_mod);
//
//        // Add force laws
//        MAKE_PTR(ImmersedBoundaryLinearMembraneForce<2>, p_boundary_force);
//        p_mod->AddImmersedBoundaryForce(p_boundary_force);
//
//        /*
//         * We first test that the perimeter elasticity forces are calculated correctly. We do this by invoking
//         * ImmersedBoundarySimulationModifier::AddImmersedBoundaryForceContributions() and checking
//         * against hand-calculated values.
//         */
//        p_mod->AddImmersedBoundaryForceContributions();
//        c_vector<double,2> force_on_node;
//
//        force_on_node = nodes[0]->rGetAppliedForce();
//        TS_ASSERT_DELTA(force_on_node[0], -0.2, 1e-10);
//        TS_ASSERT_DELTA(force_on_node[1],  0.0, 1e-10); // 0.0 due to GetVectorFromAtoB - unrealistic node spacing
//
//        force_on_node = nodes[1]->rGetAppliedForce();
//        TS_ASSERT_DELTA(force_on_node[0],  0.2, 1e-10);
//        TS_ASSERT_DELTA(force_on_node[1],  0.1, 1e-10);
//
//        force_on_node = nodes[2]->rGetAppliedForce();
//        TS_ASSERT_DELTA(force_on_node[0],  0.0, 1e-10);
//        TS_ASSERT_DELTA(force_on_node[1], -0.4, 1e-10);
//
//        force_on_node = nodes[3]->rGetAppliedForce();
//        TS_ASSERT_DELTA(force_on_node[0],  0.2, 1e-10);
//        TS_ASSERT_DELTA(force_on_node[1],  0.3, 1e-10);
//
//        force_on_node = nodes[4]->rGetAppliedForce();
//        TS_ASSERT_DELTA(force_on_node[0], -0.2, 1e-10);
//        TS_ASSERT_DELTA(force_on_node[1],  0.0, 1e-10); // 0.0 due to GetVectorFromAtoB - unrealistic node spacing
//
//        /*
//         * Next, we check that the forces are spread properly onto the fluid grids.  First
//         * we clear forces and add known forces to particular nodes and check their propagation
//         * onto the fluid grid by comparing with hand-calculated values
//         */
//        p_mod->ClearForces();
//        force_on_node[0] = 1.0;
//        force_on_node[1] = 2.0;
//        nodes[0]->AddAppliedForceContribution(force_on_node);
//
//        // Propagate the forces to the fluid force grids and grab references to them
//        p_mod->PropagateForcesToFluidGrid();
//        std::vector<std::vector<double> > force_grid_x = p_mod->mFluidForceGridX;
//        std::vector<std::vector<double> > force_grid_y = p_mod->mFluidForceGridY;
//
//        // First check that there is no force propagated outside the sphere of influence of node 0
//        for (unsigned y = 0; y < p_mod->mNumGridPtsY; y++)
//        {
//            for (unsigned x = 0; x < p_mod->mNumGridPtsX; x++)
//            {
//                if ( x < 3 || x > 6)
//                {
//                    if ( y < 1 || y > 4)
//                    {
//                        TS_ASSERT_LESS_THAN(fabs(force_grid_x[y][x]), 1e-10);
//                        TS_ASSERT_LESS_THAN(fabs(force_grid_x[y][x]), 1e-10);
//                    }
//                }
//            }
//        }
//
//        // As the system is symmetric, we only need two distances for the 1D delta function
//        double spacing = 1.0 / 10.0;
//
//        double small_delta = (0.25 * (1.0 + cos(M_PI * 0.05 / (2 * spacing)))) / spacing;
//        double large_delta = (0.25 * (1.0 + cos(M_PI * 0.15 / (2 * spacing)))) / spacing;
//
//        // Test 1D delta function
//        TS_ASSERT_DELTA(small_delta, p_mod->Delta1D(0.05, 0.1), 1e-10);
//        TS_ASSERT_DELTA(large_delta, p_mod->Delta1D(0.15, 0.1), 1e-10);
//
//        // Test first grid point
//        TS_ASSERT_DELTA(force_grid_x[1][3], 1.0 * large_delta * large_delta, 1e-10);
//        TS_ASSERT_DELTA(force_grid_y[1][3], 2.0 * large_delta * large_delta, 1e-10);
//
//        // Test second grid point
//        TS_ASSERT_DELTA(force_grid_x[3][3], 1.0 * large_delta * small_delta, 1e-10);
//        TS_ASSERT_DELTA(force_grid_y[3][3], 2.0 * large_delta * small_delta, 1e-10);
//
//        // Test third grid point
//        TS_ASSERT_DELTA(force_grid_x[2][4], 1.0 * small_delta * small_delta, 1e-10);
//        TS_ASSERT_DELTA(force_grid_y[2][4], 2.0 * small_delta * small_delta, 1e-10);
    }

    void xTestUpwindSchemeImplementation()
    {
//        MAKE_PTR(ImmersedBoundarySimulationModifier < 2 > , p_mod);
//
//        /*
//         *  Two small matrices of random numbers represent the fluid velocity field.  These are read in,
//         *  along with two matrices representing the upwind difference of this field, calculated by hand.
//         *  The by-hand and c++ implementations of the upwind difference scheme are compared.
//         */
//        unsigned num_gridpts_y = 2;
//        unsigned num_gridpts_x = 3;
//
//        // Use private method to set necessary member variables so we don't have to make a mesh and cell population
//        p_mod->SetMemberVariablesForTesting(num_gridpts_y, num_gridpts_x);
//
//        // Set up space for the two matrices representing the velocity fields
//        std::vector<std::vector<double> > vel_x; p_mod->SetupGrid(vel_x);
//        std::vector<std::vector<double> > vel_y; p_mod->SetupGrid(vel_y);
//
//        // Set up space for the two matrices representing the upwind difference of these fields, calculated by hand
//        std::vector<std::vector<double> > hand_upwind_x; p_mod->SetupGrid(hand_upwind_x);
//        std::vector<std::vector<double> > hand_upwind_y; p_mod->SetupGrid(hand_upwind_y);
//
//        // Read in vel_x
//        ifstream f_vel_x;
//        f_vel_x.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/vel_x.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_vel_x >> vel_x[y][x];
//            }
//        }
//        f_vel_x.close();
//
//        // Read in vel_y
//        ifstream f_vel_y;
//        f_vel_y.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/vel_y.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_vel_y >> vel_y[y][x];
//            }
//        }
//        f_vel_y.close();
//
//        // Read in hand_upwind_x
//        ifstream f_hand_upwind_x;
//        f_hand_upwind_x.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/upwind_x.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_hand_upwind_x >> hand_upwind_x[y][x];
//            }
//        }
//        f_hand_upwind_x.close();
//
//        // Read in hand_upwind_y
//        ifstream f_hand_upwind_y;
//        f_hand_upwind_y.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/upwind_y.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_hand_upwind_y >> hand_upwind_y[y][x];
//            }
//        }
//        f_hand_upwind_y.close();
//
//        /*
//         * Now, we calculate the upwind difference using the C++ implementation, and compare it to
//         * the matrices which have been calculated by hand
//         */
//        std::vector<std::vector<double> > cpp_upwind_x;
//        std::vector<std::vector<double> > cpp_upwind_y;
//
//        // Calculate the upwind difference using code in ImmersedBoundarySimulationModifier::UpwindScheme()
//        p_mod->UpwindScheme(vel_x, vel_y, cpp_upwind_x, cpp_upwind_y);
//
//        // Loop through the hand and c++ calculated grids to check the values agree
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                TS_ASSERT_LESS_THAN(fabs(hand_upwind_x[y][x] - cpp_upwind_x[y][x]), 1e-10);
//                TS_ASSERT_LESS_THAN(fabs(hand_upwind_y[y][x] - cpp_upwind_y[y][x]), 1e-10);
//            }
//        }
    }

    void xTestFourierTransformMethods()
    {
//        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_mod);
//
//        /*
//         * Testing for fft uses 32(y) * 16(x) matrices randomly generated using randn() in Matlab,
//         * as well as the DFT and inverse DFT using fft2() and ifft2() in Matlab.  Relevant matrices
//         * are exported from Matlab and read in during this test to compare against the fftw implementation.
//         */
//        unsigned num_gridpts_y = 32;
//        unsigned num_gridpts_x = 16;
//
//        // Use private method to set necessary member variables so we don't have to make a mesh and cell population
//        p_mod->SetMemberVariablesForTesting(num_gridpts_y, num_gridpts_x);
//
//        // Set up space for the two matrices randomly generated in Matlab
//        std::vector<std::vector<double> > randn_mat_a; p_mod->SetupGrid(randn_mat_a);
//        std::vector<std::vector<double> > randn_mat_b; p_mod->SetupGrid(randn_mat_b);
//
//        // Set up space for the matrix obtained by fft2(randn_mat_a) in Matlab
//        std::vector<std::vector<std::complex<double> > > matlab_fft2_mat_a; p_mod->SetupGrid(matlab_fft2_mat_a);
//
//        // Set up space for the matrix obtained by real(ifft2(randn_mat_a + i * randn_mat_b)) in Matlab
//        std::vector<std::vector<double> > matlab_ifft2_aplusib; p_mod->SetupGrid(matlab_ifft2_aplusib);
//
//        // Read in matrix_a
//        ifstream f_randn_mat_a;
//        f_randn_mat_a.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/matrix_a.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_randn_mat_a >> randn_mat_a[y][x];
//            }
//        }
//        f_randn_mat_a.close();
//
//        // Read in matrix_b
//        ifstream f_randn_mat_b;
//        f_randn_mat_b.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/matrix_b.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_randn_mat_b >> randn_mat_b[y][x];
//            }
//        }
//        f_randn_mat_b.close();
//
//        // Read in the fft of matrix_a
//        ifstream f_real, f_imag;
//        f_real.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/ffta_real.dat");
//        f_imag.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/ffta_imag.dat");
//        double real_temp;
//        double imag_temp;
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_real >> real_temp;
//                f_imag >> imag_temp;
//                matlab_fft2_mat_a[y][x] = real_temp + (p_mod->mI) * imag_temp;
//            }
//        }
//        f_real.close();
//        f_imag.close();
//
//        // Read in real part of the inverse transform of mat_a + i * mat_b;
//        ifstream f_ifft;
//        f_ifft.open("cell_based/test/data/TestImmersedBoundaryPdeSolveMethods/ifftaplusib_real.dat");
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                f_ifft >> matlab_ifft2_aplusib[y][x];
//            }
//        }
//        f_ifft.close();
//
//        // Create matrix_a + i * matrix_b
//        std::vector<std::vector<std::complex<double> > > mat_aplusib; p_mod->SetupGrid(mat_aplusib);
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                mat_aplusib[y][x] = randn_mat_a[y][x] + (p_mod->mI) * randn_mat_b[y][x];
//            }
//        }
//
//        /*
//         * Next, we set up space for two matrices which will store results of the C++ implementation of the DFT
//         * (using FFTW), which will be tested against the Matlab DFT results which were read in above
//         */
//        std::vector<std::vector<std::complex<double> > > cpp_fftw_mat_a;
//        std::vector<std::vector<double> > cpp_ifftw_aplusib;
//
//        // Perform FFT using FFTW code in ImmersedBoundarySimulationModifier::Fft2DForwardRealToComplex()
//        p_mod->Fft2DForwardRealToComplex(randn_mat_a, cpp_fftw_mat_a);
//
//        // Perform FFT using FFTW code in ImmersedBoundarySimulationModifier::Fft2DInverseComplexToReal()
//        p_mod->Fft2DInverseComplexToReal(mat_aplusib, cpp_ifftw_aplusib);
//
//        /*
//         * Loop through the arrays to check values agree
//         */
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                TS_ASSERT_LESS_THAN(abs(matlab_fft2_mat_a[y][x] - cpp_fftw_mat_a[y][x]), 1e-10);
//                TS_ASSERT_LESS_THAN(fabs(matlab_ifft2_aplusib[y][x] - cpp_ifftw_aplusib[y][x]), 1e-10);
//            }
//        }
//
//        /*
//         * Next, we check that an FFT followed by inverse FFT gives back the original matrix
//         */
//        std::vector<std::vector<std::complex<double> > > cpp_fftw_mat_b;
//        std::vector<std::vector<double> > cpp_ifftw_fftw_mat_b;
//
//        p_mod->Fft2DForwardRealToComplex(randn_mat_b, cpp_fftw_mat_b);
//        p_mod->Fft2DInverseComplexToReal(cpp_fftw_mat_b, cpp_ifftw_fftw_mat_b);
//        for (unsigned y = 0; y < num_gridpts_y; y++)
//        {
//            for (unsigned x = 0; x < num_gridpts_x; x++)
//            {
//                TS_ASSERT_LESS_THAN(fabs(randn_mat_b[y][x] - cpp_ifftw_fftw_mat_b[y][x]), 1e-10);
//            }
//        }
    }

    void TestBoostMultiarray()
    {
//        std::string filename = "./projects/ImmersedBoundary/src/fftw.wisdom";
//        int wisdom_flag = fftw_import_wisdom_from_filename(filename.c_str());
//
//        // 1 means it's read correctly, 0 indicates a failure
//        TS_ASSERT_EQUALS(wisdom_flag, 1);
//
//        // Create a 3D array that is 64 x 64 x 64
//        typedef boost::multi_array<std::complex<double>, 2> complex_array;
//        typedef boost::multi_array<double, 2> real_array;
//
//        complex_array complex_input(boost::extents[4096][4096]);
//        complex_array complex_output_1(boost::extents[4096][4096]);
//        complex_array complex_output_2(boost::extents[4096][4096]);
//        real_array real_input(boost::extents[4096][4096]);
//
//
//        std::complex<double>* c_it = complex_input.origin();
//        double*               r_it = real_input.origin();
//        for (; c_it < (complex_input.origin() + complex_input.num_elements()); ++c_it, ++r_it)
//        {
//            *c_it = RandomNumberGenerator::Instance()->ranf() + 1i * RandomNumberGenerator::Instance()->ranf();
//            *r_it = RandomNumberGenerator::Instance()->ranf();
//        }
//
//        Timer timer;
//        timer.Reset();
//
//        c_it = complex_input.origin();
//        r_it = real_input.origin();
//        std::complex<double>* out_it = complex_output_1.origin();
//        for (; out_it < (complex_output_1.origin() + complex_output_1.num_elements()); ++c_it, ++r_it, ++out_it)
//        {
//            *out_it = (*r_it) * (*c_it);
//        }
//
//        double it_time = timer.GetElapsedTime();
//
//        for (unsigned i=0; i < 4096; i++)
//        {
//            for (unsigned j=0; j < 4096; j++)
//            {
//                complex_output_2[i][j] = real_input[i][j] * complex_input[i][j];
//            }
//        }
//
//        double loop_time = timer.GetElapsedTime() - it_time;
//
//
//        for (unsigned i=0; i < 4096; i++)
//        {
//            for (unsigned j=0; j < 4096; j++)
//            {
//                TS_ASSERT_DELTA(complex_output_1[i][j].real(), complex_output_2[i][j].real(),  1e-10);
//                TS_ASSERT_DELTA(complex_output_1[i][j].imag(), complex_output_2[i][j].imag(),  1e-10);
//            }
//        }
//
//        PRINT_2_VARIABLES(it_time, loop_time);
    }

    void xTestFluidSolve()
    {
//        // Create a vector of nodes forming a rectangle in (0,1)x(0,1)
//        double irrational = 0.01 * sqrt(2); // ensure our point of interest isn't on a gird point
//        std::vector<Node<2>*> nodes;
//        nodes.push_back(new Node<2>(0, true, 0.4, 0.5));
//        nodes.push_back(new Node<2>(1, true, 0.5 + irrational, 0.5 + irrational));
//        nodes.push_back(new Node<2>(2, true, 0.5, 0.6));
//        nodes.push_back(new Node<2>(3, true, 0.4, 0.6));
//
//        // Create a vector of immersed boundary elements and create an element with the nodes above
//        std::vector<ImmersedBoundaryElement<2,2>*> ib_element;
//        ib_element.push_back(new ImmersedBoundaryElement<2,2>(0, nodes));
//
//        // Create a mesh with the nodes and elements vectors
//        ImmersedBoundaryMesh<2,2>* p_mesh = new ImmersedBoundaryMesh<2,2>(nodes, ib_element, 10, 10);
//
//        // Set element parameters
//        ImmersedBoundaryElement<2,2>* p_elem = p_mesh->GetElement(0u);
//        p_elem->SetMembraneSpringConstant(1.0);
//        p_elem->SetMembraneRestLength(0.4);
//
//        p_mesh->SetNumGridPtsX(16);
//        p_mesh->SetNumGridPtsY(16);
//
//        std::vector<CellPtr> cells;
//        MAKE_PTR(DifferentiatedCellProliferativeType, p_diff_type);
//        CellsGenerator<UniformCellCycleModel, 2> cells_generator;
//        cells_generator.GenerateBasicRandom(cells, p_mesh->GetNumElements(), p_diff_type);
//
//        ImmersedBoundaryCellPopulation<2> cell_population(*p_mesh, cells);
//
//        OffLatticeSimulation<2> simulator(cell_population);
//
//        // Add main immersed boundary simulation modifier
//        MAKE_PTR(ImmersedBoundarySimulationModifier<2>, p_mod);
//        simulator.AddSimulationModifier(p_mod);
//        simulator.SetDt(0.01);
//
//        // Calculate the perimeter elastic forces and ensure they match with by-hand calculations
//        p_mod->SetupConstantMemberVariables(cell_population);
//        p_mod->ClearForces();
//
//        c_vector<double, 2> force_on_node;
//        force_on_node[0] = 1.0;
//        force_on_node[1] = 0.0;
//
//        nodes[1]->AddAppliedForceContribution(force_on_node);
//
//        p_mod->PropagateForcesToFluidGrid();
//        std::vector<std::vector<double> > force_grid_x = p_mod->mFluidForceGridX;
//        std::vector<std::vector<double> > force_grid_y = p_mod->mFluidForceGridY;
//
//        // Check that there is no force propagated in the y-direction
//        for (unsigned y = 0; y < p_mod->mNumGridPtsY; y++)
//        {
//            for (unsigned x = 0; x < p_mod->mNumGridPtsX; x++)
//            {
//                TS_ASSERT_LESS_THAN(fabs(force_grid_y[y][x]), 1e-15);
//            }
//        }
//
//        SimulationTime::Instance()->SetEndTimeAndNumberOfTimeSteps(0.01,1);
//
//        p_mod->SolveNavierStokesSpectral();
//
////        const std::vector<std::vector<double> > vel_grid_x = p_mesh->rGetFluidVelocityGridX();
////        const std::vector<std::vector<double> > vel_grid_y = p_mesh->rGetFluidVelocityGridY();
//
////        // All velocity should be in the positive x direction
////        for (unsigned y = 0; y < p_mod->mNumGridPtsY; y++)
////        {
////            for (unsigned x = 0; x < p_mod->mNumGridPtsX; x++)
////            {
////                TS_ASSERT(vel_grid_x[y][x] >= 0.0);
////            }
////        }
//
//        p_mod->ClearForces();
//        for (unsigned y = 0; y < p_mod->mNumGridPtsY; y++)
//        {
//            for (unsigned x = 0; x < p_mod->mNumGridPtsX; x++)
//            {
//                p_mod->mFluidForceGridX[y][x] = 1.0;
//            }
//        }
    }
};

#endif /*TESTIMMERSEDBOUNDARYPDESOLVEMETHODS_HPP_*/
