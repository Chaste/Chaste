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


#ifndef _TESTMOREMECHANICS_HPP_
#define _TESTMOREMECHANICS_HPP_

#include <cxxtest/TestSuite.h>
#include <boost/shared_ptr.hpp>
#include <boost/assign.hpp>

#include "UblasCustomFunctions.hpp"
#include "IncompressibleNonlinearElasticitySolver.hpp"
#include "CompressibleNonlinearElasticitySolver.hpp"
#include "MooneyRivlinMaterialLaw.hpp"
#include "CompressibleMooneyRivlinMaterialLaw.hpp"
#include "ExponentialMaterialLaw.hpp"
#include "CompressibleExponentialLaw.hpp"
#include "NonlinearElasticityTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "VtkNonlinearElasticitySolutionWriter.hpp"


// This test tests some some useful but developmental mechanics code
//  1. SetUsePetscDirectSolve() which can be used to do a direct linear solve. Seems very useful (big gains in nonlinear solve convergence)
//     if there is enough memory. Petsc LU factorisation can only be in sequential though. Hopefully be replaced by MUMPs in the future.
//  2. SetTakeFullFirstNewtonStep() [The following odd behaviour has been observed: for some problems the solver
//     will fail in the first Newton iteration, with the residual not decreasing in the direction of the Newton
//     update, but if you take a full Newton step anyway (increasing the residual-norm), the solver then converges
//     perfectly. This method allows the user to choose this option
//  3. InterpolateMechanicsSolutionToNewMesh() is a useful method (defined in this test) for using the solution on a coarse
//     mesh as a guess solution on a fine mesh, may go to src later..

/**
 *  Interpolate the mechanics solution vector from a coarse mesh mesh to a fine mesh. Idea is that you
 *  solve first on a coarse mesh, interpolate the solution onto a fine mesh, then use this interpolated
 *  solution as the initial guess for the fine mesh problem.
 *
 *  The method takes a mechanics solution, which will be of the form
 *    [ U1 V1 W1 U2 V2 W2 .. Un Vn Wn ]
 *  for compressible problems, or
 *    [ U1 V1 W1 p1 U2 V2 W2 p2 .. Un Vn Wn pn]
 *  incompressible problems, and interpolates onto a finer mesh, returning an initial guess of the
 *  same structure.
 */
template<unsigned DIM>
void InterpolateMechanicsSolutionToNewMesh(QuadraticMesh<3>& rCoarseMesh, std::vector<double>& rCoarseSolution,
                                           QuadraticMesh<3>& rFineMesh, std::vector<double>& rFineSolution,
                                           CompressibilityType compressibilityType)
{
    unsigned NUM_UNKNOWNS = (compressibilityType==INCOMPRESSIBLE ? DIM+1 : DIM);

    if (rCoarseSolution.size() != rCoarseMesh.GetNumNodes()*NUM_UNKNOWNS)
    {
        EXCEPTION("rCoarseSolution not correct size");
    }
    if (rFineSolution.size() != rFineMesh.GetNumNodes()*NUM_UNKNOWNS)
    {
        rFineSolution.resize(rFineMesh.GetNumNodes()*NUM_UNKNOWNS);
    }

    c_vector<double, (DIM+1)*(DIM+2)/2> quad_basis;

    for (unsigned i=0; i<rFineMesh.GetNumNodes(); i++)
    {
        // find containing elements and weights in coarse mesh
        ChastePoint<DIM> point = rFineMesh.GetNode(i)->GetPoint();
        unsigned coarse_element_index = rCoarseMesh.GetContainingElementIndex(point,
                                                                              false);
        Element<DIM,DIM>* p_coarse_element = rCoarseMesh.GetElement(coarse_element_index);
        c_vector<double,DIM+1> weight = p_coarse_element->CalculateInterpolationWeights(point);

        c_vector<double,DIM> xi;
        xi(0) = weight(1);
        xi(1) = weight(2);
        if (DIM==3)
        {
            xi(2) = weight(3);
        }

        QuadraticBasisFunction<DIM>::ComputeBasisFunctions(xi, quad_basis);

        // interpolate (u,p) (don't do anything for p if compressible)
        c_vector<double,DIM+1> fine_solution = zero_vector<double>(DIM+1);

        for (unsigned elem_node_index=0; elem_node_index<(DIM+1)*(DIM+2)/2; elem_node_index++)
        {
            unsigned coarse_node = p_coarse_element->GetNodeGlobalIndex(elem_node_index);
            c_vector<double,DIM+1> coarse_solution_at_node;

            for (unsigned j=0; j<DIM; j++)
            {
                coarse_solution_at_node(j) = rCoarseSolution[NUM_UNKNOWNS*coarse_node + j];
            }
            if (compressibilityType==INCOMPRESSIBLE)
            {
                coarse_solution_at_node(DIM) = rCoarseSolution[NUM_UNKNOWNS*coarse_node + DIM];
            }

            fine_solution += coarse_solution_at_node*quad_basis(elem_node_index);
        }

        for (unsigned j=0; j<DIM; j++)
        {
            rFineSolution[NUM_UNKNOWNS*i + j] = fine_solution(j);
        }

        if (compressibilityType==INCOMPRESSIBLE)
        {
            rFineSolution[NUM_UNKNOWNS*i + DIM] = fine_solution(DIM);

            // Whilst the returned p from a solve is defined properly at all nodes, during the solve linear basis functions are
            // used for p and therefore p not computed explicitly at internal nodes, and the solver solves for p=0 at these internal
            // nodes. (After the solve, p is interpolated from vertices to internal nodes)
            if (rFineMesh.GetNode(i)->IsInternal())
            {
                rFineSolution[NUM_UNKNOWNS*i + DIM] = 0.0;
            }
        }
    }
}


// Solve a problem where some pressure is applied in the outward direction on the bottom (Z=0) surface
//
// Returns the number of newton iterations taken
unsigned SolvePressureOnUnderside(QuadraticMesh<3>& rMesh, std::string outputDirectory,
                                  std::vector<double>& rSolution,
                                  bool useSolutionAsGuess,
                                  double scaleFactor = 1.0)
{
    MooneyRivlinMaterialLaw<3> law(1.0, 0.0);

    std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(rMesh, 0, 0.0);

    std::vector<BoundaryElement<2,3>*> boundary_elems;

    double pressure = 0.001;

    for (TetrahedralMesh<3,3>::BoundaryElementIterator iter
           = rMesh.GetBoundaryElementIteratorBegin();
         iter != rMesh.GetBoundaryElementIteratorEnd();
         ++iter)
    {
        BoundaryElement<2,3>* p_element = *iter;
        double Z = p_element->CalculateCentroid()[2];
        if (fabs(Z)<1e-6)
        {
            boundary_elems.push_back(p_element);
        }
    }

    SolidMechanicsProblemDefinition<3> problem_defn(rMesh);
    problem_defn.SetMaterialLaw(INCOMPRESSIBLE,&law);
    problem_defn.SetZeroDisplacementNodes(fixed_nodes);

    problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure*scaleFactor);
    problem_defn.SetVerboseDuringSolve();

    IncompressibleNonlinearElasticitySolver<3> solver(rMesh,problem_defn,outputDirectory);
    solver.SetWriteOutputEachNewtonIteration();

    // cover the SetTakeFullFirstNewtonStep() method, and the SetUsePetscDirectSolve() method
    solver.SetTakeFullFirstNewtonStep();
    solver.SetUsePetscDirectSolve();

    if (useSolutionAsGuess)
    {
        if (solver.rGetCurrentSolution().size()!=rSolution.size())
        {
            EXCEPTION("Badly-sized input");
        }
        for (unsigned i=0; i<rSolution.size(); i++)
        {
            solver.rGetCurrentSolution()[i] = rSolution[i];
        }
    }

    if (scaleFactor < 1.0)
    {
        try
        {
            solver.Solve();
        }
        catch (Exception& e)
        {
            // not final Solve, so don't quit
            WARNING(e.GetMessage());
        }
    }
    else
    {
        solver.Solve();
        solver.CreateCmguiOutput();

        VtkNonlinearElasticitySolutionWriter<3> vtk_writer(solver);
        vtk_writer.SetWriteElementWiseStrains(DEFORMATION_TENSOR_C);
        vtk_writer.Write();
    }

    rSolution.clear();
    rSolution.resize(solver.rGetCurrentSolution().size());
    for (unsigned i=0; i<rSolution.size(); i++)
    {
        rSolution[i] = solver.rGetCurrentSolution()[i];
    }

    return solver.GetNumNewtonIterations();
}


// Similar to above but uses a compressible material
unsigned SolvePressureOnUndersideCompressible(QuadraticMesh<3>& rMesh, std::string outputDirectory,
                                              std::vector<double>& rSolution,
                                              bool useSolutionAsGuess,
                                              double scaleFactor = 1.0)
{
    CompressibleExponentialLaw<3> law;
    std::vector<unsigned> fixed_nodes = NonlinearElasticityTools<3>::GetNodesByComponentValue(rMesh, 0, 0.0);

    std::vector<BoundaryElement<2,3>*> boundary_elems;

    double pressure = 0.001;

    for (TetrahedralMesh<3,3>::BoundaryElementIterator iter
           = rMesh.GetBoundaryElementIteratorBegin();
         iter != rMesh.GetBoundaryElementIteratorEnd();
         ++iter)
    {
        BoundaryElement<2,3>* p_element = *iter;
        double Z = p_element->CalculateCentroid()[2];
        if (fabs(Z)<1e-6)
        {
            boundary_elems.push_back(p_element);
        }
    }

    SolidMechanicsProblemDefinition<3> problem_defn(rMesh);
    problem_defn.SetMaterialLaw(COMPRESSIBLE,&law);
    problem_defn.SetZeroDisplacementNodes(fixed_nodes);

    problem_defn.SetApplyNormalPressureOnDeformedSurface(boundary_elems, pressure*scaleFactor);


    problem_defn.SetVerboseDuringSolve();

    CompressibleNonlinearElasticitySolver<3> solver(rMesh,problem_defn,outputDirectory);
    solver.SetWriteOutputEachNewtonIteration();

    // cover the SetTakeFullFirstNewtonStep() method, and the SetUsePetscDirectSolve() method
    solver.SetTakeFullFirstNewtonStep();
    solver.SetUsePetscDirectSolve();

    if (useSolutionAsGuess)
    {
        if (solver.rGetCurrentSolution().size()!=rSolution.size())
        {
            EXCEPTION("Badly-sized input");
        }
        for (unsigned i=0; i<rSolution.size(); i++)
        {
            solver.rGetCurrentSolution()[i] = rSolution[i];
        }
    }

    if (scaleFactor < 1.0)
    {
        try
        {
            solver.Solve();
        }
        catch (Exception& e)
        {
            // not final Solve, so don't quit
            WARNING(e.GetMessage());
        }
    }
    else
    {
        solver.Solve();
        solver.CreateCmguiOutput();

        VtkNonlinearElasticitySolutionWriter<3> vtk_writer(solver);
        vtk_writer.SetWriteElementWiseStrains(DEFORMATION_TENSOR_C);
        vtk_writer.Write();
    }

    rSolution.clear();
    rSolution.resize(solver.rGetCurrentSolution().size());
    for (unsigned i=0; i<rSolution.size(); i++)
    {
        rSolution[i] = solver.rGetCurrentSolution()[i];
    }

    return solver.GetNumNewtonIterations();
}

class TestMoreMechanics : public CxxTest::TestSuite
{
public:
    void TestBarPressureOnUnderside()
    {
        EXIT_IF_PARALLEL; // petsc's direct solve only runs one 1 proc - MUMPS may be the answer for direct solves in parallel

        std::string base = "BarPressureOnUnderside";
        std::stringstream output_dir;


        /////////////////////////////////////
        // Solve on coarse mesh..
        /////////////////////////////////////
        QuadraticMesh<3> mesh;
        double h = 1.0;
        mesh.ConstructRegularSlabMesh(h, 10.0, 1.0, 1.0);
        output_dir << base << "h_" << h;
        std::vector<double> solution;
        unsigned num_iters = SolvePressureOnUnderside(mesh, output_dir.str(), solution, false, true);
        TS_ASSERT_EQUALS(num_iters, 6u);

        /////////////////////////////////////
        // Solve on slightly finer mesh..
        /////////////////////////////////////
        h = 0.5;
        QuadraticMesh<3> mesh2;
        mesh2.ConstructRegularSlabMesh(h, 10.0, 1.0, 1.0);
        std::vector<double> solution2;

        InterpolateMechanicsSolutionToNewMesh<3>(mesh, solution, mesh2, solution2, INCOMPRESSIBLE);

        output_dir.str("");
        output_dir << base << "h_" << h;
        num_iters = SolvePressureOnUnderside(mesh2, output_dir.str(), solution2, true, true);
        TS_ASSERT_EQUALS(num_iters, 2u);

        ///////////////////////////////////////
        // Could continue onto finer meshes..
        ///////////////////////////////////////
    }

    void TestBarPressureOnUndersideCompressible()
    {
        EXIT_IF_PARALLEL; // petsc's direct solve only runs one 1 proc - MUMPS may be the answer for direct solves in parallel

        std::string base = "BarPressureOnUndersideCompressible";
        std::stringstream output_dir;

        /////////////////////////////////////
        // Solve on coarse mesh..
        /////////////////////////////////////
        QuadraticMesh<3> mesh;
        double h = 1.0;
        mesh.ConstructRegularSlabMesh(h, 10.0, 1.0, 1.0);
        output_dir << base << "h_" << h;
        std::vector<double> solution;
        unsigned num_iters = SolvePressureOnUndersideCompressible(mesh, output_dir.str(), solution, false, true);
        TS_ASSERT_EQUALS(num_iters, 5u);
    }
};

#endif //_TESTMOREMECHANICS_HPP_
