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

#ifndef TESTSOLVELAPLACIANWITHQUADRATICS_HPP_
#define TESTSOLVELAPLACIANWITHQUADRATICS_HPP_

#include <cxxtest/TestSuite.h>
#include "TetrahedralMesh.hpp"
#include <petsc.h>
#include <vector>
#include <cmath>
#include "BoundaryConditionsContainer.hpp"
#include "TrianglesMeshReader.hpp"
#include "PetscTools.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "QuadraticBasisFunction.hpp"
#include "QuadraticMesh.hpp"
#include "LinearSystem.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "GaussianQuadratureRule.hpp"
#include "EllipticPdeWithLinearSource.hpp"
#include "SimpleLinearEllipticSolver.hpp"

/*
 * Assemblers which use quadratic bases (ie the elasticity solvers) are outside the
 * other assembler hierarchy. Therefore, to test the quadratic bases/structure we
 * first wrote this QuadraticLaplacianAssemblerSolver, and then used it when writing
 * the elasticity assemblers.
 */
template<unsigned DIM>
class QuadraticLaplacianAssemblerSolver
{
private:
    double mCoeffOfU;
    double mConstant;

    LinearSystem* mpLinearSystem;
    QuadraticMesh<DIM>* mpQuadMesh;
    BoundaryConditionsContainer<DIM,DIM,1>* mpBoundaryConditions;
    GaussianQuadratureRule<DIM>* mpQuadRule;

    static const unsigned NUM_BASES_PER_ELEMENT = (DIM+1)*(DIM+2)/2;
    static const unsigned STENCIL_SIZE = NUM_BASES_PER_ELEMENT; // multiplied by PROBLEM_DIM

    virtual void AssembleOnElement( Element<DIM, DIM>& rElement,
                                    c_matrix<double, STENCIL_SIZE, STENCIL_SIZE >& rAElem,
                                    c_vector<double, STENCIL_SIZE>& rBElem,
                                    bool assembleVector,
                                    bool assembleMatrix)
    {
        /**
         * \todo This assumes that the Jacobian is constant on an element.
         * This is true for linear basis functions, but not for any other type of
         * basis function. <--- not true, as we don't have curvilinear elements!
         */
        c_matrix<double, DIM, DIM> jacobian, inverse_jacobian;
        double jacobian_determinant;
        mpQuadMesh->GetInverseJacobianForElement(rElement.GetIndex(), jacobian, jacobian_determinant, inverse_jacobian);

        // Initialise element contributions to zero
        if (assembleMatrix)
        {
            rAElem.clear();
        }

        if (assembleVector)
        {
            rBElem.clear();
        }

        // Allocate memory for the basis functions values and derivative values
        c_vector<double, NUM_BASES_PER_ELEMENT> phi;
        c_matrix<double, DIM, NUM_BASES_PER_ELEMENT> grad_phi;

        // Loop over Gauss points
        for (unsigned quad_index=0; quad_index < mpQuadRule->GetNumQuadPoints(); quad_index++)
        {
            const ChastePoint<DIM>& quad_point = mpQuadRule->rGetQuadPoint(quad_index);

            QuadraticBasisFunction<DIM>::ComputeBasisFunctions(quad_point, phi);

            if (assembleMatrix)
            {
                QuadraticBasisFunction<DIM>::ComputeTransformedBasisFunctionDerivatives(quad_point, inverse_jacobian, grad_phi);
            }

            // Location of the Gauss point in the original element will be stored in x
            // Where applicable, u will be set to the value of the current solution at x
            ChastePoint<DIM> x(0,0,0);

            c_vector<double,1> u = zero_vector<double>(1);
            c_matrix<double,1,DIM> grad_u = zero_matrix<double>(1,DIM);
            double wJ = jacobian_determinant * mpQuadRule->GetWeight(quad_index);

            // Create rAElem and rBElem
            if (assembleMatrix)
            {
                noalias(rAElem) += ComputeMatrixTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
            }
            if (assembleVector)
            {
                noalias(rBElem) += ComputeVectorTerm(phi, grad_phi, x, u, grad_u, &rElement) * wJ;
            }
        }
    }

    virtual c_matrix<double,STENCIL_SIZE,STENCIL_SIZE> ComputeMatrixTerm(
        c_vector<double, NUM_BASES_PER_ELEMENT>& rPhi,
        c_matrix<double, DIM, NUM_BASES_PER_ELEMENT>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
            return   prod( trans(rGradPhi), rGradPhi )
                   - mCoeffOfU*outer_prod(rPhi, rPhi);
    }

    virtual c_vector<double,STENCIL_SIZE> ComputeVectorTerm(
        c_vector<double, NUM_BASES_PER_ELEMENT>& rPhi,
        c_matrix<double, DIM, NUM_BASES_PER_ELEMENT>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        return mConstant * rPhi;
    }

    virtual void AssembleSystem(bool assembleVector, bool assembleMatrix
                                /*Vec currentSolutionOrGuess=NULL, double currentTime=0.0*/)
    {
        // Check we've actually been asked to do something!
        assert(assembleVector || assembleMatrix);

        // Zero the matrix/vector if it is to be assembled
        if (assembleVector)
        {
            mpLinearSystem->ZeroRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->ZeroLhsMatrix();
        }

        c_matrix<double, STENCIL_SIZE, STENCIL_SIZE> a_elem;
        c_vector<double, STENCIL_SIZE> b_elem;

        // Loop over elements
        for (typename AbstractTetrahedralMesh<DIM, DIM>::ElementIterator iter = mpQuadMesh->GetElementIteratorBegin();
             iter != mpQuadMesh->GetElementIteratorEnd();
             ++iter)
        {
            Element<DIM,DIM>& element = *iter;

            if (element.GetOwnership() == true)
            {
                AssembleOnElement(element, a_elem, b_elem, assembleVector, assembleMatrix);

                //unsigned p_indices[STENCIL_SIZE];
                //element.GetStiffnessMatrixGlobalIndices(PROBLEM_DIM, p_indices);

                unsigned p_indices[STENCIL_SIZE];
                for (unsigned i=0; i<STENCIL_SIZE; i++)
                {
                    p_indices[i] = element.GetNodeGlobalIndex(i);
                }


                if (assembleMatrix)
                {
                    mpLinearSystem->AddLhsMultipleValues(p_indices, a_elem);
                }

                if (assembleVector)
                {
                    mpLinearSystem->AddRhsMultipleValues(p_indices, b_elem);
                }
            }
        }

        if (assembleVector)
        {
            mpLinearSystem->FinaliseRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->SwitchWriteModeLhsMatrix();
        }

        // Apply dirichlet boundary conditions
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*mpLinearSystem, assembleMatrix);

        if (assembleVector)
        {
            mpLinearSystem->FinaliseRhsVector();
        }
        if (assembleMatrix)
        {
            mpLinearSystem->FinaliseLhsMatrix();
        }
    }

public:
    QuadraticLaplacianAssemblerSolver(QuadraticMesh<DIM>* pMesh,
                                      BoundaryConditionsContainer<DIM,DIM,1>* pBcc)
        : mpQuadMesh(pMesh),
          mpBoundaryConditions(pBcc)
    {
        assert(pMesh);
        assert(pBcc);

        unsigned preallocation=mpQuadMesh->CalculateMaximumNodeConnectivityPerProcess(); //The midpoint, as given above

        mpLinearSystem = new LinearSystem(mpQuadMesh->GetNumNodes(), preallocation);

        mpQuadRule = new GaussianQuadratureRule<DIM>(2);

        mCoeffOfU = 0.0;
        mConstant = 1.0;
    }

    virtual ~QuadraticLaplacianAssemblerSolver()
    {
        delete mpLinearSystem;
        delete mpQuadRule;
    }

    Vec Solve()
    {
        AssembleSystem(true, true);
        return mpLinearSystem->Solve();
    }

    void SetPdeConstants(double coeffOfU, double constant)
    {
        mCoeffOfU = coeffOfU;
        mConstant = constant;
    }
};

class TestSolveLaplacianWithQuadratics : public CxxTest::TestSuite
{
public:
    // solve u'' + 1 = 0
    void TestSolveLaplacianWithQuadratics1d()
    {
        QuadraticMesh<1> quad_mesh;
        TrianglesMeshReader<1,1> mesh_reader("mesh/test/data/1D_0_to_1_10_elements_quadratic",2,1,false);
        quad_mesh.ConstructFromMeshReader(mesh_reader);

        BoundaryConditionsContainer<1,1,1> bcc;
        bcc.DefineZeroDirichletOnMeshBoundary(&quad_mesh);

        QuadraticLaplacianAssemblerSolver<1> solver(&quad_mesh, &bcc);
        solver.SetPdeConstants(0.0, 1.0);

        Vec solution = solver.Solve();
        ReplicatableVector sol_repl(solution);

        for (unsigned i=0; i<quad_mesh.GetNumNodes(); i++)
        {
            double x = quad_mesh.GetNode(i)->rGetLocation()[0];
            double u = sol_repl[i];
            double u_correct = 0.5*x*(1-x);

            TS_ASSERT_DELTA(u, u_correct, 1e-10);
        }

        PetscTools::Destroy(solution);
    }

    void TestSolveLaplacianWithQuadratics2d()
    {
        // Solve using quadratics..
        QuadraticMesh<2> quad_mesh;
        TrianglesMeshReader<2,2> mesh_reader1("mesh/test/data/square_128_elements_quadratic",2,1,false);
        quad_mesh.ConstructFromMeshReader(mesh_reader1);

        BoundaryConditionsContainer<2,2,1> bcc_quads;
        bcc_quads.DefineZeroDirichletOnMeshBoundary(&quad_mesh);

        QuadraticLaplacianAssemblerSolver<2> solver_quads(&quad_mesh, &bcc_quads);
        solver_quads.SetPdeConstants(1.0, 1.0);

        Vec solution_quads = solver_quads.Solve();
        ReplicatableVector sol_quads_repl(solution_quads);

        // Solve using linears
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        EllipticPdeWithLinearSource<2> pde(1.0, 1.0);

        BoundaryConditionsContainer<2,2,1> bcc_lin;
        bcc_lin.DefineZeroDirichletOnMeshBoundary(&mesh);

        SimpleLinearEllipticSolver<2,2> solver_lin(&mesh,&pde,&bcc_lin);

        Vec solution_lin = solver_lin.Solve();
        ReplicatableVector sol_lin_repl(solution_lin);

        /*
         * Compare results - the following assumes the vertex nodes in the
         * quad mesh are nodes 0-63, i.e. they come before all the internal
         * nodes.
         */
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), quad_mesh.GetNumVertices());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double lin_x = mesh.GetNode(i)->rGetLocation()[0];
            double lin_y = mesh.GetNode(i)->rGetLocation()[1];
            double quad_x = quad_mesh.GetNode(i)->rGetLocation()[0];
            double quad_y = quad_mesh.GetNode(i)->rGetLocation()[1];
            TS_ASSERT_DELTA(lin_x, quad_x, 1e-8);
            TS_ASSERT_DELTA(lin_y, quad_y, 1e-8);

            double u_1 = sol_lin_repl[i];
            double u_2 = sol_quads_repl[i];
            // max value of the solution is about 0.08, choose a tolerance of
            // 5% of that (wouldn't expect them to be exactly the same).
            TS_ASSERT_DELTA(u_1, u_2, 0.08*5e-2);
        }

        PetscTools::Destroy(solution_lin);
        PetscTools::Destroy(solution_quads);
    }

    void TestSolveLaplacianWithQuadratics2dReordered()
    {
        // Solve using quadratics..
        QuadraticMesh<2> quad_mesh;
        TrianglesMeshReader<2,2> mesh_reader1("mesh/test/data/square_128_elements_quadratic_reordered",2,1,false);
        quad_mesh.ConstructFromMeshReader(mesh_reader1);

        BoundaryConditionsContainer<2,2,1> bcc_quads;
        bcc_quads.DefineZeroDirichletOnMeshBoundary(&quad_mesh);

        QuadraticLaplacianAssemblerSolver<2> solver_quads(&quad_mesh, &bcc_quads);
        solver_quads.SetPdeConstants(1.0, 1.0);

        Vec solution_quads = solver_quads.Solve();
        ReplicatableVector sol_quads_repl(solution_quads);

        // Solve using linears
        TrianglesMeshReader<2,2> mesh_reader("mesh/test/data/square_128_elements");

        TetrahedralMesh<2,2> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        EllipticPdeWithLinearSource<2> pde(1.0, 1.0);

        BoundaryConditionsContainer<2,2,1> bcc_lin;
        bcc_lin.DefineZeroDirichletOnMeshBoundary(&mesh);

        SimpleLinearEllipticSolver<2,2> solver_lin(&mesh,&pde,&bcc_lin);

        Vec solution_lin = solver_lin.Solve();
        ReplicatableVector sol_lin_repl(solution_lin);

        /*
         * Compare results - the following assumes the vertex nodes in the
         * quad mesh are nodes 0-63, i.e. they come before all the internal
         * nodes.
         */
        TS_ASSERT_EQUALS(mesh.GetNumNodes(), quad_mesh.GetNumVertices());
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            unsigned quad_index=i;
            //Quad mesh has a minor permutation: vertex node 4 now appears at index 81
            if (i==4)
            {
                quad_index=81;
            }

            double lin_x = mesh.GetNode(i)->rGetLocation()[0];
            double lin_y = mesh.GetNode(i)->rGetLocation()[1];
            double quad_x = quad_mesh.GetNode(quad_index)->rGetLocation()[0];
            double quad_y = quad_mesh.GetNode(quad_index)->rGetLocation()[1];
            TS_ASSERT_DELTA(lin_x, quad_x, 1e-8);
            TS_ASSERT_DELTA(lin_y, quad_y, 1e-8);

            double u_1 = sol_lin_repl[i];
            double u_2 = sol_quads_repl[quad_index];
            // max value of the solution is about 0.08, choose a tolerance of
            // 5% of that (wouldn't expect them to be exactly the same).
            TS_ASSERT_DELTA(u_1, u_2, 0.08*5e-2);
        }

        PetscTools::Destroy(solution_lin);
        PetscTools::Destroy(solution_quads);
    }


    void TestSolveLaplacianWithQuadratics3d()
    {
        // Solve using quadratics
        QuadraticMesh<3> quad_mesh;
        TrianglesMeshReader<3,3> mesh_reader1("mesh/test/data/cube_1626_elements_quadratic",2,1,false);
        quad_mesh.ConstructFromMeshReader(mesh_reader1);

        BoundaryConditionsContainer<3,3,1> bcc_quads;
        bcc_quads.DefineZeroDirichletOnMeshBoundary(&quad_mesh);

        QuadraticLaplacianAssemblerSolver<3> solver_quads(&quad_mesh, &bcc_quads);
        solver_quads.SetPdeConstants(1.0, 1.0);

        Vec solution_quads = solver_quads.Solve();
        ReplicatableVector sol_quads_repl(solution_quads);

        // Solve using linears
        TrianglesMeshReader<3,3> mesh_reader("mesh/test/data/cube_1626_elements");
        TetrahedralMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(mesh_reader);

        EllipticPdeWithLinearSource<3> pde(1.0, 1.0);

        BoundaryConditionsContainer<3,3,1> bcc_lin;
        bcc_lin.DefineZeroDirichletOnMeshBoundary(&mesh);

        SimpleLinearEllipticSolver<3,3> solver_lin(&mesh,&pde,&bcc_lin);

        Vec solution_lin = solver_lin.Solve();
        ReplicatableVector sol_lin_repl(solution_lin);

        // compare results - the following assumes the vertex nodes in the
        // quad mesh come before all the internal nodes
        for (unsigned i=0; i<mesh.GetNumNodes(); i++)
        {
            double u_1 = sol_lin_repl[i];
            double u_2 = sol_quads_repl[i];

            // max value of the solution is about 0.059, choose a tolerance of
            // 5% of that (wouldn't expect them to be exactly the same).
            TS_ASSERT_DELTA(u_1, u_2, 0.059*5e-2);

            //double x = mesh.GetNode(i)->rGetLocation()[0];
            //double y = mesh.GetNode(i)->rGetLocation()[1];
            //double z = mesh.GetNode(i)->rGetLocation()[2];
            //
            //std::cout << x << " " << y << " " << z << " " << u_1 << " "
            //          <<  u_2 << " " << fabs(u_1-u_2) << "\n";
        }

        PetscTools::Destroy(solution_lin);
        PetscTools::Destroy(solution_quads);
    }
};

#endif /*TESTSOLVELAPLACIANWITHQUADRATICS_HPP_*/
