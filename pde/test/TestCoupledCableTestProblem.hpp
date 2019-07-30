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

#ifndef TESTCOUPLEDCABLETESTPROBLEM_HPP_
#define TESTCOUPLEDCABLETESTPROBLEM_HPP_

#include <cxxtest/TestSuite.h>


#include "MixedDimensionMesh.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"
#include "StiffnessMatrixAssembler.hpp"
#include "PetscSetupAndFinalize.hpp"
#include "AbstractStaticLinearPdeSolver.hpp"
#include "AbstractFeCableIntegralAssembler.hpp"
#include "ConstBoundaryCondition.hpp"
#include "NaturalNeumannSurfaceTermAssembler.hpp"
#include "OutputFileHandler.hpp"


////////////////////////////////////////////////////////////////
//
//
//    The test solves the test problem given in the pdf
//    attached to ticket #1835.
//
//
////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////
// Assembler for integrals over cable elements
//////////////////////////////////////////////////////////
template<unsigned DIM>
class CoupledCableTestProblemCableComponentAssembler
    :  public AbstractFeCableIntegralAssembler<DIM,DIM,2,false,true,NORMAL>
{
private:
    static const unsigned PROBLEM_DIM=2;
    double mBeta;

    c_matrix<double,PROBLEM_DIM*2,PROBLEM_DIM*2> ComputeCableMatrixTerm(
        c_vector<double, 2>& rPhi,
        c_matrix<double, DIM, 2>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double,PROBLEM_DIM, DIM>& rGradU,
        Element<1,DIM>* pElement)
    {
        c_matrix<double,PROBLEM_DIM*2, PROBLEM_DIM*2> ret = zero_matrix<double>(PROBLEM_DIM*2, PROBLEM_DIM*2);

        double sigma_i = 1+rX[2] + 0.5*rX[2]*rX[2];

        for (unsigned i=0; i<2; i++)
        {
            for (unsigned j=0; j<2; j++)
            {
                ret(2*i,  2*j)   =  mBeta*rPhi(i)*rPhi(j);
                ret(2*i+1,2*j)   = -mBeta*rPhi(i)*rPhi(j);
                ret(2*i,  2*j+1) = -mBeta*rPhi(i)*rPhi(j);
                ret(2*i+1,2*j+1) =  mBeta*rPhi(i)*rPhi(j);

                for (unsigned dim=0; dim<DIM; dim++)
                {
                    ret(2*i+1,2*j+1) += sigma_i*rGradPhi(dim,i)*rGradPhi(dim,j);
                }
            }
        }

        return ret;
    }

public:
    CoupledCableTestProblemCableComponentAssembler(MixedDimensionMesh<DIM,DIM>* pMesh, double beta)
        : AbstractFeCableIntegralAssembler<DIM,DIM,2,false,true,NORMAL>(pMesh),
          mBeta(beta)
    {
    }
};


//////////////////////////////////////////////////////////
// Assembler for integrals over volume elements
/////////////////////////////////////////////////////////
template<unsigned DIM>
class CoupledCableTestProblemVolumeComponentAssembler
    : public AbstractFeVolumeIntegralAssembler<DIM,DIM,2,true,true,NORMAL>
{
private:
   static const unsigned PROBLEM_DIM=2;

   c_matrix<double,PROBLEM_DIM*(DIM+1),PROBLEM_DIM*(DIM+1)> ComputeMatrixTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,PROBLEM_DIM>& rU,
        c_matrix<double,PROBLEM_DIM, DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        c_matrix<double,PROBLEM_DIM*(DIM+1),PROBLEM_DIM*(DIM+1)> ret = zero_matrix<double>(PROBLEM_DIM*(DIM+1),PROBLEM_DIM*(DIM+1));

        for (unsigned i=0; i<DIM+1; i++)
        {
            for (unsigned j=0; j<DIM+1; j++)
            {
                // stiffness matrix in 'top-left block'
                for (unsigned dim=0; dim<DIM; dim++)
                {
                    ret(2*i,2*j) += rGradPhi(dim,i)*rGradPhi(dim,j);
                }
            }
        }
        return ret;
    }

    c_vector<double,2*(DIM+1)> ComputeVectorTerm(
        c_vector<double, DIM+1>& rPhi,
        c_matrix<double, DIM, DIM+1>& rGradPhi,
        ChastePoint<DIM>& rX,
        c_vector<double,2>& rU,
        c_matrix<double,2,DIM>& rGradU,
        Element<DIM,DIM>* pElement)
    {
        return zero_vector<double>(2*(DIM+1));
    }

public:
    CoupledCableTestProblemVolumeComponentAssembler(AbstractTetrahedralMesh<DIM,DIM>* pMesh)
        : AbstractFeVolumeIntegralAssembler<DIM,DIM,2,true,true,NORMAL>(pMesh)
    {
    }
};

//////////////////////////////////////////////////////////////
// Solver which uses assemblers to create the linear system
//////////////////////////////////////////////////////////////

template<unsigned DIM>
class CoupledCableTestProblemSolver: public AbstractStaticLinearPdeSolver<DIM,DIM,2>
{
private:
    CoupledCableTestProblemVolumeComponentAssembler<DIM>* mpVolumeIntegralsAssembler;
    CoupledCableTestProblemCableComponentAssembler<DIM>* mpCableIntegralsAssembler;
    NaturalNeumannSurfaceTermAssembler<DIM,DIM,2>* mpNeumannSurfaceTermsAssembler;

    BoundaryConditionsContainer<DIM,DIM,2>* mpBoundaryConditions;

    void SetupLinearSystem(Vec currentSolution, bool computeMatrix)
    {
        // assemble the volume integral and Neumann surface part
        mpVolumeIntegralsAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix(), true);
        mpVolumeIntegralsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), true);
        mpVolumeIntegralsAssembler->Assemble();

        mpNeumannSurfaceTermsAssembler->SetVectorToAssemble(this->mpLinearSystem->rGetRhsVector(), false);
        mpNeumannSurfaceTermsAssembler->Assemble();

        // assemble the cable integral part
        mpCableIntegralsAssembler->SetMatrixToAssemble(this->mpLinearSystem->rGetLhsMatrix(), false /*don't zero matrix!*/);
        mpCableIntegralsAssembler->Assemble();

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->SwitchWriteModeLhsMatrix();

        // apply the Dirichlet boundary conditions
        mpBoundaryConditions->ApplyDirichletToLinearProblem(*(this->mpLinearSystem), computeMatrix);

        if (computeMatrix)
        {
            ApplyIdentityBlock();
        }

        this->mpLinearSystem->FinaliseRhsVector();
        this->mpLinearSystem->FinaliseLhsMatrix();
    }

    // Hardcoded for the expected mesh
    void ApplyIdentityBlock()
    {
        for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator
                iter=this->mpMesh->GetNodeIteratorBegin();
             iter != this->mpMesh->GetNodeIteratorEnd();
             ++iter)
        {
            Node<DIM>& r_node = *iter;
            double x = r_node.rGetLocation()[0];
            double y = r_node.rGetLocation()[1];
            double r = sqrt(x*x+y*y);

            if (fabs(r) > 1e-6) // if r>0, ie not a cable node
            {
                // Put 1.0 on the diagonal - the rest of the row must already be zero
                PetscInt index = 2*iter->GetIndex()+1;
                PetscMatTools::SetElement(this->mpLinearSystem->rGetLhsMatrix(), index, index, 1.0);
            }
        }
    }


public:
    CoupledCableTestProblemSolver(MixedDimensionMesh<DIM,DIM>* pMesh,
                                  BoundaryConditionsContainer<DIM,DIM,2>* pBoundaryConditions,
                                  double beta)
         : AbstractStaticLinearPdeSolver<DIM,DIM,2>(pMesh),
           mpBoundaryConditions(pBoundaryConditions)
    {
        mpVolumeIntegralsAssembler = new CoupledCableTestProblemVolumeComponentAssembler<DIM>(pMesh);
        mpCableIntegralsAssembler = new CoupledCableTestProblemCableComponentAssembler<DIM>(pMesh,beta);
        mpNeumannSurfaceTermsAssembler = new NaturalNeumannSurfaceTermAssembler<DIM,DIM,2>(pMesh,pBoundaryConditions);
    }

    ~CoupledCableTestProblemSolver()
    {
        delete mpCableIntegralsAssembler;
        delete mpVolumeIntegralsAssembler;
        delete mpNeumannSurfaceTermsAssembler;
    }
};


class MyBoundaryCondition : public AbstractBoundaryCondition<3>
{
private:
    double mFactor;

public:
    double GetValue(const ChastePoint<3>& rX) const
    {
        double r = sqrt(rX[0]*rX[0] + rX[1]*rX[1]);
        return mFactor*log(r);
    }

    MyBoundaryCondition(double factor)
        : AbstractBoundaryCondition<3>(),
          mFactor(factor)
    {
    }
};


class TestCoupledCableTestProblem : public CxxTest::TestSuite
{
public:
    void TestSolvingTestProblem()
    {
        std::string mesh_base("mesh/test/data/mixed_dimension_meshes/cylinder_refined");
        TrianglesMeshReader<3,3> reader(mesh_base);
        MixedDimensionMesh<3,3> mesh;
        mesh.ConstructFromMeshReader(reader);

        BoundaryConditionsContainer<3,3,2> bcc;

        //////////////////////////////////////////
        // Dirichlet BCs to phi_e and phi_i
        //////////////////////////////////////////

        for (MixedDimensionMesh<3,3>::BoundaryNodeIterator iter
                = mesh.GetBoundaryNodeIteratorBegin();
             iter != mesh.GetBoundaryNodeIteratorEnd();
             iter++)
        {
            const Node<3>* p_node = &(**iter); //Get pointer to the current node from the iterator
            double x = p_node->rGetLocation()[0];
            double y = p_node->rGetLocation()[1];
            double z = p_node->rGetLocation()[2];
            double r = sqrt(x*x+y*y);

            if (fabs(r-1)<1e-3)
            {
                // apply BC phi_e = 0
                ConstBoundaryCondition<3>* p_bc = new ConstBoundaryCondition<3>(0.0);
                bcc.AddDirichletBoundaryCondition(p_node, p_bc, 0);
            }

            if (fabs(r)<1e-6 && fabs(z)<1e-6) // r=0, z=0, ie bottom cable node
            {
                // apply BC phi_i = 1
                ConstBoundaryCondition<3>* p_bc = new ConstBoundaryCondition<3>(1.0);
                bcc.AddDirichletBoundaryCondition(p_node, p_bc, 1);
            }

            if (fabs(r)<1e-6 && fabs(z-1)<1e-6) // r=0, z=1, ie top cable node
            {
                // apply BC phi_i = 2
                ConstBoundaryCondition<3>* p_bc = new ConstBoundaryCondition<3>(2.0);
                bcc.AddDirichletBoundaryCondition(p_node, p_bc, 1);
            }
        }

        //////////////////////////////////////////
        // Neumann BCs to phi_e
        //////////////////////////////////////////

        MyBoundaryCondition* p_top_neumann_bc = new MyBoundaryCondition(-1.0/(2*M_PI));
        MyBoundaryCondition* p_bottom_neumann_bc = new MyBoundaryCondition(1.0/(2*M_PI));

        for (MixedDimensionMesh<3,3>::BoundaryElementIterator iter
             = mesh.GetBoundaryElementIteratorBegin();
            iter != mesh.GetBoundaryElementIteratorEnd();
            ++iter)
        {
            c_vector<double,3> centroid = (*iter)->CalculateCentroid();
            double z = centroid(2);
            if (fabs(z-1.0) < 1e-5)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_top_neumann_bc, 0);
            }
            if (fabs(z) < 1e-5)
            {
                bcc.AddNeumannBoundaryCondition(*iter, p_bottom_neumann_bc, 0);
            }
        }

        // This is the radius of the true fibre.
        // For R >> 0, the approximation of the fibre as 1D manifold breaks down and the analytical solution is invalid
        double R = 0.015;
        double beta = 4*M_PI/(2*log(R)-1+4*M_PI);

        // Solve
        CoupledCableTestProblemSolver<3> cable_solver(&mesh,&bcc,beta);
        Vec result = cable_solver.Solve();

        ReplicatableVector result_repl(result);
        for (AbstractTetrahedralMesh<3,3>::NodeIterator current_node = mesh.GetNodeIteratorBegin();
             current_node != mesh.GetNodeIteratorEnd();
             ++current_node)
        {
            double x = current_node->GetPoint()[0];
            double y = current_node->GetPoint()[1];
            double z = current_node->GetPoint()[2];
            double r = sqrt(x*x+y*y);

            unsigned index = current_node->GetIndex();
            double phi_e = result_repl[2*index];
            double phi_i = result_repl[2*index+1];
            if (fabs(r)<1e-6)
            {
                double phi_i_exact = 1+z;
                // Tolerance is quite high, as the mesh is fairly coarse for this problem (node spacing ~2mm)
                TS_ASSERT_DELTA(phi_i, phi_i_exact, 3e-2);
            }
            else
            {
                double phi_e_exact = -(1+z)*log(r)/(2*M_PI);

                // Tolerance is quite high, as the mesh is fairly coarse for this problem (node spacing ~2mm)
                TS_ASSERT_DELTA(phi_e, phi_e_exact, 2e-1);

                // Errors should be higher closer to the fibre and lower further away.
                // The error weighted by r should be more uniform across the mesh.
                ///\todo: Ideally, this would be checked using the H^1_\alpha norm. Dependent on #1868.
                TS_ASSERT_DELTA(r*phi_e, r*phi_e_exact, 3e-2);

                // check dummy phi_i variable is correctly set to zero
                TS_ASSERT_DELTA(phi_i, 0.0, 1e-12);
            }
        }

        PetscTools::Destroy(result);
    }
};

#endif // TESTCOUPLEDCABLETESTPROBLEM_HPP_
