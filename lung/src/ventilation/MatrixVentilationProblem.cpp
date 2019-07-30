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

#include "MatrixVentilationProblem.hpp"
#include "ReplicatableVector.hpp"
#include "Warnings.hpp"

MatrixVentilationProblem::MatrixVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex)
    : AbstractVentilationProblem(rMeshDirFilePath, rootIndex),
      mpLinearSystem(nullptr),
      mSolution(nullptr)
{

    // We solve for flux at every edge and for pressure at each node/bifurcation
    // Note pipe flow equation has 3 variables and flux balance has 3 variables (at a bifurcation)
    // preallocating 5 non-zeros allows for 4-way branching
    mSolution = PetscTools::CreateVec(mMesh.GetNumNodes()+mMesh.GetNumElements());
    mpLinearSystem = new LinearSystem(mSolution, 5u);
    mpLinearSystem->SetAbsoluteTolerance(1e-5);

#ifdef LUNG_USE_UMFPACK
// LCOV_EXCL_START
    mpLinearSystem->SetPcType("lu");
    PetscTools::SetOption("-pc_factor_mat_solver_package", "umfpack");
    mpLinearSystem->SetKspType("gmres");
// LCOV_EXCL_STOP
#elif LUNG_USE_KLU
// LCOV_EXCL_START
    mpLinearSystem->SetPcType("lu");
    PetscTools::SetOption("-pc_factor_mat_solver_package", "klu");
    mpLinearSystem->SetKspType("gmres");
// LCOV_EXCL_STOP
#else
    mpLinearSystem->SetPcType("jacobi");
    PetscTools::SetOption("-ksp_diagonal_scale","");
    PetscTools::SetOption("-ksp_diagonal_scale_fix","");
    mpLinearSystem->SetKspType("fgmres");
#endif

    mFluxScaling = 1;//mViscosity;
}

MatrixVentilationProblem::~MatrixVentilationProblem()
{
    if (mpLinearSystem)
    {
        delete mpLinearSystem;
    }
    if (mSolution)
    {
        PetscTools::Destroy(mSolution);
    }
}

void MatrixVentilationProblem::SetMeshInMilliMetres()
{
    AbstractVentilationProblem::SetMeshInMilliMetres();

    /* In SI units, the flux is measured in m^3 (rather small amounts) but the
     * pressure is measured in Pascal.  These differ by orders of magnitude at the top of the tree and
     * even more at the bottom of the tree where resistance is greater.  Here we scale the flux up to get
     * them closer in magnitude to the pressures
     */
    mFluxScaling = 1e10; // 1e9 would be enough to put the fluxes in mm^3/s rather than m^3/s we may scale further to account for small airways
}

void MatrixVentilationProblem::SetPressureAtBoundaryNode(const Node<3>& rNode, double pressure)
{
    if (rNode.IsBoundaryNode() == false)
    {
        EXCEPTION("Boundary conditions cannot be set at internal nodes");
    }
    unsigned pressure_index =  mMesh.GetNumElements() +  rNode.GetIndex();

    mpLinearSystem->SetMatrixElement(pressure_index, pressure_index,  1.0);
    mpLinearSystem->SetRhsVectorElement(pressure_index, pressure);
    PetscVecTools::SetElement(mSolution, pressure_index, pressure); // Make a good guess
}

void MatrixVentilationProblem::SetFluxAtBoundaryNode(const Node<3>& rNode, double flux)
{
    if (rNode.IsBoundaryNode() == false)
    {
        EXCEPTION("Boundary conditions cannot be set at internal nodes");
    }

    // In a <1,3> mesh a boundary node will be associated with exactly one edge.
    // Flux boundary conditions are set in the system matrix using
    // the node index for the row and the edge index for the column.
    // The row associated with the leaf node is used so that the edge's row
    // can still be used to solve for flux/pressure.
    unsigned edge_index = *( rNode.ContainingElementsBegin() );
    unsigned pressure_index =  mMesh.GetNumElements() +  rNode.GetIndex();

    mpLinearSystem->SetMatrixElement(pressure_index, edge_index,  1.0);
    mpLinearSystem->SetRhsVectorElement(pressure_index, flux*mFluxScaling);
    PetscVecTools::SetElement(mSolution, edge_index, flux*mFluxScaling); // Make a good guess
}

double MatrixVentilationProblem::GetFluxAtOutflow()
{
    if (PetscTools::IsSequential())
    {
        return PetscVecTools::GetElement(mSolution, mOutletNodeIndex) / mFluxScaling;
    }
    else
    {
        ///\todo #2300 Massively inefficient:
        std::vector<double> flux, pressure;
        GetSolutionAsFluxesAndPressures(flux, pressure);
        return flux[mOutletNodeIndex];
    }
}

void MatrixVentilationProblem::Assemble(bool dynamicReassemble)
{
    PetscInt lo, hi;
    mpLinearSystem->GetOwnershipRange(lo, hi);

    if (dynamicReassemble)
    {
        // Sanity checks
        assert(mDynamicResistance);
        // Check that mSolution is valid?
    }

    // Assemble the Poiseuille flow pipe equations
    // Poiseuille flow at each edge
    for (AbstractTetrahedralMesh<1,3>::ElementIterator iter = mMesh.GetElementIteratorBegin();
         iter != mMesh.GetElementIteratorEnd();
         ++iter)
    {
        unsigned element_index = iter->GetIndex();
        // Only assemble if the row for this element is locally owned
        if ((unsigned) lo <=  element_index && element_index < (unsigned) hi)
        {
            /* Poiseuille flow gives:
             *  pressure_node_1 - pressure_node_2 - resistance * flux = 0
             */
            // Note that we can only do this if mSolution is valid AND we own the local part
            // "flux" is only needed in the dynamic case
            double flux = PetscVecTools::GetElement(mSolution, element_index)/mFluxScaling;
            double resistance = CalculateResistance(*(iter), dynamicReassemble, flux);
            // The solution vector has flux on elements first, then pressure at nodes.
            unsigned pressure_index_0 =  mMesh.GetNumElements() +  iter->GetNodeGlobalIndex(0);
            unsigned pressure_index_1 =  mMesh.GetNumElements() +  iter->GetNodeGlobalIndex(1);

            mpLinearSystem->SetMatrixElement(element_index, element_index, -resistance/mFluxScaling);
            mpLinearSystem->SetMatrixElement(element_index, pressure_index_0,  1.0);
            mpLinearSystem->SetMatrixElement(element_index, pressure_index_1, -1.0);
        }
    }

    if (dynamicReassemble)
    {
        // In the reassemble case, we don't have to repeat the flux balance equations
        // because the matrix components will be the same.
        return;
    }
    // Assemble the flux-balance equations
    for (AbstractTetrahedralMesh<1,3>::NodeIterator iter =mMesh.GetNodeIteratorBegin();
          iter != mMesh.GetNodeIteratorEnd();
          ++iter)
    {
        if (!(iter->IsBoundaryNode()) )
        {
            unsigned pressure_index =  mMesh.GetNumElements() +  iter->GetIndex();
            /* Flux balance at each internal node (only one internal node in our case)
            * flux_in - flux_out_left - flux_out_right = 0
            */
            for (Node<3>::ContainingElementIterator element_iterator = iter->ContainingElementsBegin();
                    element_iterator != iter->ContainingElementsEnd();
                    ++element_iterator)
            {
                unsigned el_index = *element_iterator;
                //We regard flux as coming in if this node is listed second.
                double flux_out = 1.0;
                if (mMesh.GetElement(el_index)->GetNodeGlobalIndex(1) == iter->GetIndex())
                {
                    flux_out = -1.0;
                }
                mpLinearSystem->SetMatrixElement(pressure_index, el_index, flux_out);
            }
        }
    }


}

void MatrixVentilationProblem::Solve()
{
    Assemble();
    mpLinearSystem->AssembleFinalLinearSystem();
    PetscVecTools::Finalise(mSolution);
    //mpLinearSystem->DisplayMatrix();
    //mpLinearSystem->DisplayRhs();

    //This call is when the solution vector may have been used before (at the previous timestep for example)
    Vec last_solution = mSolution;
    mSolution = mpLinearSystem->Solve(last_solution);
    PetscTools::Destroy(last_solution);
    if (mDynamicResistance)
    {
        double relative_diff = DBL_MAX;
        Vec difference;
        VecDuplicate(mSolution, &difference);
        do
        {
            Assemble(true);
            mpLinearSystem->AssembleFinalLinearSystem();
            Vec old_solution = mSolution;
            mSolution = mpLinearSystem->Solve(old_solution); // However, seeding with the old solution doesn't appear to change the convergence of the KSP solver
            PetscVecTools::WAXPY(difference, -1.0, mSolution, old_solution);
            double l_inf_diff;
            VecNorm(difference, NORM_INFINITY, &l_inf_diff);
            PetscTools::Destroy(old_solution);

            double l_inf_soln;
            VecNorm(mSolution, NORM_INFINITY, &l_inf_soln);
            relative_diff = l_inf_diff/l_inf_soln;
        }
        while (relative_diff > DBL_EPSILON);
        PetscTools::Destroy(difference);
    }
    //PetscVecTools::Display(mSolution);
}

void MatrixVentilationProblem::GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges,
                                                         std::vector<double>& rPressuresOnNodes)
{
    ReplicatableVector solution_vector_repl( mSolution );
    unsigned num_elem = mMesh.GetNumElements();
//    double max_flux = 0.0;
//    double max_scaled_flux = 0.0;
//    double max_pressure = 0.0;

    rFluxesOnEdges.resize(num_elem);
    for (unsigned i=0; i<num_elem; i++)
    {
        rFluxesOnEdges[i] = solution_vector_repl[i]/mFluxScaling;
//        if (fabs(solution_vector_repl[i]) > max_scaled_flux)
//        {
//            max_scaled_flux = fabs(solution_vector_repl[i]);
//            max_flux = rFluxesOnEdges[i];
//        }
    }

    rPressuresOnNodes.resize(mMesh.GetNumNodes());
    for (unsigned i=0; i<mMesh.GetNumNodes(); i++)
    {
        rPressuresOnNodes[i] = solution_vector_repl[i+num_elem];
//        if (fabs(rPressuresOnNodes[i])>max_pressure)
//        {
//            max_pressure = fabs(rPressuresOnNodes[i]);
//        }
    }
//    PRINT_5_VARIABLES(max_flux, max_scaled_flux, max_pressure, max_scaled_flux/max_pressure, max_flux/max_pressure);
}


