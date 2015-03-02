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

#include "MatrixVentilationProblem.hpp"
#include "TrianglesMeshReader.hpp"
#include "ReplicatableVector.hpp"
#include "Warnings.hpp"
//#include "Debug.hpp"

MatrixVentilationProblem::MatrixVentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex)
    : mOutletNodeIndex(rootIndex),
      mpLinearSystem(NULL),
      mDynamicResistance(false),
      mRadiusOnEdge(false),
      mViscosity(1.92e-5),
      mDensity(1.15),
      mSolution(NULL),
      mLengthScaling(1)
{
    TrianglesMeshReader<1,3> mesh_reader(rMeshDirFilePath);
    mMesh.ConstructFromMeshReader(mesh_reader);

    if (mMesh.GetNode(mOutletNodeIndex)->IsBoundaryNode() == false)
    {
        EXCEPTION("Outlet node is not a boundary node");
    }

    // We solve for flux at every edge and for pressure at each node/bifurcation
    // Note pipe flow equation has 3 variables and flux balance has 3 variables (at a bifurcation)
    // preallocating 5 non-zeros allows for 4-way branching
    mSolution = PetscTools::CreateVec(mMesh.GetNumNodes()+mMesh.GetNumElements());
    mpLinearSystem = new LinearSystem(mSolution, 5u);
    mpLinearSystem->SetAbsoluteTolerance(1e-5);

    #ifdef LUNG_USE_UMFPACK
    #define COVERAGE_IGNORE
        mpLinearSystem->SetPcType("lu");
        PetscOptionsSetValue("-pc_factor_mat_solver_package", "umfpack");
        mpLinearSystem->SetKspType("gmres");
    #undef COVERAGE_IGNORE
    #else
        mpLinearSystem->SetPcType("jacobi");
        PetscOptionsSetValue("-ksp_diagonal_scale","");
        PetscOptionsSetValue("-ksp_diagonal_scale_fix","");
        mpLinearSystem->SetKspType("fgmres");
    #endif

    /*
     * Set up the Acinar units at the terminals
     *
     * \todo Currently this is hard coded using a Swan acinar unit
     * (with functional residual capacity appropriate initial values). Long term this
     * should be factored out into a factory to allow setup of other acinar units and
     * other initial conditions
     */
    //unsigned num_acinar = mMesh.GetNumBoundaryNodes() - 1;
    //double acinus_volume = 1.2e9/num_acinar; //Assumes a residual capacity of 1.2l (x10^9 in mm^3)

    mFluxScaling = 1;//mViscosity;
//    mFluxScaling = DBL_MAX;  // Adjusted to the least resistant of the terminal branches
//    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mMesh.GetBoundaryNodeIteratorBegin();
//          iter != mMesh.GetBoundaryNodeIteratorEnd();
//          ++iter)
//    {
//        //Calculates the resistance of the terminal bronchiole.
//        //This should be updated dynamically during the simulation
//        c_vector<double, 3> dummy;
//        double length;
//        unsigned edge_index = *( (*iter)->ContainingElementsBegin() );
//        mMesh.GetWeightedDirectionForElement(edge_index, dummy, length);
//
//        double radius = (*iter)->rGetNodeAttributes()[0];
//
//        double resistance = 8.0*mViscosity*length/(M_PI*SmallPow(radius, 4));
//        if ((*iter)->GetIndex() == mOutletNodeIndex)
//        {
//            // Use the top of the tree to get an indication of resistances
//            mFluxScaling = resistance;
//        }
//        else
//        {
//            ///\todo We need to do some parallel here in order to load balance
//            Swan2012AcinarUnit* p_acinus = new Swan2012AcinarUnit;
//
//            p_acinus->SetStretchRatio(1.26); //Stretch ratio appropriate for a lung at functional residual capacity
//            p_acinus->SetUndeformedVolume(acinus_volume);
//            p_acinus->SetPleuralPressure(-0.49); //Pleural pressure at FRC in kPa
//            p_acinus->SetAirwayPressure(-0.49);
//            p_acinus->SetTerminalBronchioleResistance(resistance);
//
//            mAcinarUnits.push_back(p_acinus);
//            ///\todo make a map for node->acinar (rather than a vector)?
//        }
//    }
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
//    for (unsigned i=0; i<mAcinarUnits.size(); i++)
//    {
//        delete mAcinarUnits[i];
//    }
}

void MatrixVentilationProblem::SetMeshInMilliMetres()
{
    mLengthScaling = 1e-3;
}


void MatrixVentilationProblem::SetOutflowPressure(double pressure)
{
    SetPressureAtBoundaryNode(*(mMesh.GetNode(mOutletNodeIndex)), pressure);
}

void MatrixVentilationProblem::SetOutflowFlux(double flux)
{
    SetFluxAtBoundaryNode(*(mMesh.GetNode(mOutletNodeIndex)), flux);
}

void MatrixVentilationProblem::SetConstantInflowPressures(double pressure)
{
    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter =mMesh.GetBoundaryNodeIteratorBegin();
          iter != mMesh.GetBoundaryNodeIteratorEnd();
          ++iter)
     {
         if ((*iter)->GetIndex() != mOutletNodeIndex)
         {
             //Boundary conditions at each boundary/leaf node
             SetPressureAtBoundaryNode(*(*iter), pressure);
         }
     }
}

void MatrixVentilationProblem::SetConstantInflowFluxes(double flux)
{
    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter =mMesh.GetBoundaryNodeIteratorBegin();
          iter != mMesh.GetBoundaryNodeIteratorEnd();
          ++iter)
     {
         if ((*iter)->GetIndex() != mOutletNodeIndex)
         {
             SetFluxAtBoundaryNode(*(*iter), flux);
         }
     }
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
        //Only assemble if the row for this element is locally owned
        if ( (unsigned) lo <=  element_index && element_index < (unsigned) hi )
        {
            /* Poiseuille flow gives:
             *  pressure_node_1 - pressure_node_2 - resistance * flux = 0
             */
            //Resistance is based on radius, length and viscosity
            double radius = 0.0;
            if (mRadiusOnEdge)
            {
                radius = iter->GetAttribute();
            }
            else
            {
                radius = ( iter->GetNode(0)->rGetNodeAttributes()[0] + iter->GetNode(1)->rGetNodeAttributes()[0]) / 2.0;
            }

            c_vector<double, 3> dummy;
            double length;
            mMesh.GetWeightedDirectionForElement(element_index, dummy, length);

            radius *= mLengthScaling;
            length *= mLengthScaling;

            double resistance = 8.0*mViscosity*length/(M_PI*SmallPow(radius, 4));
            if ( dynamicReassemble )
            {
                /* Pedley et al. 1970
                 * http://dx.doi.org/10.1016/0034-5687(70)90094-0
                 * also Swan et al. 2012. 10.1016/j.jtbi.2012.01.042 (page 224)
                 * Standard Poiseuille equation is similar to Pedley's modified Eq1. and matches Swan Eq5.
                 * R_p = 128*mu*L/(pi*d^4) = 8*mu*L/(pi*r^4)
                 *
                 * Pedley Eq 2 and Swan Eq4:
                 *  Z = C/(4*sqrt(2)) * sqrt(Re*d/l) = (C/4)*sqrt(Re*r/l)
                 * Pedley suggests that C = 1.85
                 * R_r = Z*R_p
                 *
                 * Reynold's number in a pipe is
                 * Re = rho*v*d/mu (where d is a characteristic length scale - diameter of pipe)
                 * since flux = v*area
                 * Re = Q * d/(mu*area) = 2*rho*Q/(mu*pi*r) ... (see Swan p 224)
                 *
                 *
                 * The upshot of this calculation is that the resistance is scaled with sqrt(Q)
                 */
                // Note that we can only do this if mSolution is valid AND we own the local part
                double flux = PetscVecTools::GetElement(mSolution, element_index)/mFluxScaling;
                double reynolds_number = fabs( 2.0 * mDensity * flux / (mViscosity * M_PI * radius) );
                double c = 1.85;
                double z = (c/4.0) * sqrt(reynolds_number * radius / length);

                // Pedley's method will only increase the resistance
                if (z > 1.0)
                {
                    resistance *= z;
                }
            }
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
}

Vec MatrixVentilationProblem::GetSolution()
{
    return mSolution;
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


void MatrixVentilationProblem::SetRadiusOnEdge(bool isOnEdges)
{
    mRadiusOnEdge = isOnEdges;
}

TetrahedralMesh<1, 3>& MatrixVentilationProblem::rGetMesh()
{
    return mMesh;
}


#ifdef CHASTE_VTK

void MatrixVentilationProblem::WriteVtk(const std::string& rDirName, const std::string& rFileBaseName)
{
    VtkMeshWriter<1, 3> vtk_writer(rDirName, rFileBaseName, false);
    AddDataToVtk(vtk_writer, "");
    vtk_writer.WriteFilesUsingMesh(mMesh);

}

void MatrixVentilationProblem::AddDataToVtk(VtkMeshWriter<1, 3>& rVtkWriter,
        const std::string& rSuffix)
{
    std::vector<double> pressures;
    std::vector<double> fluxes;
    GetSolutionAsFluxesAndPressures(fluxes, pressures);
    rVtkWriter.AddCellData("Flux"+rSuffix, fluxes);
    rVtkWriter.AddPointData("Pressure"+rSuffix, pressures);
}

#endif // CHASTE_VTK


void MatrixVentilationProblem::Solve(TimeStepper& rTimeStepper,
        void (*pBoundaryConditionFunction)(MatrixVentilationProblem*, double, const Node<3>&),
        const std::string& rDirName, const std::string& rFileBaseName)
{
#ifdef CHASTE_VTK
    VtkMeshWriter<1, 3> vtk_writer(rDirName, rFileBaseName, false);
#endif

    bool first_step=true;
    while (!rTimeStepper.IsTimeAtEnd())
    {
        if (first_step)
        {
            // Do a solve at t=0 before advancing time.
            first_step = false;
        }
        else
        {
            rTimeStepper.AdvanceOneTimeStep();
        }
        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mMesh.GetBoundaryNodeIteratorBegin();
                 iter != mMesh.GetBoundaryNodeIteratorEnd();
                 ++iter )
        {
            if ((*iter)->GetIndex() != mOutletNodeIndex)
            {
                //Boundary conditions at each boundary/leaf node
                pBoundaryConditionFunction(this, rTimeStepper.GetTime(), *(*iter));
            }
        }

        // Regular solve
        Solve();

        std::ostringstream suffix_name;
        suffix_name <<  "_" << std::setw(6) << std::setfill('0') << rTimeStepper.GetTotalTimeStepsTaken();
#ifdef CHASTE_VTK
        AddDataToVtk(vtk_writer, suffix_name.str());
#endif
    }

#ifdef CHASTE_VTK
    vtk_writer.WriteFilesUsingMesh(mMesh);
#endif
}

void MatrixVentilationProblem::SolveProblemFromFile(const std::string& rInFilePath, const std::string& rOutFileDir, const std::string& rOutFileName)
{
    std::ifstream file(FileFinder(rInFilePath).GetAbsolutePath().c_str(), std::ios::binary);
    if (!file.is_open())
    {
        EXCEPTION("Could not open file "+rInFilePath);
    }
    std::string key, unit;
    double value;
    while (!file.eof())
    {
        file >> key >> value >> unit;
        if (file.fail())
        {
            break;
        }
        if (key == "RHO_AIR")
        {
            SetDensity(value);
        }
        else if (key == "MU_AIR")
        {
            SetViscosity(value);
        }
        else if (key == "PRESSURE_OUT")
        {
            SetOutflowPressure(value);
        }
        else if (key == "PRESSURE_IN")
        {
            SetConstantInflowPressures(value);
        }
        else
        {
            WARNING("The key "+ key+ " is not recognised yet");
        }
    }
    Solve();
#ifdef CHASTE_VTK
    WriteVtk(rOutFileDir, rOutFileName);
#endif
}






