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

#include "VentilationProblem.hpp"
#include "TrianglesMeshReader.hpp"
#include "Warnings.hpp"
//#include "Debug.hpp"
//#include "Timer.hpp"
VentilationProblem::VentilationProblem(const std::string& rMeshDirFilePath, unsigned rootIndex)
    : AbstractVentilationProblem(rMeshDirFilePath, rootIndex),
      mDynamicResistance(false),
      mRadiusOnEdge(false),
      mDensity(1.15e-6),
      mFluxGivenAtInflow(false),
      mTerminalInteractionMatrix(NULL),
      mNumNonZeroesPerRow(25u), //See note in header definition
      mTerminalFluxChangeVector(NULL),
      mTerminalPressureChangeVector(NULL),
      mTerminalKspSolver(NULL)
{
    Initialise(rMeshDirFilePath);
}


VentilationProblem::VentilationProblem(AbstractAcinarUnitFactory* pAcinarUnitFactory,
                                       const std::string& rMeshDirFilePath,
                                       unsigned rootIndex) :
                                                     AbstractVentilationProblem(rMeshDirFilePath, rootIndex),
                                                     mDynamicResistance(false),
                                               mRadiusOnEdge(false),
                                               mDensity(1.15e-6),
                                               mFluxGivenAtInflow(false),
                                               mTerminalInteractionMatrix(NULL),
                                               mTerminalFluxChangeVector(NULL),
                                               mTerminalPressureChangeVector(NULL),
                                               mTerminalKspSolver(NULL)
{
    Initialise(rMeshDirFilePath);
    pAcinarUnitFactory->SetMesh(&mMesh);

    ///\todo We need to do some parallel here in order to load balance

    //Set up acinar units using the factory
    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mMesh.GetBoundaryNodeIteratorBegin();
          iter != mMesh.GetBoundaryNodeIteratorEnd();
          ++iter)
    {
        if ((*iter)->GetIndex() != mOutletNodeIndex)
        {

            AbstractAcinarUnit* p_acinus = pAcinarUnitFactory->CreateAcinarUnitForNode(*iter);

            //Sets the terminal bronchiole resistance of the acinar unit.
            c_vector<double, 3> dummy;
            double length;
            unsigned edge_index = *( (*iter)->ContainingElementsBegin() );
            mMesh.GetWeightedDirectionForElement(edge_index, dummy, length);

            double radius = (*iter)->rGetNodeAttributes()[0];

            double resistance = 8.0*mViscosity*length/(M_PI*SmallPow(radius, 4));
            p_acinus->SetTerminalBronchioleResistance(resistance);

            mAcinarUnits[(*iter)->GetIndex()] = p_acinus;
        }
    }
}


void VentilationProblem::Initialise(const std::string& rMeshDirFilePath)
{
    mFlux.resize(mMesh.GetNumElements());
    mPressure.resize(mMesh.GetNumNodes());
}


VentilationProblem::~VentilationProblem()
{
    if(mAcinarUnits.size() > 0)
    {
        for (unsigned i=0; i<mAcinarUnits.size(); i++)
        {
            delete mAcinarUnits[i];
        }
    }

    /* Remove the PETSc context used in the iterative solver */
    if (mTerminalInteractionMatrix)
    {
        PetscTools::Destroy(mTerminalInteractionMatrix);
        PetscTools::Destroy(mTerminalFluxChangeVector);
        PetscTools::Destroy(mTerminalPressureChangeVector);
        KSPDestroy(PETSC_DESTROY_PARAM(mTerminalKspSolver));
    }
}


void VentilationProblem::SolveDirectFromFlux()
{
    /* Work back through the node iterator looking for internal nodes: bifurcations or joints
     *
     * Each parent flux is equal to the sum of its children
     * Note that we can't iterate all the way back to the root node, but that's okay
     * because the root is not a bifurcation.  Also note that we can't use a NodeIterator
     * because it does not have operator--
     */
    for (unsigned node_index = mMesh.GetNumNodes() - 1; node_index > 0; --node_index)
    {
        Node<3>* p_node = mMesh.GetNode(node_index);
        if (p_node->IsBoundaryNode() == false)
        {
            Node<3>::ContainingElementIterator element_iterator = p_node->ContainingElementsBegin();
            // This assertion will trip if a node is not actually in the airway tree and will prevent operator++ from hanging
            assert(element_iterator != p_node->ContainingElementsEnd());
            unsigned parent_index = *element_iterator;
            ++element_iterator;

            for (mFlux[parent_index]=0.0; element_iterator != p_node->ContainingElementsEnd(); ++element_iterator)
            {
                mFlux[parent_index] += mFlux[*element_iterator];
            }
        }
    }
    // Poiseuille flow at each edge
    for (AbstractTetrahedralMesh<1,3>::ElementIterator iter = mMesh.GetElementIteratorBegin();
         iter != mMesh.GetElementIteratorEnd();
         ++iter)
    {
        /* Poiseuille flow gives:
         *  pressure_node_1 - pressure_node_2 - resistance * flux = 0
         */
        double flux = mFlux[iter->GetIndex()];
        double resistance = CalculateResistance(*iter, mDynamicResistance, flux);
        unsigned pressure_index_parent =  iter->GetNodeGlobalIndex(0);
        unsigned pressure_index_child  =  iter->GetNodeGlobalIndex(1);
        mPressure[pressure_index_child] = mPressure[pressure_index_parent] - resistance*flux;
    }
}

void VentilationProblem::SetupIterativeSolver()
{
    //double start = Timer::GetElapsedTime();
    mNumNonZeroesPerRow = std::min(mNumNonZeroesPerRow, mMesh.GetNumBoundaryNodes()-1);
    MatCreateSeqAIJ(PETSC_COMM_SELF, mMesh.GetNumBoundaryNodes()-1, mMesh.GetNumBoundaryNodes()-1, mNumNonZeroesPerRow, NULL, &mTerminalInteractionMatrix);
    PetscMatTools::SetOption(mTerminalInteractionMatrix, MAT_SYMMETRIC);
    PetscMatTools::SetOption(mTerminalInteractionMatrix, MAT_SYMMETRY_ETERNAL);

    /* Map each edge to its terminal descendants so that we can keep track
     * of which terminals can affect each other via a particular edge.
     */
    mEdgeDescendantNodes.resize(mMesh.GetNumElements());
    unsigned terminal_index=0;
    //First set up all the boundary nodes
    for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mMesh.GetBoundaryNodeIteratorBegin();
              iter != mMesh.GetBoundaryNodeIteratorEnd();
              ++iter)
    {
        unsigned node_index = (*iter)->GetIndex();
        if (node_index != mOutletNodeIndex)
        {
            unsigned parent_index =  *((*iter)->ContainingElementsBegin());
            mTerminalToNodeIndex[terminal_index] = node_index;
            mTerminalToEdgeIndex[terminal_index] = parent_index;
            mEdgeDescendantNodes[parent_index].insert(terminal_index++);
        }
    }
    // The outer loop here is for special cases where we find an internal node before its descendants
    // In this case we need to scan the tree more than once (up to log(N)) before the sets are correctly propagated
    while (mEdgeDescendantNodes[mOutletNodeIndex].size() != terminal_index)
    {
        //Work back up the tree making the unions of the sets of descendants
        for (unsigned node_index = mMesh.GetNumNodes() - 1; node_index > 0; --node_index)
        {
            Node<3>* p_node = mMesh.GetNode(node_index);
            if (p_node->IsBoundaryNode() == false)
            {
                Node<3>::ContainingElementIterator element_iterator = p_node->ContainingElementsBegin();
                unsigned parent_index = *element_iterator;
                ++element_iterator;
                for (;element_iterator != p_node->ContainingElementsEnd(); ++element_iterator)
                {
                    mEdgeDescendantNodes[parent_index].insert(mEdgeDescendantNodes[*element_iterator].begin(),mEdgeDescendantNodes[*element_iterator].end());
                }
            }
        }
    }

    FillInteractionMatrix(false);

    assert( terminal_index == mMesh.GetNumBoundaryNodes()-1);
    VecCreateSeq(PETSC_COMM_SELF, terminal_index, &mTerminalFluxChangeVector);
    VecCreateSeq(PETSC_COMM_SELF, terminal_index, &mTerminalPressureChangeVector);

    KSPCreate(PETSC_COMM_SELF, &mTerminalKspSolver);
#if ( PETSC_VERSION_MAJOR==3 && PETSC_VERSION_MINOR>=5 )
    KSPSetOperators(mTerminalKspSolver, mTerminalInteractionMatrix, mTerminalInteractionMatrix);
#else
    KSPSetOperators(mTerminalKspSolver, mTerminalInteractionMatrix, mTerminalInteractionMatrix, SAME_PRECONDITIONER);
#endif
    KSPSetFromOptions(mTerminalKspSolver);
    KSPSetUp(mTerminalKspSolver);
//    PRINT_VARIABLE(Timer::GetElapsedTime() - start);
}

void VentilationProblem::FillInteractionMatrix(bool redoExisting)
{
    assert(!redoExisting);
    if (redoExisting)
    {
//        MatZeroEntries(mTerminalInteractionMatrix);
    }
    // Use the descendant sets to build the mTerminalInteractionMatrix structure of resistances
    for (AbstractTetrahedralMesh<1,3>::ElementIterator iter = mMesh.GetElementIteratorBegin();
            iter != mMesh.GetElementIteratorEnd();
            ++iter)
    {
        unsigned parent_index = iter->GetIndex();
        double parent_resistance=0.0;
        if (redoExisting)
        {
 //           assert(mDynamicResistance);
 //           parent_resistance=CalculateResistance(*iter);
        }
        else
        {
            parent_resistance=CalculateResistance(*iter, true, mFlux[parent_index]);
        }
        std::vector<PetscInt> indices( mEdgeDescendantNodes[parent_index].begin(), mEdgeDescendantNodes[parent_index].end() );
        if (mEdgeDescendantNodes[parent_index].size() <= mNumNonZeroesPerRow)
        {
            std::vector<double> resistance_to_add(indices.size()*indices.size(), parent_resistance);

            MatSetValues(mTerminalInteractionMatrix,
                         indices.size(), (PetscInt*) &indices[0],
                         indices.size(), (PetscInt*) &indices[0], &resistance_to_add[0], ADD_VALUES);
        }
        else
        {
            ///\todo Does this make so much difference?
            for (unsigned i=0; i<indices.size(); i++)
            {
                MatSetValue(mTerminalInteractionMatrix, indices[i], indices[i], parent_resistance, ADD_VALUES);
            }
        }
    }
    PetscMatTools::Finalise(mTerminalInteractionMatrix);
}

void VentilationProblem::SolveIterativelyFromPressure()
{
    if (mTerminalInteractionMatrix == NULL)
    {
        SetupIterativeSolver();
    }
//    double start = Timer::GetElapsedTime();
    /* Now use the pressure boundary conditions to determine suitable flux boundary conditions
     * and iteratively update them until we are done
     */
    assert(mPressure[mOutletNodeIndex] == mPressureCondition[mOutletNodeIndex]);

    unsigned max_iterations=1000;
    unsigned num_terminals = mMesh.GetNumBoundaryNodes()-1u;
    double pressure_tolerance = 1e-4;
    bool converged=false;
    double last_norm_pressure_change;
    Vec old_terminal_pressure_change;
    VecDuplicate(mTerminalPressureChangeVector, &old_terminal_pressure_change);
    for (unsigned iteration = 0; iteration < max_iterations && converged==false; iteration++)
    {
        for (unsigned terminal=0; terminal<num_terminals; terminal++)
        {
            unsigned node_index = mTerminalToNodeIndex[terminal];
            // How far we are away from matching this boundary condition.
            double delta_pressure = mPressure[node_index] - mPressureCondition[node_index];
            // Offset the first iteration
            if (iteration == 0)
            {
                delta_pressure += mPressureCondition[mOutletNodeIndex];
            }
            VecSetValue(mTerminalPressureChangeVector, terminal, delta_pressure, INSERT_VALUES);
        }
        double norm_pressure_change;
        VecNorm(mTerminalPressureChangeVector, NORM_2, &last_norm_pressure_change);

        if (last_norm_pressure_change < pressure_tolerance)
        {
            converged = true;
            break;
        }

        VecCopy(mTerminalPressureChangeVector, old_terminal_pressure_change);
        KSPSolve(mTerminalKspSolver, mTerminalPressureChangeVector, mTerminalFluxChangeVector);
        double* p_mTerminalFluxChangeVector;
        VecGetArray(mTerminalFluxChangeVector, &p_mTerminalFluxChangeVector);



        for (unsigned terminal=0; terminal<num_terminals; terminal++)
        {
            double estimated_mTerminalFluxChangeVector=p_mTerminalFluxChangeVector[terminal];
            unsigned edge_index = mTerminalToEdgeIndex[terminal];
            mFlux[edge_index] +=  estimated_mTerminalFluxChangeVector;
        }
        SolveDirectFromFlux();
        /* Look at the magnitude of the response */

        for (unsigned terminal=0; terminal<num_terminals; terminal++)
        {
            unsigned node_index = mTerminalToNodeIndex[terminal];
            // How far we are away from matching this boundary condition.
            double delta_pressure = mPressure[node_index] - mPressureCondition[node_index];
            // Offset the first iteration
            if (iteration == 0)
            {
                delta_pressure += mPressureCondition[mOutletNodeIndex];
            }
            VecSetValue(mTerminalPressureChangeVector, terminal, delta_pressure, INSERT_VALUES);
        }
        VecNorm(mTerminalPressureChangeVector, NORM_2, &norm_pressure_change);
        if (norm_pressure_change < pressure_tolerance)
        {
            converged = true;
            break;
        }
        double pressure_change_dot_product;
        VecDot(mTerminalPressureChangeVector, old_terminal_pressure_change, &pressure_change_dot_product);
        if (pressure_change_dot_product < 0.0)
        {
            /* The pressure correction has changed sign
             *  * so we have overshot the root
             *  * back up by setting a correction factor
             */
            double terminal_flux_correction =  last_norm_pressure_change / (last_norm_pressure_change + norm_pressure_change) - 1.0;
///\todo some sqrt response?
//            if (mDynamicResistance)
//            {
//                terminal_flux_correction *= 0.99;
//            }
            for (unsigned terminal=0; terminal<num_terminals; terminal++)
            {
                double estimated_mTerminalFluxChangeVector=p_mTerminalFluxChangeVector[terminal];
                unsigned edge_index = mTerminalToEdgeIndex[terminal];
                mFlux[edge_index] +=  terminal_flux_correction*estimated_mTerminalFluxChangeVector;
            }
            SolveDirectFromFlux();
        }
    }
    if(!converged)
    {
        NEVER_REACHED;
    }

    PetscTools::Destroy(old_terminal_pressure_change);
//    PRINT_2_VARIABLES(Timer::GetElapsedTime() - start, mDynamicResistance);
}

double VentilationProblem::CalculateResistance(Element<1,3>& rElement, bool usePedley, double flux)
{
    //Resistance is based on radius, length and viscosity
    double radius = 0.0;
    if (mRadiusOnEdge)
    {
        radius = rElement.GetAttribute();
    }
    else
    {
        radius = ( rElement.GetNode(0)->rGetNodeAttributes()[0] + rElement.GetNode(1)->rGetNodeAttributes()[0]) / 2.0;
    }

    c_vector<double, 3> dummy;
    double length;
    mMesh.GetWeightedDirectionForElement(rElement.GetIndex(), dummy, length);

    double resistance = 8.0*mViscosity*length/(M_PI*SmallPow(radius, 4));
    if ( usePedley )
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
        double reynolds_number = fabs( 2.0 * mDensity * flux / (mViscosity * M_PI * radius) );
        double c = 1.85;
        double z = (c/4.0) * sqrt(reynolds_number * radius / length);
        // Pedley's method will only increase the resistance
        if (z > 1.0)
        {
///\todo #2300 A small reduction in density (1.51 -> 1.15) has stopped unit test from triggering
/// Pedley correction.
#define COVERAGE_IGNORE            
            resistance *= z;
#undef COVERAGE_IGNORE
        }
    }
    return resistance;
}

void VentilationProblem::SetOutflowPressure(double pressure)
{
    SetPressureAtBoundaryNode(*(mMesh.GetNode(mOutletNodeIndex)), pressure);
    mPressure[mOutletNodeIndex] = pressure;
}

void VentilationProblem::SetConstantInflowPressures(double pressure)
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

void VentilationProblem::SetConstantInflowFluxes(double flux)
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

void VentilationProblem::SetPressureAtBoundaryNode(const Node<3>& rNode, double pressure)
{
    if (rNode.IsBoundaryNode() == false)
    {
        EXCEPTION("Boundary conditions cannot be set at internal nodes");
    }
    assert(mFluxGivenAtInflow == false);

    // Store the requirement in a map for the direct solver
    mPressureCondition[rNode.GetIndex()] = pressure;
}

double VentilationProblem::GetPressureAtBoundaryNode(const Node<3>& rNode)
{
// //This doesn't actually present any problems -- we could call it GetPressureAtNode
//    if (rNode.IsBoundaryNode() == false)
//    {
//        EXCEPTION("Boundary conditions cannot be got at internal nodes");
//    }
    return mPressure[rNode.GetIndex()];
}

double VentilationProblem::GetFluxAtOutflow()
{
    return mFlux[mOutletNodeIndex];
}

void VentilationProblem::SetFluxAtBoundaryNode(const Node<3>& rNode, double flux)
{
    if (rNode.IsBoundaryNode() == false)
    {
        EXCEPTION("Boundary conditions cannot be set at internal nodes");
    }
    mFluxGivenAtInflow = true;

    // In a <1,3> mesh a boundary node will be associated with exactly one edge.
    unsigned edge_index = *( rNode.ContainingElementsBegin() );

    // Seed the information for a direct solver
    mFlux[edge_index] = flux;
}



void VentilationProblem::Solve()
{
    if (mFluxGivenAtInflow)
    {
        SolveDirectFromFlux();
    }
    else
    {
        SolveIterativelyFromPressure();
    }

}


void VentilationProblem::GetSolutionAsFluxesAndPressures(std::vector<double>& rFluxesOnEdges,
                                                         std::vector<double>& rPressuresOnNodes)
{
    rFluxesOnEdges = mFlux;
    rPressuresOnNodes = mPressure;
}


void VentilationProblem::SetRadiusOnEdge(bool isOnEdges)
{
    mRadiusOnEdge = isOnEdges;
}

TetrahedralMesh<1, 3>& VentilationProblem::rGetMesh()
{
    return mMesh;
}


#ifdef CHASTE_VTK

void VentilationProblem::WriteVtk(const std::string& rDirName, const std::string& rFileBaseName)
{
    VtkMeshWriter<1, 3> vtk_writer(rDirName, rFileBaseName, false);
    AddDataToVtk(vtk_writer, "");
    vtk_writer.WriteFilesUsingMesh(mMesh);

}

void VentilationProblem::AddDataToVtk(VtkMeshWriter<1, 3>& rVtkWriter,
                                      const std::string& rSuffix)
{
    std::vector<double> pressures;
    std::vector<double> fluxes;
    GetSolutionAsFluxesAndPressures(fluxes, pressures);
    rVtkWriter.AddCellData("Flux"+rSuffix, fluxes);
    rVtkWriter.AddPointData("Pressure"+rSuffix, pressures);

    if(mAcinarUnits.size() > 0)
    {
        std::vector<double> volumes(mMesh.GetNumNodes());
        std::vector<double> stretch_ratios(mMesh.GetNumNodes());

        for (AbstractTetrahedralMesh<1,3>::BoundaryNodeIterator iter = mMesh.GetBoundaryNodeIteratorBegin();
                 iter != mMesh.GetBoundaryNodeIteratorEnd();
                 ++iter)
        {
            if( (*iter)->GetIndex() != mOutletNodeIndex)
            {
                volumes[(*iter)->GetIndex()] = mAcinarUnits[(*iter)->GetIndex()]->GetVolume();
                stretch_ratios[(*iter)->GetIndex()] = mAcinarUnits[(*iter)->GetIndex()]->GetStretchRatio();
            }
        }
        rVtkWriter.AddPointData("Volume"+rSuffix, volumes);
        rVtkWriter.AddPointData("Stretch"+rSuffix, stretch_ratios);
    }
}

#endif // CHASTE_VTK


void VentilationProblem::Solve(TimeStepper& rTimeStepper,
        void (*pBoundaryConditionFunction)(VentilationProblem*, TimeStepper& rTimeStepper, const Node<3>&),
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
                pBoundaryConditionFunction(this, rTimeStepper, *(*iter));
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

void VentilationProblem::SolveProblemFromFile(const std::string& rInFilePath, const std::string& rOutFileDir, const std::string& rOutFileName)
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






