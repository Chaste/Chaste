/*

Copyright (c) 2005-2016, University of Oxford.
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

#include "OffLatticeSimulation.hpp"

#include <cmath>
#include <boost/make_shared.hpp>

#include "VertexBasedCellPopulation.hpp"
#include "T2SwapCellKiller.hpp"
#include "Cylindrical2dMesh.hpp"
#include "Cylindrical2dVertexMesh.hpp"
#include "AbstractTwoBodyInteractionForce.hpp"
#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "StepSizeException.hpp"
#include "SmartPointers.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OffLatticeSimulation(AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells)
    : AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells)
{
    if (!dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OffLatticeSimulations require a subclass of AbstractOffLatticeCellPopulation.");
    }

    if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation)))
    {
        // For VertexBasedCellPopulations we automatically add a T2SwapCellKiller. In order to inhibit T2 swaps
        // the user needs to set the threshold for T2 swaps in the mesh to 0.
        VertexBasedCellPopulation<SPACE_DIM>* p_vertex_based_cell_population = dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&rCellPopulation);
        MAKE_PTR_ARGS(T2SwapCellKiller<SPACE_DIM>, p_t2_swap_cell_killer, (p_vertex_based_cell_population));
        this->AddCellKiller(p_t2_swap_cell_killer);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddForce(boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > pForce)
{
    mForceCollection.push_back(pForce);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllForces()
{
    mForceCollection.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::AddCellPopulationBoundaryCondition(boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > pBoundaryCondition)
{
    mBoundaryConditions.push_back(pBoundaryCondition);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RemoveAllCellPopulationBoundaryConditions()
{
    mBoundaryConditions.clear();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetNumericalMethod(boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > pNumericalMethod)
{
    mpNumericalMethod = pNumericalMethod;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const boost::shared_ptr<AbstractNumericalMethod<ELEMENT_DIM, SPACE_DIM> > OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::GetNumericalMethod() const
{
    return mpNumericalMethod;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::UpdateCellLocationsAndTopology()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);

    double time_advanced_so_far = 0;
    double target_time_step  = this->mDt;
    double present_time_step = this->mDt;

    while (time_advanced_so_far < target_time_step)
    {
        // Store the initial node positions (these may be needed when applying boundary conditions)    
        std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > old_node_locations;

        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
            node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
            ++node_iter)
        {
            old_node_locations[&(*node_iter)] = (node_iter)->rGetLocation();
        }

        // Try to update node positions according to the numerical method 
        try
        {
            mpNumericalMethod->UpdateAllNodePositions(present_time_step);
            ApplyBoundaries(old_node_locations);

            // Successful time step! Update time_advanced_so_far
            time_advanced_so_far += present_time_step;

            // If using adaptive timestep, then increase the present_time_step (by 1% for now)
            if (mpNumericalMethod->HasAdaptiveTimestep())
            {
                ///\todo #2087 Make this a settable member variable
                double timestep_increase = 0.01;
                present_time_step = std::min((1+timestep_increase)*present_time_step, target_time_step - time_advanced_so_far);
            }

        }
        catch (StepSizeException& e)
        {
            // Detects if a node has travelled too far in a single time step
            if (mpNumericalMethod->HasAdaptiveTimestep())
            {
                // If adaptivity is switched on, revert node locations and choose a suitably smaller time step
                RevertToOldLocations(old_node_locations);
                present_time_step = std::min(e.GetSuggestedNewStep(), target_time_step - time_advanced_so_far);
            }
            else
            {
                // If adaptivity is switched off, terminate with an error
                EXCEPTION(e.what());
            }
        }
    }

    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::RevertToOldLocations(std::map<Node<SPACE_DIM>*, c_vector<double, SPACE_DIM> > oldNodeLoctions)
{
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
        node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
        ++node_iter)
    {
        (node_iter)->rGetModifiableLocation() = oldNodeLoctions[&(*node_iter)];
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::ApplyBoundaries(std::map<Node<SPACE_DIM>*,c_vector<double, SPACE_DIM> > oldNodeLoctions)
{
    // Apply any boundary conditions
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        (*bcs_iter)->ImposeBoundaryCondition(oldNodeLoctions);
    }

    // Verify that each boundary condition is now satisfied
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator bcs_iter = mBoundaryConditions.begin();
         bcs_iter != mBoundaryConditions.end();
         ++bcs_iter)
    {
        if (!((*bcs_iter)->VerifyBoundaryCondition()))
        {
            EXCEPTION("The cell population boundary conditions are incompatible.");
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::CalculateCellDivisionVector(CellPtr pParentCell)
{
    ///\todo #2800 refactor division rules to follow more closely the vertex case

    // This static cast is always valid as the constructor tests that the cell population is AbstractOffLattice
    return static_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))->CalculateCellDivisionVector(pParentCell);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::WriteVisualizerSetupFile()
{
    if (PetscTools::AmMaster())
    {
        for (unsigned i=0; i<this->mForceCollection.size(); i++)
        {
            // This may cause compilation problems, probably due to AbstractTwoBodyInteractionForce not having two template parameters
            ///\todo Check whether this comment is still valid

            boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > p_force = this->mForceCollection[i];
            if (boost::dynamic_pointer_cast<AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM> >(p_force))
            {
                double cutoff = (boost::static_pointer_cast<AbstractTwoBodyInteractionForce<ELEMENT_DIM,SPACE_DIM> >(p_force))->GetCutOffLength();
                *(this->mpVizSetupFile) << "Cutoff\t" << cutoff << "\n";
            }
        }

        // This is a quick and dirty check to see if the mesh is periodic
        if (bool(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&this->mrCellPopulation)))
        {
           if (bool(dynamic_cast<Cylindrical2dMesh*>(&(dynamic_cast<MeshBasedCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation))->rGetMesh()))))
           {
               *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
           }
        }
        else if (bool(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&this->mrCellPopulation)))
        {
           if (bool(dynamic_cast<Cylindrical2dVertexMesh*>(&(dynamic_cast<VertexBasedCellPopulation<SPACE_DIM>*>(&(this->mrCellPopulation))->rGetMesh()))))
           {
               *this->mpVizSetupFile << "MeshWidth\t" << this->mrCellPopulation.GetWidth(0) << "\n";
           }
        }
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::SetupSolve()
{
    // Clear all forces
    for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = this->mrCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != this->mrCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        node_iter->ClearAppliedForce();
    }

    // Use a forward Euler method by default, unless a numerical method has been specified already
    if (mpNumericalMethod == NULL)
    {
        mpNumericalMethod = boost::make_shared<ForwardEulerNumericalMethod<ELEMENT_DIM, SPACE_DIM> >();
    }
    mpNumericalMethod->SetCellPopulation(dynamic_cast<AbstractOffLatticeCellPopulation<ELEMENT_DIM,SPACE_DIM>*>(&(this->mrCellPopulation)));
    mpNumericalMethod->SetForceCollection(&mForceCollection);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over forces
    *rParamsFile << "\n\t<Forces>\n";
    for (typename std::vector<boost::shared_ptr<AbstractForce<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mForceCollection.begin();
         iter != mForceCollection.end();
         ++iter)
    {
        // Output force details
        (*iter)->OutputForceInfo(rParamsFile);
    }
    *rParamsFile << "\t</Forces>\n";

    // Loop over cell population boundary conditions
    *rParamsFile << "\n\t<CellPopulationBoundaryConditions>\n";
    for (typename std::vector<boost::shared_ptr<AbstractCellPopulationBoundaryCondition<ELEMENT_DIM,SPACE_DIM> > >::iterator iter = mBoundaryConditions.begin();
         iter != mBoundaryConditions.end();
         ++iter)
    {
        // Output cell Boundary condition details
        (*iter)->OutputCellPopulationBoundaryConditionInfo(rParamsFile);
    }
    *rParamsFile << "\t</CellPopulationBoundaryConditions>\n";

    // Output numerical method details
    *rParamsFile << "\n\t<NumericalMethod>\n";
    mpNumericalMethod->OutputNumericalMethodInfo(rParamsFile);
    *rParamsFile << "\t</NumericalMethod>\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void OffLatticeSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    // No new parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulation<ELEMENT_DIM,SPACE_DIM>::OutputSimulationParameters(rParamsFile);
}

// Explicit instantiation
template class OffLatticeSimulation<1,1>;
template class OffLatticeSimulation<1,2>;
template class OffLatticeSimulation<2,2>;
template class OffLatticeSimulation<1,3>;
template class OffLatticeSimulation<2,3>;
template class OffLatticeSimulation<3,3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(OffLatticeSimulation)
