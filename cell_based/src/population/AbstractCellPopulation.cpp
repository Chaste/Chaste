/*

Copyright (c) 2005-2013, University of Oxford.
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

#include "AbstractCellPopulation.hpp"
#include "AbstractOdeBasedCellCycleModel.hpp"
#include "Exception.hpp"
#include "SmartPointers.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AbstractCellPopulation( AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                                    std::vector<CellPtr>& rCells,
                                    const std::vector<unsigned> locationIndices)
    : mrMesh(rMesh),
      mCells(rCells.begin(), rCells.end()),
      mCentroid(zero_vector<double>(SPACE_DIM)),
      mpCellPropertyRegistry(CellPropertyRegistry::Instance()->TakeOwnership()),
      mOutputResultsForChasteVisualizer(true),
      mOutputCellIdData(false),
      mOutputCellMutationStates(false),
      mOutputCellAncestors(false),
      mOutputCellProliferativeTypes(false),
      mOutputCellVariables(false),
      mOutputCellCyclePhases(false),
      mOutputCellAges(false),
      mOutputCellVolumes(false)
{
    // To avoid double-counting problems, clear the passed-in cells vector
    rCells.clear();

    if (!locationIndices.empty())
    {
        // There must be a one-one correspondence between cells and location indices
        if (mCells.size() != locationIndices.size())
        {
            EXCEPTION("There is not a one-one correspondence between cells and location indices");
        }
    }

    // Set up the map between location indices and cells
    mLocationCellMap.clear();
    mCellLocationMap.clear();

    std::list<CellPtr>::iterator it = mCells.begin();
    for (unsigned i=0; it != mCells.end(); ++it, ++i)
    {
        // Give each cell a pointer to the property registry (we have taken ownership in this constructor).
        (*it)->rGetCellPropertyCollection().SetCellPropertyRegistry(mpCellPropertyRegistry.get());
    }

    /*
     * Initialise cell counts to zero.
     *
     * Note: In its current form the code requires each cell-cycle model
     * to comprise four phases (G1, S, G2, M). This is reflected in the
     * explicit use of the variable NUM_CELL_CYCLE_PHASES below.
     */
    mCellCyclePhaseCount = std::vector<unsigned>(NUM_CELL_CYCLE_PHASES);
    for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AbstractCellPopulation(AbstractMesh<ELEMENT_DIM, SPACE_DIM>& rMesh)
    : mrMesh(rMesh)
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::~AbstractCellPopulation()
{
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::InitialiseCells()
{
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->InitialiseCellCycleModel();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetDataOnAllCells(const std::string& dataName, double dataValue)
{
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->GetCellData()->SetItem(dataName, dataValue);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
AbstractMesh<ELEMENT_DIM, SPACE_DIM>& AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::list<CellPtr>& AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::rGetCells()
{
    return mCells;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetCellAncestorsToLocationIndices()
{
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (mCellLocationMap[(*cell_iter).get()]));
        cell_iter->SetAncestor(p_cell_ancestor);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<unsigned> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellAncestors()
{
    std::set<unsigned> remaining_ancestors;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        remaining_ancestors.insert(cell_iter->GetAncestor());
    }
    return remaining_ancestors;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellMutationStateCount()
{
    if (!mOutputCellMutationStates)
    {
        EXCEPTION("Call SetOutputCellMutationStates(true) before using this function");
    }

    // An ordering must be specified for cell mutation states and cell proliferative types
    SetDefaultCellMutationStateAndProliferativeTypeOrdering();

    const std::vector<boost::shared_ptr<AbstractCellProperty> >& r_cell_properties =
        mpCellPropertyRegistry->rGetAllCellProperties();

    std::vector<unsigned> cell_mutation_state_count;
    for (unsigned i=0; i<r_cell_properties.size(); i++)
    {
        if (r_cell_properties[i]->IsSubType<AbstractCellMutationState>())
        {
            cell_mutation_state_count.push_back(r_cell_properties[i]->GetCellCount());
        }
    }

    return cell_mutation_state_count;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellProliferativeTypeCount()
{
    if (!mOutputCellProliferativeTypes)
    {
        EXCEPTION("Call SetOutputCellProliferativeTypes(true) before using this function");
    }

    // An ordering must be specified for cell mutation states and cell proliferative types
    SetDefaultCellMutationStateAndProliferativeTypeOrdering();

    const std::vector<boost::shared_ptr<AbstractCellProperty> >& r_cell_properties =
        mpCellPropertyRegistry->rGetAllCellProperties();

    std::vector<unsigned> cell_proliferative_type_count;
    for (unsigned i=0; i<r_cell_properties.size(); i++)
    {
        if (r_cell_properties[i]->IsSubType<AbstractCellProliferativeType>())
        {
            cell_proliferative_type_count.push_back(r_cell_properties[i]->GetCellCount());
        }
    }

    return cell_proliferative_type_count;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
const std::vector<unsigned>& AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::rGetCellCyclePhaseCount() const
{
    if (!mOutputCellCyclePhases)
    {
        EXCEPTION("Call SetOutputCellCyclePhases(true) before using this function");
    }
    return mCellCyclePhaseCount;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
CellPtr AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellUsingLocationIndex(unsigned index)
{
    // Get the set of pointers to cells corresponding to this location index
    std::set<CellPtr> cells = mLocationCellMap[index];

    // If there is only one cell attached return the cell. Note currently only one cell per index.
    if (cells.size()==1)
    {
        return *(cells.begin());
    }
    if (cells.size()==0)
    {
        EXCEPTION("Location index input argument does not correspond to a Cell");
    }
    else
    {
        EXCEPTION("Multiple cells are attached to a single location index.");
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::set<CellPtr> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellsUsingLocationIndex(unsigned index)
{
    // Return the set of pointers to cells corresponding to this location index, note the set may be empty.
    return mLocationCellMap[index];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::IsCellAttachedToLocationIndex(unsigned index)
{
    // Get the set of pointers to cells corresponding to this location index
    std::set<CellPtr> cells = mLocationCellMap[index];

    // check if there is a cell attached to the location index.
    if (cells.size()==0)
    {
        return false;
    }
    else
    {
        return true;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    // Clear the maps
    mLocationCellMap[index].clear();
       mCellLocationMap.erase(pCell.get());

    // Replace with new cell
    mLocationCellMap[index].insert(pCell);

    // Do other half of the map
    mCellLocationMap[pCell.get()] = index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    mLocationCellMap[index].insert(pCell);
    mCellLocationMap[pCell.get()] = index;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::RemoveCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    std::set<CellPtr>::iterator cell_iter = mLocationCellMap[index].find(pCell);

    if (cell_iter == mLocationCellMap[index].end())
    {
        EXCEPTION("Tried to remove a cell which is not attached to the given location index");
    }
    else
    {
        mLocationCellMap[index].erase(cell_iter);
        mCellLocationMap.erase(pCell.get());
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::MoveCellInLocationMap(CellPtr pCell, unsigned old_index, unsigned new_index)
{
    // Remove the cell from its current location
    RemoveCellUsingLocationIndex(old_index, pCell);

    // Add it to the new location
    AddCellUsingLocationIndex(new_index, pCell);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
unsigned AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetLocationIndexUsingCell(CellPtr pCell)
{
    // Check the cell is in the map
    assert(this->mCellLocationMap.find(pCell.get()) != this->mCellLocationMap.end());

    return mCellLocationMap[pCell.get()];
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
boost::shared_ptr<CellPropertyRegistry> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellPropertyRegistry()
{
    return mpCellPropertyRegistry;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetDefaultCellMutationStateAndProliferativeTypeOrdering()
{
    boost::shared_ptr<CellPropertyRegistry> p_registry = GetCellPropertyRegistry();
    if (!p_registry->HasOrderingBeenSpecified())
    {
        std::vector<boost::shared_ptr<AbstractCellProperty> > mutations_and_proliferative_types;
        mutations_and_proliferative_types.push_back(p_registry->Get<WildTypeCellMutationState>());
        mutations_and_proliferative_types.push_back(p_registry->Get<ApcOneHitCellMutationState>());
        mutations_and_proliferative_types.push_back(p_registry->Get<ApcTwoHitCellMutationState>());
        mutations_and_proliferative_types.push_back(p_registry->Get<BetaCateninOneHitCellMutationState>());
        mutations_and_proliferative_types.push_back(p_registry->Get<StemCellProliferativeType>());
        mutations_and_proliferative_types.push_back(p_registry->Get<TransitCellProliferativeType>());
        mutations_and_proliferative_types.push_back(p_registry->Get<DifferentiatedCellProliferativeType>());
        p_registry->SpecifyOrdering(mutations_and_proliferative_types);
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCentroidOfCellPopulation()
{
    mCentroid = zero_vector<double>(SPACE_DIM);
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        mCentroid += GetLocationOfCellCentre(*cell_iter);
    }
    mCentroid /= this->GetNumRealCells();

    return mCentroid;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::UpdateCellProcessLocation()
{
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);

    if (PetscTools::AmMaster())
    {
        if (mOutputResultsForChasteVisualizer)
        {
            mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes");
            mpVizBoundaryNodesFile = output_file_handler.OpenOutputFile("results.vizboundarynodes");
            mpVizCellProliferativeTypesFile = output_file_handler.OpenOutputFile("results.vizcelltypes");
        }
        if (mOutputCellAncestors)
        {
            mpVizCellAncestorsFile = output_file_handler.OpenOutputFile("results.vizancestors");
        }
        if (mOutputCellMutationStates)
        {
            // An ordering must be specified for cell mutation states and cell proliferative types
            SetDefaultCellMutationStateAndProliferativeTypeOrdering();

            mpCellMutationStatesFile = output_file_handler.OpenOutputFile("cellmutationstates.dat");

            *mpCellMutationStatesFile << "Time\t ";

            const std::vector<boost::shared_ptr<AbstractCellProperty> >& r_cell_properties =
                mpCellPropertyRegistry->rGetAllCellProperties();

            std::vector<unsigned> cell_mutation_state_count;
            for (unsigned i=0; i<r_cell_properties.size(); i++)
            {
                if (r_cell_properties[i]->IsSubType<AbstractCellMutationState>())
                {
                    *mpCellMutationStatesFile << r_cell_properties[i]->GetIdentifier() << "\t ";
                }
            }
            *mpCellMutationStatesFile << "\n";
        }
        if (mOutputCellProliferativeTypes)
        {
            mpCellProliferativeTypesFile = output_file_handler.OpenOutputFile("celltypes.dat");
        }
        if (mOutputCellVariables)
        {
            mpCellVariablesFile = output_file_handler.OpenOutputFile("cellvariables.dat");
        }
        if (mOutputCellCyclePhases)
        {
            mpCellCyclePhasesFile = output_file_handler.OpenOutputFile("cellcyclephases.dat");
            mpVizCellProliferativePhasesFile = output_file_handler.OpenOutputFile("results.vizcellphases");
        }
        if (mOutputCellAges)
        {
            mpCellAgesFile = output_file_handler.OpenOutputFile("cellages.dat");
        }
        if (mOutputCellIdData)
        {
            mpCellIdFile = output_file_handler.OpenOutputFile("loggedcell.dat");
        }
        if (this->mOutputCellVolumes)
        {
            mpCellVolumesFile = output_file_handler.OpenOutputFile("cellareas.dat");
        }
    }

    mDirPath = rDirectory;
#ifdef CHASTE_VTK
    mpVtkMetaFile = output_file_handler.OpenOutputFile("results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::CloseOutputFiles()
{
    // In parallel all files are closed after writing
    if (PetscTools::IsSequential())
    {
        if (mOutputResultsForChasteVisualizer)
        {
            mpVizNodesFile->close();
            mpVizBoundaryNodesFile->close();
            mpVizCellProliferativeTypesFile->close();
        }
        if (mOutputCellMutationStates)
        {
            mpCellMutationStatesFile->close();
        }
        if (mOutputCellProliferativeTypes)
        {
            mpCellProliferativeTypesFile->close();
        }
        if (mOutputCellVariables)
        {
            mpCellVariablesFile->close();
        }
        if (mOutputCellCyclePhases)
        {
            mpCellCyclePhasesFile->close();
            mpVizCellProliferativePhasesFile->close();
        }
        if (mOutputCellAncestors)
        {
            mpVizCellAncestorsFile->close();
        }
        if (mOutputCellAges)
        {
            mpCellAgesFile->close();
        }
        if (mOutputCellIdData)
        {
            mpCellIdFile->close();
        }
        if (this->mOutputCellVolumes)
        {
            mpCellVolumesFile->close();
        }
    }
#ifdef CHASTE_VTK
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::ResetCellCounters()
{
    for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GenerateCellResults(CellPtr pCell)
{
    unsigned location_index = this->GetLocationIndexUsingCell(pCell);

    unsigned colour = STEM_COLOUR;

    if (mOutputCellCyclePhases)
    {
        // Update mCellCyclePhaseCount
        switch (pCell->GetCellCycleModel()->GetCurrentCellCyclePhase())
        {
            case G_ZERO_PHASE:
                mCellCyclePhaseCount[0]++;
                break;
            case G_ONE_PHASE:
                mCellCyclePhaseCount[1]++;
                break;
            case S_PHASE:
                mCellCyclePhaseCount[2]++;
                break;
            case G_TWO_PHASE:
                mCellCyclePhaseCount[3]++;
                break;
             case M_PHASE:
                 mCellCyclePhaseCount[4]++;
                break;
            default:
                NEVER_REACHED;
        }
        *mpVizCellProliferativePhasesFile << pCell->GetCellCycleModel()->GetCurrentCellCyclePhase() << " ";
    }

    if (mOutputCellAncestors)
    {
        // Set colour dependent on cell ancestor and write to file
        colour = pCell->GetAncestor();
        if (colour == UNSIGNED_UNSET)
        {
            // Set the file to -1 to mark this case.
            colour = 1;
            *mpVizCellAncestorsFile << "-";
        }
        *mpVizCellAncestorsFile << colour << " ";
    }

    // Set colour dependent on cell proliferative type
    colour = pCell->GetCellProliferativeType()->GetColour();

    if (mOutputCellMutationStates)
    {
        // Set colour dependent on cell mutation state
        if (!pCell->GetMutationState()->IsType<WildTypeCellMutationState>())
        {
            colour = pCell->GetMutationState()->GetColour();
        }
        if (pCell->HasCellProperty<CellLabel>())
        {
            CellPropertyCollection collection = pCell->rGetCellPropertyCollection().GetProperties<CellLabel>();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
            colour = p_label->GetColour();
        }
    }

    if (pCell->HasCellProperty<ApoptoticCellProperty>() || pCell->HasApoptosisBegun())
    {
        // For any type of cell set the colour to this if it is undergoing apoptosis
        colour = APOPTOSIS_COLOUR;
    }

    // Write cell variable data to file if required
    if (mOutputCellVariables && dynamic_cast<AbstractOdeBasedCellCycleModel*>(pCell->GetCellCycleModel()) )
    {
        // Write location index corresponding to cell
        *mpCellVariablesFile << location_index << " ";

        // Write cell variables
        std::vector<double> proteins = (static_cast<AbstractOdeBasedCellCycleModel*>(pCell->GetCellCycleModel()))->GetProteinConcentrations();
        for (unsigned i=0; i<proteins.size(); i++)
        {
            *mpCellVariablesFile << proteins[i] << " ";
        }
    }

    // Write cell age data to file if required
    if (mOutputCellAges)
    {
        // Write location index corresponding to cell
        *mpCellAgesFile << location_index << " ";

        // Write cell location
        c_vector<double, SPACE_DIM> cell_location = GetLocationOfCellCentre(pCell);

        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *mpCellAgesFile << cell_location[i] << " ";
        }

        // Write cell age
        *mpCellAgesFile << pCell->GetAge() << " ";
    }
    if (mOutputResultsForChasteVisualizer)
    {
        *mpVizCellProliferativeTypesFile << colour << " ";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::WriteCellResultsToFiles()
{
    if (mOutputResultsForChasteVisualizer)
    {
        *mpVizCellProliferativeTypesFile << "\n";
    }
    if (mOutputCellAncestors)
    {
        *mpVizCellAncestorsFile << "\n";
    }

    // Write cell mutation state data to file if required
    if (mOutputCellMutationStates)
    {
        std::vector<unsigned> mutation_state_count = GetCellMutationStateCount();

        for (unsigned i=0; i<mutation_state_count.size(); i++)
        {
            *mpCellMutationStatesFile << mutation_state_count[i] << "\t";
        }
        *mpCellMutationStatesFile << "\n";
    }

    // Write cell proliferative type data to file if required
    if (mOutputCellProliferativeTypes)
    {
        std::vector<unsigned> proliferative_type_count = GetCellProliferativeTypeCount();

        for (unsigned i=0; i<proliferative_type_count.size(); i++)
        {
            *mpCellProliferativeTypesFile << proliferative_type_count[i] << "\t";
        }
        *mpCellProliferativeTypesFile << "\n";
    }

    if (mOutputCellVariables)
    {
        *mpCellVariablesFile << "\n";
    }

    // Write cell cycle phase data to file if required
    if (mOutputCellCyclePhases)
    {
        for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
        {
            *mpCellCyclePhasesFile << mCellCyclePhaseCount[i] << "\t";
        }
        *mpCellCyclePhasesFile << "\n";

        // The data for this is output in GenerateCellResults()
        *mpVizCellProliferativePhasesFile << "\n";
    }

    // Write cell age data to file if required
    if (mOutputCellAges)
    {
        *mpCellAgesFile << "\n";
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::WriteTimeAndNodeResultsToFiles()
{
    OutputFileHandler output_file_handler(mDirPath, false);

    PetscTools::BeginRoundRobin();
    {
        if (!PetscTools::AmMaster() || SimulationTime::Instance()->IsEndTimeAndNumberOfTimeStepsSetUp())
        {
            mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes", std::ios::app);
            mpVizBoundaryNodesFile = output_file_handler.OpenOutputFile("results.vizboundarynodes", std::ios::app);
        }
        if (PetscTools::AmMaster())
        {
            double time = SimulationTime::Instance()->GetTime();

            *mpVizNodesFile << time << "\t";
            *mpVizBoundaryNodesFile << time << "\t";
        }
        // Write node data to file
        for (typename AbstractMesh<ELEMENT_DIM, SPACE_DIM>::NodeIterator node_iter = mrMesh.GetNodeIteratorBegin();
                node_iter != mrMesh.GetNodeIteratorEnd();
                ++node_iter)
        {
            if (!node_iter->IsDeleted())
            {
                const c_vector<double,SPACE_DIM>& position = node_iter->rGetLocation();

                for (unsigned i=0; i<SPACE_DIM; i++)
                {
                    *mpVizNodesFile << position[i] << " ";
                }
                *mpVizBoundaryNodesFile << node_iter->IsBoundaryNode() << " ";
            }
        }
        if (PetscTools::AmTopMost())
        {
            *mpVizNodesFile << "\n";
            *mpVizBoundaryNodesFile << "\n";
        }

        mpVizNodesFile->close();
        mpVizBoundaryNodesFile->close();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::WriteResultsToFiles()
{
    if (mOutputResultsForChasteVisualizer)
    {
        WriteTimeAndNodeResultsToFiles();
    }

    double time = SimulationTime::Instance()->GetTime();

    if (mOutputResultsForChasteVisualizer)
    {
        *mpVizCellProliferativeTypesFile << time << "\t";
    }
    if (mOutputCellAncestors)
    {
        *mpVizCellAncestorsFile << time << "\t";
    }
    if (mOutputCellMutationStates)
    {
        *mpCellMutationStatesFile << time << "\t";
    }
    if (mOutputCellProliferativeTypes)
    {
        *mpCellProliferativeTypesFile << time << "\t";
    }
    if (mOutputCellVariables)
    {
        *mpCellVariablesFile << time << "\t";
    }
    if (mOutputCellCyclePhases)
    {
        *mpCellCyclePhasesFile << time << "\t";
        *mpVizCellProliferativePhasesFile << time << "\t";
    }
    if (mOutputCellAges)
    {
        *mpCellAgesFile << time << "\t";
    }
    if (this->mOutputCellVolumes)
    {
        WriteCellVolumeResultsToFile();
    }

    ResetCellCounters();

    GenerateCellResultsAndWriteToFiles();

    // Write logged cell data if required
    if (mOutputCellIdData)
    {
        WriteCellIdDataToFile();
    }

    // VTK can only be written in 2 or 3 dimensions
    if (SPACE_DIM > 1)
    {
         WriteVtkResultsToFile();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::WriteCellIdDataToFile()
{
    // Write time to file
    *mpCellIdFile << SimulationTime::Instance()->GetTime();

    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = Begin();
         cell_iter != End();
         ++cell_iter)
    {
        unsigned cell_id = cell_iter->GetCellId();
        unsigned location_index = mCellLocationMap[(*cell_iter).get()];
        *mpCellIdFile << " " << cell_id << " " << location_index;

        c_vector<double, SPACE_DIM> coords = GetLocationOfCellCentre(*cell_iter);
        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            *mpCellIdFile << " " << coords[i];
        }
    }
    *mpCellIdFile << "\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationInfo(out_stream& rParamsFile)
{
    std::string cell_population_type = GetIdentifier();

    *rParamsFile << "\t<" << cell_population_type << ">\n";
    OutputCellPopulationParameters(rParamsFile);
    *rParamsFile << "\t</" << cell_population_type << ">\n";
    *rParamsFile << "\n";
    *rParamsFile << "\t<CellCycleModels>\n";

    /**
     * Loop over cells and generate a set of cell-cycle model classes
     * that are present in the population.
     *
     * \todo this currently ignores different parameter regimes (#1453)
     */
    std::set<std::string> unique_cell_cycle_models;
    std::vector<CellPtr> first_cell_with_unique_CCM;
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        std::string identifier = cell_iter->GetCellCycleModel()->GetIdentifier();
        if (unique_cell_cycle_models.count(identifier) == 0)
        {
            unique_cell_cycle_models.insert(identifier);
            first_cell_with_unique_CCM.push_back((*cell_iter));
        }
    }

    // Loop over unique cell-cycle models
    for (unsigned i=0; i<first_cell_with_unique_CCM.size(); i++)
    {
        // Output cell-cycle model details
        first_cell_with_unique_CCM[i]->GetCellCycleModel()->OutputCellCycleModelInfo(rParamsFile);
    }

    *rParamsFile << "\t</CellCycleModels>\n";
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputResultsForChasteVisualizer>" << mOutputResultsForChasteVisualizer << "</OutputResultsForChasteVisualizer>\n";
    *rParamsFile << "\t\t<OutputCellIdData>" << mOutputCellIdData << "</OutputCellIdData>\n";
    *rParamsFile << "\t\t<OutputCellMutationStates>" << mOutputCellMutationStates << "</OutputCellMutationStates>\n";
    *rParamsFile << "\t\t<OutputCellAncestors>" << mOutputCellAncestors << "</OutputCellAncestors>\n";
    *rParamsFile << "\t\t<OutputCellProliferativeTypes>" << mOutputCellProliferativeTypes << "</OutputCellProliferativeTypes>\n";
    *rParamsFile << "\t\t<OutputCellVariables>" << mOutputCellVariables << "</OutputCellVariables>\n";
    *rParamsFile << "\t\t<OutputCellCyclePhases>" << mOutputCellCyclePhases << "</OutputCellCyclePhases>\n";
    *rParamsFile << "\t\t<OutputCellAges>" << mOutputCellAges << "</OutputCellAges>\n";
    *rParamsFile << "\t\t<OutputCellVolumes>" << mOutputCellVolumes << "</OutputCellVolumes>\n";
}

///////////////////////////////////////////////////////////////////////
// Getter methods
///////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputResultsForChasteVisualizer()
{
    return mOutputResultsForChasteVisualizer;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellIdData()
{
    return mOutputCellIdData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellMutationStates()
{
    return mOutputCellMutationStates;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellAncestors()
{
    return mOutputCellAncestors;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellProliferativeTypes()
{
    return mOutputCellProliferativeTypes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellVariables()
{
    return mOutputCellVariables;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellCyclePhases()
{
    return mOutputCellCyclePhases;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellAges()
{
    return mOutputCellAges;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetOutputCellVolumes()
{
    return mOutputCellVolumes;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddPopulationWriter(AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>* pWriter)
{
    //\ todo Check that we don't already have a writer of type pWriter.
    mCellPopulationWriters.push_back(pWriter);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::AddCellWriter(AbstractCellWriter<ELEMENT_DIM, SPACE_DIM>* pWriter)
{
    //\ todo Check that we don't already have a writer of type pWriter.
    mCellWriters.push_back(pWriter);
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputResultsForChasteVisualizer(bool outputResultsForChasteVisualizer)
{
    mOutputResultsForChasteVisualizer = outputResultsForChasteVisualizer;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellIdData(bool writeCellIdData)
{
    mOutputCellIdData = writeCellIdData;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellMutationStates(bool outputCellMutationStates)
{
    mOutputCellMutationStates = outputCellMutationStates;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellAncestors(bool outputCellAncestors)
{
    mOutputCellAncestors = outputCellAncestors;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellProliferativeTypes(bool outputCellProliferativeTypes)
{
    mOutputCellProliferativeTypes = outputCellProliferativeTypes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellVariables(bool outputCellVariables)
{
    mOutputCellVariables = outputCellVariables;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellCyclePhases(bool outputCellCyclePhases)
{
    mOutputCellCyclePhases = outputCellCyclePhases;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellAges(bool outputCellAges)
{
    mOutputCellAges = outputCellAges;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::SetOutputCellVolumes(bool outputCellVolumes)
{
    mOutputCellVolumes = outputCellVolumes;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
bool AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::IsRoomToDivide(CellPtr pCell)
{
    return true;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,SPACE_DIM> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetSizeOfCellPopulation()
{
    // Compute the centre of mass of the cell population
    c_vector<double,SPACE_DIM> centre = GetCentroidOfCellPopulation();

    // Loop over cells and find the maximum distance from the centre of mass in each dimension
    c_vector<double,SPACE_DIM> max_distance_from_centre = zero_vector<double>(SPACE_DIM);
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        c_vector<double,SPACE_DIM> cell_location = GetLocationOfCellCentre(*cell_iter);
        c_vector<double,SPACE_DIM> displacement;
        displacement = centre - cell_location;

        for (unsigned i=0; i<SPACE_DIM; i++)
        {
            if (displacement[i] > max_distance_from_centre[i])
            {
                max_distance_from_centre[i] = displacement[i];
            }
        }
    }

    return max_distance_from_centre;
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class AbstractCellPopulation<1,1>;
template class AbstractCellPopulation<1,2>;
template class AbstractCellPopulation<2,2>;
template class AbstractCellPopulation<1,3>;
template class AbstractCellPopulation<2,3>;
template class AbstractCellPopulation<3,3>;
