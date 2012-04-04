/*

Copyright (c) 2005-2012, University of Oxford.
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

template<unsigned DIM>
AbstractCellPopulation<DIM>::AbstractCellPopulation( AbstractMesh<DIM, DIM>& rMesh,
									std::vector<CellPtr>& rCells,
                                    const std::vector<unsigned> locationIndices)
    : mrMesh(rMesh),
      mCells(rCells.begin(), rCells.end()),
      mCentroid(zero_vector<double>(DIM)),
      mpCellPropertyRegistry(CellPropertyRegistry::Instance()->TakeOwnership()),
      mOutputCellIdData(false),
      mOutputCellMutationStates(false),
      mOutputCellAncestors(false),
      mOutputCellProliferativeTypes(false),
      mOutputCellVariables(false),
      mOutputCellCyclePhases(false),
      mOutputCellAges(false),
      mOutputCellVolumes(false)
{
    // There must be at least one cell
    assert(!mCells.empty());

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
    std::list<CellPtr>::iterator it = mCells.begin();
    for (unsigned i=0; it != mCells.end(); ++it, ++i)
    {
        unsigned index = locationIndices.empty() ? i : locationIndices[i]; // assume that the ordering matches
        SetCellUsingLocationIndex(index,*it);
        mCellLocationMap[(*it).get()] = index;

        // Give each cell a pointer to the property registry (we have taken ownership in this constructor).
        (*it)->rGetCellPropertyCollection().SetCellPropertyRegistry(mpCellPropertyRegistry.get());
    }

    /*
     * Initialise cell counts to zero.
     *
     * Note: In its current form the code requires each cell-cycle model
     * to comprise four phases (G1, S, G2, M) and requires a cell to have
     * one of three possible proliferative types (STEM, TRANSIT and
     * DIFFERENTIATED). This is reflected in the explicit use of the
     * variables NUM_CELL_PROLIFERATIVE_TYPES and NUM_CELL_CYCLE_PHASES
     * below.
     */
    mCellProliferativeTypeCount = std::vector<unsigned>(NUM_CELL_PROLIFERATIVE_TYPES);
    for (unsigned i=0; i<mCellProliferativeTypeCount.size(); i++)
    {
        mCellProliferativeTypeCount[i] = 0;
    }

    mCellCyclePhaseCount = std::vector<unsigned>(NUM_CELL_CYCLE_PHASES);
    for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
}

template<unsigned DIM>
AbstractCellPopulation<DIM>::AbstractCellPopulation(AbstractMesh<DIM, DIM>& rMesh)
			: mrMesh(rMesh)
{
}

template<unsigned DIM>
AbstractCellPopulation<DIM>::~AbstractCellPopulation()
{
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::InitialiseCells()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        cell_iter->InitialiseCellCycleModel();
    }
}

template<unsigned DIM>
AbstractMesh<DIM, DIM>& AbstractCellPopulation<DIM>::rGetMesh()
{
    return mrMesh;
}

template<unsigned DIM>
std::list<CellPtr>& AbstractCellPopulation<DIM>::rGetCells()
{
    return mCells;
}

template<unsigned DIM>
unsigned AbstractCellPopulation<DIM>::GetNumRealCells()
{
    unsigned counter = 0;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        counter++;
    }
    return counter;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetCellAncestorsToLocationIndices()
{
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
    	MAKE_PTR_ARGS(CellAncestor, p_cell_ancestor, (mCellLocationMap[(*cell_iter).get()]));
        cell_iter->SetAncestor(p_cell_ancestor);
    }
}

template<unsigned DIM>
std::set<unsigned> AbstractCellPopulation<DIM>::GetCellAncestors()
{
    std::set<unsigned> remaining_ancestors;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter=this->Begin(); cell_iter!=this->End(); ++cell_iter)
    {
        remaining_ancestors.insert(cell_iter->GetAncestor());
    }
    return remaining_ancestors;
}

template<unsigned DIM>
std::vector<unsigned> AbstractCellPopulation<DIM>::GetCellMutationStateCount()
{
    if (!mOutputCellMutationStates)
    {
        EXCEPTION("Call SetOutputCellMutationStates(true) before using this function");
    }

    // An ordering must have been specified for cell mutation states
    SetDefaultMutationStateOrdering();

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

template<unsigned DIM>
const std::vector<unsigned>& AbstractCellPopulation<DIM>::rGetCellProliferativeTypeCount() const
{
    if (!mOutputCellProliferativeTypes)
    {
        EXCEPTION("Call SetOutputCellProliferativeTypes(true) before using this function");
    }
    return mCellProliferativeTypeCount;
}

template<unsigned DIM>
const std::vector<unsigned>& AbstractCellPopulation<DIM>::rGetCellCyclePhaseCount() const
{
    if (!mOutputCellCyclePhases)
    {
        EXCEPTION("Call SetOutputCellCyclePhases(true) before using this function");
    }
    return mCellCyclePhaseCount;
}

template<unsigned DIM>
CellPtr AbstractCellPopulation<DIM>::GetCellUsingLocationIndex(unsigned index)
{
    // Get a pointer to the cell corresponding to this location index
    CellPtr p_cell = mLocationCellMap[index];

    // Unless this pointer is null, return the cell
    if (p_cell)
    {
        return p_cell;
    }
    else
    {
        EXCEPTION("Location index input argument does not correspond to a Cell");
    }
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::IsCellAttachedToLocationIndex(unsigned index)
{
    // Check if pointer to cell is NULL
    if (mLocationCellMap[index])
    {
        return true;
    }
    else
    {
        return false;
    }
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetCellUsingLocationIndex(unsigned index, CellPtr pCell)
{
    // This currently assumes a one to one relationship between cells
    mLocationCellMap[index]= pCell;
}

template<unsigned DIM>
unsigned AbstractCellPopulation<DIM>::GetLocationIndexUsingCell(CellPtr pCell)
{
    return mCellLocationMap[pCell.get()];
}

template<unsigned DIM>
boost::shared_ptr<CellPropertyRegistry> AbstractCellPopulation<DIM>::GetCellPropertyRegistry()
{
    return mpCellPropertyRegistry;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetDefaultMutationStateOrdering()
{
    boost::shared_ptr<CellPropertyRegistry> p_registry = GetCellPropertyRegistry();
    if (!p_registry->HasOrderingBeenSpecified())
    {
        std::vector<boost::shared_ptr<AbstractCellProperty> > mutations;
        mutations.push_back(p_registry->Get<WildTypeCellMutationState>());
        mutations.push_back(p_registry->Get<ApcOneHitCellMutationState>());
        mutations.push_back(p_registry->Get<ApcTwoHitCellMutationState>());
        mutations.push_back(p_registry->Get<BetaCateninOneHitCellMutationState>());
        p_registry->SpecifyOrdering(mutations);
    }
}

template<unsigned DIM>
c_vector<double, DIM> AbstractCellPopulation<DIM>::GetCentroidOfCellPopulation()
{
    mCentroid = zero_vector<double>(DIM);
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        mCentroid += GetLocationOfCellCentre(*cell_iter);
    }
    mCentroid /= this->GetNumRealCells();

    return mCentroid;
}

//////////////////////////////////////////////////////////////////////////////
//                             Output methods                               //
//////////////////////////////////////////////////////////////////////////////


template<unsigned DIM>
void AbstractCellPopulation<DIM>::CreateOutputFiles(const std::string& rDirectory, bool cleanOutputDirectory)
{
    OutputFileHandler output_file_handler(rDirectory, cleanOutputDirectory);

    if(PetscTools::AmMaster())
    {
		mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes");
		mpVizBoundaryNodesFile = output_file_handler.OpenOutputFile("results.vizboundarynodes");
		mpVizCellProliferativeTypesFile = output_file_handler.OpenOutputFile("results.vizcelltypes");

		if (mOutputCellAncestors)
		{
			mpVizCellAncestorsFile = output_file_handler.OpenOutputFile("results.vizancestors");
		}
		if (mOutputCellMutationStates)
		{
			// An ordering must be specified for cell mutation states
			SetDefaultMutationStateOrdering();

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

template<unsigned DIM>
void AbstractCellPopulation<DIM>::CloseOutputFiles()
{
	// In parallel all files are closed after writing
	if(PetscTools::IsSequential())
	{
		mpVizNodesFile->close();
		mpVizBoundaryNodesFile->close();
		mpVizCellProliferativeTypesFile->close();


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

template<unsigned DIM>
void AbstractCellPopulation<DIM>::GenerateCellResults(unsigned locationIndex,
                                              std::vector<unsigned>& rCellProliferativeTypeCounter,
                                              std::vector<unsigned>& rCellCyclePhaseCounter)
{
    unsigned colour = STEM_COLOUR;

    CellPtr p_cell = GetCellUsingLocationIndex(locationIndex);

    if (mOutputCellCyclePhases)
    {
        // Update rCellCyclePhaseCounter
        switch (p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase())
        {
            case G_ZERO_PHASE:
                rCellCyclePhaseCounter[0]++;
                break;
            case G_ONE_PHASE:
                rCellCyclePhaseCounter[1]++;
                break;
            case S_PHASE:
                rCellCyclePhaseCounter[2]++;
                break;
            case G_TWO_PHASE:
                rCellCyclePhaseCounter[3]++;
                break;
             case M_PHASE:
                rCellCyclePhaseCounter[4]++;
                break;
            default:
                NEVER_REACHED;
        }
        *mpVizCellProliferativePhasesFile << p_cell->GetCellCycleModel()->GetCurrentCellCyclePhase() << " ";
    }

    if (mOutputCellAncestors)
    {
        // Set colour dependent on cell ancestor and write to file
        colour = p_cell->GetAncestor();
        if (colour == UNSIGNED_UNSET)
        {
            // Set the file to -1 to mark this case.
            colour = 1;
            *mpVizCellAncestorsFile << "-";
        }
        *mpVizCellAncestorsFile << colour << " ";
    }

    // Set colour dependent on cell type
    switch (p_cell->GetCellCycleModel()->GetCellProliferativeType())
    {
        case STEM:
            colour = STEM_COLOUR;
            if (mOutputCellProliferativeTypes)
            {
                rCellProliferativeTypeCounter[0]++;
            }
            break;
        case TRANSIT:
            colour = TRANSIT_COLOUR;
            if (mOutputCellProliferativeTypes)
            {
                rCellProliferativeTypeCounter[1]++;
            }
            break;
        case DIFFERENTIATED:
            colour = DIFFERENTIATED_COLOUR;
            if (mOutputCellProliferativeTypes)
            {
                rCellProliferativeTypeCounter[2]++;
            }
            break;
        default:
            NEVER_REACHED;
    }

    if (mOutputCellMutationStates)
    {
        // Set colour dependent on cell mutation state
        if (!p_cell->GetMutationState()->IsType<WildTypeCellMutationState>())
        {
            colour = p_cell->GetMutationState()->GetColour();
        }
        if (p_cell->HasCellProperty<CellLabel>())
        {
            CellPropertyCollection collection = p_cell->rGetCellPropertyCollection().GetProperties<CellLabel>();
            boost::shared_ptr<CellLabel> p_label = boost::static_pointer_cast<CellLabel>(collection.GetProperty());
            colour = p_label->GetColour();
        }
    }

    if (p_cell->HasCellProperty<ApoptoticCellProperty>() || p_cell->HasApoptosisBegun())
    {
        // For any type of cell set the colour to this if it is undergoing apoptosis
        colour = APOPTOSIS_COLOUR;
    }

    // Write cell variable data to file if required
    if (mOutputCellVariables && dynamic_cast<AbstractOdeBasedCellCycleModel*>(p_cell->GetCellCycleModel()) )
    {
        // Write location index corresponding to cell
        *mpCellVariablesFile << locationIndex << " ";

        // Write cell variables
        std::vector<double> proteins = (static_cast<AbstractOdeBasedCellCycleModel*>(p_cell->GetCellCycleModel()))->GetProteinConcentrations();
        for (unsigned i=0; i<proteins.size(); i++)
        {
            *mpCellVariablesFile << proteins[i] << " ";
        }
    }

    // Write cell age data to file if required
    if (mOutputCellAges)
    {
        // Write location index corresponding to cell
        *mpCellAgesFile << locationIndex << " ";

        // Write cell location
        c_vector<double, DIM> cell_location = GetLocationOfCellCentre(p_cell);

        for (unsigned i=0; i<DIM; i++)
        {
            *mpCellAgesFile << cell_location[i] << " ";
        }

        // Write cell age
        *mpCellAgesFile << p_cell->GetAge() << " ";
    }

    *mpVizCellProliferativeTypesFile << colour << " ";

}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::WriteCellResultsToFiles(std::vector<unsigned>& rCellProliferativeTypeCounter,
                                                  std::vector<unsigned>& rCellCyclePhaseCounter)
{
    *mpVizCellProliferativeTypesFile << "\n";

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

    // Write cell type data to file if required
    if (mOutputCellProliferativeTypes)
    {
        for (unsigned i=0; i<mCellProliferativeTypeCount.size(); i++)
        {
            mCellProliferativeTypeCount[i] = rCellProliferativeTypeCounter[i];
            *mpCellProliferativeTypesFile << rCellProliferativeTypeCounter[i] << "\t";
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
            mCellCyclePhaseCount[i] = rCellCyclePhaseCounter[i];
            *mpCellCyclePhasesFile << rCellCyclePhaseCounter[i] << "\t";
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

template<unsigned DIM>
void AbstractCellPopulation<DIM>::WriteTimeAndNodeResultsToFiles()
{
    OutputFileHandler output_file_handler(mDirPath, false);

    PetscTools::BeginRoundRobin();
    {
		if(!PetscTools::AmMaster() || SimulationTime::Instance()->IsEndTimeAndNumberOfTimeStepsSetUp())
		{
			mpVizNodesFile = output_file_handler.OpenOutputFile("results.viznodes", std::ios::app);
			mpVizBoundaryNodesFile = output_file_handler.OpenOutputFile("results.vizboundarynodes", std::ios::app);
		}
		if(PetscTools::AmMaster())
		{
			double time = SimulationTime::Instance()->GetTime();

			*mpVizNodesFile << time << "\t";
			*mpVizBoundaryNodesFile << time << "\t";
		}
		// Write node data to file
		for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = mrMesh.GetNodeIteratorBegin();
				node_iter != mrMesh.GetNodeIteratorEnd();
				++node_iter)
		{
			if (!node_iter->IsDeleted())
			{
				const c_vector<double,DIM>& position = node_iter->rGetLocation();

				for (unsigned i=0; i<DIM; i++)
				{
					*mpVizNodesFile << position[i] << " ";
				}
				*mpVizBoundaryNodesFile << node_iter->IsBoundaryNode() << " ";
			}
		}
		if(PetscTools::AmTopMost())
		{
			*mpVizNodesFile << "\n";
			*mpVizBoundaryNodesFile << "\n";
		}

		mpVizNodesFile->close();
		mpVizBoundaryNodesFile->close();
    }
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::WriteResultsToFiles()
{
    WriteTimeAndNodeResultsToFiles();
    double time = SimulationTime::Instance()->GetTime();

    *mpVizCellProliferativeTypesFile << time << "\t";

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
    GenerateCellResultsAndWriteToFiles();
    // Write logged cell data if required
    if (mOutputCellIdData)
    {
        WriteCellIdDataToFile();
    }

    // VTK can only be written in 2 or 3 dimensions
    if (DIM > 1)
    {
         WriteVtkResultsToFile();
    }
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::WriteCellIdDataToFile()
{
    // Write time to file
    *mpCellIdFile << SimulationTime::Instance()->GetTime();

    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = Begin();
         cell_iter != End();
         ++cell_iter)
    {
        unsigned cell_id = cell_iter->GetCellId();
        unsigned location_index = mCellLocationMap[(*cell_iter).get()];
        *mpCellIdFile << " " << cell_id << " " << location_index;

        c_vector<double, DIM> coords = GetLocationOfCellCentre(*cell_iter);
        for (unsigned i=0; i<DIM; i++)
        {
            *mpCellIdFile << " " << coords[i];
        }
    }
    *mpCellIdFile << "\n";
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::OutputCellPopulationInfo(out_stream& rParamsFile)
{
    std::string cell_population_type = GetIdentifier();

    *rParamsFile << "\t<" << cell_population_type << ">\n";
    OutputCellPopulationParameters(rParamsFile);
    *rParamsFile << "\t</" << cell_population_type << ">\n";
    *rParamsFile << "\n";
    *rParamsFile << "\t<CellCycleModels>\n";

    /*
     * Loop over cells and generate a set of cell-cycle model classes
     * that are present in the population.
     *
     * \todo this currently ignores different parameter regimes (#1453)
     */
    std::set<std::string> unique_cell_cycle_models;
    std::vector<CellPtr> first_cell_with_unique_CCM;
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
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

template<unsigned DIM>
void AbstractCellPopulation<DIM>::OutputCellPopulationParameters(out_stream& rParamsFile)
{
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

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellIdData()
{
    return mOutputCellIdData;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellMutationStates()
{
    return mOutputCellMutationStates;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellAncestors()
{
    return mOutputCellAncestors;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellProliferativeTypes()
{
    return mOutputCellProliferativeTypes;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellVariables()
{
    return mOutputCellVariables;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellCyclePhases()
{
    return mOutputCellCyclePhases;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellAges()
{
    return mOutputCellAges;
}

template<unsigned DIM>
bool AbstractCellPopulation<DIM>::GetOutputCellVolumes()
{
    return mOutputCellVolumes;
}

///////////////////////////////////////////////////////////////////////
// Setter methods
///////////////////////////////////////////////////////////////////////

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellIdData(bool writeCellIdData)
{
    mOutputCellIdData = writeCellIdData;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellMutationStates(bool outputCellMutationStates)
{
    mOutputCellMutationStates = outputCellMutationStates;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellAncestors(bool outputCellAncestors)
{
    mOutputCellAncestors = outputCellAncestors;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellProliferativeTypes(bool outputCellProliferativeTypes)
{
    mOutputCellProliferativeTypes = outputCellProliferativeTypes;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellVariables(bool outputCellVariables)
{
    mOutputCellVariables = outputCellVariables;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellCyclePhases(bool outputCellCyclePhases)
{
    mOutputCellCyclePhases = outputCellCyclePhases;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellAges(bool outputCellAges)
{
    mOutputCellAges = outputCellAges;
}

template<unsigned DIM>
void AbstractCellPopulation<DIM>::SetOutputCellVolumes(bool outputCellVolumes)
{
    mOutputCellVolumes = outputCellVolumes;
}

template<unsigned DIM>
c_vector<double,DIM> AbstractCellPopulation<DIM>::GetSizeOfCellPopulation()
{
    // Compute the centre of mass of the cell population
    c_vector<double,DIM> centre = GetCentroidOfCellPopulation();

    // Loop over cells and find the maximum distance from the centre of mass in each dimension
    c_vector<double,DIM> max_distance_from_centre = zero_vector<double>(DIM);
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = this->Begin();
         cell_iter != this->End();
         ++cell_iter)
    {
        c_vector<double,DIM> cell_location = GetLocationOfCellCentre(*cell_iter);
        c_vector<double,DIM> displacement = centre - cell_location;

        for (unsigned i=0; i<DIM; i++)
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

template class AbstractCellPopulation<1>;
template class AbstractCellPopulation<2>;
template class AbstractCellPopulation<3>;
