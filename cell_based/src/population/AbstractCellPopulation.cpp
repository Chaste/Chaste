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
#include "PetscTools.hpp"

// Writers
#include "NodeLocationWriter.hpp"
#include "BoundaryNodeWriter.hpp"
#include "CellPopulationElementWriter.hpp"
#include "CellMutationStatesWriter.hpp"
#include "CellProliferativeTypesCountWriter.hpp"
#include "CellProliferativePhasesCountWriter.hpp"
#include "VoronoiDataWriter.hpp"
#include "VertexSwapWriters.hpp"
#include "CellWriters.hpp"

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
    DeleteWriters();
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::DeleteWriters()
{
    for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
    {
        delete mCellPopulationWriters[i];
    }
    for (unsigned i=0; i<mCellWriters.size(); i++)
    {
        delete mCellWriters[i];
    }
    mCellPopulationWriters.clear();
    mCellWriters.clear();
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
    return mCellMutationStateCount;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::vector<unsigned> AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GetCellProliferativeTypeCount()
{
    if (!mOutputCellProliferativeTypes)
    {
        EXCEPTION("Call SetOutputCellProliferativeTypes(true) before using this function");
    }
    return mCellProliferativeTypesCount;
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
    mDirPath = rDirectory;

    // Destroy any writers that are already open.
    DeleteWriters();

    if (mOutputResultsForChasteVisualizer)
    {
        AddPopulationWriter(new NodeLocationWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
        AddPopulationWriter(new BoundaryNodeWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
        AddCellWriter(new CellProliferativeTypesWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellAncestors)
    {
        AddCellWriter(new CellAncestorWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellMutationStates)
    {
        AddPopulationWriter(new CellMutationStatesWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellProliferativeTypes)
    {
        AddPopulationWriter(new CellProliferativeTypesCountWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellVariables)
    {
        AddCellWriter(new CellVariablesWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellCyclePhases)
    {
        AddPopulationWriter(new CellProliferativePhasesCountWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
        AddCellWriter(new CellProliferativePhasesWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellAges)
    {
        AddCellWriter(new CellAgesWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (mOutputCellIdData)
    {
        AddCellWriter(new CellIdWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }
    if (this->mOutputCellVolumes)
    {
        AddCellWriter(new CellVolumesWriter<ELEMENT_DIM, SPACE_DIM>(rDirectory));
    }

    if (cleanOutputDirectory)
    {
        // Use output handler to clean the directory.
        OutputFileHandler cleaner_handler(mDirPath, true);
    }

#ifdef CHASTE_VTK
    OutputFileHandler output_file_handler(mDirPath, cleanOutputDirectory);
    mpVtkMetaFile = output_file_handler.OpenOutputFile("results.pvd");
    *mpVtkMetaFile << "<?xml version=\"1.0\"?>\n";
    *mpVtkMetaFile << "<VTKFile type=\"Collection\" version=\"0.1\" byte_order=\"LittleEndian\" compressor=\"vtkZLibDataCompressor\">\n";
    *mpVtkMetaFile << "    <Collection>\n";
#endif //CHASTE_VTK
}


template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::CloseOutputFiles()
{
    for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
    {
        mCellPopulationWriters[i]->CloseFile();
    }
    for (unsigned i=0; i<mCellWriters.size(); i++)
    {
        mCellWriters[i]->CloseFile();
    }

#ifdef CHASTE_VTK
    *mpVtkMetaFile << "    </Collection>\n";
    *mpVtkMetaFile << "</VTKFile>\n";
    mpVtkMetaFile->close();
#endif //CHASTE_VTK
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::OpenWritersFiles()
{
    for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
    {
        mCellPopulationWriters[i]->OpenOutputFile();
        mCellPopulationWriters[i]->WriteHeader(this);
    }
    for (unsigned i=0; i<mCellWriters.size(); i++)
    {
        mCellWriters[i]->OpenOutputFile();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::OpenWritersFilesForAppend()
{
    for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
    {
        mCellPopulationWriters[i]->OpenOutputFileForAppend();
    }
    for (unsigned i=0; i<mCellWriters.size(); i++)
    {
        mCellWriters[i]->OpenOutputFileForAppend();
    }
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::ResetCellCounters()
{
    for (unsigned i=0; i<mCellCyclePhaseCount.size(); i++)
    {
        mCellCyclePhaseCount[i] = 0;
    }
    mCellProliferativeTypesCount.clear();
    mCellMutationStateCount.clear();

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::GenerateCellResults()
{
    for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->Begin();
            cell_iter != this->End();
            ++cell_iter)
    {
        if (mOutputCellCyclePhases)
        {
            // Update mCellCyclePhaseCount
            switch ((*cell_iter)->GetCellCycleModel()->GetCurrentCellCyclePhase())
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
        }
    }

    // Reduce results onto all processes.
    if (PetscTools::IsParallel() && mOutputCellCyclePhases)
    {
        std::vector<unsigned> phase_counts(mCellCyclePhaseCount.size(), 0u);

        for (unsigned i=0; i<phase_counts.size(); i++)
        {
            MPI_Allreduce(&mCellCyclePhaseCount[i], &phase_counts[i], 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        }

        mCellCyclePhaseCount = phase_counts;
    }

    /*
     * Calculate proliferative types count.
     */

    // An ordering must be specified for cell mutation states and cell proliferative types
    SetDefaultCellMutationStateAndProliferativeTypeOrdering();

    const std::vector<boost::shared_ptr<AbstractCellProperty> >& r_cell_properties =
        mpCellPropertyRegistry->rGetAllCellProperties();

    for (unsigned i=0; i<r_cell_properties.size(); i++)
    {
        if (r_cell_properties[i]->IsSubType<AbstractCellProliferativeType>())
        {
            mCellProliferativeTypesCount.push_back(r_cell_properties[i]->GetCellCount());
        }
    }

    // Reduce results onto all processes.
    if (PetscTools::IsParallel())
    {
        // Make sure the vector on each process has the same size
        unsigned local_size = mCellProliferativeTypesCount.size();
        unsigned global_size;

        MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED, MPI_MAX, PetscTools::GetWorld());
        mCellProliferativeTypesCount.resize(global_size, 0u);

        std::vector<unsigned> types_counts(mCellProliferativeTypesCount.size(), 0u);

        for (unsigned i=0; i<types_counts.size(); i++)
        {
            MPI_Allreduce(&mCellProliferativeTypesCount[i], &types_counts[i], 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        }

        mCellProliferativeTypesCount = types_counts;
    }

    /**
     * Calculate mutationstates count.
     */
    for (unsigned i=0; i<r_cell_properties.size(); i++)
    {
        if (r_cell_properties[i]->IsSubType<AbstractCellMutationState>())
        {
            mCellMutationStateCount.push_back(r_cell_properties[i]->GetCellCount());
        }
    }

    // Reduce results onto all processes.
    if (PetscTools::IsParallel())
    {
        // Make sure the vector on each process has the same size
        unsigned local_size = mCellMutationStateCount.size();
        unsigned global_size;

        MPI_Allreduce(&local_size, &global_size, 1, MPI_UNSIGNED, MPI_MAX, PetscTools::GetWorld());
        mCellMutationStateCount.resize(global_size, 0u);

        std::vector<unsigned> mutation_counts(mCellMutationStateCount.size(), 0u);

        for (unsigned i=0; i<mutation_counts.size(); i++)
        {
            MPI_Allreduce(&mCellMutationStateCount[i], &mutation_counts[i], 1, MPI_UNSIGNED, MPI_SUM, PetscTools::GetWorld());
        }

        mCellMutationStateCount = mutation_counts;
    }

}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::WriteResultsToFiles()
{
    if (!(mCellWriters.empty() && mCellPopulationWriters.empty()))
    {
        ResetCellCounters();

        GenerateCellResults();

        PetscTools::BeginRoundRobin();
        {
            OpenWritersFilesForAppend();

            // The master writes time stamps.
            if (PetscTools::AmMaster())
            {
                for (unsigned i=0; i<mCellWriters.size(); i++)
                {
                    mCellWriters[i]->WriteTimeStamp();
                }
                for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
                {
                    mCellPopulationWriters[i]->WriteTimeStamp();
                }
            }

            // Every process writes to file.
            for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
            {
                AcceptPopulationWriter(mCellPopulationWriters[i]);
            }
            for (typename AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>::Iterator cell_iter = this->Begin();
                    cell_iter != this->End();
                    ++cell_iter)
            {
                for (unsigned i=0; i<mCellWriters.size(); i++)
                {
                    AcceptCellWriter(mCellWriters[i], *cell_iter);
                }
            }

            // The top-most adds a newline
            if (PetscTools::AmTopMost())
            {
                for (unsigned i=0; i<mCellWriters.size(); i++)
                {
                    mCellWriters[i]->WriteNewline();
                }
                for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
                {
                    mCellPopulationWriters[i]->WriteNewline();
                }
            }

            // Then the files are closed.
            for (unsigned i=0; i<mCellPopulationWriters.size(); i++)
            {
                mCellPopulationWriters[i]->CloseFile();
            }
            for (unsigned i=0; i<mCellWriters.size(); i++)
            {
                mCellWriters[i]->CloseFile();
            }
        }
        PetscTools::EndRoundRobin();
    }

    // VTK can only be written in 2 or 3 dimensions.
    if (SPACE_DIM > 1)
    {
         WriteVtkResultsToFile();
    }
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

        // Note that we define this vector before setting it as otherwise the profiling build will break (see #2367)
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

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
std::pair<unsigned,unsigned> AbstractCellPopulation<ELEMENT_DIM,SPACE_DIM>::CreateOrderedPair(unsigned index1, unsigned index2)
{
    assert(index1 != index2);

    std::pair<unsigned, unsigned> ordered_pair;

    if (index1 < index2)
    {
        ordered_pair.first = index1;
        ordered_pair.second = index2;
    }
    else
    {
        ordered_pair.first = index2;
        ordered_pair.second = index1;
    }
    return ordered_pair;
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
