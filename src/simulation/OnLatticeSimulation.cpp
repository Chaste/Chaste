/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

*/

#include "OnLatticeSimulation.hpp"

#include "AbstractCentreBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"

#include "CellBasedEventHandler.hpp"
#include "LogFile.hpp"
#include "Version.hpp"
#include "ExecutableSupport.hpp"


template<unsigned DIM>
OnLatticeSimulation<DIM>::OnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                                bool deleteCellPopulationInDestructor,
                                                bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation, deleteCellPopulationInDestructor, initialiseCells),
      mOutputCellVelocities(false)
{
    if (!dynamic_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OnLatticeSimulations require a PottsBasedCellPopulation.");
    }

    mpStaticCastCellPopulation = static_cast<PottsBasedCellPopulation<DIM>*>(&rCellPopulation);

    this->mDt = 1.0/120.0; // 30 seconds
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    std::vector<c_vector<double, DIM> > old_cell_locations;

    if (mOutputCellVelocities)
    {
        unsigned num_cells = this->mrCellPopulation.GetNumRealCells();

        // Get the previous cell centre locations to calculate velocity
        old_cell_locations.reserve(num_cells);

        for (unsigned cell_index=0; cell_index<num_cells; cell_index++)
        {
            old_cell_locations[cell_index] = mpStaticCastCellPopulation->rGetMesh().GetCentroidOfElement(cell_index);
        }
    }

    std::vector<c_vector<double, DIM> > forces(1, zero_vector<double>(DIM));

    // Update node positions
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    this->mrCellPopulation.UpdateNodeLocations(forces, this->mDt);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);

    // Write cell velocities to file if required
    if (mOutputCellVelocities)
    {
        if (SimulationTime::Instance()->GetTimeStepsElapsed()%this->mSamplingTimestepMultiple == 0)
        {
            unsigned num_cells = this->mrCellPopulation.GetNumRealCells();

            *mpCellVelocitiesFile << SimulationTime::Instance()->GetTime() << "\t";
            for (unsigned cell_index=0; cell_index<num_cells; cell_index++)
            {
                const c_vector<double,DIM>& position = mpStaticCastCellPopulation->rGetMesh().GetCentroidOfElement(cell_index);
                c_vector<double, DIM> velocity = (position - old_cell_locations[cell_index])/this->mDt;
                *mpCellVelocitiesFile << cell_index  << " ";
                for (unsigned i=0; i<DIM; i++)
                {
                    *mpCellVelocitiesFile << position[i] << " ";
                }
                for (unsigned i=0; i<DIM; i++)
                {
                    *mpCellVelocitiesFile << velocity[i] << " ";
                }
            }
            *mpCellVelocitiesFile << "\n";
        }
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::SetupSolve()
{
    if (mOutputCellVelocities)
    {
        OutputFileHandler output_file_handler2(this->mSimulationOutputDirectory+"/", false);
        mpCellVelocitiesFile = output_file_handler2.OpenOutputFile("cellvelocities.dat");
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AfterSolve()
{
    if (mOutputCellVelocities)
    {
        mpCellVelocitiesFile->close();
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > > update_rule_collection = mpStaticCastCellPopulation->rGetUpdateRuleCollection();

    // Loop over the collection of update rules and output info for each
    *rParamsFile << "\n\t<UpdateRules>\n";
    for (typename std::vector<boost::shared_ptr<AbstractPottsUpdateRule<DIM> > >::iterator iter = update_rule_collection.begin();
         iter != update_rule_collection.end();
         ++iter)
    {
        (*iter)->OutputUpdateRuleInfo(rParamsFile);
    }
    *rParamsFile << "\t</UpdateRules>\n";
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractPottsUpdateRule<DIM> > pUpdateRule)
{
    mpStaticCastCellPopulation->AddUpdateRule(pUpdateRule);
}

template<unsigned DIM>
bool OnLatticeSimulation<DIM>::GetOutputCellVelocities()
{
    return mOutputCellVelocities;
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::SetOutputCellVelocities(bool outputCellVelocities)
{
    mOutputCellVelocities = outputCellVelocities;
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t<OutputCellVelocities>" << mOutputCellVelocities << "</OutputCellVelocities>\n";

    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class OnLatticeSimulation<1>;
template class OnLatticeSimulation<2>;
template class OnLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulation)
