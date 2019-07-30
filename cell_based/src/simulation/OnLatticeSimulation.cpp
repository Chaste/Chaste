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

#include "OnLatticeSimulation.hpp"
#include "AbstractOnLatticeCellPopulation.hpp"
#include "CellBasedEventHandler.hpp"

template<unsigned DIM>
OnLatticeSimulation<DIM>::OnLatticeSimulation(AbstractCellPopulation<DIM>& rCellPopulation,
                                              bool deleteCellPopulationInDestructor,
                                              bool initialiseCells)
    : AbstractCellBasedSimulation<DIM>(rCellPopulation,
                                       deleteCellPopulationInDestructor,
                                       initialiseCells)
{
    if (!dynamic_cast<AbstractOnLatticeCellPopulation<DIM>*>(&rCellPopulation))
    {
        EXCEPTION("OnLatticeSimulations require a subclass of AbstractOnLatticeCellPopulation.");
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::AddUpdateRule(boost::shared_ptr<AbstractUpdateRule<DIM> > pUpdateRule)
{
    // This static_cast is fine, since otherwise an exception would have been thrown in the constructor
    static_cast<AbstractOnLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->AddUpdateRule(pUpdateRule);
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::RemoveAllUpdateRules()
{
    // This static_cast is fine, since otherwise an exception would have been thrown in the constructor
    static_cast<AbstractOnLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->RemoveAllUpdateRules();
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellLocationsAndTopology()
{
    CellBasedEventHandler::BeginEvent(CellBasedEventHandler::POSITION);
    static_cast<AbstractOnLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->UpdateCellLocations(this->mDt);
    CellBasedEventHandler::EndEvent(CellBasedEventHandler::POSITION);
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::UpdateCellPopulation()
{
    bool update_cell_population_this_timestep = (this->mInitialiseCells) || (SimulationTime::Instance()->GetTimeStepsElapsed() != 0);
    if (update_cell_population_this_timestep)
    {
        AbstractCellBasedSimulation<DIM>::UpdateCellPopulation();
    }
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputAdditionalSimulationSetup(out_stream& rParamsFile)
{
    // Loop over the collection of update rules and output info for each
    *rParamsFile << "\n\t<UpdateRules>\n";

    // This static_cast is fine, since otherwise an exception would have been thrown in the constructor
    std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > > collection =
            static_cast<AbstractOnLatticeCellPopulation<DIM>*>(&(this->mrCellPopulation))->GetUpdateRuleCollection();

    for (typename std::vector<boost::shared_ptr<AbstractUpdateRule<DIM> > >::iterator iter = collection.begin();
         iter != collection.end();
         ++iter)
    {
        (*iter)->OutputUpdateRuleInfo(rParamsFile);
    }

    *rParamsFile << "\t</UpdateRules>\n";
}

template<unsigned DIM>
void OnLatticeSimulation<DIM>::OutputSimulationParameters(out_stream& rParamsFile)
{
    // Call method on direct parent class
    AbstractCellBasedSimulation<DIM>::OutputSimulationParameters(rParamsFile);
}

// Explicit instantiation
template class OnLatticeSimulation<1>;
template class OnLatticeSimulation<2>;
template class OnLatticeSimulation<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(OnLatticeSimulation)
