/*

Copyright (C) University of Oxford, 2005-2010

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
#include "SimpleWntCellCycleModelCellsGeneratorForVertex.hpp"
#include "WildTypeCellMutationState.hpp"

template<unsigned DIM>
AbstractCellCycleModel* SimpleWntCellCycleModelCellsGeneratorForVertex<DIM>::CreateCellCycleModel()
{
	SimpleWntCellCycleModel* p_model = new SimpleWntCellCycleModel();
	p_model->SetDimension(DIM);
    return p_model;
}


template<unsigned DIM>
double SimpleWntCellCycleModelCellsGeneratorForVertex<DIM>::GetTypicalTransitCellCycleTime()
{
    return TissueConfig::Instance()->GetTransitCellG1Duration()
            + TissueConfig::Instance()->GetSG2MDuration();
}


template<unsigned DIM>
double SimpleWntCellCycleModelCellsGeneratorForVertex<DIM>::GetTypicalStemCellCycleTime()
{
    return TissueConfig::Instance()->GetStemCellG1Duration()
            + TissueConfig::Instance()->GetSG2MDuration();
}

// \TODO merge this with the GenreateForCrypt method in AbsractCellsGenerator
template<unsigned DIM>
void SimpleWntCellCycleModelCellsGeneratorForVertex<DIM>::GenerateForVertexCrypt(std::vector<TissueCell>& rCells,
                                 VertexMesh<2,2>& rMesh,
                                 bool randomBirthTimes,
                                 bool initialiseCells)
{
    #define COVERAGE_IGNORE
    assert(DIM==2);
    #undef COVERAGE_IGNORE

    RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

    unsigned num_cells = rMesh.GetNumElements();


    AbstractCellCycleModel* p_cell_cycle_model = NULL;
    double typical_transit_cycle_time;
    double typical_stem_cycle_time;

    rCells.clear();
    rCells.reserve(num_cells);
    boost::shared_ptr<AbstractCellMutationState> p_state(new WildTypeCellMutationState);
    for (unsigned i=0; i<rMesh.GetNumElements(); i++)
    {
        CellProliferativeType cell_type;

        p_cell_cycle_model = CreateCellCycleModel();
        typical_transit_cycle_time = this->GetTypicalTransitCellCycleTime();
        typical_stem_cycle_time = GetTypicalStemCellCycleTime();

        double birth_time = 0.0;
        if (randomBirthTimes)
        {
            birth_time = -p_random_num_gen->ranf();
        }

        cell_type = TRANSIT;
        birth_time *= typical_transit_cycle_time; // hours

        TissueCell cell(cell_type, p_state, p_cell_cycle_model);
        if (initialiseCells)
        {
            cell.InitialiseCellCycleModel();
        }

        cell.SetBirthTime(birth_time);

        rCells.push_back(cell);
    }
}



/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

//template class SimpleWntCellCycleModelCellsGenerator<1>;
template class SimpleWntCellCycleModelCellsGeneratorForVertex<2>;
//template class SimpleWntCellCycleModelCellsGenerator<3>;
