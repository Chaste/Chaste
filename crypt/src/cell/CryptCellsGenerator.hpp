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

#ifndef CRYPTCELLSGENERATOR_HPP_
#define CRYPTCELLSGENERATOR_HPP_

#include <boost/mpl/integral_c.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/if.hpp>

#include "CellsGenerator.hpp"

#include "CellPropertyRegistry.hpp"
#include "TetrahedralMesh.hpp"
#include "VertexMesh.hpp"

#include "StochasticDurationGenerationBasedCellCycleModel.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "WntCellCycleModel.hpp"
#include "SimpleWntCellCycleModel.hpp"
#include "StochasticWntCellCycleModel.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisOne.hpp"
#include "VanLeeuwen2009WntSwatCellCycleModelHypothesisTwo.hpp"
#include "Exception.hpp"


/**
 *  Small helper method, which returns whether the two classes given as the template parameters
 *  are identical or not
 */
template<class T1, class T2>
bool ClassesAreSame()
{
    using namespace boost::mpl;
    using namespace boost;
    typedef typename if_< is_same<T1, T2>, integral_c<unsigned, 1>, integral_c<unsigned, 0> >::type selector_t;
    return  (selector_t()==1);
}


/**
 * A subclass of CellsGenerator that generates cells for crypt simulations.
 *
 * It is templated over types of cell-cycle model.
 */
template<class CELL_CYCLE_MODEL>
class CryptCellsGenerator : public CellsGenerator<CELL_CYCLE_MODEL,2>
{
public:

    /**
     * Generates cells of a specified cell cycle type under the correct
     * crypt conditions and gives random ages if required,
     * or gives them an age of 0.0 - creates least work for solver startup.
     *
     * @param rCells  An empty cells vector for this function to fill up
     * @param pMesh  The crypt mesh (can be cylindrical)
     * @param locationIndices the node indices corresponding to real cells
     * @param randomBirthTimes  Whether to assign the cells random birth times
     *    (this can be expensive computationally with ODE models)
     * @param y0  below this line cells are generation 0 (defaults to 0.3)
     * @param y1  below this line cells are generation 1 (defaults to 2.0)
     * @param y2  below this line cells are generation 2 (defaults to 3.0)
     * @param y3  below this line cells are generation 3 (defaults to 4.0)
     * @param initialiseCells  whether to initialise the cell-cycle models as each
     *   cell is created
     */
    void Generate(std::vector<CellPtr>& rCells,
                  AbstractMesh<2,2>* pMesh,
                  const std::vector<unsigned> locationIndices,
                  bool randomBirthTimes,
                  double y0 = 0.3,
                  double y1 = 2.0,
                  double y2 = 3.0,
                  double y3 = 4.0,
                  bool initialiseCells = false);
};


template<class CELL_CYCLE_MODEL>
void CryptCellsGenerator<CELL_CYCLE_MODEL>::Generate(
                                      std::vector<CellPtr>& rCells,
                                      AbstractMesh<2,2>* pMesh,
                                      const std::vector<unsigned> locationIndices,
                                      bool randomBirthTimes,
                                      double y0,
                                      double y1,
                                      double y2,
                                      double y3,
                                      bool initialiseCells)
{
    CellPropertyRegistry::Instance()->Clear();

    RandomNumberGenerator* p_random_num_gen = RandomNumberGenerator::Instance();

    rCells.clear();

    unsigned mesh_size;
    if (dynamic_cast<TetrahedralMesh<2,2>*>(pMesh))
    {
        mesh_size = pMesh->GetNumNodes();
        unsigned num_cells = locationIndices.empty() ? pMesh->GetNumNodes() : locationIndices.size();
        rCells.reserve(num_cells);
    }
    else
    {
        /*
         * We cannot directly assert this dynamic cast. This is because the assert macro
         * doesn't understand the <2,2> and thinks that it is being passed two arguments.
         */
        bool is_vertex_mesh = (dynamic_cast<VertexMesh<2,2>*>(pMesh));
        if (!is_vertex_mesh)
        {
            NEVER_REACHED;
        }
        mesh_size = static_cast<VertexMesh<2,2>*>(pMesh)->GetNumElements();
        rCells.reserve(mesh_size);
    }

    for (unsigned i=0; i<mesh_size; i++)
    {
        CellProliferativeType cell_type;
        unsigned generation;

        double y = 0.0;

        if (dynamic_cast<TetrahedralMesh<2,2>*>(pMesh))
        {
            if (locationIndices.empty())
            {
                y = pMesh->GetNode(i)->GetPoint().rGetLocation()[1];

            }
            else if (std::find(locationIndices.begin(), locationIndices.end(), i) != locationIndices.end())
            {
                y = pMesh->GetNode(i)->GetPoint().rGetLocation()[1];
            }
        }
        else
        {
            /*
             * We cannot directly assert this dynamic cast. This is because the assert macro
             * doesn't understand the <2,2> and thinks that it is being passed two arguments.
             */
            bool is_vertex_mesh = (dynamic_cast<VertexMesh<2,2>*>(pMesh));
            if (!is_vertex_mesh)
            {
                NEVER_REACHED;
            }
            y = dynamic_cast<VertexMesh<2,2>*>(pMesh)->GetCentroidOfElement(i)[1];
        }

        CELL_CYCLE_MODEL* p_cell_cycle_model = new CELL_CYCLE_MODEL;
        p_cell_cycle_model->SetDimension(2);

        double typical_transit_cycle_time = p_cell_cycle_model->GetAverageTransitCellCycleTime();
        double typical_stem_cycle_time = p_cell_cycle_model->GetAverageStemCellCycleTime();

        double birth_time = 0.0;
        if (randomBirthTimes)
        {
            birth_time = -p_random_num_gen->ranf();
        }

        if (y <= y0)
        {
            cell_type = STEM;
            generation = 0;
            birth_time *= typical_stem_cycle_time; // hours
        }
        else if (y < y1)
        {
            cell_type = TRANSIT;
            generation = 1;
            birth_time *= typical_transit_cycle_time; // hours
        }
        else if (y < y2)
        {
            cell_type = TRANSIT;
            generation = 2;
            birth_time *= typical_transit_cycle_time; // hours
        }
        else if (y < y3)
        {
            cell_type = TRANSIT;
            generation = 3;
            birth_time *= typical_transit_cycle_time; // hours
        }
        else
        {
            cell_type = p_cell_cycle_model->CanCellTerminallyDifferentiate() ? DIFFERENTIATED : TRANSIT;
            generation = 4;
            birth_time *= typical_transit_cycle_time; // hours
        }

        if (dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(p_cell_cycle_model))
        {
            dynamic_cast<AbstractSimpleGenerationBasedCellCycleModel*>(p_cell_cycle_model)->SetGeneration(generation);
        }
        p_cell_cycle_model->SetCellProliferativeType(cell_type);

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));

        if (initialiseCells)
        {
            p_cell->InitialiseCellCycleModel();
        }

        p_cell->SetBirthTime(birth_time);

        if (locationIndices.empty())
        {
            rCells.push_back(p_cell);
        }
        else if (std::find(locationIndices.begin(), locationIndices.end(), i) != locationIndices.end())
        {
            rCells.push_back(p_cell);
        }
    }
}

#endif /* CRYPTCELLSGENERATOR_HPP_ */
