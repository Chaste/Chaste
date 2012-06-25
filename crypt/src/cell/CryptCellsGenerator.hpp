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

        boost::shared_ptr<AbstractCellProperty> p_state(CellPropertyRegistry::Instance()->Get<WildTypeCellMutationState>());

        CellPtr p_cell(new Cell(p_state, p_cell_cycle_model));
        p_cell->SetCellProliferativeType(cell_type);

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
