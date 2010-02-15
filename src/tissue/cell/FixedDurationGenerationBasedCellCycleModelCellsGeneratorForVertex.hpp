/*

Copyright (C) University of Oxford, 2005-2009

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
#ifndef FIXEDDURATIONGENERATIONBASSEDCELLCYCLEMODELCELLSGENERATORFORVERTEX_HPP_
#define FIXEDDURATIONGENERATIONBASSEDCELLCYCLEMODELCELLSGENERATORFORVERTEX_HPP_

#include "AbstractCellsGenerator.hpp"
#include "FixedDurationGenerationBasedCellCycleModel.hpp"
#include "VertexMesh.hpp"

/**
 * A helper class for generating a vector of cells with
 * SimpleWntCellCycleModels for a given vertex mesh.
 * 
 * \todo this should be added to the AbstractCellsGenerator class when 
 * vertex models are added to the main code.
 */
template<unsigned DIM>
class FixedDurationGenerationBasedCellCycleModelCellsGeneratorForVertex : public AbstractCellsGenerator<DIM>
{
public:

    /**
     * @return a pointer to a new StochasticDurationGenerationBasedCellCycleModel.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * @return default cell cycle time for a transit cell.
     */
    double GetTypicalTransitCellCycleTime();

    /**
     * @return default cell cycle time for a transit cell.
     */
    double GetTypicalStemCellCycleTime();

    /**
     * @return true (cells can always differentiate).
     */
    virtual bool CellsCanDifferentiate();

    /**
     * Generates cells of a specified cell cycle type under the correct
     * crypt conditions and gives random ages if required,
     * or gives them an age of 0.0 - creates least work for solver startup.
     *
     * @param rCells  An empty cells vector for this function to fill up
     * @param rMesh  The crypt mesh (A VERTEX MESH)
     * @param locationIndices the node indices corresponding to real cells
     * @param randomBirthTimes  Whether to assign the cells random birth times
     *    (this can be expensive computationally with ODE models)
     * @param y0  below this line cells are generation 0 (defaults to 0.8)
     * @param y1  below this line cells are generation 1 (defaults to 2.0)
     * @param y2  below this line cells are generation 2 (defaults to 3.0)
     * @param y3  below this line cells are generation 3 (defaults to 4.0)
     * @param initialiseCells  whether to initialise the cell cycle models as each
     *   cell is created
     * 
     * \todo this should be refactored with GenerateForCrypt in AbstractCellsGenerator when Vertex code is released
     */
    void GenerateForVertexCrypt(std::vector<TissueCell>& rCells,
                                  VertexMesh<2,2>& rMesh,
                                  const std::vector<unsigned> locationIndices,
                                  bool randomBirthTimes,
                                  double y0 = 0.8,
                                  double y1 = 2.0,
                                  double y2 = 3.0,
                                  double y3 = 4.0,
                                  bool initialiseCells = false);
};

#endif /*FIXEDDURATIONGENERATIONBASSEDCELLCYCLEMODELCELLSGENERATORFORVERTEX_HPP_*/
