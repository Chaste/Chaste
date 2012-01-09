/*

Copyright (C) University of Oxford, 2005-2012

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

#ifndef CELLCYCLEPHASES_HPP_
#define CELLCYCLEPHASES_HPP_

/**
 * Possible phases of the cell cycle.
 *
 * When our cells 'divide' they are actually entering M phase,
 * so a cell progresses round the cell cycle in the following sequence from its birth time
 * Divide-> M -> G0/G1 -> S -> G2 -> Divide.
 *
 * G0 is a cell which stays in the G1 phase and is not going to divide. (i.e. quiescent or differentiated.)
 */
typedef enum CellCyclePhase_
{
    G_ZERO_PHASE,
    G_ONE_PHASE,
    S_PHASE,
    G_TWO_PHASE,
    M_PHASE
} CellCyclePhase;

static const unsigned NUM_CELL_CYCLE_PHASES=5;

#endif /*CELLCYCLEPHASES_HPP_*/
