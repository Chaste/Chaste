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

#ifndef NASHHUNTERPOLEZEROLAW_HPP_
#define NASHHUNTERPOLEZEROLAW_HPP_

#include "PoleZeroMaterialLaw.hpp"

/**
 *  The pole-zero law but with parameters set using values from Chapter 41 of "Cardiac Mechano-Electric
 *  Feedback and Arrhythmias: from Pipette to Patient" (eds Franz, Kohl, Sachs), Remme, Nash and Hunter, 2005.
 *
 *  The stiffness are in KILOPASCALS
 */
template<unsigned DIM>
class NashHunterPoleZeroLaw : public PoleZeroMaterialLaw<DIM>
{
friend class TestMaterialLaws;

public:

    /**
     * Constructor.
     */
    NashHunterPoleZeroLaw();
};

#endif /*NASHHUNTERPOLEZEROLAW_HPP_*/
