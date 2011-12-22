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
#ifndef MECHANICSEVENTHANDLER_HPP_
#define MECHANICSEVENTHANDLER_HPP_

#include "GenericEventHandler.hpp"

/**
 * An event handler class with event types suitable for cardiac electromechanics
 * simulations.
 *
 * It also contains events suitable to most generic PDE solvers too.
 */
class MechanicsEventHandler : public GenericEventHandler<7,MechanicsEventHandler>
{
public:

    /** Definition of mechanics event types. */
    typedef enum
    {
        ASSEMBLE=0,
        SOLVE,
        UPDATE,
        ALL_MECH,
        NON_MECH,
        OUTPUT,
        ALL
    } MechanicsEventType;

    /** Character array holding mechanics event names. There are seven mechanics events. */
    static const char* EventName[7];
};

#endif /*MECHANICSEVENTHANDLER_HPP_*/
