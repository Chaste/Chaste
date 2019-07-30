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

#ifndef HEARTEVENTHANDLER_HPP_
#define HEARTEVENTHANDLER_HPP_

#include "GenericEventHandler.hpp"

/**
 * An event handler class with event types suitable for cardiac electrophysiological
 * simulations.
 *
 * It also contains events suitable to most generic PDE solvers too.
 */
class HeartEventHandler : public GenericEventHandler<16, HeartEventHandler>
{
public:

    /** Character array holding heart event names. There are eleven heart events. */
    static const char* EventName[16];

    /** Definition of heart event types. */
    typedef enum
    {
        READ_MESH=0,
        INITIALISE,
        ASSEMBLE_SYSTEM,
        SOLVE_ODES,
        COMMUNICATION,
        ASSEMBLE_RHS,
        NEUMANN_BCS,
        DIRICHLET_BCS,
        SOLVE_LINEAR_SYSTEM,
        WRITE_OUTPUT,
        DATA_CONVERSION,
        POST_PROC,
        USER1,
        USER2,
        USER3,
        EVERYTHING
    } EventType;
};

#endif /*HEARTEVENTHANDLER_HPP_*/
