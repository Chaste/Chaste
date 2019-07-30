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

#ifndef _FUNCTIONALBOUNDARYCONDITION_HPP_
#define _FUNCTIONALBOUNDARYCONDITION_HPP_

#include "AbstractBoundaryCondition.hpp"

/**
 * A boundary condition that takes a function pointer in its constructor, and
 * evaluates the function to determine the value of the condition at a given
 * point.
 */
template<unsigned SPACE_DIM>
class FunctionalBoundaryCondition : public AbstractBoundaryCondition<SPACE_DIM>
{
private:

    /**
     * The function pointer used to determine the value of the boundary condition at a given point.
     *
     * @param rX a point in space
     */
    double (*mFunction)(const ChastePoint<SPACE_DIM>& rX);

public:

    /**
     * Typical use:
     *  pBoundaryCondition = new FunctionalBoundaryCondition(&function_name);
     *
     * @param func Pointer to a function to be used for evaluating this boundary
     *     condition
     */
    FunctionalBoundaryCondition(double (*func)(const ChastePoint<SPACE_DIM>& rX));

    /**
     * @return the value of the boundary condition at a given point.
     *
     * @param rX a point in space
     */
    double GetValue(const ChastePoint<SPACE_DIM>& rX) const;
};

#endif //_FUNCTIONALBOUNDARYCONDITION_HPP_
