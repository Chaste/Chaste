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

#ifndef CELLWISEODESYSTEMINFORMATION_HPP_
#define CELLWISEODESYSTEMINFORMATION_HPP_

#include "AbstractOdeSystemInformation.hpp"

/**
 * Concrete implementation of AbstractOdeSystemInformation designed for use where
 * some of the system information *can* vary across instances.  As with
 * OdeSystemInformation it is templated by ODE system class, to aid developers of
 * ODE system classes - only a specialisation of the Initialise method is required.
 *
 * Note: unexpected behaviour can occur if ODE system objects are copied
 * (via copy constructor or operator=).  The AbstractOdeSystem maintains a smart
 * pointer (boost::shared_ptr) to the system information object.  Hence both the
 * original and the copy will share the same information object.
 */
template<class ODE_SYSTEM>
class CellwiseOdeSystemInformation : public AbstractOdeSystemInformation
{
public:

    /**
     * Default constructor; calls Initialise.
     *
     * Designed to be used as follows by ODE system classes in their constructors:
     *   mpSystemInfo.reset(new CellwiseOdeSystemInformation<CLASS>);
     */
    CellwiseOdeSystemInformation();

protected:

    /**
     * Generic implementation of Initialise, which does nothing.
     *
     * Developers should specialise this method to their ODE system.  For example,

        template<>
        void CellwiseOdeSystemInformation<MyNewOdeSystem>::Initialise()
        {
            this->mVariableNames.push_back("Variable_1");
            this->mVariableUnits.push_back("Units_1");
            this->mInitialConditions.push_back(0.0);

            this->mInitialised = true;
        }
     */
    void Initialise();
};

template<class ODE_SYSTEM>
CellwiseOdeSystemInformation<ODE_SYSTEM>::CellwiseOdeSystemInformation()
{
    CellwiseOdeSystemInformation<ODE_SYSTEM>::Initialise();
}

template<class ODE_SYSTEM>
void CellwiseOdeSystemInformation<ODE_SYSTEM>::Initialise()
{
}

#endif /*CELLWISEODESYSTEMINFORMATION_HPP_*/
