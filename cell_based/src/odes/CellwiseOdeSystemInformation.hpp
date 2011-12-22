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
