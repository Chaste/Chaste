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

#ifndef ABSTRACTMODIFIER_HPP_
#define ABSTRACTMODIFIER_HPP_

/**
 * This family of classes are used to add simple functions into cell models to modify
 * particular quantities on the fly.  Rather than the model using the quantity directly
 * in computing its right-hand side, it calls calc() with the current value and uses
 * the result of that instead.
 *
 * Clearly for this to work the cell model must be modified to include calls to instances
 * of these classes.  PyCml has some experimental support for this, generating subclasses
 * of AbstractCardiacCellWithModifiers.
 */
class AbstractModifier
{
  public:
    /**
     * Default constructor.
     */
    AbstractModifier(void)
    {
    }

    /**
     * Default destructor.
     */
    virtual ~AbstractModifier()
    {
    }

    /**
     * Pure virtual function which must be overriden in subclasses to actually
     * perform the modification.
     *
     * @param param  the current value of the quantity which is being modified
     * @param time  the current simulation time
     */
    virtual double Calc(double param, double time) = 0;
};


#endif  //ABSTRACTMODIFIER_HPP_

