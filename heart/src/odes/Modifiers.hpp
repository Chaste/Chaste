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

#ifndef MODIFIERS_HPP_
#define MODIFIERS_HPP_

#include "AbstractModifier.hpp"
#include <cmath>

/**
 * This class allows modification of parameters by a scale factor.
 */
class FactorModifier : public AbstractModifier
{
private:
    /** Factor to multiply parameter of interest by. */
    double mFactor;

public:
    /**
     * Constructor
     * @param factor  scale factor to use, defaults to 1 (i.e. no effect)
     */
    FactorModifier(double factor=1)
        : mFactor(factor)
    {
    }

    /**
     * Perform the modification.
     *
     * @param param  the current value of the quantity which is being modified
     * @param time  the current simulation time
     */
    virtual double Calc(double param, double time)
    {
        return (param * mFactor);
    }
};

/**
 * This is just an example class to show how you might specify a custom
 * modifier to change a parameter through time. In this case it implements
 * a sin(time)*default_parameter factor modifier.
 */
class TimeModifier : public AbstractModifier
{
public:
    /**
     * Constructor
     */
    TimeModifier()
    {
    }

    /**
     * Perform the modification.
     *
     * @param param  the current value of the quantity which is being modified
     * @param time  the current simulation time
     */
    virtual double Calc(double param, double time)
    {
        return param * sin(time);
    }
};

/**
 * This class just returns a fixed value, regardless of the parameter's default or the time.
 */
class FixedModifier : public AbstractModifier
{
private:
    /** Fixed value to clamp parameter at */
    double mValue;

public:
    /**
     * Constructor
     * @param value  The fixed value to use.
     */
    FixedModifier(double value)
        : mValue(value)
    {
    }

    /**
     * Perform the modification.
     *
     * @param param  the current value of the quantity which is being modified
     * @param time  the current simulation time
     * @return  the fixed value (ignores inputs)
     */
    virtual double Calc(double param, double time)
    {
        return mValue;
    }
};

/**
 * This class returns the parameter's default value and does not modify it.
 */
class DummyModifier : public AbstractModifier
{
private:

public:
    /**
     * Default constructor
     */
    DummyModifier()
    {
    }

    /**
     * This calculate does nothing and returns the parameter 'unharmed'.
     *
     * @param param  the current value of the quantity which is being modified
     * @param time  the current simulation time
     */
    virtual double Calc(double param, double time)
    {
        return param;
    }
};


#endif  //MODIFIERS_HPP_
