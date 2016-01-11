/*

Copyright (c) 2005-2016, University of Oxford.
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

#ifndef _FLUIDSOURCE_HPP_
#define _FLUIDSOURCE_HPP_

#include "ChastePoint.hpp"

/**
 * A fluid source in an Immersed Boundary Mesh, used in ImmersedBoundary simulations.
 */
template<unsigned SPACE_DIM>
class FluidSource
{
private:

    /** The index of this fluid source within the mesh. */
    unsigned mIndex;

    /** The location of this fluid source within the mesh. */
    c_vector<double, SPACE_DIM> mLocation;

    /** The strength of the fluid source. */
    double mStrength;

    /** Whether the fluid source is associated with an Immersed Boundary Element. */
    bool mIsSourceAssociatedWithElement;

    /** Index of the Immersed Boundary Element associated with this fluid source. */
    unsigned mAssociatedElementIndex;

public:

    /**
     * Constructor that takes in the fluid source's location as a ChastePoint.
     *
     * @param index  the index of the fluid source in the mesh
     * @param point  the location of the fluid source in the mesh
     */
    FluidSource(unsigned index, ChastePoint<SPACE_DIM> point);

    /**
     * Constructor that takes in the fluid source's location as a c_vector.
     *
     * @param index  the index of the fluid source in the mesh
     * @param location  the location of the fluid source in the mesh
     */
    FluidSource(unsigned index, c_vector<double, SPACE_DIM> location);

    /**
     * Constructor that takes the coordinates of the fluid source's location as separate input arguments.
     *
     * @param index  the index of the fluid source in the mesh
     * @param v1 the x-coordinate of the fluid source in the mesh (defaults to 0)
     * @param v2 the y-coordinate of the fluid source in the mesh (defaults to 0)
     * @param v3 the z-coordinate of the fluid source in the mesh (defaults to 0)
     */
    FluidSource(unsigned index, double v1=0.0, double v2=0.0, double v3=0.0);


    /**
     * Explicit destructor.
     */
    ~FluidSource();

    /**
     * @return the index of this fluid source in the mesh.
     */
    unsigned GetIndex() const;

    /**
     * Set the index of the fluid source in the mesh.
     *
     * @param index the index of the fluid source
     */
    void SetIndex(unsigned index);

    /**
     * @return the fluid source's location as a ChastePoint.
     */
    ChastePoint<SPACE_DIM> GetPoint() const;

    /**
     * @return the fluid source's location as a c_vector.
     *
     * The returned location may not be modified; if you want that functionality use
     * rGetModifiableLocation instead.
     */
    const c_vector<double, SPACE_DIM>& rGetLocation() const;

    /**
     * @return the fluid source's location as a c_vector.
     *
     * Don't forget to assign the result of this call to a reference!
     */
    c_vector<double, SPACE_DIM>& rGetModifiableLocation();

    /**
     * @return the strength of this fluid source.
     */
    double GetStrength() const;

    /**
     * Set the new strength of the fluid source.
     *
     * @param strength of the fluid source
     */
    void SetStrength(double strength);

    /**
     * Set whether the fluid source is associated with an element.
     *
     * @param whether there is an association
     */
    void SetIfSourceIsAssociatedWithElement(bool associated);

    /**
     * @return whether the fluid source is associated with an element.
     */
    bool IsSourceAssociatedWithElement();

    /**
     * @return the index of the element associated to this fluid source.
     */
    unsigned GetAssociatedElementIndex() const;

    /**
     * Set the index of the element associated to this fluid source.
     *
     * @param index of the element associated to this fluid source.
     */
    void SetAssociatedElementIndex(unsigned associatedElementIndex);
};

#endif //_FLUIDSOURCE_HPP_
