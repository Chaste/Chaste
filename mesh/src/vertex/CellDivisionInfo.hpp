/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef CELLDIVISIONINFO_HPP_
#define CELLDIVISIONINFO_HPP_

#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "ChasteSerialization.hpp"
#include "UblasIncludes.hpp"

/**
 * Records information about cell division event
 */
template <unsigned int SPACE_DIM>
struct CellDivisionInfo
{
    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & mLocation;
        archive & mDaughterLocation1;
        archive & mDaughterLongAxis1;
        archive & mDaughterLocation2;
        archive & mDaughterLongAxis2;
        archive & mDivisionAxis;
    }

    /** The centroid of the mother cell */
    c_vector<double, SPACE_DIM> mLocation;

    /** The centroid of the first daughter cell */
    c_vector<double, SPACE_DIM> mDaughterLocation1;

    /** The orientation of the first daughter cell */
    c_vector<double, SPACE_DIM> mDaughterLongAxis1;

    /** The centroid of the second daughter cell */
    c_vector<double, SPACE_DIM> mDaughterLocation2;

    /** The orientation of the second daughter cell */
    c_vector<double, SPACE_DIM> mDaughterLongAxis2;

    /** The orientation of the division axis */
    c_vector<double, SPACE_DIM> mDivisionAxis;
};

#endif /* CELLDIVISIONINFO_HPP_ */
