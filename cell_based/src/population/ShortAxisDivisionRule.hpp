/*

Copyright (c) 2005-2014, University of Oxford.
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

#ifndef SHORTAXISDIVISIONRULE_HPP_
#define SHORTAXISDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "VertexBasedCellPopulation.hpp"

template<unsigned SPACE_DIM> class VertexBasedCellPopulation;

/**
 * A class for Vertex-based cell populations to use to generate the short axis of
 * a vertex cell, for use in cell division. This is the default rule that is used
 * in most of the vertex simulations.
 */
template <unsigned SPACE_DIM>
class ShortAxisDivisionRule
{
public:
    /**
     * Default Constructor.
     */
    ShortAxisDivisionRule(){};

    /**
     * Empty destructor.
     */
    virtual ~ShortAxisDivisionRule(){};

    /**
     * Calculate the vector that will divide the two halves of the existing Vertex cell
     * to form the boundary between parent and daughter cell.
     *
     * @param pParentCell  The existing vertex cell
     * @param rCellPopulation  The Vertex cell population
     * @return the division vector.
     */
    c_vector<double, SPACE_DIM> CalculateCellDivisionVector(CellPtr pParentCell,
                                                            VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation);

private:
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * Note that serialization of the mesh and cells is handled by load/save_construct_data.
     *
     * Note also that member data related to writers is not saved - output must
     * be set up again by the caller after a restart.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // When this is in an inheritance tree it will want to call a method like this:
        //archive & boost::serialization::base_object<AbstractVertexDivisionRule<SPACE_DIM> >(*this);
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ShortAxisDivisionRule)

#endif // SHORTAXISDIVISIONRULE_HPP_
