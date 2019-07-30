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

#ifndef FIXEDVERTEXBASEDDIVISIONRULE_HPP_
#define FIXEDVERTEXBASEDDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "AbstractVertexBasedDivisionRule.hpp"
#include "VertexBasedCellPopulation.hpp"

// Forward declaration prevents circular include chain
template<unsigned SPACE_DIM> class VertexBasedCellPopulation;
template<unsigned SPACE_DIM> class AbstractVertexBasedDivisionRule;

/**
 * A class to generate a division vector of unit length specified in
 * the class constructor.
 *
 * This helper class is used in TestVertexhBasedCellPopulation.hpp.
 */
template<unsigned SPACE_DIM>
class FixedVertexBasedDivisionRule : public AbstractVertexBasedDivisionRule<SPACE_DIM>
{
private:

    /**
     * The specified division vector.
     * Initialized in the constructor.
     */
    c_vector<double, SPACE_DIM> mDivisionVector;

    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractVertexBasedDivisionRule<SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param rDivisionVector the specified division vector
     */
    FixedVertexBasedDivisionRule(c_vector<double, SPACE_DIM>& rDivisionVector);

    /**
     * Empty destructor.
     */
    virtual ~FixedVertexBasedDivisionRule()
    {
    }

    /**
     * @return mDivisionVector.
     */
    const c_vector<double, SPACE_DIM>& rGetDivisionVector() const;

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The vertex-based cell population
     * @return mDivisionVector.
     */
    virtual c_vector<double, SPACE_DIM> CalculateCellDivisionVector(CellPtr pParentCell,
        VertexBasedCellPopulation<SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FixedVertexBasedDivisionRule)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FixedVertexBasedDivisionRule.
 */
template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const FixedVertexBasedDivisionRule<SPACE_DIM>* t, const unsigned int file_version)
{
    // Archive c_vector one component at a time
    c_vector<double, SPACE_DIM> vector = t->rGetDivisionVector();
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar << vector[i];
    }
}

/**
 * De-serialize constructor parameters and initialize a FixedVertexBasedDivisionRule.
 */
template<class Archive, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, FixedVertexBasedDivisionRule<SPACE_DIM>* t, const unsigned int file_version)
{
    // Archive c_vector one component at a time
    c_vector<double, SPACE_DIM> vector;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar >> vector[i];
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)FixedVertexBasedDivisionRule<SPACE_DIM>(vector);
}
}
} // namespace ...

#endif // FIXEDVERTEXBASEDDIVISIONRULE_HPP_
