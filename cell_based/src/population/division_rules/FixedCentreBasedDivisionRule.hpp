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

#ifndef FIXEDCENTREBASEDDIVISIONRULE_HPP_
#define FIXEDCENTREBASEDDIVISIONRULE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "AbstractCentreBasedDivisionRule.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

// Forward declaration prevents circular include chain
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedCellPopulation;
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM> class AbstractCentreBasedDivisionRule;

/**
 * A class to generate two daughter cell positions, one given by the
 * position of the dividing cell and the other specified in the
 * class constructor.
 *
 * This helper class is used in TestMeshBasedCellPopulation.hpp and
 * TestNodeBasedCellPopulation.hpp.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM=ELEMENT_DIM>
class FixedCentreBasedDivisionRule : public AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>
{
private:

    /**
     * The specified location of the new daughter cell.
     * Initialized in the constructor.
     */
    c_vector<double, SPACE_DIM> mDaughterLocation;

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
        archive & boost::serialization::base_object<AbstractCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

public:

    /**
     * Default constructor.
     *
     * @param rDaughterLocation the specified location of the daughter cell
     */
    FixedCentreBasedDivisionRule(c_vector<double, SPACE_DIM>& rDaughterLocation);

    /**
     * Empty destructor.
     */
    virtual ~FixedCentreBasedDivisionRule()
    {
    }

    /**
     * @return mDaughterLocation.
     */
    const c_vector<double, SPACE_DIM>& rGetDaughterLocation() const;

    /**
     * Overridden CalculateCellDivisionVector() method.
     *
     * @param pParentCell  The cell to divide
     * @param rCellPopulation  The centre-based cell population
     *
     * @return the two daughter cell positions.
     */
    virtual std::pair<c_vector<double, SPACE_DIM>, c_vector<double, SPACE_DIM> > CalculateCellDivisionVector(CellPtr pParentCell,
        AbstractCentreBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>& rCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(FixedCentreBasedDivisionRule)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a FixedCentreBasedDivisionRule.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const FixedCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Archive c_vector one component at a time
    c_vector<double, SPACE_DIM> location = t->rGetDaughterLocation();
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar << location[i];
    }
}

/**
 * De-serialize constructor parameters and initialize a FixedCentreBasedDivisionRule.
 */
template<class Archive, unsigned ELEMENT_DIM, unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, FixedCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>* t, const unsigned int file_version)
{
    // Archive c_vector one component at a time
    c_vector<double, SPACE_DIM> location;
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        ar >> location[i];
    }

    // Invoke inplace constructor to initialise instance
    ::new(t)FixedCentreBasedDivisionRule<ELEMENT_DIM, SPACE_DIM>(location);
}
}
} // namespace ...

#endif // FIXEDCENTREBASEDDIVISIONRULE_HPP_
