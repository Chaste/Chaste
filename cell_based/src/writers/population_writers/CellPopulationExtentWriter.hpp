/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef CELLPOPULATIONEXTENTWRITER_HPP_
#define CELLPOPULATIONEXTENTWRITER_HPP_

#include <cmath>

#include "AbstractCellPopulationWriter.hpp"
#include "UblasVectorInclude.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * A writer that records how far the cell population extends along each axis at every output
 * timestep, measured as the component of its radius of gyration along that axis: the
 * root-mean-square displacement of its nodes from their centroid, taken one axis at a time.
 *
 * One line is written per output timestep, containing the simulation time followed by the value
 * for each axis in turn.
 *
 * A radius of gyration is used in preference to the bounding box because a bounding box is fixed by
 * just two nodes, whichever happen to lie furthest apart, which makes it both noisy and sensitive to
 * how a discrete arrangement of nodes is oriented. A population built on a lattice can report quite
 * different bounding boxes along different axes while being genuinely isotropic, simply because the
 * lattice reaches further along one axis than another. Averaging over every node removes both
 * problems.
 *
 * Note that the value written is an RMS displacement, not a width: for a uniformly filled sphere of
 * radius R it is R/sqrt(5) on each axis rather than 2R. It is proportional to size, so ratios and
 * strains are unaffected, but it should not be compared directly against a bounding box.
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CellPopulationExtentWriter : public AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM>
{
private:

    /** Needed for serialization. */
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
        archive & boost::serialization::base_object<AbstractCellPopulationWriter<ELEMENT_DIM, SPACE_DIM> >(*this);
    }

    /**
     * Write the per-axis radius of gyration of the given population. This depends only on where the
     * nodes are, so every Visit() method below defers to it.
     *
     * This is templated over the population type rather than taking an
     * AbstractCellPopulation<ELEMENT_DIM, SPACE_DIM>: every population type except the mesh-based
     * one is templated over SPACE_DIM alone, so for a writer instantiated with ELEMENT_DIM
     * different from SPACE_DIM they share no common base with it.
     *
     * @param pCellPopulation a pointer to the population to write
     */
    template<class POPULATION_TYPE>
    void WriteExtent(POPULATION_TYPE* pCellPopulation)
    {
        c_vector<double, SPACE_DIM> centroid = zero_vector<double>(SPACE_DIM);
        unsigned num_nodes = 0u;

        for (auto node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            centroid += node_iter->rGetLocation();
            num_nodes++;
        }

        // A population with no nodes has no centroid to measure from
        if (num_nodes == 0u)
        {
            for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
            {
                *this->mpOutStream << 0.0 << " ";
            }
            return;
        }

        centroid /= static_cast<double>(num_nodes);

        c_vector<double, SPACE_DIM> sum_of_squares = zero_vector<double>(SPACE_DIM);
        for (auto node_iter = pCellPopulation->rGetMesh().GetNodeIteratorBegin();
             node_iter != pCellPopulation->rGetMesh().GetNodeIteratorEnd();
             ++node_iter)
        {
            for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
            {
                const double displacement = node_iter->rGetLocation()[dim] - centroid[dim];
                sum_of_squares[dim] += displacement * displacement;
            }
        }

        for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
        {
            *this->mpOutStream << sqrt(sum_of_squares[dim] / static_cast<double>(num_nodes)) << " ";
        }
    }

public:

    /**
     * Default constructor.
     */
    CellPopulationExtentWriter();

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(MeshBasedCellPopulation<ELEMENT_DIM, SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(CaBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(NodeBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(PottsBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(VertexBasedCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(ImmersedBoundaryCellPopulation<SPACE_DIM>* pCellPopulation);

    /**
     * Visit the population and write its extent.
     *
     * @param pCellPopulation a pointer to the population to visit
     */
    virtual void Visit(SemBasedCellPopulation<SPACE_DIM>* pCellPopulation);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_ALL_DIMS(CellPopulationExtentWriter)

#endif /*CELLPOPULATIONEXTENTWRITER_HPP_*/
