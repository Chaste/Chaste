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

#ifndef ABSTRACTSEMRANDOMFORCE_HPP_
#define ABSTRACTSEMRANDOMFORCE_HPP_

#include <vector>

#include "AbstractForce.hpp"
#include "ClassIsAbstract.hpp"
#include "SemBasedCellPopulation.hpp"
#include "UblasVectorInclude.hpp"

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

/**
 * Base class for random forces in the subcellular element method.
 *
 * This class implements the overdamped Langevin force scaling used in the SEM
 * paper:
 *
 *   F = eta * sqrt(2 D / dt) * z
 *
 * where eta is the node damping constant, D is the diffusion constant, and z is
 * a dimensionless unit random variable supplied by subclasses.
 */
template<unsigned DIM>
class AbstractSemRandomForce : public AbstractForce<DIM>
{
private:

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mDiffusionConstant;
    }

protected:

    /** Diffusion constant used in the SEM Langevin noise term. */
    double mDiffusionConstant;

    /**
     * Generate dimensionless unit random noise for each supplied node.
     *
     * @param rNodes the SEM nodes to sample noise for
     *
     * @return a vector of DIM-dimensional unit noise vectors
     */
    virtual std::vector<c_vector<double, DIM> > GenerateUnitNoise(const std::vector<Node<DIM>*>& rNodes)=0;

public:

    /**
     * Constructor.
     */
    AbstractSemRandomForce();

    /**
     * Destructor.
     */
    virtual ~AbstractSemRandomForce() = default;

    /**
     * @return the diffusion constant
     */
    double GetDiffusionConstant() const;

    /**
     * @param diffusionConstant the new diffusion constant
     */
    void SetDiffusionConstant(double diffusionConstant);

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation) override;

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;
};

TEMPLATED_CLASS_IS_ABSTRACT_1_UNSIGNED(AbstractSemRandomForce)

#endif /*ABSTRACTSEMRANDOMFORCE_HPP_*/
