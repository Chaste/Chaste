/*

Copyright (c) 2005-2013, University of Oxford.
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

#ifndef GAMMADISTRIBUTEDSTOCHASTICDURATIONCELLCYCLEMODEL_HPP_
#define GAMMADISTRIBUTEDSTOCHASTICDURATIONCELLCYCLEMODEL_HPP_

#include "AbstractSimpleCellCycleModel.hpp"
#include "RandomNumberGenerator.hpp"

/**
 * A stochastic cell-cycle model where cells keep dividing with a stochastic G1 duration
 * drawn from a gamma distribution with specified shape and scale parameters.
 */
class GammaDistributedStochasticDurationCellCycleModel : public AbstractSimpleCellCycleModel
{
    friend class TestSimpleCellCycleModels;

private:

    /** The shape parameter of the gamma distribution. This must be a positive real number. */
    double mShape;

    /** The scale parameter of the gamma distribution. This must be a positive real number. */
    double mScale;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the cell-cycle model and random number generator, never used directly - boost uses this.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSimpleCellCycleModel>(*this);

        // Make sure the RandomNumberGenerator singleton gets saved too
        SerializableSingleton<RandomNumberGenerator>* p_wrapper = RandomNumberGenerator::Instance()->GetSerializationWrapper();
        archive & p_wrapper;

        archive & mShape;
        archive & mScale;
    }

public:

    /**
     * Constructor.
     */
    GammaDistributedStochasticDurationCellCycleModel();

    /**
     * Overridden SetG1Duration().
     */
    void SetG1Duration();

    /**
     * Overridden builder method to create new copies of this cell-cycle model.
     * @return a pointer to the GammaDistributedStochasticDurationCellCycleModel created.
     */
    AbstractCellCycleModel* CreateCellCycleModel();

    /**
     * Set mShape.
     *
     * @param shape the value of the shape parameter
     */
    void SetShape(double shape);

    /**
     * @return mScale.
     *
     * @param scale the value of the scale parameter
     */
    void SetScale(double scale);

    /**
     * @return mShape.
     */
    double GetShape();

    /**
     * @return mScale.
     */
    double GetScale();

    /**
     * Overridden OutputCellCycleModelParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputCellCycleModelParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(GammaDistributedStochasticDurationCellCycleModel)

#endif /* GAMMADISTRIBUTEDSTOCHASTICDURATIONCELLCYCLEMODEL_HPP_ */
