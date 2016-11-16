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

#ifndef GENERALMONOLAYERVERTEXMESHFORCE_HPP_
#define GENERALMONOLAYERVERTEXMESHFORCE_HPP_

#include "AbstractForce.hpp"

/**
 * \todo Define this class and write down energy expression (#2850)
 */
class GeneralMonolayerVertexMeshForce : public AbstractForce<3>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<3> >(*this);
        archive & mTargetApicalArea;
        archive & mApicalAreaParameter;
        archive & mApicalEdgeParameter;
        archive & mTargetBasalArea;
        archive & mBasalAreaParameter;
        archive & mBasalEdgeParameter;
        archive & mLateralEdgeParameter;
        archive & mTargetVolume;
        archive & mVolumeParameter;
    }

protected:

    ///\todo #2850 Consider target area/volume growth in 3D

    /**
     * Target area for each cell's apical surface.
     * Initialised to 0 in the constructor.
     */
    double mTargetApicalArea;

    /**
     * Strength of each cell's apical surface area term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mApicalAreaParameter;

    /**
     * Strength of each apical cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mApicalEdgeParameter;

    /**
     * Target area for each cell's basal surface.
     * Initialised to 0 in the constructor.
     */
    double mTargetBasalArea;

    /**
     * Strength of each cell's basal surface area term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mBasalAreaParameter;

    /**
     * Strength of each basal cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mBasalEdgeParameter;

    /**
     * Strength of each lateral (apico-basal) cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mLateralEdgeParameter;

    /**
     * Target volume for each cell.
     * Initialised to 0 in the constructor.
     */
    double mTargetVolume;

    /**
     * Strength of each cell's volume term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mVolumeParameter;

public:

    /**
     * Constuctor.
     */
    GeneralMonolayerVertexMeshForce();

    /**
     * Destructor.
     */
    virtual ~GeneralMonolayerVertexMeshForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculate the force on each node in the vertex-based cell population based on the energy expression.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation);

    /**
     * Set mApicalEdgeParameter, mApicalAreaParameter and mTargetApicalArea.
     *
     * @param lineParameter the new value of mApicalEdgeParameter
     * @param areaParameter the new value of mApicalAreaParameter (defaults to 0)
     * @param targetArea the new value of mTargetApicalArea (defaults to 0)
     */
    void SetApicalParameters(const double lineParameter, const double areaParameter=0, const double targetArea=0);

    /**
     * Set mBasalEdgeParameter, mBasalAreaParameter and mTargetBasalArea.
     *
     * @param lineParameter the new value of mBasalEdgeParameter
     * @param areaParameter the new value of mBasalAreaParameter (defaults to 0)
     * @param targetArea the new value of mTargetBasalArea (defaults to 0)
     */
    void SetBasalParameters(const double lineParameter, const double areaParameter=0, const double targetArea=0);

    /**
     * Set mLateralEdgeParameter.
     *
     * @param parameter the new value of mLateralEdgeParameter
     */
    void SetLateralParameter(const double parameter);

    /**
     * Set mVolumeParameter and mTargetVolume.
     *
     * @param volumeParameter the new value of mVolumeParameter
     * @param targetVolume the new value of mTargetVolume (defaults to 0)
     */
    void SetVolumeParameters(const double volumeParameter, const double targetVolume=0);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(GeneralMonolayerVertexMeshForce)

#endif /*GENERALMONOLAYERVERTEXMESHFORCE_HPP_*/
