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

#ifndef GENERALMONOLAYERVERTEXMESHFORCE_HPP_
#define GENERALMONOLAYERVERTEXMESHFORCE_HPP_

#include "AbstractForce.hpp"

/**
 * A force for 3D "monolayer" vertex-based cell populations, in which each cell is a
 * prism-like VertexElement<3,3> with a distinguished apical face, basal face and a set
 * of lateral faces (see MonolayerVertexMeshGenerator).
 *
 * The force on each vertex is minus the gradient, with respect to that vertex's position,
 * of a free energy that is a sum over all cells of three kinds of contribution:
 *
 *  - a volume-elasticity term penalising the deviation of each cell's volume V from a
 *    target volume V_0, with modulus set by SetVolumeParameters();
 *  - an area (surface-tension) term for the apical, basal and lateral faces, proportional
 *    to face area, with independent moduli for apical, basal and lateral faces;
 *  - an edge (line-tension) term for the apical, basal and lateral edges, proportional to
 *    edge length, with independent moduli for apical, basal and lateral edges.
 *
 * The apical, basal, lateral and volume parameters are set independently via
 * SetApicalParameters(), SetBasalParameters(), SetLateralParameter() and
 * SetVolumeParameters(). The contributions are assembled by the helper methods
 * AddVolumeContribution(), AddAreaContribution() and AddEdgeContribution().
 */
class GeneralMonolayerVertexMeshForce : public AbstractForce<3>
{
private:
    friend class TestForces;

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<AbstractForce<3> >(*this);
        archive& mTargetApicalArea;
        archive& mApicalAreaParameter;
        archive& mApicalEdgeParameter;
        archive& mTargetBasalArea;
        archive& mBasalAreaParameter;
        archive& mBasalEdgeParameter;
        archive& mLateralAreaParameter;
        archive& mLateralEdgeParameter;
        archive& mTargetVolume;
        archive& mVolumeParameter;
    }

protected:
    ///\todo #480 Consider target area/volume growth in 3D

    /**
     * Target area for each cell's apical surface.
     * Initialised to 0 in the constructor.
     */
    double mTargetApicalArea = 0.0;

    /**
     * Strength of each cell's apical surface area term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mApicalAreaParameter = 0.0;

    /**
     * Strength of each apical cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mApicalEdgeParameter = 0.0;

    /**
     * Target area for each cell's basal surface.
     * Initialised to 0 in the constructor.
     */
    double mTargetBasalArea = 0.0;

    /**
     * Strength of each cell's basal surface area term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mBasalAreaParameter = 0.0;

    /**
     * Strength of each basal cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mBasalEdgeParameter = 0.0;

    /**
     * Strength of each cell's lateral surface area term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mLateralAreaParameter = 0.0;

    /**
     * Strength of each lateral (apico-basal) cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mLateralEdgeParameter = 0.0;

    /**
     * Target volume for each cell.
     * Initialised to 0 in the constructor.
     */
    double mTargetVolume = 0.0;

    /**
     * Strength of each cell's volume term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mVolumeParameter = 0.0;

    /**
     * Helper function which is called by AddForceContribution.
     * @param pCellPopulation pointer to the vertex based cell population
     */
    virtual void AddVolumeContribution(VertexBasedCellPopulation<3>* pCellPopulation);

    /**
     * Helper function which is called by AddForceContribution.
     * @param pCellPopulation pointer to the vertex based cell population
     */
    virtual void AddAreaContribution(VertexBasedCellPopulation<3>* pCellPopulation);

    /**
     * Helper function which is called by AddForceContribution.
     * @param pCellPopulation pointer to the vertex based cell population
     */
    virtual void AddEdgeContribution(VertexBasedCellPopulation<3>* pCellPopulation);

public:
    /**
     * Constructor.
     */
    GeneralMonolayerVertexMeshForce() = default;

    /**
     * Destructor.
     */
    ~GeneralMonolayerVertexMeshForce() override = default;

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculate the force on each node in the vertex-based cell population based on the energy expression.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation) override;

    /**
     * Set mApicalEdgeParameter, mApicalAreaParameter and mTargetApicalArea.
     *
     * @param lineParameter the new value of mApicalEdgeParameter
     * @param areaParameter the new value of mApicalAreaParameter (defaults to 0)
     * @param targetArea the new value of mTargetApicalArea (defaults to 0)
     */
    void SetApicalParameters(const double lineParameter, const double areaParameter = 0, const double targetArea = 0);

    /**
     * Set mBasalEdgeParameter, mBasalAreaParameter and mTargetBasalArea.
     *
     * @param lineParameter the new value of mBasalEdgeParameter
     * @param areaParameter the new value of mBasalAreaParameter (defaults to 0)
     * @param targetArea the new value of mTargetBasalArea (defaults to 0)
     */
    void SetBasalParameters(const double lineParameter, const double areaParameter = 0, const double targetArea = 0);

    /**
     * Set mLateralEdgeParameter.
     *
     * @param lineParameter the new value of mLateralEdgeParameter
     * @param areaParameter the new value of mLateralAreaParameter
     */
    void SetLateralParameter(const double lineParameter, const double areaParameter = 0);

    /**
     * Set mVolumeParameter and mTargetVolume.
     *
     * @param volumeParameter the new value of mVolumeParameter
     * @param targetVolume the new value of mTargetVolume (defaults to 0)
     */
    void SetVolumeParameters(const double volumeParameter, const double targetVolume = 0);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(GeneralMonolayerVertexMeshForce)

#endif /*GENERALMONOLAYERVERTEXMESHFORCE_HPP_*/
