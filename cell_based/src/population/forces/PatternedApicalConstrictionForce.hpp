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

#include "GeneralMonolayerVertexMeshForce.hpp"

///\todo document class (#2850)
class PatternedApicalConstrictionForce : public GeneralMonolayerVertexMeshForce
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<GeneralMonolayerVertexMeshForce >(*this);
    }

    /**
     * Target area for each patterned cell's apical surface.
     * Initialised to 0 in the constructor.
     */
    double mPatternedTargetApicalArea;

    /**
     * Strength of each patterned cell's apical surface area term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mPatternedApicalAreaParameter;

    /**
     * Strength of each patterned apical cell-cell interface length term in the energy expression.
     * Initialised to 0 in the constructor.
     */
    double mPatternedApicalEdgeParameter;

public:

    /**
     * Constuctor.
     */
    PatternedApicalConstrictionForce();

    /**
     * Destructor.
     */
    virtual ~PatternedApicalConstrictionForce();

    /**
     * Set mPatternedTargetApicalArea, mPatternedApicalAreaParameter and mPatternedApicalEdgeParameter.
     *
     * @param patternedApicalEdgeParameter the new value of mPatternedTargetApicalArea
     * @param patternedApicalAreaParameter the new value of mPatternedApicalAreaParameter
     * @param patternedTargetApicalArea the new value of mPatternedApicalEdgeParameter
     */
    void SetPatternedApicalParameter(const double patternedApicalEdgeParameter,
                                     const double patternedApicalAreaParameter,
                                     const double patternedTargetApicalArea);

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculate the force on each node in the vertex-based cell population based on the energy expression.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    virtual void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(MisraForce)
