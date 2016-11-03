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

template<unsigned DIM>
class GeneralMonolayerVertexMeshForce : public AbstractForce<DIM>
{
private:

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mTargetApicalArea;
        archive & mApicalareaParameter;
        archive & mApicalEdgeParameter;
        archive & mTargetBasalArea;
        archive & mBasalareaParameter;
        archive & mBasalEdgeParameter;
        archive & mLateralEdgeParameter;
        archive & mTargetVolume;
        archive & mVolumeParameter;
    }

protected:

    double mTargetApicalArea;
    double mApicalareaParameter;
    double mApicalEdgeParameter;

    double mTargetBasalArea;
    double mBasalareaParameter;
    double mBasalEdgeParameter;

    double mLateralEdgeParameter;

    double mTargetVolume;
    double mVolumeParameter;

public:
    GeneralMonolayerVertexMeshForce();

    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    void SetApicalParameter(const double lineParameter, const double areaParameter=0, const double targetArea=0);

    void SetBasalParameter(const double lineParameter, const double areaParameter=0, const double targetArea=0);

    void SetLateralParameter(const double parameter);

    void SetVolumeParameter(const double volumeParameter, const double targetVolume=0);

    void OutputForceParameters(out_stream& rParamsFile);
};
#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS1(GeneralMonolayerVertexMeshForce,3)

#endif /*GENERALMONOLAYERVERTEXMESHFORCE_HPP_*/
