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

#ifndef BIELMEIERFORCE_HPP_
#define BIELMEIERFORCE_HPP_

#include "GeneralMonolayerVertexMeshForce.hpp"

/**
 * A force class for use in vertex-based model simulations. This force is based on the
 * energy function proposed by Bielmeier et al in the following paper:
 *
 * Christina Bielmeier, Silvanus Alt, Vanessa Weichselberger, Marco La Fortezza,
 * Hartmann Harz, Frank Jülicher, Guillaume Salbreux, Anne-Kathrin Classen.
 * Interface Contractility between Differently Fated Cells Drives Cell Elimination and Cyst Formation.
 * Current Biology, 26(5), pp.563-574.
 * http://dx.doi.org/10.1016/j.cub.2015.12.063
 */
class BielmeierForce : public GeneralMonolayerVertexMeshForce
{
private:
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<GeneralMonolayerVertexMeshForce>(*this);
        archive& mEcmSpringConstant;
        archive& mExternalSurfaceTensionParameter;
    }

protected:
    /**
     * Spring modulus of elastic bond attaching the tissue to the external extracellular matrix (ECM)
     */
    double mEcmSpringConstant;

    /**
     * External in-plane surface tension parameter which constrains the area of tissue.
     */
    double mExternalSurfaceTensionParameter;

public:
    /**
     * Default constructor.
     *
     * @param springConstant initial value of mEcmSpringConstant
     * @param ExternalSurfaceTensionParameter initial value of mExternalSurfaceTensionParameter
     */
    BielmeierForce(const double springConstant = 5, const double ExternalSurfaceTensionParameter = -4.2);

    /**
     * Destructor.
     */
    ~BielmeierForce() override = default;

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculate the force on each node in the vertex-based cell population based on the energy function in Bielmeier et al's paper.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation) override;

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile) override;
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BielmeierForce)

#endif /*BIELMEIERFORCE_HPP_*/
