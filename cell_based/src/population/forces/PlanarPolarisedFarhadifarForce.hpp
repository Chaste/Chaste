/*

Copyright (c) 2005-2023, University of Oxford.
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

#ifndef PLANARPOLARISEDFARHADIFARFORCE_HPP_
#define PLANARPOLARISEDFARHADIFARFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "FarhadifarForce.hpp"


/**
 * A force class for use in vertex-based simulations. This force is based on the
 * energy function proposed by Farhadifar et al in  Curr. Biol., 2007, 17, 2095-2104, 
 * but with a planar polarised line tension parameter, similar to that proposed by 
 * Rauzi et al in Nat. Cell Biol., 2008, 10, 1401-1410.
 */
template<unsigned DIM>
class PlanarPolarisedFarhadifarForce : public FarhadifarForce<DIM>
{
friend class TestForces;

private:

    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<FarhadifarForce<DIM> >(*this);
        archive & mPlanarPolarisedLineTensionMultiplier;
    }

protected:

    /**
     * A scalar that multiplies the strength of the line tension term in the model for edges 
     * whose angle relative to the x axis are between 45 degrees and 135 degrees.
     */
    double mPlanarPolarisedLineTensionMultiplier;


public:

    /**
     * Constructor.
     */
    PlanarPolarisedFarhadifarForce();

    /**
     * Destructor.
     */
    virtual ~PlanarPolarisedFarhadifarForce();

    /**
     * Overridden GetLineTensionParameter() method.
     * 
     * Get the line tension parameter for the edge between two given nodes.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     * @param rVertexCellPopulation reference to the cell population
     *
     * @return the line tension parameter for this edge.
     */
    virtual double GetLineTensionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB, VertexBasedCellPopulation<DIM>& rVertexCellPopulation);

    /**
     * @return mPlanarPolarisedLineTensionMultiplier
     */
    double GetPlanarPolarisedLineTensionMultiplier();

    /**
     * Set mPlanarPolarisedLineTensionMultiplier.
     *
     * @param planarPolarisedLineTensionMultiplier the new value of mPlanarPolarisedLineTensionMultiplier
     */
    void SetPlanarPolarisedLineTensionMultiplier(double planarPolarisedLineTensionMultiplier);

    /**
     * @return mLineTensionParameter
     */
    double GetLineTensionParameter();

    /**
     * @return mBoundaryLineTensionParameter
     */
    double GetBoundaryLineTensionParameter();

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PlanarPolarisedFarhadifarForce)

#endif /*PLANARPOLARISEDFARHADIFARFORCE_HPP_*/
