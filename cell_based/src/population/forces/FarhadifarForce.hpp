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

#ifndef FARHADIFARFORCE_HPP_
#define FARHADIFARFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include "Exception.hpp"

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

#include <iostream>

/**
 * A force class for use in Vertex-based simulations. This force is based on the
 * Energy function proposed by Farhadifar et al in  Curr. Biol., 2007, 17, 2095-2104.
 */


template<unsigned DIM>
class FarhadifarForce : public AbstractForce<DIM>
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
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mAreaElasticityParameter;
        archive & mPerimeterContractilityParameter;
        archive & mLineTensionParameter;
        archive & mBoundaryLineTensionParameter;
    }

protected:

    /**
     * The strength of the area term in the model. Corresponds to K_alpha in Farhadifar's paper.
     */
    double mAreaElasticityParameter;

    /**
     * The strength of the perimeter term in the model. Corresponds to Gamma_alpha in Farhadifar's paper.
     */
    double mPerimeterContractilityParameter;

    /**
     * The strength of the line tension term in the model. Lambda_{i,j} in Farhadifar's paper.
     */
    double mLineTensionParameter;

    /**
     * The strength of the line tension at the boundary. This term does correspond to Lambda_{i,j} in Farhadifar's paper.
     */
    double mBoundaryLineTensionParameter;


public:

    /**
     * Constructor.
     */
    FarhadifarForce();

    /**
     * Destructor.
     */
    virtual ~FarhadifarForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the energy function
     * Farhadifar's model.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
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
     * @return mAreaElasticityParameter
     */
    double GetAreaElasticityParameter();

    /**
     * @return mPerimeterContractilityParameter
     */
    double GetPerimeterContractilityParameter();

    /**
     * @return mLineTensionParameter
     */
    double GetLineTensionParameter();

    /**
     * @return mBoundaryLineTensionParameter
     */
    double GetBoundaryLineTensionParameter();

    /**
     * Set mAreaElasticityParameter.
     *
     * @param areaElasticityParameter the new value of mAreaElasticityParameter
     */
    void SetAreaElasticityParameter(double areaElasticityParameter);

    /**
     * Set mPerimeterContractilityParameter.
     *
     * @param perimeterContractilityParameter the new value of perimterContractilityParameter
     */
    void SetPerimeterContractilityParameter(double perimeterContractilityParameter);

    /**
     * Set mLineTensionParameter.
     *
     * @param lineTensionParameter the new value of mLineTensionParameter
     */
    void SetLineTensionParameter(double lineTensionParameter);

    /**
     * Set mBoundaryLineTensionParameter.
     *
     * @param boundaryLineTensionParameter the new value of mBoundaryLineTensionParameter
     */
    void SetBoundaryLineTensionParameter(double boundaryLineTensionParameter);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(FarhadifarForce)

#endif /*FARHADIFARFORCE_HPP_*/
