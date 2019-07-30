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

#ifndef WELIKYOSTERFORCE_HPP_
#define WELIKYOSTERFORCE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "VertexBasedCellPopulation.hpp"

/**
 * A force class for use in vertex-based simulations, based on a mechanical
 * model proposed by M. Weliky and G. Oster ("The mechanical basis of cell rearrangement.
 * I. Epithelial morphogenesis during Fundulus epiboly", Development 109:373-386).
 *
 * The default values for the two model parameter member variables are our own best
 * estimates, since they are not given in the Weliky & Oster paper.
 */
template<unsigned DIM>
class WelikyOsterForce  : public AbstractForce<DIM>
{
friend class TestForces;

private:

    /**
     * Area parameter. Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     */
    double mWelikyOsterAreaParameter;

    /**
     * Perimeter parameter. Has units of kg s^-2 (cell size at equilibrium rest length)^-1.
     */
    double mWelikyOsterPerimeterParameter;

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
        archive & mWelikyOsterAreaParameter;
        archive & mWelikyOsterPerimeterParameter;
    }

public:

    /**
     * Constructor.
     */
    WelikyOsterForce();

    /**
     * Destructor.
     */
    ~WelikyOsterForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based cell population based on the
     * Weliky Oster model.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * @return mWelikyOsterAreaParameter.
     */
    double GetWelikyOsterAreaParameter();

    /**
     * @return mWelikyOsterPerimeterParameter.
     */
    double GetWelikyOsterPerimeterParameter();

    /**
     * Set mWelikyOsterAreaParameter.
     *
     * @param welikyOsterAreaParameter the new value of mWelikyOsterAreaParameter
     */
    void SetWelikyOsterAreaParameter(double welikyOsterAreaParameter);

    /**
     * Set mWelikyOsterPerimeterParameter.
     *
     * @param welikyOsterPerimeterParameter the new value of mWlikyOsterPerimeterParameter
     */
    void SetWelikyOsterPerimeterParameter(double welikyOsterPerimeterParameter);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(WelikyOsterForce)

#endif /*WELIKYOSTERFORCE_HPP_*/
