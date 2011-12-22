/*

Copyright (C) University of Oxford, 2005-2011

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Chaste is free software: you can redistribute it and/or modify it
under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 2.1 of the License, or
(at your option) any later version.

Chaste is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
License for more details. The offer of Chaste under the terms of the
License is subject to the License being interpreted in accordance with
English Law and subject to any action against the University of Oxford
being under the jurisdiction of the English Courts.

You should have received a copy of the GNU Lesser General Public License
along with Chaste. If not, see <http://www.gnu.org/licenses/>.

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
     * @param rForces reference to vector of forces on nodes
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces,
                              AbstractCellPopulation<DIM>& rCellPopulation);

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
