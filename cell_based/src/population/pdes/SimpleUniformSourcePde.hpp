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

#ifndef SIMPLEUNIFORMSOURCEPDE_HPP_
#define SIMPLEUNIFORMSOURCEPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearEllipticPde.hpp"

/**
 *  A simple nutrient PDE which is not directly coupled to the cell population.
 */
template<unsigned DIM>
class SimpleUniformSourcePde : public AbstractLinearEllipticPde<DIM,DIM>
{
    friend class TestCellBasedPdes;

private:

    /** Needed for serialization.*/
    friend class boost::serialization::access;
    /**
     * Serialize the PDE and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
       archive & boost::serialization::base_object<AbstractLinearEllipticPde<DIM, DIM> >(*this);
       archive & mCoefficient;
    }

    /** Coefficient of consumption of nutrient by cells. */
    double mCoefficient;

public:

    /**
     * Constructor.
     *
     * @param coefficient the coefficient of consumption of nutrient by cells (defaults to 0.0)
     */
    SimpleUniformSourcePde(double coefficient=0.0);

    /**
     * @return mCoefficient
     */
    double GetCoefficient() const;

    /**
     * Overridden ComputeConstantInUSourceTerm() method.
     *
     * @param rX The point in space
     * @param pElement The element
     *
     * @return the constant in u part of the source term, i.e g(x) in
     *  Div(D Grad u)  +  f(x)u + g(x) = 0.
     */
    double ComputeConstantInUSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement);

    /**
     * Overridden ComputeLinearInUCoeffInSourceTerm() method.
     *
     * @param rX The point in space
     * @param pElement the element
     *
     * @return the coefficient of u in the linear part of the source term, i.e f(x) in
     *  Div(D Grad u)  +  f(x)u + g(x) = 0.
     */
    double ComputeLinearInUCoeffInSourceTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @param rX The point in space at which the diffusion term is computed
     *
     * @return a matrix.
     */
    c_matrix<double,DIM,DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleUniformSourcePde)

#endif /*SIMPLEUNIFORMSOURCEPDE_HPP_*/
