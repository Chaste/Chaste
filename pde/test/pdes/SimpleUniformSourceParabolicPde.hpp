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

#ifndef SIMPLEUNIFORMSOURCEPARABOLICPDE_HPP_
#define SIMPLEUNIFORMSOURCEPARABOLICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearParabolicPde.hpp"
#include "ChastePoint.hpp"
#include "Element.hpp"

/**
 * A simple parabolic PDE used in tests.
 */
template <unsigned DIM>
class SimpleUniformSourceParabolicPde : public AbstractLinearParabolicPde<DIM,DIM>
{
    friend class TestCellBasedParabolicPdes;

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
       archive & boost::serialization::base_object<AbstractLinearParabolicPde<DIM,DIM> >(*this);
       archive & mCoefficient;
    }

    /** Constant source term (rate of production) of the dependent variable. */
    double mCoefficient;

public:

    /**
     * Constructor.
     *
     * @param coefficient the source term (rate of production) of the dependent variable (defaults to 0.0)
     */
    SimpleUniformSourceParabolicPde(double coefficient=0.0);

    /**
     * @return mCoefficient
     */
    double GetCoefficient() const;

    /**
     * Overridden ComputeSourceTerm() method.
     *
     * @param rX the point in space at which the nonlinear source term is computed
     * @param u the value of the dependent variable at the point
     * @param pElement The element
     *
     * @return the the source term.
     */
    double ComputeSourceTerm(const ChastePoint<DIM>& rX, double u, Element<DIM,DIM>* pElement=NULL);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @param rX the point in space at which the diffusion term is computed
     * @param pElement the mesh element that x is contained in (optional; defaults to NULL).
     *
     * @return a matrix.
     */
    c_matrix<double, DIM, DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX, Element<DIM,DIM>* pElement=NULL);

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * @return the function c(x) in "c(x) du/dt = Grad.(DiffusionTerm(x)*Grad(u))+LinearSourceTerm(x)+NonlinearSourceTerm(x, u)"
     *
     * @param rX the point in space at which the function c is computed
     */
    double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(SimpleUniformSourceParabolicPde)

#endif //SIMPLEUNIFORMSOURCEPARABOLICPDE_HPP_
