/*

Copyright (c) 2005-2017, University of Oxford.
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

#ifndef UNIFORMSOURCEPARABOLICPDE_HPP_
#define UNIFORMSOURCEPARABOLICPDE_HPP_

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractLinearParabolicPdeSystem.hpp"

/**
 * An elliptic PDE to be solved numerically using the finite element method, for
 * coupling to a cell-based simulation.
 *
 * The PDE takes the form
 *
 * du/dt = Grad.(Grad(u)) + k,
 *
 * where the scalar k is specified by the member mSourceCoefficient, whose value
 * must be set in the constructor.
 *
 * Thus, there is no direct coupling between the cell-based simulation and the
 * terms of the PDE; here, the cell population just defines the spatial domain
 * on which to solve the PDE.
 */
template <unsigned DIM>
class UniformSourceParabolicPde : public AbstractLinearParabolicPdeSystem<DIM, DIM>
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
       archive & boost::serialization::base_object<AbstractLinearParabolicPdeSystem<DIM,DIM> >(*this);
       archive & mSourceCoefficient;
    }

    /** Constant source term (rate of production) of the dependent variable. */
    double mSourceCoefficient;

public:

    /**
     * Constructor.
     *
     * @param sourceCoefficient the source term coefficient, k (defaults to 0.0)
     */
    UniformSourceParabolicPde(double sourceCoefficient=0.0);

    /**
     * @return mSourceCoefficient
     */
    double GetCoefficient() const;

    /**
     * Overridden ComputeSourceTerm() method.
     *
     * @return the source term coefficient, k
     *
     * @param rX a point in space (unused)
     * @param rU the dependent variable, u, at the point x as a vector (unused)
     * @param pdeIndex the index of the PDE (unused)
     * @param pElement The mesh element that x is contained in (optional; unused)
     */
    virtual double ComputeSourceTerm(const ChastePoint<DIM>& rX,
                                     c_vector<double,1>& rU,
                                     unsigned pdeIndex,
                                     Element<DIM, DIM>* pElement=nullptr);

    /**
     * Overridden ComputeDiffusionTerm() method.
     *
     * @return the identity matrix
     *
     * @param rX a point in space (unused)
     * @param pdeIndex the index of the PDE (unused)
     * @param pElement The mesh element that x is contained in (optional; unused)
     */
    virtual c_matrix<double, DIM, DIM> ComputeDiffusionTerm(const ChastePoint<DIM>& rX,
                                                            unsigned pdeIndex,
                                                            Element<DIM, DIM>* pElement=nullptr);

    /**
     * Overridden ComputeDuDtCoefficientFunction() method.
     *
     * @return 1.0
     *
     * @param rX a point in space (unused)
     * @param pdeIndex the index of the PDE (unused)
     */
    double ComputeDuDtCoefficientFunction(const ChastePoint<DIM>& rX, unsigned pdeIndex);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(UniformSourceParabolicPde)

#endif //UNIFORMSOURCEPARABOLICPDE_HPP_
