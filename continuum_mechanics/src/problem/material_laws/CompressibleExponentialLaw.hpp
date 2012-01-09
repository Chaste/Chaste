/*

Copyright (C) University of Oxford, 2005-2012

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
#ifndef COMPRESSIBLEEXPONENTIALLAW_HPP_
#define COMPRESSIBLEEXPONENTIALLAW_HPP_


#include "AbstractCompressibleMaterialLaw.hpp"

/**
 *  The compressible exponential material law implemented in
 *
 *  Uysk, Effect of Laminar Orthotropic Myofiber Architecture on Regional Stress and Strain in the Canine Left Ventricle,
 *  Journal of Elasticity, 2000.
 *
 *  W = a[exp(Q)-1]/2 + c (J ln(J) - J + 1)
 *  where Q = sum b_{MN} E_{MN}^2
 *
 *  The exponential term is the same form as in the SchmidCosta law, although the parameters
 *  here are those given in the paper cited above, not the same as in the SchmidCosta class.
 *
 *  Note, by default, the fibre direction is assumed to be THE X-DIRECTION,
 *  and the sheet direction the Y-DIRECTION (ie sheets in the XY plane). Call
 *  SetChangeOfBasisMatrix() before ComputeStressAndStressDerivative(), with
 *  the matrix P = [fibre_vec, sheet_vec, normal_vec] if this is not the case.
 */
template<unsigned DIM>
class CompressibleExponentialLaw : public AbstractCompressibleMaterialLaw<DIM>
{
private:

    /** Parameter a. (kPa) */
    double mA;

    /** Matrix of parameters b (dimensionless). */
    std::vector<std::vector<double> > mB;

    /** Compressibility parameter, c in W = a[exp(Q)-1]/2 + c (J ln(J) - J + 1) */
    double mCompressibilityParam;

    /** identity matrix. */
    c_matrix<double,DIM,DIM> mIdentity;

public:

    /** Constructor. */
    CompressibleExponentialLaw();

    /**
     *  Compute the (2nd Piola Kirchoff) stress T and the stress derivative dT/dE for
     *  a given strain.
     *
     *  NOTE: the strain E is not expected to be passed in, instead the Lagrangian
     *  deformation tensor C is required (recall, E = 0.5(C-I))
     *
     *  dT/dE is a fourth-order tensor, where dT/dE[M][N][P][Q] = dT^{MN}/dE_{PQ}
     *
     *  @param rC The Lagrangian deformation tensor (F^T F)
     *  @param rInvC The inverse of C. Should be computed by the user. (Change this?)
     *  @param pressure the current pressure
     *  @param rT the stress will be returned in this parameter
     *  @param rDTdE the stress derivative will be returned in this parameter, assuming
     *    the final parameter is true
     *  @param computeDTdE a boolean flag saying whether the stress derivative is
     *    required or not.
     */
    void ComputeStressAndStressDerivative(c_matrix<double,DIM,DIM>& rC,
                                          c_matrix<double,DIM,DIM >& rInvC,
                                          double pressure,
                                          c_matrix<double,DIM,DIM>& rT,
                                          FourthOrderTensor<DIM,DIM,DIM,DIM>& rDTdE,
                                          bool computeDTdE);

    /** Get method for the parameter a */
    double GetA()
    {
        return mA;
    }

    /**  Get method for the parameter b (the values which multiply the strains in Q) */
    std::vector<std::vector<double> > GetB()
    {
        return mB;
    }

    /** Get method for compressibility parameter */
    double GetCompressibilityParam()
    {
        return mCompressibilityParam;
    }
};

#endif /* COMPRESSIBLEEXPONENTIALLAW_HPP_ */
