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
friend class TestCompressibleLawTransverselyIsotropic; // note: this law is not transversely isotropic with the default parameters

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

    /** @return the parameter a */
    double GetA()
    {
        return mA;
    }

    /**  @return the parameter b (the values which multiply the strains in Q) */
    std::vector<std::vector<double> > GetB()
    {
        return mB;
    }

    /** @return compressibility parameter */
    double GetCompressibilityParam()
    {
        return mCompressibilityParam;
    }
};

#endif /* COMPRESSIBLEEXPONENTIALLAW_HPP_ */
