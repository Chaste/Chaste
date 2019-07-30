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


#ifndef POLEZEROMATERIALLAW_HPP_
#define POLEZEROMATERIALLAW_HPP_

#include "AbstractIncompressibleMaterialLaw.hpp"
#include "Exception.hpp"

/**
 *  Pole-zero material law, as stated in: "Computational mechanics of the heart: from
 *  tissue structure to ventricular function" Nash, Hunter, J. Elasticity, 2000; or
 *  in Chapter 41 of "Cardiac Mechano-Electric Feedback and Arrhythmias: from
 *  Pipette to Patient" (eds Franz, Kohl, Sachs), Remme, Nash and Hunter, 2005.
 *
 *  W = Sum_{M,N=1..3} k_{MN}  E_{MN}^2 / (a_{MN} - E_{MN})^b_{MN}
 *
 *  This class doesn't set parameter values, see NashHunterPoleZeroLaw for
 *  a derived class which sets cardiac parameter values.
 *
 *  Not isotropic, so inherits directly from AbstractIncompressibleMaterialLaw
 *
 *  Note, by default, the fibre direction is assumed to be THE X-DIRECTION,
 *  and the sheet direction the Y-DIRECTION (ie sheets in the XY plane). Call
 *  SetChangeOfBasisMatrix() before ComputeStressAndStressDerivative(), with
 *  the matrix P = [fibre_vec, sheet_vec, normal_vec] if this is not the case.
 *
 */
template<unsigned DIM>
class PoleZeroMaterialLaw : public AbstractIncompressibleMaterialLaw<DIM>
{
friend class TestMaterialLaws;

private:

    /** Matrix of parameters k. */
    std::vector<std::vector<double> > mK;

    /** Matrix of parameters a. */
    std::vector<std::vector<double> > mA;

    /** Matrix of parameters b. */
    std::vector<std::vector<double> > mB;

    /** Identity matrix. */
    c_matrix<double,DIM,DIM> mIdentity;

protected:

    /**
     * Protected default constructor doing nothing. Just saw inherited classes
     * can be instantiated and THEN set up the parameters
     */
    PoleZeroMaterialLaw();

    /**
     * Set k, a, and b. To be called by the constuctor or a child class
     * Set comments for constructor.
     *
     * @param k  the parameter k
     * @param a  the parameter a
     * @param b  the parameter b
     */
    void SetParameters(std::vector<std::vector<double> > k,
                       std::vector<std::vector<double> > a,
                       std::vector<std::vector<double> > b);

public:

    /**
     * Constructor, taking in parameters k_i, a_i, b_i as matrices.
     * These matrices must be of size DIM-by-DIM and must be symmetric
     *
     * Note: using the k_1..k_6 convention,  k_4 = 2*k[0][1] = 2*k[1][0], etc
     *
     * @param k  the parameter k
     * @param a  the parameter a
     * @param b  the parameter b
     */
     PoleZeroMaterialLaw(std::vector<std::vector<double> > k,
                         std::vector<std::vector<double> > a,
                         std::vector<std::vector<double> > b);

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
                                          c_matrix<double,DIM,DIM>& rInvC,
                                          double                    pressure,
                                          c_matrix<double,DIM,DIM>& rT,
                                          FourthOrderTensor<DIM,DIM,DIM,DIM>&   rDTdE,
                                          bool                      computeDTdE);

    /**
     * @return the pressure corresponding to E=0, ie C=identity.
     */
    double GetZeroStrainPressure();

    /**
     * Scale the dimensional material parameters (ie the K's).
     *
     * @param scaleFactor
     */
    void ScaleMaterialParameters(double scaleFactor);
};

#endif /*POLEZEROMATERIALLAW_HPP_*/
