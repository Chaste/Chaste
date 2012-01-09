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
     * Get the pressure corresponding to E=0, ie C=identity.
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
