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

#include "SimpleNonlinearEllipticSolver.hpp"

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_matrix<double,1*(ELEMENT_DIM+1),1*(ELEMENT_DIM+1)> SimpleNonlinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>::ComputeMatrixTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> ret;

    c_matrix<double, SPACE_DIM, SPACE_DIM> f_of_u = mpNonlinearEllipticPde->ComputeDiffusionTerm(rX, rU(0));
    c_matrix<double, SPACE_DIM, SPACE_DIM> f_of_u_prime = mpNonlinearEllipticPde->ComputeDiffusionTermPrime(rX, rU(0));

    // LinearSourceTerm(x) not needed as it is a constant wrt u
    double forcing_term_prime = mpNonlinearEllipticPde->ComputeNonlinearSourceTermPrime(rX, rU(0));

    // Note rGradU is a 1 by SPACE_DIM matrix, the 1 representing the dimension of
    // u (ie in this problem the unknown is a scalar). r_gradu_0 is rGradU as a vector
    matrix_row< c_matrix<double, 1, SPACE_DIM> > r_gradu_0(rGradU, 0);
    c_vector<double, SPACE_DIM> temp1 = prod(f_of_u_prime, r_gradu_0);
    c_vector<double, ELEMENT_DIM+1> temp1a = prod(temp1, rGradPhi);

    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values1 = outer_prod(temp1a, rPhi);
    c_matrix<double, SPACE_DIM, ELEMENT_DIM+1> temp2 = prod(f_of_u, rGradPhi);
    c_matrix<double, ELEMENT_DIM+1, ELEMENT_DIM+1> integrand_values2 = prod(trans(rGradPhi), temp2);
    c_vector<double, ELEMENT_DIM+1> integrand_values3 = forcing_term_prime * rPhi;

    ret = integrand_values1 + integrand_values2 - outer_prod( scalar_vector<double>(ELEMENT_DIM+1), integrand_values3);

    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
c_vector<double,1*(ELEMENT_DIM+1)> SimpleNonlinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>::ComputeVectorTerm(
        c_vector<double, ELEMENT_DIM+1>& rPhi,
        c_matrix<double, SPACE_DIM, ELEMENT_DIM+1>& rGradPhi,
        ChastePoint<SPACE_DIM>& rX,
        c_vector<double,1>& rU,
        c_matrix<double,1,SPACE_DIM>& rGradU,
        Element<ELEMENT_DIM,SPACE_DIM>* pElement)
{
    c_vector<double, 1*(ELEMENT_DIM+1)> ret;

    // For solving an AbstractNonlinearEllipticEquation
    // d/dx [f(U,x) du/dx ] = -g
    // where g(x,U) is the forcing term
    double forcing_term = mpNonlinearEllipticPde->ComputeLinearSourceTerm(rX);
    forcing_term += mpNonlinearEllipticPde->ComputeNonlinearSourceTerm(rX, rU(0));

    c_matrix<double, ELEMENT_DIM, ELEMENT_DIM> FOfU = mpNonlinearEllipticPde->ComputeDiffusionTerm(rX, rU(0));

    // Note rGradU is a 1 by SPACE_DIM matrix, the 1 representing the dimension of
    // u (ie in this problem the unknown is a scalar). rGradU0 is rGradU as a vector.
    matrix_row< c_matrix<double, 1, SPACE_DIM> > rGradU0(rGradU, 0);
    c_vector<double, ELEMENT_DIM+1> integrand_values1 =
        prod(c_vector<double, ELEMENT_DIM>(prod(rGradU0, FOfU)), rGradPhi);

    ret = integrand_values1 - (forcing_term * rPhi);
    return ret;
}

template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
SimpleNonlinearEllipticSolver<ELEMENT_DIM,SPACE_DIM>::SimpleNonlinearEllipticSolver(
                              AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>* pMesh,
                              AbstractNonlinearEllipticPde<SPACE_DIM>* pPde,
                              BoundaryConditionsContainer<ELEMENT_DIM, SPACE_DIM, 1>* pBoundaryConditions)
    :  AbstractNonlinearAssemblerSolverHybrid<ELEMENT_DIM,SPACE_DIM,1>(pMesh,pBoundaryConditions),
       mpNonlinearEllipticPde(pPde)
{
    assert(pPde!=nullptr);
}

// Explicit instantiation
template class SimpleNonlinearEllipticSolver<1,1>;
template class SimpleNonlinearEllipticSolver<2,2>;
template class SimpleNonlinearEllipticSolver<3,3>;
