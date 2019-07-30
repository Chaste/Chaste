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

#include "UblasCustomFunctions.hpp"

c_vector<double, 1> Create_c_vector(double x)
{
    c_vector<double, 1> v;
    v[0] = x;
    return v;
}

c_vector<double, 2> Create_c_vector(double x, double y)
{
    c_vector<double, 2> v;
    v[0] = x;
    v[1] = y;
    return v;
}

c_vector<double, 3> Create_c_vector(double x, double y, double z)
{
    c_vector<double, 3> v;
    v[0] = x;
    v[1] = y;
    v[2] = z;
    return v;
}

c_vector<double, 3> CalculateEigenvectorForSmallestNonzeroEigenvalue(c_matrix<double, 3, 3>& rA)
{
    //Check for symmetry
    if (norm_inf(rA - trans(rA)) > 10 * DBL_EPSILON)
    {
        EXCEPTION("Matrix should be symmetric");
    }

    // Find the eigenvector by brute-force using the power method.
    // We can't use the inverse method, because the matrix might be singular

    c_matrix<double, 3, 3> copy_A(rA);
    //Eigenvalue 1
    c_vector<double, 3> eigenvec1 = scalar_vector<double>(3, 1.0);

    double eigen1 = CalculateMaxEigenpair(copy_A, eigenvec1);

    // Take out maximum eigenpair
    c_matrix<double, 3, 3> wielandt_reduce_first_vector = identity_matrix<double>(3, 3);
    wielandt_reduce_first_vector -= outer_prod(eigenvec1, eigenvec1);
    copy_A = prod(wielandt_reduce_first_vector, copy_A);

    c_vector<double, 3> eigenvec2 = scalar_vector<double>(3, 1.0);
    double eigen2 = CalculateMaxEigenpair(copy_A, eigenvec2);

    // Take out maximum (second) eigenpair
    c_matrix<double, 3, 3> wielandt_reduce_second_vector = identity_matrix<double>(3, 3);
    wielandt_reduce_second_vector -= outer_prod(eigenvec2, eigenvec2);
    copy_A = prod(wielandt_reduce_second_vector, copy_A);

    c_vector<double, 3> eigenvec3 = scalar_vector<double>(3, 1.0);
    double eigen3 = CalculateMaxEigenpair(copy_A, eigenvec3);

    //Look backwards through the eigenvalues, checking that they are non-zero
    if (eigen3 >= DBL_EPSILON)
    {
        return eigenvec3;
    }
    if (eigen2 >= DBL_EPSILON)
    {
        return eigenvec2;
    }
    UNUSED_OPT(eigen1);
    assert(eigen1 > DBL_EPSILON);
    return eigenvec1;
}

double CalculateMaxEigenpair(c_matrix<double, 3, 3>& rA, c_vector<double, 3>& rEigenvector)
{
    double norm = 0.0;
    double step = DBL_MAX;
    while (step > DBL_EPSILON) //Machine precision
    {
        c_vector<double, 3> old_value(rEigenvector);
        rEigenvector = prod(rA, rEigenvector);
        norm = norm_2(rEigenvector);
        rEigenvector /= norm;
        if (norm < DBL_EPSILON)
        {
            //We don't care about a zero eigenvector, so don't polish it
            break;
        }
        step = norm_inf(rEigenvector - old_value);
    }
    return norm;
}
