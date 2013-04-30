/*

Copyright (c) 2005-2013, University of Oxford.
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

c_vector<double,3> CalculateEigenvectorForSmallestNonzeroEigenvalue(c_matrix<double, 3, 3>& rA)
{
    //Check for symmetry
    if (norm_inf( rA - trans(rA)) > 10*DBL_EPSILON)
    {
        EXCEPTION("Matrix should be symmetric");
    }
    PetscBLASInt info;
    c_vector<PetscReal, 3> eigenvalues_real_part;
    c_vector<PetscReal, 3> eigenvalues_imaginary_part;
    c_vector<PetscScalar, 4*3 > workspace;
    c_matrix<PetscScalar, 3, 3> right_eigenvalues;

    char dont_compute_left_evectors = 'N';
    char compute_right_evectors = 'V';

    PetscBLASInt matrix_size = 3;
    PetscBLASInt matrix_ld = matrix_size;
    PetscBLASInt workspace_size = 4*matrix_size;

    c_matrix<PetscScalar, 3, 3> a_transpose;
    noalias(a_transpose) = trans(rA);

    // PETSc alias for dgeev or dgeev_
    LAPACKgeev_(&dont_compute_left_evectors, &compute_right_evectors,
                &matrix_size, a_transpose.data(), &matrix_ld,
                eigenvalues_real_part.data(), eigenvalues_imaginary_part.data(),
                NULL, &matrix_ld,
                right_eigenvalues.data(), &matrix_ld,
                workspace.data(), &workspace_size,
                &info);
    assert(info==0);

    // If this fails a complex eigenvalue was found
    assert(norm_2(eigenvalues_imaginary_part) < DBL_EPSILON);

    unsigned index_of_smallest = UINT_MAX;
    double min_eigenvalue = DBL_MAX;

    for (unsigned i=0; i<3; i++)
    {
        double eigen_magnitude = fabs(eigenvalues_real_part(i));
        if (eigen_magnitude < min_eigenvalue && eigen_magnitude >= DBL_EPSILON)
        {
            // A zero eigenvalue is ignored
            min_eigenvalue = eigen_magnitude;
            index_of_smallest = i;
        }
        //Check for positive semi-definite
        assert(eigenvalues_real_part(i) > -DBL_EPSILON);
    }
    assert (min_eigenvalue != DBL_MAX);
    assert (index_of_smallest != UINT_MAX);
    assert (min_eigenvalue >= DBL_EPSILON);

    c_vector<double, 3> output;
    output(0) = right_eigenvalues(index_of_smallest, 0);
    output(1) = right_eigenvalues(index_of_smallest, 1);
    output(2) = right_eigenvalues(index_of_smallest, 2);

    //--- AT THIS POINT WE COULD RETURN THE EIGENVECTOR. ---
    // Find the eigenvector by brute-force power method.
    // We can't use the inverse method, because the matrix might be singular

    c_matrix<double,3,3> copy_A(rA);
    //Eigenvalue 1
    c_vector<double, 3> eigenvec1 = scalar_vector<double>(3, 1.0);
    double eigen, norm;
    //Eigenvector
    eigen=DBL_MAX;
    norm=0.0;
    while (fabs(eigen - norm) >DBL_EPSILON) //Machine precision
    {
        eigen = norm;
        eigenvec1 = prod(copy_A, eigenvec1);
        norm = norm_2(eigenvec1);
        eigenvec1 /= norm;
    }
    double eigen1 = eigen;

    // Take out maximum eigenpair
    c_matrix<double, 3, 3> wielandt_reduce_first_vector = identity_matrix<double>(3,3);
    wielandt_reduce_first_vector -= outer_prod(eigenvec1, eigenvec1);
    copy_A = prod(wielandt_reduce_first_vector, copy_A);

    //Eigenvalue 2
    //Eigenvector
    eigen=DBL_MAX;
    norm=0.0;
    c_vector<double, 3> eigenvec2 = scalar_vector<double>(3, 1.0);
    while (fabs(eigen - norm) >DBL_EPSILON) //Machine precision
    {
        eigen = norm;
        eigenvec2 = prod(copy_A, eigenvec2);
        norm = norm_2(eigenvec2);
        eigenvec2 /= norm;
    }
    double eigen2 = eigen;

    // Take out maximum eigenpair
    c_matrix<double, 3, 3> wielandt_reduce_second_vector = identity_matrix<double>(3,3);
    wielandt_reduce_second_vector -= outer_prod(eigenvec2, eigenvec2);
    copy_A = prod(wielandt_reduce_second_vector, copy_A);
    //Eigenvalue 3
    //Eigenvector
    eigen=DBL_MAX;
    norm=0.0;
    c_vector<double, 3> eigenvec3 = scalar_vector<double>(3, 1.0);
    while (fabs(eigen - norm) >DBL_EPSILON) //Machine precision
    {
        eigen = norm;
        eigenvec3 = prod(copy_A, eigenvec3);
        norm = norm_2(eigenvec3);
        eigenvec3 /= norm;
    }
    //Check that we can do as well as LAPACK geev
    if (eigen >= DBL_EPSILON)
    {
        assert(CompareDoubles::WithinAbsoluteTolerance(eigen, min_eigenvalue, 4*DBL_EPSILON));
    }
    else if (eigen2>=DBL_EPSILON)
    {
        assert(CompareDoubles::WithinAbsoluteTolerance(eigen2, min_eigenvalue, 4*DBL_EPSILON));
    }
    else
    {
        assert(CompareDoubles::WithinAbsoluteTolerance(eigen1, min_eigenvalue, 4*DBL_EPSILON));
    }

    return output;
}
