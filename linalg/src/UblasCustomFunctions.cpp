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
    int info;
    c_vector<double, 3> eigenvalues_real_part;
    c_vector<double, 3> eigenvalues_imaginary_part;
    c_vector<double, 4*3 > workspace;
    c_matrix<double, 3, 3> right_eigenvalues;

    char dont_compute_left_evectors = 'N';
    char compute_right_evectors = 'V';

    int matrix_size = 3;
    int matrix_ld = matrix_size;
    int workspace_size = 4*matrix_size;

    c_matrix<double, 3, 3> a_transpose;
    noalias(a_transpose) = trans(rA);

    // PETSc alias for dgeev or dgeev_
    LAPACKgeev_(&dont_compute_left_evectors, &compute_right_evectors,
           &matrix_size, a_transpose.data(),&matrix_ld,
           eigenvalues_real_part.data(), eigenvalues_imaginary_part.data(),
           NULL, &matrix_ld,
           right_eigenvalues.data(),&matrix_ld,
           workspace.data(),&workspace_size,
           &info);
    assert(info==0);

    // If this fails a complex eigenvalue was found
    assert(norm_2(eigenvalues_imaginary_part) < DBL_EPSILON);

    unsigned index_of_smallest=UINT_MAX;
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
    }
    assert (min_eigenvalue != DBL_MAX);
    assert (index_of_smallest != UINT_MAX);
    assert (min_eigenvalue >= DBL_EPSILON);

    c_vector<double, 3> output;
    output(0) = right_eigenvalues(index_of_smallest, 0);
    output(1) = right_eigenvalues(index_of_smallest, 1);
    output(2) = right_eigenvalues(index_of_smallest, 2);

    return output;
}
