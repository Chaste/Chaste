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

#ifndef _BASISFUNCTIONSCHECKERS_HPP_
#define _BASISFUNCTIONSCHECKERS_HPP_

#include <cxxtest/TestSuite.h>
#include "LinearBasisFunction.hpp"
#include "QuadraticBasisFunction.hpp"
#include <vector>

/*
 * Test whether the basis functions are correctly equal to 1 or 0 at the
 * nodes of the canonical element.
 *
 * NOTE: This is a badly written test, because it assumes evaluationPoints
 * are the nodes in the canonical element, in the correct order.
 */
template <unsigned ELEMENT_DIM>
class BasisFunctionsCheckers
{
public:

    void checkLinearBasisFunctions(std::vector<ChastePoint<ELEMENT_DIM>*> evaluationPoints)
    {
        unsigned size = evaluationPoints.size(); // number of evalutation points and basis functions too

        assert(size==ELEMENT_DIM+1);

        std::vector<double> basis_function_vector; // store results of evalutation

        double expected_evaluation;

        for (unsigned point_index=0; point_index<size; point_index++)
        {
            c_vector<double, ELEMENT_DIM+1> basis_function_vector;
            LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunctions(*(evaluationPoints[point_index]), basis_function_vector);

            for (unsigned func_index=0; func_index<size; func_index ++)
            {
                if (func_index == point_index)
                {
                    expected_evaluation = 1.0;
                }
                else
                {
                    expected_evaluation = 0.0;
                }
                TS_ASSERT_DELTA(basis_function_vector(func_index), expected_evaluation, 1e-12);

                TS_ASSERT_DELTA(LinearBasisFunction<ELEMENT_DIM>::ComputeBasisFunction(*(evaluationPoints[point_index]),func_index),
                                expected_evaluation,
                                1e-12);
            }
        }
    }

    void checkQuadraticBasisFunctions(std::vector<ChastePoint<ELEMENT_DIM>*> evaluationPoints)
    {
        unsigned size = evaluationPoints.size(); // number of evalutation points and basis functions too

        assert(size==(ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2);

        std::vector<double> basis_function_vector; // store results of evalutation

        double expected_evaluation;

        for (unsigned point_index=0; point_index<size; point_index++)
        {
            c_vector<double, (ELEMENT_DIM+1)*(ELEMENT_DIM+2)/2> basis_function_vector;
            QuadraticBasisFunction<ELEMENT_DIM>::ComputeBasisFunctions(*(evaluationPoints[point_index]), basis_function_vector);

            for (unsigned func_index=0; func_index<size; func_index ++)
            {
                if (func_index == point_index)
                {
                    expected_evaluation = 1.0;
                }
                else
                {
                    expected_evaluation = 0.0;
                }
                TS_ASSERT_DELTA(basis_function_vector(func_index), expected_evaluation, 1e-12);

                TS_ASSERT_DELTA(QuadraticBasisFunction<ELEMENT_DIM>::ComputeBasisFunction(*(evaluationPoints[point_index]),func_index),
                                expected_evaluation,
                                1e-12);
            }
        }
    }
};

#endif //_BASISFUNCTIONSCHECKERS_HPP_
