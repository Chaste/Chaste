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
#ifndef CARDIACNEWTONSOLVER_HPP_
#define CARDIACNEWTONSOLVER_HPP_

#include <cmath>
#include "IsNan.hpp"
#include "UblasCustomFunctions.hpp"
#include "AbstractBackwardEulerCardiacCell.hpp"
#include "Warnings.hpp"

/**
 * Specialised Newton solver for solving the nonlinear systems arising when
 * simulating a cardiac cell using Backward Euler.
 *
 * The class is templated by the size of the nonlinear system, and uses the
 * singleton pattern to ensure only 1 solver for any given system size is created.
 * This allows us to be both computationally and memory efficient.
 *
 * It would be nice to have a test of this class directly, but you need a cardiac
 * cell in order to test it.  So all tests occur when testing particular cardiac
 * cells, e.g. the LuoRudy1991BackwardEuler.
 */
template<unsigned SIZE, typename CELLTYPE>
class CardiacNewtonSolver
{
public:
    /**
     * Call this method to obtain a solver instance.
     *
     * @return a single instance of the class
     */
    static CardiacNewtonSolver<SIZE, CELLTYPE>* Instance()
    {
        static CardiacNewtonSolver<SIZE, CELLTYPE> inst;
        return &inst;
    }

    /**
     * Use Newton's method to solve the given cell for the next timestep.
     *
     * @param rCell  the cell to solve
     * @param time  the current time
     * @param rCurrentGuess  the current guess at a solution.  Will be updated on exit.
     */
    void Solve(CELLTYPE &rCell,
               double time,
               double rCurrentGuess[SIZE])
    {
        unsigned counter = 0;
        const double eps = 1e-6; // JonW tolerance

        // check that the initial guess that was given gives a valid residual
        rCell.ComputeResidual(time, rCurrentGuess, mResidual.data());
        double norm_of_residual = norm_inf(mResidual);
        assert(!std::isnan(norm_of_residual));
        double norm_of_update = 0.0; //Properly initialised in the loop
        do
        {
            // Calculate Jacobian for current guess
            rCell.ComputeJacobian(time, rCurrentGuess, mJacobian);

            // Solve Newton linear system for mUpdate, given mJacobian and mResidual
            SolveLinearSystem();

            // Update norm (JonW style)
            norm_of_update = norm_inf(mUpdate);

            // Update current guess and recalculate residual
            for (unsigned i=0; i<SIZE; i++)
            {
                rCurrentGuess[i] -= mUpdate[i];
            }
            double norm_of_previous_residual = norm_of_residual;
            rCell.ComputeResidual(time, rCurrentGuess, mResidual.data());
            norm_of_residual = norm_inf(mResidual);
            if (norm_of_residual > norm_of_previous_residual && norm_of_update > eps)
            {
                //Second part of guard:
                //Note that if norm_of_update < eps (converged) then it's
                //likely that both the residual and the previous residual were
                //close to the root.

                //Work out where the biggest change in the guess has happened.
                double relative_change_max = 0.0;
                unsigned relative_change_direction = 0;
                for (unsigned i=0; i<SIZE; i++)
                {
                    double relative_change = fabs(mUpdate[i]/rCurrentGuess[i]);
                    if (relative_change > relative_change_max)
                    {
                        relative_change_max = relative_change;
                        relative_change_direction = i;
                    }
                }

                if (relative_change_max > 1.0)
                {
                    //Only walk 0.2 of the way in that direction (put back 0.8)
                    rCurrentGuess[relative_change_direction] += 0.8*mUpdate[relative_change_direction];
                    rCell.ComputeResidual(time, rCurrentGuess, mResidual.data());
                    norm_of_residual = norm_inf(mResidual);
                    WARNING("Residual increasing and one direction changing radically - back tracking in that direction");
                }
            }
            counter++;

            // avoid infinite loops
            if (counter > 15)
            {
#define COVERAGE_IGNORE
                EXCEPTION("Newton method diverged in CardiacNewtonSolver::Solve()");
#undef COVERAGE_IGNORE
            }
        }
        while (norm_of_update > eps);

#define COVERAGE_IGNORE
#ifndef NDEBUG
        if (norm_of_residual > 2e-10)
        { //This line is for correlation - in case we use norm_of_residual as convergence criterion
            WARN_ONCE_ONLY("Newton iteration terminated because update vector norm is small, but residual norm is not small.");
        }
#endif // NDEBUG
#undef COVERAGE_IGNORE
    }

protected:
    /** Singleton pattern - protected default constructor. */
    CardiacNewtonSolver()
    {}
    /** Singleton pattern - protected copy constructor.  Not implemented. */
    CardiacNewtonSolver(const CardiacNewtonSolver<SIZE, CELLTYPE>&);
    /** Singleton pattern - protected assignment operator.  Not implemented. */
    CardiacNewtonSolver<SIZE, CELLTYPE>& operator= (const CardiacNewtonSolver<SIZE, CELLTYPE>&);

    /**
     * Solve a linear system to calculate the Newton update step
     *
     * This is solving
     *  Jacbian . update = residual
     * for update given values of the Jacobian matrix and residual
     *
     * The implementation does Gaussian elimination with no pivotting and no underflow checking
     */
    void SolveLinearSystem()
    {
        for (unsigned i=0; i<SIZE; i++)
        {
            for (unsigned ii=i+1; ii<SIZE; ii++)
            {
                double fact = mJacobian[ii][i]/mJacobian[i][i];
                for (unsigned j=i; j<SIZE; j++)
                {
                    mJacobian[ii][j] -= fact*mJacobian[i][j];
                }
                mResidual[ii] -= fact*mResidual[i];
            }
        }
        for (unsigned i=SIZE; i-- > 0; )
        {
            mUpdate[i] = mResidual[i];
            for (unsigned j=i+1; j<SIZE; j++)
            {
                mUpdate[i] -= mJacobian[i][j]*mUpdate[j];
            }
            mUpdate[i] /= mJacobian[i][i];
        }
    }

private:
    /** Working memory : residual vector */
    c_vector<double, SIZE> mResidual;
    /** Working memory : Jacobian matrix */
    double mJacobian[SIZE][SIZE];
    /** Working memory : update vector */
    c_vector<double, SIZE> mUpdate;
};

#endif /*CARDIACNEWTONSOLVER_HPP_*/
