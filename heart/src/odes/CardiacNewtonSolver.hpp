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
// LCOV_EXCL_START
                EXCEPTION("Newton method diverged in CardiacNewtonSolver::Solve()");
// LCOV_EXCL_STOP
            }
        }
        while (norm_of_update > eps);

// LCOV_EXCL_START
#ifndef NDEBUG
        if (norm_of_residual > 2e-10)
        { //This line is for correlation - in case we use norm_of_residual as convergence criterion
            WARN_ONCE_ONLY("Newton iteration terminated because update vector norm is small, but residual norm is not small.");
        }
#endif // NDEBUG
// LCOV_EXCL_STOP
    }

protected:
    /** Singleton pattern - protected default constructor. */
    CardiacNewtonSolver()
    {}
    /** Singleton pattern - protected copy constructor.  Not implemented. */
    CardiacNewtonSolver(const CardiacNewtonSolver<SIZE, CELLTYPE>&);
    /** Singleton pattern - protected assignment operator.  Not implemented. @return (would be reference)*/
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
