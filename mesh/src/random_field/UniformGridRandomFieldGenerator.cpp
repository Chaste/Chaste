/*

Copyright (c) 2005-2018, University of Oxford.
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

#include "UniformGridRandomFieldGenerator.hpp"

#include <iomanip>
#include <fstream>
#include <numeric>

// Spectra includes (for the eigenvalue and eigenvector calculations)
//#include <SymEigsSolver.h>
//#include <MatOp/SparseGenMatProd.h>

#include "Exception.hpp"
#include "FileFinder.hpp"
#include "Node.hpp"
#include "OutputFileHandler.hpp"
#include "RandomNumberGenerator.hpp"
#include "Warnings.hpp"

#include "Debug.hpp"

template <unsigned SPACE_DIM>
UniformGridRandomFieldGenerator<SPACE_DIM>::UniformGridRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                                                            std::array<double, SPACE_DIM> upperCorner,
                                                                            std::array<unsigned, SPACE_DIM> numGridPts,
                                                                            std::array<bool, SPACE_DIM> periodicity,
                                                                            double traceProportion,
                                                                            double lengthScale)
        : mLowerCorner(lowerCorner),
          mUpperCorner(upperCorner),
          mNumGridPts(numGridPts),
          mPeriodicity(periodicity),
          mTraceProportion(traceProportion),
          mLengthScale(lengthScale),
          mCacheDir("CachedRandomFields/")
{
    // Calculate the total number of grid points
    mNumTotalGridPts = std::accumulate(mNumGridPts.begin(), mNumGridPts.end(), 1u, std::multiplies<unsigned>());

    // Check parameters are sensible
    {
        for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
        {
            assert(mLowerCorner[dim] < mUpperCorner[dim]);
            assert(mNumGridPts[dim] > 0);
        }

        assert(mTraceProportion > 0.0);
        assert(mTraceProportion < 1.0);
        assert(mLengthScale > 0.0);
    }

    // Calculate the grid spacings
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        mGridSpacing[dim] = (mUpperCorner[dim] - mLowerCorner[dim]) / mNumGridPts[dim];
        mOneOverGridSpacing[dim] = mNumGridPts[dim] / (mUpperCorner[dim] - mLowerCorner[dim]);
    }
}


template <unsigned SPACE_DIM>
void UniformGridRandomFieldGenerator<SPACE_DIM>::CalculateEigenDecomposition()
{
    // Eigen::SparseMatrix<double> cov_matrix = CalculateCovarianceMatrix();

    // // mTraceProportion * mNumTotalGridPts is an upper bound for the number of eigenvalues to calculate
    // const double trace_threshold = mTraceProportion * mNumTotalGridPts;

    // // Calculate incremental proportions of the eigenvalues
    // // \todo: refine how this is done. There should be a way to provide a good residual vec from previous calculation
    // for (const double& proportion : {0.2, 0.4, 0.6, 0.8, 1.0})
    // {
    //     const long num_evals_to_compute = std::min(std::lround(proportion * trace_threshold) + 1l, mNumTotalGridPts - 1l);

    //     Spectra::SparseGenMatProd<double> op(cov_matrix);
    //     Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double>> eigs(
    //             &op, num_evals_to_compute, mNumTotalGridPts);

    //     eigs.init();
    //     eigs.compute();

    //     if (eigs.info() != Spectra::SUCCESSFUL)
    //     {
    //         EXCEPTION("Spectra decomposition was not successful.");
    //     }

    //     PRINT_2_VARIABLES(eigs.eigenvalues().array().sum(), trace_threshold);

    //     // Go straight back to the top of the loop at this point if we didn't calculate enough eigenvalues
    //     // but only continue if we haven't already calculated the maximum number of eigenvalues
    //     if (eigs.eigenvalues().array().sum() < trace_threshold && num_evals_to_compute < mNumTotalGridPts - 1l)
    //     {
    //         continue;
    //     }

    //     /*
    //      * The .max(0.0) here ensure eigenvalues are not negative when the square root is taken. This can only happen
    //      * due to rounding error in the calculation, so any negative eigenvalues will have negligible magnitude, so we
    //      * don't mind the slight error.
    //      */
    //     const Eigen::ArrayXd& evals = eigs.eigenvalues().array().max(0.0);

    //     // Decide how many eigenvalues to keep
    //     double cumulative_sum = 0.0;
    //     mNumEigenvals = num_evals_to_compute;
    //     for (unsigned e_val_idx = 0; e_val_idx < evals.size(); ++e_val_idx)
    //     {
    //         cumulative_sum += evals.coeffRef(e_val_idx);

    //         if (cumulative_sum > trace_threshold && e_val_idx < mNumEigenvals)
    //         {
    //             mNumEigenvals = e_val_idx + 1;
    //             break;
    //         }
    //     }

    //     // Store only those eigenvectors that we need
    //     const auto e_vectors = eigs.eigenvectors();
    //     mScaledEigenvecs.resize(mNumTotalGridPts, mNumEigenvals);
    //     for (unsigned e_val = 0; e_val < mNumEigenvals; ++e_val)
    //     {
    //         mScaledEigenvecs.col(e_val) = std::sqrt(evals.coeffRef(e_val)) * e_vectors.col(e_val);
    //     }

    //     PRINT_VARIABLE(mNumEigenvals);

    //     // Break out of the loop to prevent any unnecessary eigenvalue calculations
    //     break;
    // }
}

template <unsigned SPACE_DIM>
Eigen::SparseMatrix<double> UniformGridRandomFieldGenerator<SPACE_DIM>::CalculateCovarianceMatrix() const noexcept
{
    // Create a node at each grid location
    std::vector<Node<SPACE_DIM>> nodes;

    // \todo: can probably remove this if I think about it...
    switch(SPACE_DIM)
    {
        case 1:
        {
            unsigned idx = 0u;
            for (unsigned x = 0; x < mNumGridPts[0]; ++x)
            {
                const double x_pos = mGridSpacing[0] * x;
                nodes.emplace_back(Node<SPACE_DIM>(idx, false, x_pos));
                idx++;
            }
            break;
        }
        case 2:
        {
            unsigned idx = 0u;
            for (unsigned y = 0; y < mNumGridPts[1]; ++y)
            {
                const double y_pos = mGridSpacing[1] * y;
                for (unsigned x = 0; x < mNumGridPts[0]; ++x)
                {
                    const double x_pos = mGridSpacing[0] * x;
                    nodes.emplace_back(Node<SPACE_DIM>(idx, false, x_pos, y_pos));
                    idx++;
                }
            }
            break;
        }
        case 3:
        {
            unsigned idx = 0u;
            for (unsigned z = 0; z < mNumGridPts[2]; ++z)
            {
                const double z_pos = mGridSpacing[2] * z;
                for (unsigned y = 0; y < mNumGridPts[1]; ++y)
                {
                    const double y_pos = mGridSpacing[1] * y;
                    for (unsigned x = 0; x < mNumGridPts[0]; ++x)
                    {
                        const double x_pos = mGridSpacing[0] * x;
                        nodes.emplace_back(Node<SPACE_DIM>(idx, false, x_pos, y_pos, z_pos));
                        idx++;
                    }
                }
            }
            break;
        }
        default:
        NEVER_REACHED;
    }

    // Create vector of eigen triplets representing the pairwise covariance using the Gaussian covariance function
    const double minus_one_over_len_sq = -1.0 / (mLengthScale * mLengthScale);
    std::vector<Eigen::Triplet<double>> triplets;
    for (unsigned x = 0; x < mNumTotalGridPts; ++x)
    {
        for (unsigned y = 0; y < x; ++y)
        {
            const double val = std::exp(minus_one_over_len_sq * GetSquaredDistAtoB(nodes[x].rGetLocation(), nodes[y].rGetLocation()));

            if (val > 1e-9)  // \todo: remove this magic number?
            {
                // Do (x,y) and (y,x) as the distance is symmetric
                triplets.emplace_back(Eigen::Triplet<double>(x, y, val));
                triplets.emplace_back(Eigen::Triplet<double>(y, x, val));
            }
        }
    }

    // Efficiently take case of the diagonal terms
    for (unsigned x = 0; x < mNumTotalGridPts; ++x)
    {
        triplets.emplace_back(Eigen::Triplet<double>(x, x, 1.0));
    }

    // Create a sparse matrix given the grid points
    Eigen::SparseMatrix<double> cov_matrix(mNumTotalGridPts, mNumTotalGridPts);
    cov_matrix.setFromTriplets(triplets.begin(), triplets.end());

    const double occupancy_proportion = static_cast<double>(cov_matrix.nonZeros()) / (mNumTotalGridPts * mNumTotalGridPts);
    PRINT_VARIABLE(occupancy_proportion);

    return cov_matrix;
}

template <unsigned SPACE_DIM>
double UniformGridRandomFieldGenerator<SPACE_DIM>::GetSquaredDistAtoB(
    const c_vector<double, SPACE_DIM>& rLocation1,
    const c_vector<double, SPACE_DIM>& rLocation2) const noexcept
{
    double dist_squared = 0.0;

    // Loop over dimensions
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        // Get the (non-periodic) absolute difference
        double delta_this_dim = std::fabs(rLocation2[dim] - rLocation1[dim]);

        // If this dimension is periodic, recalculate the distance if necessary
        if (mPeriodicity[dim])
        {
            const double domain_width_this_dim = mUpperCorner[dim] - mLowerCorner[dim];

            delta_this_dim = std::min(delta_this_dim, domain_width_this_dim - delta_this_dim);
        }

        dist_squared += delta_this_dim * delta_this_dim;
    }

    return dist_squared;
}

template <unsigned SPACE_DIM>
std::string UniformGridRandomFieldGenerator<SPACE_DIM>::GetFilenameFromParams() const noexcept
{
    std::stringstream name;
    name << std::fixed << std::setprecision(3);

    switch(SPACE_DIM)
    {
        case 1:
        {
            name << "x_"
                 << mLowerCorner[0] << "_"
                 << mUpperCorner[0] << "_"
                 << mNumGridPts[0] << "_"
                 << mPeriodicity[0] << "_"
                 << mTraceProportion << "_"
                 << mLengthScale;
            break;
        }
        case 2:
        {
            name << "xy_"
                 << mLowerCorner[0] << "_" << mLowerCorner[1] << "_"
                 << mUpperCorner[0] << "_" << mUpperCorner[1] << "_"
                 << mNumGridPts[0] << "_" << mNumGridPts[1] << "_"
                 << mPeriodicity[0] << "_" << mPeriodicity[1] << "_"
                 << mTraceProportion << "_"
                 << mLengthScale;
            break;
        }
        case 3:
        {
            name << "xyz_"
                 << mLowerCorner[0] << "_" << mLowerCorner[1] << "_" << mLowerCorner[2] << "_"
                 << mUpperCorner[0] << "_" << mUpperCorner[1] << "_" << mUpperCorner[2] << "_"
                 << mNumGridPts[0] << "_" << mNumGridPts[1] << "_" << mNumGridPts[2] << "_"
                 << mPeriodicity[0] << "_" << mPeriodicity[1] << "_" << mPeriodicity[2] << "_"
                 << mTraceProportion << "_"
                 << mLengthScale;
            break;
        }
        default:
            NEVER_REACHED;
    }

    // Add file extension
    name << ".rfg";
    return name.str();
}

template <unsigned SPACE_DIM>
std::vector<double> UniformGridRandomFieldGenerator<SPACE_DIM>::SampleRandomField() const noexcept
{ 
    // TODO: randomise seed
    OpenSimplex2S os(123);
    auto reshape = [](const double val) {
        double distFromHalf = 2.0 * (0.5 - std::abs(0.5 - std::abs(val)));
        double strength = (1.0 - std::abs(val)) + distFromHalf * 0.2;
        return val * (1.0 - 0.18 * strength);
    };

    std::vector<double> samples(mNumTotalGridPts);
    for (unsigned int x = 0; x < mNumGridPts[0]; x++) {
        for (unsigned int y = 0; y < mNumGridPts[1]; y++) {
            // TODO: Better rng for time/continuous time
            samples[mNumGridPts[1] * y + x] = reshape(reshape(os.noise3_XYBeforeZ(x * mLengthScale, y * mLengthScale, rand())));
        }
    }
    
    return samples;
    
    /*// Generate mNumTotalGridPts normally-distributed random numbers
    std::vector<double> samples_from_n01(mNumTotalGridPts);
    for (unsigned i = 0; i < samples_from_n01.size(); ++i)
    {
        samples_from_n01[i] = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
    }

    // Generate the instance of the random field
    Eigen::VectorXd grf(mNumTotalGridPts);
    grf.setZero();
    for (unsigned j = 0; j < mScaledEigenvecs.cols(); ++j)
    {
        grf += samples_from_n01[j] * mScaledEigenvecs.col(j);
    }

    // Translate to a std::vector so that eigen objects aren't leaking out to other places in Chaste
    return std::vector<double>(grf.data(), grf.data() + grf.size());*/
}

template <unsigned SPACE_DIM>
double UniformGridRandomFieldGenerator<SPACE_DIM>::Interpolate(const std::vector<double>& rRandomField,
                                                               const c_vector<double, SPACE_DIM>& rLocation) const noexcept
{
    assert(mNumTotalGridPts == rRandomField.size());

    // Find the nearest node
    std::array<long, SPACE_DIM> lower_left;

    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        lower_left[dim] = std::floor((rLocation[dim] - mLowerCorner[dim]) / mGridSpacing[dim]);

        if (lower_left[dim] < 0)
        {
            lower_left[dim] = 0;
            WARN_ONCE_ONLY("Interpolating outside random field grid: does the random field need to be larger?");
        }
        else if (lower_left[dim] >= mNumGridPts[dim])
        {
            lower_left[dim] = mNumGridPts[dim] - 1;
            WARN_ONCE_ONLY("Interpolating outside random field grid: does the random field need to be larger?");
        }
    }

    double interpolated_value{};

    // Perform the interpolation
    switch(SPACE_DIM)
    {
        case 1:
        {
            // The value of the field at points either side of rLocation
            std::array<long, SPACE_DIM> upper_idx = lower_left;
            upper_idx[0] =  (upper_idx[0] + 1) % mNumGridPts[0];

            const double field_at_lower = rRandomField[GetLinearIndex(lower_left)];
            const double field_at_upper = rRandomField[GetLinearIndex(upper_idx)];

            // Perform a simple linear interpolation
            const std::array<double, SPACE_DIM> lower_location = GetPositionUsingGridIndex(lower_left);
            const double interpolant = (rLocation[0] - lower_location[0]) * mOneOverGridSpacing[0];

            interpolated_value = field_at_lower * (1.0 - interpolant) + field_at_upper * interpolant;
            break;
        }
        case 2:
        {
            // The value of the field at the four points of the square containing rLocation
            std::array<long, SPACE_DIM> x_upper = lower_left;
            x_upper[0] = (x_upper[0] + 1) % mNumGridPts[0];

            std::array<long, SPACE_DIM> y_upper = lower_left;
            y_upper[1] = (y_upper[1] + 1) % mNumGridPts[1];

            std::array<long, SPACE_DIM> xy_upper = x_upper;
            xy_upper[1] = y_upper[1];

            const double field_xy_lower = rRandomField[GetLinearIndex(lower_left)];
            const double field_x_upper = rRandomField[GetLinearIndex(x_upper)];
            const double field_y_upper = rRandomField[GetLinearIndex(y_upper)];
            const double field_xy_upper = rRandomField[GetLinearIndex(xy_upper)];

            // Perform a simple bilinear interpolation
            const std::array<double, SPACE_DIM> lower_location = GetPositionUsingGridIndex(lower_left);
            const double dist_x_lower = rLocation[0] - lower_location[0];
            const double dist_x_upper = mGridSpacing[0] - dist_x_lower;
            const double dist_y_lower = rLocation[1] - lower_location[1];
            const double dist_y_upper = mGridSpacing[1] - dist_y_lower;

            interpolated_value = mOneOverGridSpacing[0] * mOneOverGridSpacing[1] * (field_xy_lower * dist_x_upper * dist_y_upper +
                                                                                    field_x_upper * dist_x_lower * dist_y_upper +
                                                                                    field_y_upper * dist_x_upper * dist_y_lower +
                                                                                    field_xy_upper * dist_x_lower * dist_y_lower);

            break;
        }
        case 3:
        {
            // \todo: perform a suitable interpolation rather than just near neighbour
            interpolated_value = rRandomField[GetLinearIndex(lower_left)];

            break;
        }
        default:
            NEVER_REACHED;
    }

    return interpolated_value;
}

template <unsigned SPACE_DIM>
long UniformGridRandomFieldGenerator<SPACE_DIM>::GetLinearIndex(std::array<long, SPACE_DIM> gridIndex) const noexcept
{
    long linear_index;

    // \todo: double check I'm not calculating the transpose of each point
    switch(SPACE_DIM)
    {
        case 1:
        {
            linear_index = gridIndex[0];
            break;
        }
        case 2:
        {
            linear_index = gridIndex[0] +
                           gridIndex[1] * mNumGridPts[0];
            break;
        }
        case 3:
        {
            linear_index = gridIndex[0] +
                           gridIndex[1] * mNumGridPts[0] +
                           gridIndex[2] * mNumGridPts[0] * mNumGridPts[1];
            break;
        }
        default:
            NEVER_REACHED;
    }

    return linear_index;
}

template <unsigned SPACE_DIM>
std::array<double, SPACE_DIM> UniformGridRandomFieldGenerator<SPACE_DIM>::GetPositionUsingGridIndex(std::array<long, SPACE_DIM> gridIndex) const noexcept
{
    std::array<double, SPACE_DIM> position = {{}};

    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        position[dim] = mLowerCorner[dim] + mGridSpacing[dim] * gridIndex[dim];
    }

    return position;
}

// Explicit instantiation
template class UniformGridRandomFieldGenerator<1>;
template class UniformGridRandomFieldGenerator<2>;
template class UniformGridRandomFieldGenerator<3>;
