/*

Copyright (c) 2005-2017, University of Oxford.
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
#include <SymEigsSolver.h>
#include <MatOp/SparseGenMatProd.h>

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
                                                                            unsigned numEigenvals,
                                                                            double lengthScale)
        : mLowerCorner(lowerCorner),
          mUpperCorner(upperCorner),
          mNumGridPts(numGridPts),
          mPeriodicity(periodicity),
          mNumEigenvals(numEigenvals),
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

        assert(mNumEigenvals > 0);
        assert(mNumEigenvals < mNumTotalGridPts);
        assert(mLengthScale > 0.0);
    }

    // Calculate the grid spacings
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        mGridSpacing[dim] = (mUpperCorner[dim] - mLowerCorner[dim]) / mNumGridPts[dim];
    }

    // Check if there's a cached random field matching these parameters
    FileFinder cached_version_file(mCacheDir + GetFilenameFromParams(), RelativeTo::ChasteTestOutput);

    if (cached_version_file.Exists())
    {
        LoadFromCache(cached_version_file.GetAbsolutePath());
    }
    else
    {
        CalculateEigenDecomposition();
    }
}

template <unsigned SPACE_DIM>
UniformGridRandomFieldGenerator<SPACE_DIM>::UniformGridRandomFieldGenerator(const std::string filename)
{
    // First check the cached random field exists
    FileFinder cached_version_file(filename, RelativeTo::ChasteTestOutput);
    EXCEPT_IF_NOT(cached_version_file.Exists());

    LoadFromCache(cached_version_file.GetAbsolutePath());
}

template <unsigned SPACE_DIM>
void UniformGridRandomFieldGenerator<SPACE_DIM>::CalculateEigenDecomposition()
{
    Eigen::SparseMatrix<double> cov_matrix = CalculateCovarianceMatrix();

    Spectra::SparseGenMatProd<double> op(cov_matrix);
    Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double>> eigs(
            &op, mNumEigenvals, std::min(2 * mNumEigenvals, mNumTotalGridPts));

    eigs.init();
    eigs.compute();

    EXCEPT_IF_NOT(eigs.info() == Spectra::SUCCESSFUL);

    mEigenvals = eigs.eigenvalues();
    mEigenvecs = eigs.eigenvectors();
}

template <unsigned SPACE_DIM>
Eigen::SparseMatrix<double> UniformGridRandomFieldGenerator<SPACE_DIM>::CalculateCovarianceMatrix() const noexcept
{
    // Create a node at each grid location
    std::vector<Node<2>> nodes;

    // \todo: can probably remove this if I think about it...
    switch(SPACE_DIM)
    {
        case 1:
        {
            unsigned idx = 0u;
            for (unsigned x = 0; x < mNumGridPts[0]; ++x)
            {
                const double x_pos = mGridSpacing[0] * x;
                nodes.emplace_back(Node<2>(idx, false, x_pos));
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
                    nodes.emplace_back(Node<2>(idx, false, x_pos, y_pos));
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
                        nodes.emplace_back(Node<2>(idx, false, x_pos, y_pos, z_pos));
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
    const double length_squared = mLengthScale * mLengthScale;
    const double tol_cov = -std::log(1e-12);  // \todo: remove this magic number
    std::vector<Eigen::Triplet<double>> triplets;
    for (unsigned x = 0; x < mNumTotalGridPts; ++x)
    {
        for (unsigned y = 0; y < mNumTotalGridPts; ++y)
        {
            const double dist_squared = GetSquaredDistAtoB(nodes[x].rGetLocation(), nodes[y].rGetLocation());
            const double exponent = dist_squared / length_squared;

            if (exponent < tol_cov)
            {
                triplets.emplace_back(Eigen::Triplet<double>(x, y, std::exp(-exponent)));
            }
        }
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
                 << mNumEigenvals << "_"
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
                 << mNumEigenvals << "_"
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
                 << mNumEigenvals << "_"
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
void UniformGridRandomFieldGenerator<SPACE_DIM>::LoadFromCache(const std::string& absoluteFilePath)
{
    RandomFieldCacheHeader<SPACE_DIM> header;

    std::ifstream input_file(absoluteFilePath, std::ios::in | std::ios::binary);
    EXCEPT_IF_NOT(input_file.is_open());

    // Read the header, and populate class parameters
    input_file.read((char*) &header, sizeof(RandomFieldCacheHeader<SPACE_DIM>));
    mLowerCorner = header.mLowerCorner;
    mUpperCorner = header.mUpperCorner;
    mNumGridPts = header.mNumGridPts;
    mPeriodicity = header.mPeriodicity;
    mNumEigenvals = header.mNumEigenvals;
    mLengthScale = header.mLengthScale;

    // Calculate how many grid points there are in total, and resize the eigen data arrays accordingly
    mNumTotalGridPts = std::accumulate(mNumGridPts.begin(), mNumGridPts.end(), 1u, std::multiplies<unsigned>());
    mEigenvals.resize(mNumEigenvals);
    mEigenvecs.resize(mNumTotalGridPts, mNumEigenvals);

    // Calculate the remaining unknown: the grid spacing in each dimension, needed for interpolation
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        mGridSpacing[dim] = (mUpperCorner[dim] - mLowerCorner[dim]) / mNumGridPts[dim];
    }

    // Read the eigenvalues and eigenvectors into their respective data arrays
    input_file.read((char*) mEigenvals.data(), mNumEigenvals * sizeof(double));
    input_file.read((char*) mEigenvecs.data(), mNumTotalGridPts * mNumEigenvals * sizeof(double));

    input_file.close();
}

template <unsigned SPACE_DIM>
std::vector<double> UniformGridRandomFieldGenerator<SPACE_DIM>::SampleRandomField() const noexcept
{
    // Generate mNumTotalGridPts normally-distributed random numbers
    std::vector<double> samples_from_n01(mNumTotalGridPts);
    for (unsigned i = 0; i < samples_from_n01.size(); ++i)
    {
        samples_from_n01[i] = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
    }

    // Generate the instance of the random field
    Eigen::VectorXd grf(mNumTotalGridPts);
    grf.setZero();
    for (unsigned j = 0; j < mNumEigenvals; ++j)
    {
        grf += samples_from_n01[j] * std::sqrt(mEigenvals(j)) * mEigenvecs.col(j);
    }

    // Translate to a std::vector so that eigen objects aren't leaking out to other places in Chaste
    return std::vector<double>(grf.data(), grf.data() + grf.size());
}

template <unsigned SPACE_DIM>
void UniformGridRandomFieldGenerator<SPACE_DIM>::SaveToCache()
{
    // Get the absolute file path to the cached file, given the current parameters
    FileFinder cached_version_file(mCacheDir + GetFilenameFromParams(), RelativeTo::ChasteTestOutput);

    // If the cache does not exist, create it
    if (not cached_version_file.Exists())
    {
        OutputFileHandler results_handler(mCacheDir, false);
        auto results_file = results_handler.OpenOutputFile(GetFilenameFromParams(), std::ios::out | std::ios::binary);

        // Generate the header struct
        RandomFieldCacheHeader<SPACE_DIM> header;
        header.mLowerCorner = mLowerCorner;
        header.mUpperCorner = mUpperCorner;
        header.mNumGridPts = mNumGridPts;
        header.mPeriodicity = mPeriodicity;
        header.mNumEigenvals = mNumEigenvals;
        header.mLengthScale = mLengthScale;

        // Write the information to file
        results_file->write((char *) &header, sizeof(RandomFieldCacheHeader<SPACE_DIM>));
        results_file->write((char *) mEigenvals.data(), mEigenvals.size() * sizeof(double));
        results_file->write((char *) mEigenvecs.data(), mEigenvecs.size() * sizeof(double));
        results_file->close();
    }
}

template <unsigned SPACE_DIM>
double UniformGridRandomFieldGenerator<SPACE_DIM>::Interpolate(const std::vector<double>& rRandomField,
                                                               const c_vector<double, SPACE_DIM>& rLocation) const noexcept
{
    assert(mNumTotalGridPts == rRandomField.size());

    std::array<long, SPACE_DIM> nearest_node;

    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        nearest_node[dim] = std::floor((rLocation[dim] - mLowerCorner[dim]) / mGridSpacing[dim]);

        if (nearest_node[dim] < 0)
        {
            nearest_node[dim] = 0;
            WARNING("Interpolating outside random field grid: does the random field need to be larger?");
        }
        else if (nearest_node[dim] >= mNumGridPts[dim])
        {
            nearest_node[dim] = mNumGridPts[dim] - 1;
            WARNING("Interpolating outside random field grid: does the random field need to be larger?");
        }
    }

    long linear_index;

    // \todo: double check I'm not calculating the transpose of each point
    switch(SPACE_DIM)
    {
        case 1:
        {
            linear_index = nearest_node[0];
            break;
        }
        case 2:
        {
            linear_index = nearest_node[0] +
                           nearest_node[1] * mNumGridPts[0];
            break;
        }
        case 3:
        {
            linear_index = nearest_node[0] +
                           nearest_node[1] * mNumGridPts[0] +
                           nearest_node[2] * mNumGridPts[0] * mNumGridPts[1];
            break;
        }
        default:
            NEVER_REACHED;
    }

    return rRandomField[linear_index];
}

// Explicit instantiation
template class UniformGridRandomFieldGenerator<1>;
template class UniformGridRandomFieldGenerator<2>;
template class UniformGridRandomFieldGenerator<3>;
