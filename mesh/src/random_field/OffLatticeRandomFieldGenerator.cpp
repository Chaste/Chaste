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

#include "OffLatticeRandomFieldGenerator.hpp"

#include <algorithm>
#include <iomanip>
#include <fstream>
#include <numeric>

// Spectra includes (for the eigenvalue and eigenvector calculations)
#include <MatOp/SparseGenMatProd.h>
#include <SymEigsSolver.h>

#include "ChasteMakeUnique.hpp"
#include "Exception.hpp"
#include "RandomNumberGenerator.hpp"


template <unsigned SPACE_DIM>
OffLatticeRandomFieldGenerator<SPACE_DIM>::OffLatticeRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                                                          std::array<double, SPACE_DIM> upperCorner,
                                                                          std::array<bool, SPACE_DIM> periodicity,
                                                                          unsigned numEigenvals,
                                                                          double lengthScale,
                                                                          double boxWidth)
        : mLowerCorner(lowerCorner),
          mUpperCorner(upperCorner),
          mPeriodicity(periodicity),
          mNumEigenvals(numEigenvals),
          mLengthScale(lengthScale),
          mpBoxCollection(nullptr)
{
    // Reset the box width if the default value is being used
    if (boxWidth == DOUBLE_UNSET)
    {
        // \todo: can we remove the magic number here?
        boxWidth = std::min(4.0 * lengthScale, upperCorner[0] - lowerCorner[0]);
    }

    // Set up the box collection
    c_vector<double, 2 * SPACE_DIM> domain_size;
    for (unsigned dim = 0; dim < SPACE_DIM; ++dim)
    {
        domain_size[2 * dim] = lowerCorner[dim];
        domain_size[2 * dim + 1] = upperCorner[dim];
    }

    // \todo: Obsolete box collection should take a std::array<bool, DIM>
    bool periodic_x = periodicity[0];
    bool periodic_y = SPACE_DIM > 1 ? periodicity[1] : false;
    bool periodic_z = SPACE_DIM > 2 ? periodicity[2] : false;

    mpBoxCollection = our::make_unique<ObsoleteBoxCollection<SPACE_DIM>>(
            boxWidth,
            domain_size,
            periodic_x,
            periodic_y,
            periodic_z
    );

    mpBoxCollection->SetupLocalBoxesHalfOnly();
}

template <unsigned SPACE_DIM>
void OffLatticeRandomFieldGenerator<SPACE_DIM>::Update(const std::vector<Node<SPACE_DIM>*>& rNodes)
{
    // Special case of zero correlation length
    if (mLengthScale == 0.0)
    {
        return;
    }

    Eigen::SparseMatrix<double> cov_matrix = CalculateCovarianceMatrix(rNodes);

    Spectra::SparseGenMatProd<double> op(cov_matrix);
    Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double>> eigs(
            &op, mNumEigenvals, rNodes.size());

    eigs.init();
    eigs.compute();

    if (eigs.info() != Spectra::SUCCESSFUL)
    {
        EXCEPTION("Spectra decomposition not completed successfully.");
    }

    mSqrtEigenvals = eigs.eigenvalues().array().sqrt();
    mEigenvecs = eigs.eigenvectors();
}

template <unsigned SPACE_DIM>
void OffLatticeRandomFieldGenerator<SPACE_DIM>::TuneNumEigenvals(
        const std::vector<Node<SPACE_DIM>*>& rNodes, const double proportionOfTrace)
{
    // Special case of zero correlation length
    if (mLengthScale == 0.0)
    {
        return;
    }

    const unsigned most_eigenvals = rNodes.size() - 1;

    Eigen::SparseMatrix<double> cov_matrix = CalculateCovarianceMatrix(rNodes);

    Spectra::SparseGenMatProd<double> op(cov_matrix);
    Spectra::SymEigsSolver<double, Spectra::LARGEST_MAGN, Spectra::SparseGenMatProd<double>> eigs(
            &op, most_eigenvals, rNodes.size());

    eigs.init();
    eigs.compute();

    if (eigs.info() != Spectra::SUCCESSFUL)
    {
        EXCEPTION("Spectra decomposition not completed successfully.");
    }

    // Calculate the number of eigenvalues at which we reach the trace threshold.
    // The trace of the covariance matrix, which is also the sum of the eigenvalues, is just 1 * rNodes.size()
    const double trace_threshold = proportionOfTrace * rNodes.size();

    mNumEigenvals = UINT_MAX;

    const Eigen::ArrayXd& e_vals = eigs.eigenvalues().array();
    double cumulative_sum = 0.0;
    for (unsigned e_val_idx = 0; e_val_idx < e_vals.size(); ++e_val_idx)
    {
        cumulative_sum += e_vals.coeffRef(e_val_idx);

        if (cumulative_sum > trace_threshold)
        {
            mNumEigenvals = e_val_idx + 1;
            break;
        }
    }

    if (mNumEigenvals + 1 > rNodes.size())
    {
        EXCEPTION("Something has gone wrong tuning mNumEigenvals.  Check 0 < proportionOfTrace < 1.");
    }

    std::cout << "Random field generator: tuned to "
              << mNumEigenvals
              << " eigenvalues ("
              << (100.0 * mNumEigenvals) / rNodes.size()
              << "%) at "
              << 100.0 * proportionOfTrace
              << "% of trace.\n";
}

template <unsigned SPACE_DIM>
Eigen::SparseMatrix<double> OffLatticeRandomFieldGenerator<SPACE_DIM>::CalculateCovarianceMatrix(
        const std::vector<Node<SPACE_DIM>*>& rNodes) const noexcept
{
    // Create copies of the nodes: it is imperative to keep the oder correct, so we assign sequential indices
    std::vector<Node<SPACE_DIM>*> node_copies(rNodes.size());
    for (unsigned node_idx = 0; node_idx < rNodes.size(); ++node_idx)
    {
        node_copies[node_idx] = new Node<SPACE_DIM>(node_idx, rNodes[node_idx]->rGetLocation());
    }

    std::vector<std::pair<Node<SPACE_DIM>*, Node<SPACE_DIM>*>> node_pairs;
    mpBoxCollection->CalculateNodePairs(node_copies, node_pairs);

    // Create a vector of eigen triplets for half the interactions (not including the diagonal) representing the
    // pairwise covariance using the Gaussian covariance function
    const double tol_cov = -std::log(1e-15);  // \todo: remove this magic number
    const double length_squared = mLengthScale * mLengthScale;
    std::vector<Eigen::Triplet<double>> triplets;
    for (const auto& pair : node_pairs)
    {
        const double dist_squared = GetSquaredDistAtoB(pair.first->rGetLocation(), pair.second->rGetLocation());
        const double exponent = dist_squared / length_squared;

        if (exponent < tol_cov)
        {
            const unsigned a_idx = pair.first->GetIndex();
            const unsigned b_idx = pair.second->GetIndex();

            // Add both directions, as the box collection does not double-count the pairs
            triplets.emplace_back(Eigen::Triplet<double>(a_idx, b_idx, std::exp(-exponent)));
            triplets.emplace_back(Eigen::Triplet<double>(b_idx, a_idx, std::exp(-exponent)));
        }
    }

    // Efficiently take case of the diagonal terms
    for (unsigned node_idx = 0; node_idx < node_copies.size(); ++node_idx)
    {
        triplets.emplace_back(Eigen::Triplet<double>(node_idx, node_idx, 1.0));
    }

    // Create the sparse matrix C given the triplets, and add the identity to capture the self interactions
    Eigen::SparseMatrix<double> cov_matrix(rNodes.size(), rNodes.size());
    cov_matrix.setFromTriplets(triplets.begin(), triplets.end());

    // Clean up the temporary nodes we created
    for (const auto& node : node_copies)
    {
        delete(node);
    }

    return cov_matrix;
}

template <unsigned SPACE_DIM>
double OffLatticeRandomFieldGenerator<SPACE_DIM>::GetSquaredDistAtoB(
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
std::vector<double> OffLatticeRandomFieldGenerator<SPACE_DIM>::SampleRandomField() const noexcept
{
    // Generate a normally-distributed random number for each node in the mesh
    std::vector<double> samples_from_n01(mEigenvecs.rows());
    for (auto& sample : samples_from_n01)
    {
        sample = RandomNumberGenerator::Instance()->StandardNormalRandomDeviate();
    }

    // Special case of zero correlation length
    if (mLengthScale == 0.0)
    {
        return samples_from_n01;
    }

    // Generate the instance of the random field
    Eigen::VectorXd grf = Eigen::VectorXd::Zero(mEigenvecs.rows());
    for (unsigned j = 0; j < mEigenvecs.cols(); ++j)
    {
        const double factor = samples_from_n01[j] * mSqrtEigenvals.coeffRef(j);
        grf += factor * mEigenvecs.col(j);
    }

    // Translate to a std::vector so that eigen objects aren't leaking out to other places in Chaste
    return std::vector<double>(grf.data(), grf.data() + grf.size());
}

// Explicit instantiation
template class OffLatticeRandomFieldGenerator<1>;
template class OffLatticeRandomFieldGenerator<2>;
template class OffLatticeRandomFieldGenerator<3>;
