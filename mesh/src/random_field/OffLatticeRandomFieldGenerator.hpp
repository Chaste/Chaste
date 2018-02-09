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

#ifndef OFF_LATTICE_RANDOM_FIELD_GENERATOR_HPP_
#define OFF_LATTICE_RANDOM_FIELD_GENERATOR_HPP_

#include <array>
#include <memory>
#include <vector>

// Eigen includes
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>

#include "ChasteSerialization.hpp"
#include "Node.hpp"
#include "ObsoleteBoxCollection.hpp"
#include "UblasVectorInclude.hpp"

/**
 * An OffLatticeRandomFieldGenerator for adding spatially-correlated noise to a mesh with off-lattice nodes.
 */
template<unsigned SPACE_DIM>
class OffLatticeRandomFieldGenerator
{
private:

    friend class boost::serialization::access;
    friend class TestOffLatticeRandomFieldGenerator;

    /** Coordinates of the lower corner of the grid */
    std::array<double, SPACE_DIM> mLowerCorner;

    /** Coordinates of the upper corner of the grid */
    std::array<double, SPACE_DIM> mUpperCorner;

    /** Whether the grid is periodic in each dimension */
    std::array<bool, SPACE_DIM> mPeriodicity;

    /** Number of eigenvalues to calculate in the matrix decomposition */
    unsigned mNumEigenvals;

    /** Length scale over which the noise is to be correlated */
    double mLengthScale;

    /** Array storing the square roots of the first mNumEigenvals largest eigenvalues of the covariance matrix */
    Eigen::ArrayXd mSqrtEigenvals;

    /** Matrix storing the first mNumEigenvals eigenvectors of the covariance matrix */
    Eigen::MatrixXd mEigenvecs;

    /** An owning pointer to a box collection for efficiently keeping track of nearby nodes */
    std::unique_ptr<ObsoleteBoxCollection<SPACE_DIM>> mpBoxCollection;

    /** The number of nodes passed in at the most recent Update() */
    unsigned mNumNodesAtLastUpdate;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // Todo: how can we serialize Eigen arrays/matrices?
        archive & mLowerCorner;
        archive & mUpperCorner;
        archive & mPeriodicity;
        archive & mNumEigenvals;
        archive & mLengthScale;
    }

    /**
     * Helper method for Update().
     *
     * This matrix C is of size rNodes.size() x rNodes.size(), where C_ij is the some function of the distance
     * between locations i and j.  Here, that function is exp{-dist_squared / mLengthScale^2}.
     *
     * @param rNodes the nodes in the mesh
     * @return and Eigen::SparseMatrix<double> holding the point-point covariance data
     */
    Eigen::SparseMatrix<double> CalculateCovarianceMatrix(const std::vector<Node<SPACE_DIM>*>& rNodes) const noexcept;

    /**
     * Get the squared distance between two points, which is needed to calculate the covariance matrix.
     * This function takes into account possible periodicity in the mesh.
     *
     * @param rLocation1 the first location
     * @param rLocation2 the second location
     * @return the squared distance between rLocation1 and rLocation2
     */
    double GetSquaredDistAtoB(const c_vector<double, SPACE_DIM>& rLocation1,
                              const c_vector<double, SPACE_DIM>& rLocation2) const noexcept;

public:

    /**
     * Constructor that takes all parameters as arguments, and sets up an appropriate box collection.
     *
     * @param lowerCorner the lower corner of the rectangular grid
     * @param upperCorner the upper corner of the rectangular grid
     * @param periodicity whether the grid is periodic in each dimension
     * @param numEigenvals the number of eigenvalues to calculate
     * @param lengthScale the length scale of the correlation when calculating the covariance matrix
     * @param boxWidth the box size for the box collection. Defaults to DOUBLE_UNSET in which case it is reset appropriately.
     */
    OffLatticeRandomFieldGenerator(std::array<double, SPACE_DIM> lowerCorner,
                                   std::array<double, SPACE_DIM> upperCorner,
                                   std::array<bool, SPACE_DIM> periodicity,
                                   unsigned numEigenvals,
                                   double lengthScale,
                                   double boxWidth=DOUBLE_UNSET);


    /**
     * Default constructor
     */
    virtual ~OffLatticeRandomFieldGenerator() = default;

    /**
     * Update the members mSqrtEigenvals and mEigenvecs. This is the main workhorse of this class.
     *
     * Uses the Spectra routines to calculate the first mEigenvals eigenvalues and eigenvectors of the covariance
     * matrix. This method will throw if the Spectra eigen decomposition is unsuccessful.
     *
     * This method returns early if mLengthScale = 0 (noise is not spatially correlated at all).
     *
     * @param rNodes the nodes in the mesh
     */
    void Update(const std::vector<Node<SPACE_DIM>*>& rNodes);

    /**
     * Tune the number of eigenvalues to calculate, by calculating the majority and determining the number at which
     * the partial sum of eigenvalues exceeds proportionOfTrace * Trace(C).
     *
     * This method re-sets numEigenvals to a value different than that set in the constructor.
     *
     * This method returns early if mLengthScale = 0 (noise is not spatially correlated at all).
     *
     * @param rNodes
     * @param proportionOfTrace
     * @return the number of eigenvalues that will be calculated per update, after this tuning
     */
    unsigned TuneNumEigenvals(const std::vector<Node<SPACE_DIM>*>& rNodes, const double proportionOfTrace);

    /**
     * Sample an instance of the random field.  First, draw mNumTotalGridPts random numbers from N(0,1), and then
     * create an appropriate linear combination of the eigenvectors.
     *
     * This method returns early if mLengthScale = 0 (noise is not spatially correlated at all).
     *
     * @return A vector representing an instance of the random field.
     */
    std::vector<double> SampleRandomField() const noexcept;
};


#endif /*OFF_LATTICE_RANDOM_FIELD_GENERATOR_HPP_*/
