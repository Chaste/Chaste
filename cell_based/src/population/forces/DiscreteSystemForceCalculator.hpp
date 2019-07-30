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

#ifndef DISCRETESYSTEMFORCECALCULATOR_HPP_
#define DISCRETESYSTEMFORCECALCULATOR_HPP_

#include "AbstractTwoBodyInteractionForce.hpp"
#include "MeshBasedCellPopulation.hpp"

/**
 * A class for calculating the force and stress on each node in a mesh-based cell population.
 */
class DiscreteSystemForceCalculator
{
    friend class TestDiscreteSystemForceCalculator;

private:

    /**
     * Reference to cell population.
     */
    MeshBasedCellPopulation<2>& mrCellPopulation;

    /** The mechanics used to determine the new location of the cells. */
    std::vector<boost::shared_ptr<AbstractTwoBodyInteractionForce<2> > > mForceCollection;

    /**
     * Small parameter, used in GetSamplingAngles().
     */
    double mEpsilon;

    /** The file that the results of CalculateExtremalNormalForces. */
    out_stream mpVizStressResultsFile;

    /**
     * Given a node index and angle of intersecting line in the range (-pi,pi],
     * returns the tangential and normal forces.
     *
     * @param index the node index
     * @param theta the angle of intersection
     *
     * @return the vector of tangential and normal forces.
     */
    std::vector<double> CalculateFtAndFn(unsigned index, double theta);

    /**
     * Given a node index, returns a vector of sampling angles in the range (-pi,pi]
     * that can be used by GetExtremalAngles() to find the locations of local extrema
     * of the normal force.
     *
     * @param index
     *
     * @return the vector of sampling angles.
     */
    std::vector<double> GetSamplingAngles(unsigned index);

    /**
     * Given a node index and two sampling angles, finds the location of
     * the root of the tangential force in the interval between the two
     * angles. There is no guarantee that this will lie in (-pi,pi].
     *
     * @param index the node index
     * @param angle1 the first sampling angle
     * @param angle2 the second sampling angle
     *
     * @return the local extremum.
     */
    double GetLocalExtremum(unsigned index, double angle1, double angle2);

    /**
     * Given a vector of sampling angles in the range (-pi,pi], returns a vector
     * of extremal angles, i.e. angles at which local extrema of the normal force
     * occur, again in the range (-pi,pi].
     *
     * @param index the node index
     * @param samplingAngles the vector of sampling angles
     *
     * @return the vector of extremal angles.
     */
    std::vector<double> GetExtremalAngles(unsigned index, std::vector<double> samplingAngles);

public:

    /**
     * Constructor.
     *
     * @param rCellPopulation reference to the cell population
     * @param forceCollection vector of force laws present
     */
    DiscreteSystemForceCalculator(MeshBasedCellPopulation<2>& rCellPopulation, std::vector<boost::shared_ptr<AbstractTwoBodyInteractionForce<2> > > forceCollection);

    /**
     * @return the extremal normal forces on each node in the cell population.
     */
    std::vector< std::vector<double> > CalculateExtremalNormalForces();

    /**
     * Write results to file.
     *
     * @param simulationOutputDirectory the output directory, relative to where Chaste output is stored
     */
    void WriteResultsToFile(std::string simulationOutputDirectory);

    /**
     * Set #mEpsilon.
     *
     * @param epsilon the new value of mEpsilon
     */
    void SetEpsilon(double epsilon);
};

#endif /*DISCRETESYSTEMFORCECALCULATOR_HPP_*/
