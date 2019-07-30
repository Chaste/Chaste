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


#ifndef _PROPAGATIONPROPERTIESCALCULATOR_HPP_
#define _PROPAGATIONPROPERTIESCALCULATOR_HPP_

#include "Hdf5DataReader.hpp"
#include <string>


/**
 * Calculate physiological properties at given global mesh indices
 *  - maximum upstroke velocity at a single cell
 *  - times of upstroke at a single cell
 *  - (all) conduction velocities between two cells
 *  - (all) action potential duration at a single cell (node)
 *  - (all) action potential durations for a range of nodes
 *  - maximum transmembrane potential (maximum systolic potential) at a single cell.
 */
class PropagationPropertiesCalculator
{
private:
    /** Reader to get the data from which we use to calculate properties. */
    Hdf5DataReader* mpDataReader;
    /** Name of the variable representing the membrane potential. */
    const std::string mVoltageName;
    /** Time values */
    std::vector<double> mTimes;
    /** Which node voltages have been cached for, if any */
    unsigned mCachedNodeGlobalIndex;
    /** The cached voltages vector */
    std::vector<double> mCachedVoltages;

protected:
    /**
     * @return the voltages vector for the given node and cache it, returning a reference
     * to the cached vector.  If subsequently called with the same index, will return
     * the cached vector without re-reading from file.
     *
     * Note: will only cache the last node index used.
     *
     * @param globalNodeIndex  the index of the node to cache voltages for
     */
    std::vector<double>& rGetCachedVoltages(unsigned globalNodeIndex);

public:
    /**
     * Constructor.
     *
     * @param pDataReader  Pointer to the data reader containing the simulation.
     * @param voltageName  Optionally the name of the variable representing the
     *     membrane potential.  Defaults to "V".
     */
    PropagationPropertiesCalculator(Hdf5DataReader* pDataReader,
                                    const std::string voltageName = "V");

    /** Destructor */
    virtual ~PropagationPropertiesCalculator();

    /**
     * @return the maximum upstroke velocity at a single cell.
     * We calculate for the last upstroke found in the simulation data.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex);

     /**
     * @return the maximum upstroke velocity at a single cell.
     * We return all the max upstroke velocities for all APs.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold   The voltage threshold (we use this for marking the end of an AP)
     */
    std::vector<double> CalculateAllMaximumUpstrokeVelocities(unsigned globalNodeIndex, double threshold);

     /**
     * @return the times of upstroke at a single cell.
     * We return all the times of upstroke velocities for all APs.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold   The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<double> CalculateUpstrokeTimes(unsigned globalNodeIndex, double threshold);

    /**
     * @return the conduction velocity between two cells, i.e. the time
     * taken for an AP to propagate from one to the other. It returns
     * the value of conduction velocity of the LAST action potential
     * that reached both nodes. Throws exceptions if an AP never reached
     * one of the nodes.
     *
     * @param globalNearNodeIndex  The cell to measure from.
     * @param globalFarNodeIndex  The cell to measure to.
     * @param euclideanDistance  The distance the AP travels between the cells, along the tissue.
     */
    double CalculateConductionVelocity(unsigned globalNearNodeIndex,
                                       unsigned globalFarNodeIndex,
                                       const double euclideanDistance);

     /**
     * @return all the conduction velocities between two cells, i.e. the time
     * taken for all APs to propagate from one to the other. It returns a vector
     * containing all the conduction velocities for each of the APs that
     * reached the two nodes (only the APs that reached both nodes).
     * Throws exceptions if an AP never reached one of the nodes.
     *
     * @param globalNearNodeIndex  The cell to measure from.
     * @param globalFarNodeIndex  The cell to measure to.
     * @param euclideanDistance  The distance the AP travels between the cells, along the tissue.
     */
     std::vector<double> CalculateAllConductionVelocities(unsigned globalNearNodeIndex,
                                                          unsigned globalFarNodeIndex,
                                                          const double euclideanDistance);
    /**
     * @return the action potential duration at a single cell.
     * We calculate for the last AP found in the simulation data.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateActionPotentialDuration(const double percentage,
                                            unsigned globalNodeIndex);
    /**
     * @return the maximum transmembrane potential (maximum systolic
     * potential) at a single cell.
     * We calculate for the last AP found in the simulation data.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculatePeakMembranePotential(unsigned globalNodeIndex);

     /**
     * @return all the action potentials duration at a single cell.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold  The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<double> CalculateAllActionPotentialDurations(const double percentage,
                                                             unsigned globalNodeIndex,
                                                             double threshold);

     /**
     * @return all the action potentials duration at cells [lowerNodeIndex, upperNodeIndex-1].
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param lowerNodeIndex  First cell at which to calculate.
     * @param upperNodeIndex  One past the last cell at which to calculate.
     * @param threshold  The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<std::vector<double> > CalculateAllActionPotentialDurationsForNodeRange(const double percentage,
                                                                                       unsigned lowerNodeIndex,
                                                                                       unsigned upperNodeIndex,
                                                                                       double threshold);

     /**
      * @return all the depolarisations that occur above threshold at a single cell.
      *
      * @param globalNodeIndex the cell at which to calculate
      * @param threshold the threshold above which the depolarisations are counted
      *
      */
    std::vector<unsigned> CalculateAllAboveThresholdDepolarisations(unsigned globalNodeIndex,
                                                                    double threshold);

     /**
      * @return the depolarisations that occur above threshold at a single cell during the last recorded Ap
      *
      * @param globalNodeIndex the cell at which to calculate
      * @param threshold the threshold above which the depolarisations are counted
      *
      */
    unsigned CalculateAboveThresholdDepolarisationsForLastAp(unsigned globalNodeIndex,
                                                             double threshold);

    /**
     * Provide a new pointer to an HDF5 data reader
     *
     * @param pDataReader  An HDF5 data reader to use (needed if the existing one is deleted and a new one opened)
     */
    void SetHdf5DataReader(Hdf5DataReader* pDataReader);
};

#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_
