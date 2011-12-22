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
    Hdf5DataReader *mpDataReader;
    /** Name of the variable representing the membrane potential. */
    const std::string mVoltageName;
    /** Time values */
    std::vector<double> mTimes;
    /** Which node voltages have been cached for, if any */
    unsigned mCachedNodeGlobalIndex;
    /** The cached voltages vector */
    std::vector<double> mCachedVoltages;

    /**
     * Read the voltages vector for the given node and cache it, returning a reference
     * to the cached vector.  If subsequently called with the same index, will return
     * the cached vector without re-reading from file.
     *
     * Note: will only cache the last node index used.
     *
     * @param globalNodeIndex  the index of the node to cache voltages for
     */
    std::vector<double>& rCacheVoltages(unsigned globalNodeIndex);

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
     * Calculate the maximum upstroke velocity at a single cell.
     * We calculate for the last upstroke found in the simulation data.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex);

     /**
     * Calculate the maximum upstroke velocity at a single cell.
     * We return all the max upstroke velocities for all APs.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold   The voltage threshold (we use this for marking the end of an AP)
     */
    std::vector<double> CalculateAllMaximumUpstrokeVelocities(unsigned globalNodeIndex, double threshold);

     /**
     * Calculate the times of upstroke at a single cell.
     * We return all the times of upstroke velocities for all APs.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold   The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<double> CalculateUpstrokeTimes(unsigned globalNodeIndex, double threshold);

    /**
     * Calculate the conduction velocity between two cells, i.e. the time
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
     * Calculate all the conduction velocities between two cells, i.e. the time
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
     * Calculate the action potential duration at a single cell.
     * We calculate for the last AP found in the simulation data.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculateActionPotentialDuration(const double percentage,
                                            unsigned globalNodeIndex);
    /**
     * Calculate the maximum transmembrane potential (maximum systolic
     * potential) at a single cell.
     * We calculate for the last AP found in the simulation data.
     *
     * @param globalNodeIndex  The cell at which to calculate.
     */
    double CalculatePeakMembranePotential(unsigned globalNodeIndex);

     /**
     * Calculate all the action potentials duration at a single cell.
     *
     * @param percentage  The percentage of the amplitude to calculate for.
     * @param globalNodeIndex  The cell at which to calculate.
     * @param threshold  The voltage threshold for APD calculation (we count this as the start of an AP)
     */
    std::vector<double> CalculateAllActionPotentialDurations(const double percentage,
                                                             unsigned globalNodeIndex,
                                                             double threshold);

     /**
     * Calculate all the action potentials duration at a cells [lowerNodeIndex, upperNodeIndex-1].
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
      * Calculates all the depolarisations that occur above threshold at a single cell.
      *
      * @param globalNodeIndex the cell at which to calculate
      * @param threshold the threshold above which the depolarisations are counted
      *
      */
    std::vector<unsigned> CalculateAllAboveThresholdDepolarisations(unsigned globalNodeIndex,
                                                                    double threshold);

     /**
      * Calculates the depolarisations that occur above threshold at a single cell during the last recorded Ap
      *
      * @param globalNodeIndex the cell at which to calculate
      * @param threshold the threshold above which the depolarisations are counted
      *
      */
    unsigned CalculateAboveThresholdDepolarisationsForLastAp(unsigned globalNodeIndex,
                                                             double threshold);

};

#endif //_PROPAGATIONPROPERTIESCALCULATOR_HPP_
