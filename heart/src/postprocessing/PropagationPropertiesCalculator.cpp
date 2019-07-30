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

#include "UblasIncludes.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "CellProperties.hpp"
#include "Exception.hpp"
#include <sstream>
#include "HeartEventHandler.hpp"

PropagationPropertiesCalculator::PropagationPropertiesCalculator(Hdf5DataReader* pDataReader,
                                                                 const std::string voltageName)
    : mpDataReader(pDataReader),
      mVoltageName(voltageName),
      mTimes(mpDataReader->GetUnlimitedDimensionValues()),
      mCachedNodeGlobalIndex(UNSIGNED_UNSET)
{}

PropagationPropertiesCalculator::~PropagationPropertiesCalculator()
{
    // We don't own the data reader, so we don't destroy it.
}

double PropagationPropertiesCalculator::CalculateMaximumUpstrokeVelocity(unsigned globalNodeIndex)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes);
    return cell_props.GetLastMaxUpstrokeVelocity();
}

std::vector<double> PropagationPropertiesCalculator::CalculateAllMaximumUpstrokeVelocities(unsigned globalNodeIndex, double threshold)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes, threshold);
    return cell_props.GetMaxUpstrokeVelocities();
}

std::vector<double> PropagationPropertiesCalculator::CalculateUpstrokeTimes(unsigned globalNodeIndex, double threshold)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes, threshold);
    return cell_props.GetTimesAtMaxUpstrokeVelocity();
}

double PropagationPropertiesCalculator::CalculateActionPotentialDuration(const double percentage,
                                                                         unsigned globalNodeIndex)
{
    if (percentage < 1.0 || percentage >= 100.0)
    {
        EXCEPTION("First argument of CalculateActionPotentialDuration() is expected to be a percentage");
    }
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes);
    return cell_props.GetLastActionPotentialDuration(percentage);
}

std::vector<double> PropagationPropertiesCalculator::CalculateAllActionPotentialDurations(const double percentage,
                                                                                          unsigned globalNodeIndex,
                                                                                          double threshold)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes, threshold);
    return cell_props.GetAllActionPotentialDurations(percentage);
}

std::vector<std::vector<double> > PropagationPropertiesCalculator::CalculateAllActionPotentialDurationsForNodeRange(
        const double percentage,
        unsigned lowerNodeIndex,
        unsigned upperNodeIndex,
        double threshold)
{
    std::vector<std::vector<double> > output_data;
    output_data.reserve(upperNodeIndex-lowerNodeIndex+1);
    unsigned num_nodes_per_data_block = 100; // number of nodes
    unsigned num_complete_blocks = (upperNodeIndex-lowerNodeIndex) / num_nodes_per_data_block;
    unsigned size_last_block = (upperNodeIndex-lowerNodeIndex) % num_nodes_per_data_block;

    for (unsigned block_num=0;
         block_num<num_complete_blocks+1;
         block_num++)
    {
        unsigned num_nodes_to_read;
        if (block_num != num_complete_blocks)
        {
            num_nodes_to_read = num_nodes_per_data_block;
        }
        else
        {
            num_nodes_to_read = size_last_block;
        }

        if (num_nodes_to_read > 0)
        {
            // Read a big block of data
            unsigned low_node = lowerNodeIndex + block_num*num_nodes_per_data_block;
            unsigned high_node = low_node + num_nodes_to_read;
            std::vector<std::vector<double> > voltages = mpDataReader->GetVariableOverTimeOverMultipleNodes(mVoltageName, low_node, high_node);

            for (unsigned node_within_block=0;
                 node_within_block < num_nodes_to_read;
                 node_within_block++)
            {
                std::vector<double>& r_voltages = voltages[node_within_block];
                CellProperties cell_props(r_voltages, mTimes, threshold);
                std::vector<double> apds;
                try
                {
                    apds = cell_props.GetAllActionPotentialDurations(percentage);
                    assert(apds.size() != 0);
                }
                catch (Exception& e)
                {
                    assert(e.GetShortMessage()=="No full action potential was recorded" ||
                           e.GetShortMessage()=="AP did not occur, never exceeded threshold voltage.");
                    apds.push_back(0);
                    assert(apds.size() == 1);
                }
                output_data.push_back(apds);
            }
        }
    }
    return output_data;
}

double PropagationPropertiesCalculator::CalculatePeakMembranePotential(unsigned globalNodeIndex)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    double max = -DBL_MAX;
    for (unsigned i=0; i<r_voltages.size(); i++)
    {
        if (r_voltages[i]>max)
        {
            max = r_voltages[i];
        }
    }
    return max;
}

double PropagationPropertiesCalculator::CalculateConductionVelocity(unsigned globalNearNodeIndex,
                                                                    unsigned globalFarNodeIndex,
                                                                    const double euclideanDistance)
{
    double t_near = 0;
    double t_far = 0;
    std::vector<double>& r_near_voltages = rGetCachedVoltages(globalNearNodeIndex);
    std::vector<double> far_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalFarNodeIndex);

    CellProperties near_cell_props(r_near_voltages, mTimes);
    CellProperties far_cell_props(far_voltages, mTimes);

    //The size of each vector is the number of APs that reached that node
    unsigned aps_near_node = near_cell_props.GetMaxUpstrokeVelocities().size();
    unsigned aps_far_node = far_cell_props.GetMaxUpstrokeVelocities().size();

    //These should never be empty. If so, an exception should have been thrown in the GetMaxUpstrokeVelocities() method.
    assert(aps_near_node > 0);
    assert(aps_far_node > 0);

    //if the same number of APs reached both nodes, get the last one...
    if (aps_near_node == aps_far_node)
    {
        t_near = near_cell_props.GetTimeAtLastMaxUpstrokeVelocity();
        t_far = far_cell_props.GetTimeAtLastMaxUpstrokeVelocity();
    }
    //...otherwise get the one with the smallest value, which is the last AP to reach both nodes
    //This prevents possible calculation of negative conduction velocities
    //for repeated stimuli
    else if (aps_near_node > aps_far_node)
    {
        t_near = near_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_far_node-1];
        t_far = far_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_far_node-1];
    }
    else
    {
        t_near = near_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_near_node-1];
        t_far = far_cell_props.GetTimesAtMaxUpstrokeVelocity()[aps_near_node-1];
    }

    ///\todo remove magic number? (#1884)
    if ((globalNearNodeIndex == globalFarNodeIndex) || ( fabs(t_far - t_near) < 1e-8))
    {
        // globalNearNodeIndex and globalFarNodeIndex are the same node, preventing a 0/0
        // or
        // AP number i is happening at the same time at nodes globalNearNodeIndex and globalFarNodeIndex
        return 0.0;
    }
    else
    {
        return euclideanDistance / (t_far - t_near);
    }


}

std::vector<double> PropagationPropertiesCalculator::CalculateAllConductionVelocities(unsigned globalNearNodeIndex,
                                                                                      unsigned globalFarNodeIndex,
                                                                                      const double euclideanDistance)
{
    std::vector<double> conduction_velocities;

    std::vector<double> t_near;
    std::vector<double> t_far;
    unsigned number_of_aps = 0;

    std::vector<double>& r_near_voltages = rGetCachedVoltages(globalNearNodeIndex);
    std::vector<double> far_voltages = mpDataReader->GetVariableOverTime(mVoltageName, globalFarNodeIndex);

    CellProperties near_cell_props(r_near_voltages, mTimes);
    CellProperties far_cell_props(far_voltages, mTimes);

    t_near = near_cell_props.GetTimesAtMaxUpstrokeVelocity();
    t_far = far_cell_props.GetTimesAtMaxUpstrokeVelocity();

    //exception should have been thrown within the GetTimesAtMaxUpstrokeVelocity method if the threshold is never reached
    //and these vectors are empty
    assert(t_near.size() !=0);
    assert(t_far.size() !=0);

    //Check the node where the least number of aps is reached.
    //We will calculate only where AP reached both nodes
    if (t_near.size() > t_far.size())
    {
        number_of_aps = t_far.size();
    }
    else
    {
        number_of_aps = t_near.size();
    }
    //now fill the vector

    if (globalNearNodeIndex == globalFarNodeIndex)
    {
        // globalNearNodeIndex and globalFarNodeIndex are the same node, preventing a 0/0
        for (unsigned i = 0 ; i < number_of_aps;i++)
        {
            conduction_velocities.push_back(0.0);
        }
    }
    else
    {
        for (unsigned i = 0 ; i < number_of_aps;i++)
        {
            ///\todo remove magic number? (#1884)
            if (fabs(t_far[i] - t_near[i]) < 1e-8)
            {
                // AP number i is happening at the same time at nodes globalNearNodeIndex and globalFarNodeIndex
                conduction_velocities.push_back(0.0);
            }
            else
            {
                conduction_velocities.push_back(euclideanDistance / (t_far[i] - t_near[i]));
            }
        }
    }

    return conduction_velocities;
}


std::vector<unsigned> PropagationPropertiesCalculator::CalculateAllAboveThresholdDepolarisations(unsigned globalNodeIndex, double threshold)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes, threshold);
    return cell_props.GetNumberOfAboveThresholdDepolarisationsForAllAps();
}


unsigned PropagationPropertiesCalculator::CalculateAboveThresholdDepolarisationsForLastAp(unsigned globalNodeIndex, double threshold)
{
    std::vector<double>& r_voltages = rGetCachedVoltages(globalNodeIndex);
    CellProperties cell_props(r_voltages, mTimes, threshold);
    return cell_props.GetNumberOfAboveThresholdDepolarisationsForLastAp();
}


std::vector<double>& PropagationPropertiesCalculator::rGetCachedVoltages(unsigned globalNodeIndex)
{
    if (globalNodeIndex != mCachedNodeGlobalIndex)
    {
        mCachedVoltages = mpDataReader->GetVariableOverTime(mVoltageName, globalNodeIndex);
        mCachedNodeGlobalIndex = globalNodeIndex;
    }
    return mCachedVoltages;
}

void PropagationPropertiesCalculator::SetHdf5DataReader(Hdf5DataReader* pDataReader)
{
    mpDataReader = pDataReader;
}


