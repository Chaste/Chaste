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

#ifndef POSTPROCESSINGWRITER_HPP_
#define POSTPROCESSINGWRITER_HPP_


#include "Hdf5DataReader.hpp"
#include "PropagationPropertiesCalculator.hpp"
#include "AbstractTetrahedralMesh.hpp"
#include <string>

/**
 * Write out physiological parameters at the end of a simulation
 * - APD map
 * - Upstroke time map
 * - Upstroke Velocity map
 * - Conduction Velocity map
 *
 * N.B. You should only ever have one PostProcessingWriter around at once, as
 * multiple Hdf5Readers (a member variable of this class) seem to cause problems.
 *
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PostProcessingWriter
{
    friend class TestPostProcessingWriter;

private:
    FileFinder mDirectory; /**< The directory the HDF5 file is in */
    std::string mHdf5File; /**< The name of the HDF5 file to post-process */
    std::string mVoltageName; /**< The name of the variable representing the membrane potential */

    Hdf5DataReader* mpDataReader; /**< An HDF5 reader from which to build the PropagationPropertiesCalculator */
    PropagationPropertiesCalculator* mpCalculator; /**< PropagationPropertiesCalculator based on HDF5 data reader*/
    unsigned mLo; /**< Cache of mLo from the mesh DitributedVectorFactory */
    unsigned mHi; /**< Cache of mHi from the mesh DitributedVectorFactory */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& mrMesh;/**< A mesh used to calculate the distance map to pass to the conduction velocity calculator*/
    hsize_t mHdf5DataWriterChunkSize; /**< Chunk size parameter for Hdf5DataWriter */

public:
    /**
     * Constructor
     *
     * @param rMesh A reference to the mesh used to calculate the distance map to pass to the conduction velocity calculator.
     * @param rDirectory The directory the data is in. The output is written to \<directory\>/output
     * @param rHdf5FileName The file the data is in.
     * @param rVoltageName  (Optional) The name of the variable representing the
     *     membrane potential. It is used in the creation of the PropagationPropertiesCalculator object. Defaults to "V".
     * @param hdf5DataWriterChunkSize (Optional) Chunk size and alignment parameter to pass to Hdf5DataWriter
     */
    PostProcessingWriter(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh,
                         const FileFinder& rDirectory,
                         const std::string& rHdf5FileName,
                         const std::string& rVoltageName = "V",
                         hsize_t hdf5DataWriterChunkSize=0);

    /**
     *  Write out data files. The data that is written depends on which maps have been requested using
     *  either the XML file or HeartConfig
     */
    void WritePostProcessingFiles();

    /**
     * Destructor
     */
    ~PostProcessingWriter();

    /**
     * Method for opening a file and writing one row per node:
     * line 1: <number of upstrokes for node 0> <number of above-threshold depolarisations for node 0>
     * line 2: <number of upstrokes for node 1> <number of above-threshold depolarisations for node 1>
     * etc.
     *
     * For the nodes where the threshold isn't crossed, the 'number of upstrokes' will be 0
     * (so will the number of above-threshold depolarisations for that node)
     *
     * \todo This method ought to be private and called by the  WritePostProcessingFiles method if the user requests for it.
     *       This will be possible after modifying the schema and specifying Get and Set methods in HeartConfig
     *       to check whether the user wants this file or not
     *
     * @param  threshold - used to signify the upstroke (mV) AND to specify above which voltage value the depolarisations are counted
     */
    void WriteAboveThresholdDepolarisationFile(double threshold);

private:

    /**
     * Method that extrapolates the output variables over time at specified nodes and output all to file.
     * The use of this method is intended as follows: the user supplies a list of node indices (rNodeIndices).
     * This method outputs the time series at each node in rNodeIndices (one time series per column).
     * The node numbering is referred to the original mesh.
     *
     * Regardless of the permutation that was used in the simulation, this method will have the same output.
     *
     * There will be one file per variable contained in the hdf5 file.
     * For each file, the name will be NodalTraces_VARIABLENAME.dat.
     * So for example, in a bidomain simulation with default variable names, output files will be:
     *
     * NodalTraces_V.dat
     * NodalTraces_Phi_e.dat
     *
     * Each output file will look like (taking variable V as an example):
     *
     * Line 1: <V at time 0, Node 0> <V at time 0, Node 1> <V at time 0, Node 2> etc...
     * Line 2: <V at time 1, Node 0> <V at time 1, Node 1> <V at time 1, Node 2> etc...
     * etc...
     *
     * where Node 0, Node 1, Node 2 are the nodes included in the rIndices vector and "V" is the name of the variable.
     *
     * @param rNodeIndices the node indices (in the unpermuted mesh) that we want the output for.
     */
    void WriteVariablesOverTimeAtNodes(std::vector<unsigned>& rNodeIndices);

    /**
     * Method for opening an APD map file and writing one row per node
     * line 1: <first APD for node 0> <second APD for node 0> ...
     * line 2: <first APD for node 1> <second APD for node 1> ...
     * etc.
     *
     * Nodes where there is no APD are respresented by a single
     * 0
     *
     * @param  repolarisationPercentage eg. 90.0 for APD90
     * @param  threshold - Vm used to signify the upstroke (mV)
     */
    void WriteApdMapFile(double repolarisationPercentage, double threshold);


    /**
     * Write out times of each upstroke for each node:
     *
     * line 1: <first upstroke time for node 0> <second upstroke time for node 0> ...
     * line 2: <first upstroke time for node 1> <second upstroke time for node 1> ...
     * etc.
     *
     * If there is no upstroke then there will a blank line
     *
     * @param threshold  - Vm used to signify the upstroke (mV)
     */
    void WriteUpstrokeTimeMap(double threshold);

    /**
     * Write out velocities of each max upstroke for each node:
     *
     * line 1: <first upstroke velocity for node 0> <second upstroke velocity for node 0> ...
     * line 2: <first upstroke velocity for node 1> <second upstroke velocity for node 1> ...
     * etc.
     *
     * If there is no upstroke then there will a blank line
     *
     * @param threshold  - Vm used to signify the upstroke (mV)
     */
    void WriteMaxUpstrokeVelocityMap(double threshold);

    /**
     * Write out conduction velocity map from the given node the rest of the mesh:
     *
     * line 1: <conduction velocity for node 0 and AP 0> <conduction velocity for node 0 and AP 1> ...
     * line 2: <conduction velocity for node 1 and AP 0> <conduction velocity for node 1 and AP 1> ...
     * etc.
     *
     * Note: the line corresponding to node number originNode will contain ...
     *
     * @param originNode  - Node to compute the conduction velocity from
     * @param distancesFromOriginNode - Distance map from originNode to all the nodes in the simulation. Typically calculated with DistanceMapCalculator
     */
    void WriteConductionVelocityMap(unsigned originNode, std::vector<double> distancesFromOriginNode);

    /**
     * Method for opening a file and writing one row per node
     * line 1: <first scalar data for node 0> <second scalar data for node 0> ...
     * line 2: <first scalar data for node 1> <second scalar data for node 1> ...
     * etc.
     * @param  rDataPayload vector data for each node.  Each node's data are represented by a vector of scalars (variable length)
     * @param  rFolder subfolder for postprocessing in which to put the data.
     * @param  rFileName where to put the data.
     */
    void WriteGenericFileToMeshalyzer(std::vector<std::vector<double> >& rDataPayload, const std::string& rFolder, const std::string& rFileName);

    /**
     * Put the post-processed data into the main HDF5 results file.
     *
     * @param rDataPayload  The postprocessed quantities
     * @param rDatasetName  The name of the quantities
     * @param rDatasetUnit  The unit of the quantities
     * @param rUnlimitedVariableName  The name of the unlimited variable (defaults to "PaceNumber")
     * @param rUnlimitedVariableUnit  The unlimited variable units (defaults to "dimensionless")
     */
    void WriteOutputDataToHdf5(const std::vector<std::vector<double> >& rDataPayload,
                               const std::string& rDatasetName,
                               const std::string& rDatasetUnit,
                               const std::string& rUnlimitedVariableName = "PaceNumber",
                               const std::string& rUnlimitedVariableUnit = "dimensionless");

    /**
     * Convert a string with numbers in it into alphanumeric plus underscores.
     *
     * e.g.
     * 20 -> "_20"
     * -20 -> "_minus_20"
     * 30.2 -> "_30pt20"
     * -11.238 -> "_minus_11pt23" (always does decimals to (floor) 2d.p.)
     *
     * @param threshold  A numerical threshold which may contain minuses or a decimal point.
     * @return  A string version of the number without minuses or decimal points.
     */
    std::string ConvertToHdf5FriendlyString(double threshold);
};

#endif /*POSTPROCESSINGWRITER_HPP_*/
