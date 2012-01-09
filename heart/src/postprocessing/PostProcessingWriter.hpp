/*

Copyright (C) University of Oxford, 2005-2012

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
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class PostProcessingWriter
{
    friend class TestPostProcessingWriter;

private:
    std::string mDirectory; /**< The directory the HDF5 file is in */
    std::string mHdf5File; /**< The name of the HDF5 file to post-process */
    bool mMakeAbsolute; /**< Whether to convert #mDirectory to an absolute path */
    std::string mVoltageName; /**< The name of the variable representing the membrane potential */

    Hdf5DataReader* mpDataReader; /**< An HDF5 reader from which to build the PropagationPropertiesCalculator */
    PropagationPropertiesCalculator* mpCalculator; /**< PropagationPropertiesCalculator based on HDF5 data reader*/
    unsigned mLo; /**< Cache of mLo from the mesh DitributedVectorFactory */
    unsigned mHi; /**< Cache of mHi from the mesh DitributedVectorFactory */
    AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>& mrMesh;/**< A mesh used to calculate the distance map to pass to the conduction velocity calculator*/

public:
    /**
     * Constructor
     *
     * @param rMesh A reference to the mesh used to calculate the distance map to pass to the conduction velocity calculator.
     * @param directory The directory the data is in. The output is written to \<directory\>/output
     * @param hdf5File The file the data is in.
     * @param makeAbsolute Whether to convert the path to absolute using the OutputFileHandler (via the HdfDataReader)
     * @param voltageName  (Optional) The name of the variable representing the
     *     membrane potential. It is used in the creation of the PropagationPropertiesCalculator object. Defaults to "V".
     */
    PostProcessingWriter(AbstractTetrahedralMesh<ELEMENT_DIM, SPACE_DIM>& rMesh, std::string directory, std::string hdf5File, bool makeAbsolute, std::string voltageName = "V");

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
     *       to check whetehr the user wants this file or not
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
     * @param distancesFromOriginNode - Distance map from originNode to all the nodes in the simulation. Tipically calculated with DistanceMapCalculator
     */
    void WriteConductionVelocityMap(unsigned originNode, std::vector<double> distancesFromOriginNode);
    /**
     * Method for opening a file and writing one row per node
     * line 1: <first scalar data for node 0> <second scalar data for node 0> ...
     * line 2: <first scalar data for node 1> <second scalar data for node 1> ...
     * etc.
     * @param  rDataPayload vector data for each node.  Each node's data are represented by a vector of scalars (variable length)
     * @param  fileName where to put the data.
     */
    void WriteGenericFile(std::vector<std::vector<double> >& rDataPayload, std::string fileName);

};

#endif /*POSTPROCESSINGWRITER_HPP_*/
