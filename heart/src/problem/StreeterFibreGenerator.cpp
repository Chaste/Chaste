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

#include "UblasCustomFunctions.hpp"

#include "StreeterFibreGenerator.hpp"

#include <cmath>
#include <fstream>
#include <sstream>
#include "OutputFileHandler.hpp"
#include "Exception.hpp"
//#include "HeartRegionCodes.hpp"

// Add the citation for original Streeter paper.
#include "Citations.hpp"
static PetscBool StreeterCite = PETSC_FALSE;
const char StreeterCitation[] = "@article{streeter1969fiber,\n"
"  title={Fiber orientation in the canine left ventricle during diastole and systole},\n"
"  author={Streeter, Daniel D and Spotnitz, Henry M and Patel, Dali P and Ross, John and Sonnenblick, Edmund H},\n"
"  journal={Circulation research},\n"
"  volume={24},\n"
"  number={3},\n"
"  pages={339--347},\n"
"  year={1969},\n"
"  publisher={Am Heart Assoc}\n"
"}\n";

template<unsigned SPACE_DIM>
double StreeterFibreGenerator<SPACE_DIM>::GetAveragedThicknessLocalNode(
        const unsigned nodeIndex, const std::vector<double>& wallThickness) const
{
    if (nodeIndex < this->mpMesh->GetDistributedVectorFactory()->GetLow() ||
        nodeIndex >= this->mpMesh->GetDistributedVectorFactory()->GetHigh() )
    {
        return 0.0;  //Don't calculate this for nodes which aren't local
    }

    // Initialise the average with the value corresponding to the current node
    double average = wallThickness[nodeIndex];
    unsigned nodes_visited = 1;

    // Use a set to store visited nodes
    std::set<unsigned> visited_nodes;
    visited_nodes.insert(nodeIndex);

    Node<SPACE_DIM>* p_current_node = this->mpMesh->GetNode(nodeIndex);

    // Loop over the elements containing the given node
    for (typename Node<SPACE_DIM>::ContainingElementIterator element_iterator = p_current_node->ContainingElementsBegin();
        element_iterator != p_current_node->ContainingElementsEnd();
        ++element_iterator)
    {
        // Get a pointer to the container element
        Element<SPACE_DIM,SPACE_DIM>* p_containing_element = this->mpMesh->GetElement(*element_iterator);

       // Loop over the nodes of the element
       for (unsigned node_local_index=0;
           node_local_index<p_containing_element->GetNumNodes();
           node_local_index++)
       {
            Node<SPACE_DIM>* p_neighbour_node = p_containing_element->GetNode(node_local_index);
            unsigned neighbour_node_index = p_neighbour_node->GetIndex();

            // Check if the neighbour node has already been visited
            if (visited_nodes.find(neighbour_node_index) == visited_nodes.end())
            {
                average += wallThickness[neighbour_node_index];
                visited_nodes.insert(neighbour_node_index);
                nodes_visited++;
            }
       }
    }

    return average/nodes_visited;
}

template<unsigned SPACE_DIM>
double StreeterFibreGenerator<SPACE_DIM>::GetFibreMaxAngle(
        const c_vector<HeartRegionType, SPACE_DIM+1>& nodesRegionsForElement) const
{
    unsigned lv=0, rv=0;

    for (unsigned index=0; index<SPACE_DIM+1; index++)
    {
        switch (nodesRegionsForElement[index])
        {
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_WALL:
            case HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM:
                lv++;
                break;

            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_WALL:
            case HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM:
                rv++;
                break;

            case HeartGeometryInformation<SPACE_DIM>::UNKNOWN:
            default:
                NEVER_REACHED;
        }
    }

    // If most of the nodes are in the right ventricle
    if (rv>lv)
    {
        return M_PI/4.0;
    }

    // Anywhere else
    return M_PI/3.0;
}

template<unsigned SPACE_DIM>
StreeterFibreGenerator<SPACE_DIM>::StreeterFibreGenerator(AbstractTetrahedralMesh<SPACE_DIM,SPACE_DIM>& rMesh)
    : AbstractPerElementWriter<SPACE_DIM,SPACE_DIM,SPACE_DIM*SPACE_DIM>(&rMesh),
      mpGeometryInfo(NULL),
      mApexToBase(zero_vector<double>(SPACE_DIM)),
      mLogInfo(false)
{
    // Record a reference for the calculations performed here, can be extracted with the '-citations' flag.
    Citations::Register(StreeterCitation, &StreeterCite);

    mWallThickness.resize(rMesh.GetNumNodes());
    mAveragedWallThickness.resize(rMesh.GetNumNodes());
}

template<unsigned SPACE_DIM>
StreeterFibreGenerator<SPACE_DIM>::~StreeterFibreGenerator()
{
    delete mpGeometryInfo;
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetSurfaceFiles(
            const std::string& rEpicardiumFile,
            const std::string& rRightVentricleFile,
            const std::string& rLeftVentricleFile,
            bool indexFromZero)
{
    // Compute the distance map of each surface
     mpGeometryInfo = new HeartGeometryInformation<SPACE_DIM>(*(this->mpMesh), rEpicardiumFile, rLeftVentricleFile, rRightVentricleFile, indexFromZero);
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::WriteHeaderOnMaster()
{
    *(this->mpMasterFile) << this->mpMesh->GetNumElements();
    *(this->mpMasterFile) << std::setprecision(16);
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::PreWriteCalculations(OutputFileHandler& rOutputDirectory)
{
    assert(SPACE_DIM == 3);
    if (mpGeometryInfo == NULL)
    {
        EXCEPTION("Files defining the heart surfaces not set");
    }

    // Open files
    out_stream p_regions_file, p_thickness_file, p_ave_thickness_file;

    //Make sure that only the master process writes the log files if requested.
    bool logInfo = PetscTools::AmMaster() && mLogInfo;

    if (logInfo)
    {
        p_regions_file  = rOutputDirectory.OpenOutputFile("node_regions.data");
        p_thickness_file = rOutputDirectory.OpenOutputFile("wall_thickness.data");
        p_ave_thickness_file = rOutputDirectory.OpenOutputFile("averaged_thickness.data");
    }

    //We expect that the apex to base has been set
    if (fabs(norm_2(mApexToBase)) < DBL_EPSILON)
    {
        EXCEPTION("Apex to base vector has not been set");
    }

    // Compute wall thickness parameter
    unsigned num_nodes = this->mpMesh->GetNumNodes();
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        double dist_epi, dist_endo;

        HeartRegionType node_region = mpGeometryInfo->GetHeartRegion(node_index);

        switch(node_region)
        {
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::LEFT_VENTRICLE_WALL:
                dist_epi = mpGeometryInfo->rGetDistanceMapEpicardium()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_SURFACE:
            case HeartGeometryInformation<SPACE_DIM>::RIGHT_VENTRICLE_WALL:
                dist_epi = mpGeometryInfo->rGetDistanceMapEpicardium()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::LEFT_SEPTUM:
                dist_epi = mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::RIGHT_SEPTUM:
                dist_epi = mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index];
                dist_endo = mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index];
                break;

            case HeartGeometryInformation<SPACE_DIM>::UNKNOWN:
                // LCOV_EXCL_START
                std::cerr << "Wrong distances node: " << node_index << "\t"
                          << "Epi " << mpGeometryInfo->rGetDistanceMapEpicardium()[node_index] << "\t"
                          << "RV " << mpGeometryInfo->rGetDistanceMapRightVentricle()[node_index] << "\t"
                          << "LV " << mpGeometryInfo->rGetDistanceMapLeftVentricle()[node_index]
                          << std::endl;

                // Make wall_thickness=0 as in Martin's code
                dist_epi = 1;
                dist_endo = 0;
                break;
                // LCOV_EXCL_STOP

            default:
                NEVER_REACHED;
        }

        mWallThickness[node_index] = dist_endo / (dist_endo + dist_epi);

        if (std::isnan(mWallThickness[node_index]))
        {
            // LCOV_EXCL_START
            /*
             *  A node contained on both epicardium and lv (or rv) surfaces has wall thickness 0/0.
             *  By setting its value to 0 we consider it contained only on the lv (or rv) surface.
             */
            mWallThickness[node_index] = 0;
            // LCOV_EXCL_STOP
        }

        if (logInfo)
        {
            *p_regions_file << node_region*100 << "\n";
            *p_thickness_file << mWallThickness[node_index] << "\n";
        }
    }

    /*
     *  For each node, average its value of e with the values of all the neighbours
     */
    std::vector<double> my_averaged_wall_thickness(num_nodes);
    for (unsigned node_index=0; node_index<num_nodes; node_index++)
    {
        my_averaged_wall_thickness[node_index] = GetAveragedThicknessLocalNode(node_index, mWallThickness);
    }

    // Non-local information appear as zeros in the vector
    MPI_Allreduce(&my_averaged_wall_thickness[0], &mAveragedWallThickness[0], num_nodes,
                  MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);

    if (logInfo)
    {
        for (unsigned node_index=0; node_index<num_nodes; node_index++)
        {
             *p_ave_thickness_file << mAveragedWallThickness[node_index] << "\n";
        }
    }

    if (logInfo)
    {
        p_regions_file->close();
        p_thickness_file->close();
        p_ave_thickness_file->close();
    }
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::Visit(Element<SPACE_DIM, SPACE_DIM>* pElement,
                                              unsigned localElementIndex,
                                              c_vector<double, SPACE_DIM*SPACE_DIM>& rData)
{
    /*
     *  The gradient of the averaged thickness at the element is:
     *
     *     grad_ave_wall_thickness[element_index] = ave' * BF * inv(J)
     *
     *  being : ave, averaged thickness values of the nodes defining the element
     *          J,   the Jacobian of the element as defined in class Element.
     *                               (-1 -1 -1)
     *          BF,  basis functions ( 1  0  0)
     *                               ( 0  1  0)
     *                               ( 0  0  1)
     *
     *  Defined as u in Streeter paper.
     */
    c_vector<double, SPACE_DIM> grad_ave_wall_thickness;
    c_vector<double, SPACE_DIM+1> elem_nodes_ave_thickness;
    double element_averaged_thickness = 0.0;
    c_vector<HeartRegionType, SPACE_DIM+1> elem_nodes_region;

    for (unsigned local_node_index=0; local_node_index<SPACE_DIM+1; local_node_index++)
    {
        // Get node's global index
        unsigned global_node_index = pElement->GetNode(local_node_index)->GetIndex();

        elem_nodes_ave_thickness[local_node_index] = mAveragedWallThickness[global_node_index];
        elem_nodes_region[local_node_index] = mpGeometryInfo->GetHeartRegion(global_node_index);

        // Calculate wall thickness averaged value for the element
        element_averaged_thickness +=  mWallThickness[global_node_index];
    }

    element_averaged_thickness /= SPACE_DIM+1;

    c_matrix<double, SPACE_DIM+1, SPACE_DIM> basis_functions( zero_matrix<double>(4u,3u) );
    basis_functions(0,0) = basis_functions(0,1) = basis_functions(0,2) = -1.0;
    basis_functions(1,0) = basis_functions(2,1) = basis_functions(3,2) =  1.0;

    c_matrix<double, SPACE_DIM+1, SPACE_DIM> temp;
    c_matrix<double, SPACE_DIM, SPACE_DIM> jacobian, inverse_jacobian;
    double jacobian_det;
    unsigned element_index = pElement->GetIndex();
    this->mpMesh->GetInverseJacobianForElement(element_index, jacobian, jacobian_det, inverse_jacobian);
    noalias(temp) = prod (basis_functions, inverse_jacobian);
    noalias(grad_ave_wall_thickness) = prod(elem_nodes_ave_thickness, temp);

    grad_ave_wall_thickness /= norm_2(grad_ave_wall_thickness);

    /*
     * Normal to the gradient (v in Streeter paper) which is then the circumferential direction
     * (it will be the fibre direction after rotation)
     *
     * Computed as the cross product with the base-apex direction (originally assumed base-apex axis is x). The output vector is not normal,
     * since the angle between them may be != 90, normalise it.
     */
    c_vector<double, SPACE_DIM> fibre_direction = VectorProduct(grad_ave_wall_thickness, mApexToBase);
    fibre_direction /= norm_2(fibre_direction);

    /*
     *  Longitude direction (w in Streeter paper)
     */
    c_vector<double, SPACE_DIM> longitude_direction = VectorProduct(grad_ave_wall_thickness, fibre_direction);

    /*
     *  Compute fibre to v angle: alpha = R*(1-2e)^3
     *
     *    R is the maximum angle between the fibre and the v axis (heart region dependant)
     *    (1 - 2e)^3 scales it by a value in [-1, 1] defining the rotation of the fibre based
     *       on the position in the wall
     */
    double alpha = GetFibreMaxAngle(elem_nodes_region) * SmallPow( (1 - 2*element_averaged_thickness), 3 );

    /*
     *  Apply alpha rotation about the u axis to the orthonormal basis
     *
     *               ( u(1) v(1) w(1) )
     *   (u, v, w) = ( u(2) v(2) w(2) )
     *               ( u(3) v(3) w(3) )
     *
     *  The following matrix defines a rotation about the u axis
     *
     *                 ( 1        0           0      ) (u')
     *   R = (u, v, w) ( 0    cos(alpha) -sin(alpha) ) (v')
     *                 ( 0    sin(alpha)  cos(alpha) ) (w')
     *
     *  The rotated basis is computed like:
     *
     *                                             ( 1        0           0      )
     *  (u, v_r, w_r ) = R * (u, v, w) = (u, v, w) ( 0    cos(alpha) -sin(alpha) )
     *                                             ( 0    sin(alpha)  cos(alpha) )
     *
     *  Which simplifies to:
     *
     *   v_r =  v*cos(alpha) + w*sin(alpha)
     *   w_r = -v*sin(alpha) + w*cos(alpha)
     */
    c_vector<double, SPACE_DIM> rotated_fibre_direction = fibre_direction*cos(alpha) + longitude_direction*sin(alpha);
    c_vector<double, SPACE_DIM> rotated_longitude_direction = -fibre_direction*sin(alpha) + longitude_direction*cos(alpha);


    /*
     * Test the orthonormality of the basis
     */
    assert( fabs(norm_2(rotated_fibre_direction) - 1) < 100*DBL_EPSILON );
    assert( fabs(norm_2(grad_ave_wall_thickness) - 1) < 100*DBL_EPSILON );
    assert( fabs(norm_2(rotated_longitude_direction) - 1) < 100*DBL_EPSILON );

    assert( fabs(inner_prod(rotated_fibre_direction, grad_ave_wall_thickness)) < 100*DBL_EPSILON );
    assert( fabs(inner_prod(rotated_fibre_direction, rotated_longitude_direction)) < 100*DBL_EPSILON);
    assert( fabs(inner_prod(grad_ave_wall_thickness, rotated_longitude_direction)) < 100*DBL_EPSILON);

    /*
     *  Output the direction of the myofibre, the transverse to it in the plane
     *  of the myocite laminae and the normal to this laminae (in that order)
     *
     *  See Fig. 1 "Laminar Structure of the Heart: a mathematical model" LeGrice et al. 97
     *
     */
    rData[0] = rotated_fibre_direction[0];
    rData[1] = rotated_fibre_direction[1];
    rData[2] = rotated_fibre_direction[2];
    rData[3] = grad_ave_wall_thickness[0];
    rData[4] = grad_ave_wall_thickness[1];
    rData[5] = grad_ave_wall_thickness[2];
    rData[6] = rotated_longitude_direction[0];
    rData[7] = rotated_longitude_direction[1];
    rData[8] = rotated_longitude_direction[2];
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetApexToBase(const c_vector<double, SPACE_DIM>& apexToBase)
{
    double norm = norm_2(apexToBase);
    if (norm < DBL_EPSILON)
    {
        EXCEPTION("Apex to base vector should be non-zero");
    }
    mApexToBase = apexToBase / norm;
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetApexToBase(unsigned axis)
{
    if (axis >= SPACE_DIM)
    {
        EXCEPTION("Apex to base coordinate axis was out of range");
    }
    mApexToBase = zero_vector<double>(SPACE_DIM);
    mApexToBase[axis] = 1.0;
}

template<unsigned SPACE_DIM>
void StreeterFibreGenerator<SPACE_DIM>::SetLogInfo(bool logInfo)
{
    mLogInfo = logInfo;
}

// Explicit instantiation
// LCOV_EXCL_START
template class StreeterFibreGenerator<3>;
// LCOV_EXCL_STOP
