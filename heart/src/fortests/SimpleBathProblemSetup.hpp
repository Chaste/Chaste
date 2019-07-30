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
#ifndef SIMPLEBATHPROBLEMSETUP_HPP_
#define SIMPLEBATHPROBLEMSETUP_HPP_

/**
 * @file
 *
 * Some helper classes and functions for setting up a simple bath problem for testing.
 */

#include "AbstractCardiacCellFactory.hpp"
#include "SimpleStimulus.hpp"
#include "LuoRudy1991.hpp"
#include "HeartRegionCodes.hpp"

/**
 * A simple cell factory for bath problems, applying a SimpleStimulus for
 * 0.5ms at a single point.
 */
template<unsigned DIM, class CELLTYPE=CellLuoRudy1991FromCellML>
class BathCellFactory : public AbstractCardiacCellFactory<DIM>
{
private:
    /** The stimulus to apply */
    boost::shared_ptr<SimpleStimulus> mpStimulus;
    /** and where to apply it */
    c_vector<double,DIM> mStimulatedPoint;

public:
    /**
     * @return a newly created cell factory.
     * @param stimulusMagnitude
     * @param stimulatedPoint spatial co-ordinates of where to stimulate.
     *    Must correspond to a node location.
     */
    BathCellFactory(double stimulusMagnitude, c_vector<double,DIM> stimulatedPoint)
        : AbstractCardiacCellFactory<DIM>(),
          mpStimulus(new SimpleStimulus(stimulusMagnitude, 0.5)),
          mStimulatedPoint(stimulatedPoint)
    {
    }

    /**
     * @return a newly created new cell.
     * @param pNode  pointer to Node object for cell.
     */
    AbstractCardiacCellInterface* CreateCardiacCellForTissueNode(Node<DIM>* pNode)
    {
        // paranoia - check this is really a tissue node
        assert(HeartRegionCode::IsRegionTissue( pNode->GetRegion() ));

        // stimulate centre node normally..
        bool is_centre;

        if (DIM==1)
        {
            is_centre = (fabs(pNode->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6);
        }
        else if (DIM==2)
        {
            is_centre = (    (fabs(pNode->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(pNode->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6) );
        }
        else
        {
            is_centre = (    (fabs(pNode->GetPoint()[0]-mStimulatedPoint(0)) < 1e-6)
                          && (fabs(pNode->GetPoint()[1]-mStimulatedPoint(1)) < 1e-6)
                          && (fabs(pNode->GetPoint()[2]-mStimulatedPoint(2)) < 1e-6) );
        }

        if (is_centre)
        {
            return new CELLTYPE(this->mpSolver, mpStimulus);
        }
        else
        {
            return new CELLTYPE(this->mpSolver, this->mpZeroStimulus);
        }
    }
};

/**
 * Set everything outside a central circle in the given 2d mesh to be bath.
 *
 * @param pMesh the mesh
 * @param centreX X co-ord of tissue centre
 * @param centreY Y co-ord of tissue centre
 * @param radius radius of tissue
 */
template<class MeshType>
void SetCircularTissueIn2dMesh(MeshType* pMesh,
                               double centreX, double centreY, double radius)
{
    for (typename MeshType::ElementIterator it = pMesh->GetElementIteratorBegin();
         it != pMesh->GetElementIteratorEnd();
         ++it)
    {
        double x = it->CalculateCentroid()[0];
        double y = it->CalculateCentroid()[1];
        if ((x-centreX)*(x-centreX) + (y-centreY)*(y-centreY) > radius*radius)
        {
            it->SetAttribute(HeartRegionCode::GetValidBathId());
        }
    }
    pMesh->SetMeshHasChangedSinceLoading();
}

/**
 * Load a 2d mesh, and set everything outside a central circle to be bath.
 *
 * @param rMeshPath relative path to the mesh
 * @param centreX X co-ord of tissue centre
 * @param centreY Y co-ord of tissue centre
 * @param radius radius of tissue
 * @return the new mesh
 */
template<class MeshType>
MeshType* Load2dMeshAndSetCircularTissue(const std::string& rMeshPath,
                                         double centreX, double centreY, double radius)
{
    TrianglesMeshReader<2,2> reader(rMeshPath);
    MeshType* p_mesh = new MeshType;
    p_mesh->ConstructFromMeshReader(reader);

    SetCircularTissueIn2dMesh(p_mesh, centreX, centreY, radius);

    return p_mesh;
}

/**
 * Specialization for a parallel mesh.
 *
 * @param rMeshPath relative path to the mesh
 * @param centreX X co-ord of tissue centre
 * @param centreY Y co-ord of tissue centre
 * @param radius radius of tissue
 * @return the new mesh
 */
template<>
DistributedTetrahedralMesh<2,2>* Load2dMeshAndSetCircularTissue(const std::string& rMeshPath,
                                                                double centreX, double centreY, double radius)
{
    TrianglesMeshReader<2,2> reader(rMeshPath);
    // Force dumb partitioning so migration tests pass!
    DistributedTetrahedralMesh<2,2>* p_mesh = new DistributedTetrahedralMesh<2,2>(DistributedTetrahedralMeshPartitionType::DUMB);
    p_mesh->ConstructFromMeshReader(reader);

    SetCircularTissueIn2dMesh(p_mesh, centreX, centreY, radius);

    return p_mesh;
}

#endif /*SIMPLEBATHPROBLEMSETUP_HPP_*/
