/*

Copyright (c) 2005-2023, University of Oxford.
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

#include <map>
#include "PeriodicNodesOnlyMesh.hpp"

template<unsigned SPACE_DIM>
PeriodicNodesOnlyMesh<SPACE_DIM>::PeriodicNodesOnlyMesh(c_vector<double,SPACE_DIM> width)
        : NodesOnlyMesh<SPACE_DIM>(),
          mWidth(width)
{
    // Convert periodic widths into a boolean vector
    mIsDimPeriodic = zero_vector<bool>(SPACE_DIM);

    for (unsigned i = 0; i < SPACE_DIM; i++)
    {
        if (width[i]>0.0)
        {
            mIsDimPeriodic[i] = true;
        }
    }
}

template<unsigned SPACE_DIM>
void PeriodicNodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(double cutOffLength, c_vector<double, 2*SPACE_DIM> domainSize, int numLocalRows, c_vector<bool,SPACE_DIM> isDimPeriodic)
{

    // Check that we haven't got too many processors; the maximum number of processors is width/cutOffLength in the y (2D)/z (3D) direction
    if ( mIsDimPeriodic[SPACE_DIM-1] && ((mWidth[SPACE_DIM-1]+1e-14)/cutOffLength) < PetscTools::GetNumProcs() )
    {
        EXCEPTION("Too many processors for the periodic domain width and cut off length. NumProcs should be less than or equal to Y (in 2D) or Z (in 3D) width / Cut off length.\n");
    }
    // Container for converting the periodic dims vector to a true/false input into box collection setup
    // For each periodic dimension:

    for (unsigned i=0; i < SPACE_DIM; i++)
    {
        if (mIsDimPeriodic[i])
        {
            // Ensure that the width is a multiple of cut-off length
            if ( fmod(mWidth[i]+1e-14,cutOffLength ) > 2e-14 )
            {
                EXCEPTION("The periodic width must be a multiple of cut off length.");
            }
            // A width of two boxes gives different simulation results as some connections are considered twice.
            else if ( mWidth[i]/cutOffLength <= (2.0+1e-14) )
            {
                EXCEPTION( "The periodic domain width cannot be less than 2*CutOffLength." );
            }

            // We force the domain to the periodic widths
            domainSize[ i*2 ] = 0;
            domainSize[ i*2+1 ] = mWidth[i];
        }
    }
    NodesOnlyMesh<SPACE_DIM>::SetUpBoxCollection(cutOffLength, domainSize, PETSC_DECIDE, mIsDimPeriodic);

    this->AddNodesToBoxes();
}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PeriodicNodesOnlyMesh<SPACE_DIM>::GetPeriodicWidths() const
{
    return mWidth;
}

template<unsigned SPACE_DIM>
double PeriodicNodesOnlyMesh<SPACE_DIM>::GetWidth(const unsigned& rDimension) const
{
    double width = 0.0;
    assert(rDimension < SPACE_DIM);

    if ( !mIsDimPeriodic[rDimension] )
    {
        // If we are not a periodic dimension, we just return the width
        width = MutableMesh<SPACE_DIM,SPACE_DIM>::GetWidth(rDimension);
    }
    else
    {
        width = mWidth[rDimension];
    }
    return width;
}

template<unsigned SPACE_DIM>
void PeriodicNodesOnlyMesh<SPACE_DIM>::SetNode(unsigned nodeIndex, ChastePoint<SPACE_DIM> point, bool concreteMove)
{
    // concreteMove should always be false for NodesOnlyMesh as no elements to check
    assert(!concreteMove);

    c_vector<double, SPACE_DIM>& point_coord = point.rGetLocation();

    // Loop over the periodic dimensions and perform a periodic movement if necessary
    for ( unsigned i = 0; i < SPACE_DIM; i++ )
    {
        if ( mIsDimPeriodic[i] )
        {
            double local_point_coord = point_coord(i);
            double local_width = mWidth[i];
            // We need to make sure that the fudge factor is greater than the one used to calculate the containing
            // box in distributed box collection. Otherwise the forced domain width causes issues at the right domain edge.
            double fudge_factor = 1.0e-13;
            if ( local_point_coord >= local_width )
            {
                double new_coord = local_point_coord - local_width;

                point.SetCoordinate(i,new_coord);
            }
            else if ( local_point_coord > (local_width-fudge_factor))
            {
                // This is to ensure that the position is never equal to mWidth, which would be outside the box domain.
                // This is due to the fact that mWidth-1e-16=mWidth
                point.SetCoordinate(i,local_width-fudge_factor);
            }
            else if ( local_point_coord < 0.0 )
            {
                double new_coord = local_point_coord + local_width;

                // This is to ensure that the position is never equal to mWidth, which would be outside the box domain.
                // This is due to the fact that mWidth-1e-16=mWidth
                if ( new_coord > local_width-fudge_factor )
                {
                    new_coord = local_width-fudge_factor;
                }
                point.SetCoordinate(i, new_coord);
            }
        }
    }

    // Update the node's location
    this->GetNode(nodeIndex)->SetPoint(point);
}

template<unsigned SPACE_DIM>
unsigned PeriodicNodesOnlyMesh<SPACE_DIM>::AddNode(Node<SPACE_DIM>* pNewNode)
{
    // Call method on parent class
    unsigned node_index = NodesOnlyMesh<SPACE_DIM>::AddNode(pNewNode);

    // If necessary move it to be back on the cylinder
    ChastePoint<SPACE_DIM> new_node_point = pNewNode->GetPoint();
    SetNode(node_index, new_node_point);

    return node_index;

}

template<unsigned SPACE_DIM>
c_vector<double, SPACE_DIM> PeriodicNodesOnlyMesh<SPACE_DIM>::GetVectorFromAtoB(const c_vector<double, SPACE_DIM>& rLocation1, const c_vector<double, SPACE_DIM>& rLocation2)
{
    c_vector<double, SPACE_DIM> vector = rLocation2 - rLocation1;

    // Loop over the periodic dimensions and measure across boundaries
    // if points are more than halfway across the width apart
    for ( unsigned i = 0; i < SPACE_DIM; i++ )
    {
        if ( mIsDimPeriodic[i] )
        {
            vector(i) = fmod(vector(i),mWidth[i]);
            if ( vector(i) > 0.5*mWidth[i] )
            {
                vector(i) -= mWidth[i];
            }
            else if (vector(i) < -0.5*mWidth[i] )
            {
                vector(i) += mWidth[i];
            }
        }
    }

    return vector;
}

// Refresh mesh -> Check and then move if outside box
template<unsigned SPACE_DIM>
void PeriodicNodesOnlyMesh<SPACE_DIM>::RefreshMesh()
{
    // Loop over the periodic dimensions and if they are
    // outside the domain, get the fmod and relocate
    for ( typename std::vector<Node<SPACE_DIM> *>::iterator it_node = this->mNodes.begin();
            it_node != this->mNodes.end(); ++it_node )
    {
        c_vector<double,SPACE_DIM> & location = (*it_node)->rGetModifiableLocation();

        for ( unsigned i = 0; i < SPACE_DIM; i++ )
        {
            if ( mIsDimPeriodic[i] )
            {
                if ( location(i) < 0.0 )
                {
                    location(i) = fmod( location(i), mWidth[i] ) + mWidth[i];
                }
                else if ( location(i) >= mWidth[i] )
                {
                    location(i) = fmod( location(i), mWidth[i] );
                }
            }
        }
    }

    // Now run the base class method
    NodesOnlyMesh<SPACE_DIM>::RefreshMesh();
}

//////////////////////////////////////////////////////////////////////////
// Explicit instantiation
//////////////////////////////////////////////////////////////////////////

template class PeriodicNodesOnlyMesh<1>;
template class PeriodicNodesOnlyMesh<2>;
template class PeriodicNodesOnlyMesh<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PeriodicNodesOnlyMesh)
