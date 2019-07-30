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

#include "Electrodes.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "IsNan.hpp"
#include "HeartConfig.hpp"

template<unsigned DIM>
Electrodes<DIM>::Electrodes(AbstractTetrahedralMesh<DIM,DIM>& rMesh)
    : mpMesh(&rMesh),
      mLeftElectrodeArea(0.0),
      mRightElectrodeArea(0.0)
{
    unsigned axis_index;
    double magnitude, duration;

    HeartConfig::Instance()->GetElectrodeParameters(mGroundSecondElectrode, axis_index, magnitude, mStartTime, duration);

    assert(axis_index < DIM);
    assert(duration > 0);
    mEndTime = mStartTime + duration;
    mAreActive = false; // consider electrodes initially switched off!

    ChasteCuboid<DIM> bounding_box = mpMesh->CalculateBoundingBox();
    double global_min = bounding_box.rGetLowerCorner()[axis_index];
    double global_max = bounding_box.rGetUpperCorner()[axis_index];

    mpBoundaryConditionsContainer.reset(new BoundaryConditionsContainer<DIM,DIM,2>);


    double input_flux;
    double output_flux;

    try
    {
        ComputeElectrodesAreasAndCheckEquality(axis_index, global_min, global_max);
        input_flux = magnitude;
        output_flux = -input_flux;

    }
    catch (Exception&)
    {
        // magnitude of second electrode scaled so that left_area*magnitude_left = -right_area*magnitude_right
        input_flux = magnitude;
        output_flux = -input_flux*mLeftElectrodeArea/mRightElectrodeArea;

        // Paranoia. In case one of the areas is 0
        assert( ! std::isnan(output_flux));
        assert( output_flux != 0.0);
    }

    ConstBoundaryCondition<DIM>* p_bc_flux_in = new ConstBoundaryCondition<DIM>(input_flux);
    ConstBoundaryCondition<DIM>* p_bc_flux_out = new ConstBoundaryCondition<DIM>(output_flux);

    // loop over boundary elements and add a non-zero phi_e boundary condition (ie extracellular
    // stimulus) if (assuming axis_index=0, etc) x=lowerValue (where x is the x-value of the centroid)
    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
            = mpMesh->GetBoundaryElementIteratorBegin();
       iter != mpMesh->GetBoundaryElementIteratorEnd();
       iter++)
    {
        if (fabs((*iter)->CalculateCentroid()[axis_index] - global_min) < 1e-6)
        {
            mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_in,  1);
        }

        if (!mGroundSecondElectrode)
        {
            if (fabs((*iter)->CalculateCentroid()[axis_index] - global_max) < 1e-6)
            {
                mpBoundaryConditionsContainer->AddNeumannBoundaryCondition(*iter, p_bc_flux_out, 1);
            }
        }
    }

    // set up mGroundedNodes using opposite surface ie second electrode is
    // grounded
    if (mGroundSecondElectrode)
    {
        ConstBoundaryCondition<DIM>* p_zero_bc = new ConstBoundaryCondition<DIM>(0.0);
        // We will need to recalculate this when HasDirichletBoundaryConditions() is called.
        this->mpBoundaryConditionsContainer->ResetDirichletCommunication();

        for (typename AbstractTetrahedralMesh<DIM,DIM>::NodeIterator iter=mpMesh->GetNodeIteratorBegin();
             iter != mpMesh->GetNodeIteratorEnd();
             ++iter)
        {
            if (fabs((*iter).rGetLocation()[axis_index]-global_max) < 1e-6)
            {
                mpBoundaryConditionsContainer->AddDirichletBoundaryCondition(&(*iter), p_zero_bc, 1);
            }
        }

        //Unused boundary conditions will not be deleted by the b.c. container
        delete p_bc_flux_out;
    }
}

template<unsigned DIM>
boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,2> > Electrodes<DIM>::GetBoundaryConditionsContainer()
{
    //assert(mAreActive);
    return mpBoundaryConditionsContainer;
}

template<unsigned DIM>
bool Electrodes<DIM>::SwitchOff(double time)
{
    // This smidge has to be the same as the one below
    ///\todo remove magic number? (#1884)
    double smidge = 1e-10;
    if (mAreActive && time>mEndTime-smidge)
    {
        mAreActive = false;
        return true;
    }

    return false;
}

template<unsigned DIM>
bool Electrodes<DIM>::SwitchOn(double time)
{
    // This smidge has to be the same as the one above.
    ///\todo remove magic number? (#1884)
    double smidge = 1e-10;
    if (!mAreActive && time>=mStartTime && time<=mEndTime - smidge)
    {
        mAreActive = true;
        return true;
    }

    return false;
}

template<unsigned DIM>
void Electrodes<DIM>::ComputeElectrodesAreasAndCheckEquality(unsigned dimensionIndex, double lowerValue, double upperValue)
{
    //
    // loop over boundary elements and determine areas of the electrode faces
    //
    double local_left_area = 0.0;
    double local_right_area = 0.0;

    c_vector<double,DIM> weighted_direction;
    double jacobian_determinant;

    for (typename AbstractTetrahedralMesh<DIM,DIM>::BoundaryElementIterator iter
             = mpMesh->GetBoundaryElementIteratorBegin();
         iter != mpMesh->GetBoundaryElementIteratorEnd();
         iter++)
    {
        if (mpMesh->CalculateDesignatedOwnershipOfBoundaryElement((*iter)->GetIndex()))
        {
            if (fabs((*iter)->CalculateCentroid()[dimensionIndex] - lowerValue) < 1e-6)
            {
                mpMesh->GetWeightedDirectionForBoundaryElement((*iter)->GetIndex(), weighted_direction, jacobian_determinant);
                local_left_area += jacobian_determinant;
            }

            if (fabs((*iter)->CalculateCentroid()[dimensionIndex] - upperValue) < 1e-6)
            {
                mpMesh->GetWeightedDirectionForBoundaryElement((*iter)->GetIndex(), weighted_direction, jacobian_determinant);
                local_right_area += jacobian_determinant;
            }
        }
    }

    if (DIM == 3)
    {
        // if the dimension of the face is 1, the mapping is from the face to the canonical element [0,1], so the
        // jacobian_determinants used above will be exactly the areas of the faces
        // if the dimension of the face is 1, the mapping is from the face to the canonical element, the triangle
        // with nodes (0,0), (0,1), (1,0), which has area 0.5, so scale the jacobian_determinants by 0.5 to
        // get the face areas, ie scale the total area by 0.5:
        // (Not technically needed as it is only the ratio of these that is used but since we are calling the variables
        // 'area's we ought to be correct).
        local_left_area /= 2.0;
        local_right_area /= 2.0;
    }

    int mpi_ret = MPI_Allreduce(&local_left_area, &mLeftElectrodeArea, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    UNUSED_OPT(mpi_ret);
    assert(mpi_ret == MPI_SUCCESS);

    mpi_ret = MPI_Allreduce(&local_right_area, &mRightElectrodeArea, 1, MPI_DOUBLE, MPI_SUM, PETSC_COMM_WORLD);
    assert(mpi_ret == MPI_SUCCESS);

    if (mLeftElectrodeArea != mRightElectrodeArea)
    {
        EXCEPTION("Electrodes have different area");
    }
}

// Explicit instantiation
template class Electrodes<1>;
template class Electrodes<2>;
template class Electrodes<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Electrodes)
