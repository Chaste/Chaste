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

#ifndef ELECTRODES_HPP_
#define ELECTRODES_HPP_

#include <boost/shared_ptr.hpp>
#include "ChasteSerialization.hpp"
#include <boost/serialization/split_member.hpp>
#include <boost/serialization/shared_ptr.hpp>

#include "AbstractTetrahedralMesh.hpp"
#include "DistributedVector.hpp"
#include "BoundaryConditionsContainer.hpp"
#include "ConstBoundaryCondition.hpp"

/**
 *  A class for setting up the boundary conditions associated with electrodes.
 *  There are two modes: grounding the second electrode, in which case the first
 *  electrode has an input flux (Neumann boundary condition) of the specified
 *  magnitude, the extracellular potential is fixed on the second electrode; or
 *  not grounded, in which case the opposite electrode has an equal and opposite
 *  output flux.
 *
 *  This class assumes the given mesh is cuboid, and the electrodes are taken to be
 *  the specified opposite surfaces.
 *
 *  Note the class now includes a pointer to the mesh as a member variable, as this
 *  is required to archive (reconstruct) itself. Because boost is clever it will
 *  still only archive one copy of the mesh.
 */
template<unsigned DIM>
class Electrodes
{
friend class TestBidomainWithBathProblem;
friend class TestElectrodes;

private:
    /** Whether the second electrode is grounded */
    bool mGroundSecondElectrode;
    /** The created bcc, which BidomainProblem will use */
    boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,2> > mpBoundaryConditionsContainer;
    /** The time the electrodes are switched on */
    double mStartTime;
    /** The time the electrodes are switched off */
    double mEndTime;
    /** Whether the electrodes are currently switched on */
    bool mAreActive;
    /**
     * This is needed for archiving:
     * the boundary conditions refer to nodes and/or elements, so need
     * the mesh to be archived, but don't have a pointer to the mesh itself.
     */
    AbstractTetrahedralMesh<DIM,DIM>* mpMesh;

    /** Left electrode area*/
    double mLeftElectrodeArea;

    /** Right electrode area*/
    double mRightElectrodeArea;

    /**
     * Helper method to compute electrodes area and check they are equal. Throws if they are not.
     *
     *  @param index The value i when applying the electrodes to x_i=a and x_i=b (a<b)
     *  @param lowerValue The value a when applying the electrodes to x_i=a and x_i=b (a<b) (should
     *    be the minimum value of x_i for the given mesh)
     *  @param upperValue The value b when applying the electrodes to x_i=a and x_i=b (a<b) (should
     *    be the maximum value of x_i for the given mesh)
     */
    void ComputeElectrodesAreasAndCheckEquality(unsigned index, double lowerValue, double upperValue);

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Save the Electrodes class
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void save(Archive & archive, const unsigned int version) const
    {
        archive & mGroundSecondElectrode;
        archive & mStartTime;
        archive & mEndTime;
        archive & mAreActive;
        archive & mpMesh;
        (*ProcessSpecificArchive<Archive>::Get()) & mpBoundaryConditionsContainer;
    }
    /**
     * Load the Electrodes class
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void load(Archive & archive, const unsigned int version)
    {
        archive & mGroundSecondElectrode;
        archive & mStartTime;
        archive & mEndTime;
        archive & mAreActive;
        archive & mpMesh;

        (*ProcessSpecificArchive<Archive>::Get()) & mpBoundaryConditionsContainer;
        assert(mpBoundaryConditionsContainer);
        // Allow the new object to load itself from the archive.
        mpBoundaryConditionsContainer->LoadFromArchive(*ProcessSpecificArchive<Archive>::Get(), mpMesh);
    }
    BOOST_SERIALIZATION_SPLIT_MEMBER()

    /**
     * Private default constructor for archiving only.
     */
    Electrodes(){};

public:


    /** Constructor.
     * Needs only a reference to a mesh.
     * All other parameters from the HeartConfig class
     *
     *  @param rMesh The mesh, assumed to be a cuboid.
     */
    Electrodes(AbstractTetrahedralMesh<DIM,DIM>& rMesh); // implemented in cpp

    /**
     *  @return the boundary conditions container in which is set up the Neumann
     *  fluxes for the first electrode, and the opposite fluxes for the second
     *  electrode if the the second electrode isn't grounded
     */
    boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,2> > GetBoundaryConditionsContainer();

    /**
     *  @return whether it is time to switch off the electrodes yet. THIS ONLY RETURNS
     *  TRUE ONCE - the first appropriate time. After that the electrodes assume
     *  they have been switched off and therefore this returns false.
     *
     * @param time  the current time
     */
    bool SwitchOff(double time);

    /**
     *  @return whether it is time to switch on the electrodes yet. THIS ONLY RETURNS
     *  TRUE ONCE - the first appropriate time. After that the electrodes assume
     *  they have been switched on and therefore this returns false.
     *
     * @param time  the current time
     */
    bool SwitchOn(double time);

    /** @return the time the electrodes are switched on */
    double GetSwitchOnTime()
    {
        return mStartTime;
    }

    /** @return the time the electrodes are switched off */
    double GetSwitchOffTime()
    {
        return mEndTime;
    }

    /** @return whether the second electrode is grounded or not */
    bool HasGroundedElectrode()
    {
        return mGroundSecondElectrode;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Electrodes)

#endif /*ELECTRODES_HPP_*/
