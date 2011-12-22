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
     *  Get the boundary conditions container in which is set up the Neumann
     *  fluxes for the first electrode, and the opposite fluxes for the second
     *  electrode if the the second electrode isn't grounded
     */
    boost::shared_ptr<BoundaryConditionsContainer<DIM,DIM,2> > GetBoundaryConditionsContainer();

    /**
     *  Whether it is time to switch off the electrodes yet. THIS ONLY RETURNS
     *  TRUE ONCE - the first appropriate time. After that the electrodes assume
     *  they have been switched off and therefore this returns false.
     *
     * @param time  the current time
     */
    bool SwitchOff(double time);

    /**
     *  Whether it is time to switch on the electrodes yet. THIS ONLY RETURNS
     *  TRUE ONCE - the first appropriate time. After that the electrodes assume
     *  they have been switched on and therefore this returns false.
     *
     * @param time  the current time
     */
    bool SwitchOn(double time);

    /** Get the time the electrodes are switched on */
    double GetSwitchOnTime()
    {
        return mStartTime;
    }

    /** Get the time the electrodes are switched off */
    double GetSwitchOffTime()
    {
        return mEndTime;
    }

    /** Whether the second electrode is grounded or not */
    bool HasGroundedElectrode()
    {
        return mGroundSecondElectrode;
    }
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(Electrodes)

#endif /*ELECTRODES_HPP_*/
