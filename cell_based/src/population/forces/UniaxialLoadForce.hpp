/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef UNIAXIALLOADFORCE_HPP_
#define UNIAXIALLOADFORCE_HPP_

#include <vector>

#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>

#include "AbstractForce.hpp"

/**
 * A force that loads a specimen in uniaxial compression or tension, following the bulk rheology
 * protocol of Sandersius & Newman (2008) Phys. Biol. 5 015002, Section 3.1.
 *
 * Two slabs are identified, one at each end of the specimen along the loading axis, each
 * consisting of the nodes lying within a given thickness of that end. A total load F is then
 * shared equally between the nodes of each slab, so that a node of the upper slab feels
 * F/N_upper and a node of the lower slab feels F/N_lower, directed towards each other for
 * compression and away from each other for tension.
 *
 * Because each slab receives a total of exactly F, the two loads cancel and the net force on the
 * specimen is zero: it deforms without translating, and no boundary condition is needed to hold
 * it in place.
 *
 * Net torque is not cancelled. Two equal and opposite loads form a couple whose arm is the
 * lateral offset between the two slab centroids, so a specimen whose ends are not symmetric may
 * rotate slowly. For the meshes this was developed against the offset is either exactly zero or a
 * few hundredths of a node spacing, but it is worth checking in general.
 *
 * The loading axis is unconstrained in the other directions: the specimen is free to bulge
 * laterally as it is compressed.
 *
 * Slab membership is determined once, from the configuration in force the first time
 * AddForceContribution() is called, and is then held fixed - as a rheometer plate grips a
 * specimen. This keeps the load from silently redistributing as nodes drift across the slab
 * boundary, but it does mean the force assumes the set of nodes is stable, so it is not suitable
 * for a simulation in which cells divide or die.
 */
template<unsigned DIM>
class UniaxialLoadForce : public AbstractForce<DIM>
{
private:

    /** The total load applied to each slab. */
    double mLoad;

    /** The axis along which the specimen is loaded. */
    unsigned mLoadingAxis;

    /** The thickness of each slab, measured inwards from each end of the specimen. */
    double mSlabThickness;

    /** Whether the slabs are loaded towards each other (compression) or apart (tension). */
    bool mIsCompressive;

    /** Indices of the nodes in the slab at the upper end of the loading axis. */
    std::vector<unsigned> mUpperSlabNodeIndices;

    /** Indices of the nodes in the slab at the lower end of the loading axis. */
    std::vector<unsigned> mLowerSlabNodeIndices;

    /** Whether the slabs have been identified yet. */
    bool mSlabsIdentified;

    /**
     * Populate the two slab node index lists from the current configuration, and mark the slabs
     * as identified. Called once, from the first AddForceContribution().
     *
     * @param rCellPopulation reference to the cell population
     */
    void IdentifySlabs(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Archiving.
     */
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * The identified slabs are archived along with the parameters, so that a restarted simulation
     * continues to load the same nodes rather than re-identifying them from the deformed
     * configuration.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mLoad;
        archive & mLoadingAxis;
        archive & mSlabThickness;
        archive & mIsCompressive;
        archive & mUpperSlabNodeIndices;
        archive & mLowerSlabNodeIndices;
        archive & mSlabsIdentified;
    }

public:

    /**
     * Constructor. The load and slab thickness both default to zero and must be set before the
     * force is used; the loading axis defaults to the last one, i.e. y in 2D and z in 3D.
     */
    UniaxialLoadForce();

    /**
     * Destructor.
     */
    virtual ~UniaxialLoadForce() = default;

    /**
     * Set the total load applied to each slab.
     *
     * @param load the total load, which must be positive
     */
    void SetLoad(double load);

    /**
     * @return the total load applied to each slab
     */
    double GetLoad() const;

    /**
     * Set the axis along which the specimen is loaded.
     *
     * @param loadingAxis the axis, which must be less than DIM
     */
    void SetLoadingAxis(unsigned loadingAxis);

    /**
     * @return the axis along which the specimen is loaded
     */
    unsigned GetLoadingAxis() const;

    /**
     * Set the thickness of each slab, measured inwards from each end of the specimen. A thickness
     * of about one node spacing gives a slab one node deep.
     *
     * @param slabThickness the slab thickness, which must be positive
     */
    void SetSlabThickness(double slabThickness);

    /**
     * @return the thickness of each slab
     */
    double GetSlabThickness() const;

    /**
     * Set whether the specimen is compressed or stretched.
     *
     * @param isCompressive whether the slabs are loaded towards each other
     */
    void SetIsCompressive(bool isCompressive);

    /**
     * @return whether the slabs are loaded towards each other
     */
    bool GetIsCompressive() const;

    /**
     * The number of nodes in the slab at the upper end of the loading axis. Zero until the slabs
     * have been identified, which happens on the first call to AddForceContribution().
     *
     * @return the number of nodes in the upper slab
     */
    unsigned GetUpperSlabNodeCount() const;

    /**
     * The number of nodes in the slab at the lower end of the loading axis. Zero until the slabs
     * have been identified, which happens on the first call to AddForceContribution().
     *
     * @return the number of nodes in the lower slab
     */
    unsigned GetLowerSlabNodeCount() const;

    /**
     * Overridden AddForceContribution() method.
     *
     * @param rCellPopulation reference to the cell population
     */
    void AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

#include "SerializationExportWrapper.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(UniaxialLoadForce)

#endif /*UNIAXIALLOADFORCE_HPP_*/
