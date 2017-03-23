/*
 * BielmeierForce.hpp
 *
 *  Created on: 23 Mar 2017
 *      Author: Weijie
 */

#ifndef BIELMEIERFORCE_HPP_
#define BIELMEIERFORCE_HPP_

#include "GeneralMonolayerVertexMeshForce.hpp"

/**
 * A force class for use in vertex-based model simulations. This force is based on the
 * energy function proposed by Bielmeier et al in the following paper:
 *
 * Christina Bielmeier, Silvanus Alt, Vanessa Weichselberger, Marco La Fortezza,
 * Hartmann Harz, Frank Jülicher, Guillaume Salbreux, Anne-Kathrin Classen.
 * Interface Contractility between Differently Fated Cells Drives Cell Elimination and Cyst Formation.
 * Current Biology, 26(5), pp.563-574.
 * http://dx.doi.org/10.1016/j.cub.2015.12.063
 */
class BielmeierForce : public GeneralMonolayerVertexMeshForce
{
private:
    friend class boost::serialization::access;
    /**
     * Boost Serialization method for archiving/checkpointing.
     * Archives the object and its member variables.
     *
     * @param archive  The boost archive.
     * @param version  The current version of this class.
     */
    template <class Archive>
    void serialize(Archive& archive, const unsigned int version)
    {
        archive& boost::serialization::base_object<GeneralMonolayerVertexMeshForce>(*this);
        archive& mEcmSpringConstant;
        archive& mExternalSurfaceTensionParameter;
    }

protected:
    /**
     * Spring modulus of elastic bond attaching the tissue to the external extracellular matrix (ECM)
     */
    double mEcmSpringConstant;

    /**
     * External in-plane surface tension parameter which constrains the area of tissue.
     */
    double mExternalSurfaceTensionParameter;

public:
    /**
     * Default constructor.
     */
    BielmeierForce(const double springConstant = 5, const double ExternalSurfaceTensionParameter = -4.2);

    /**
     * Destructor.
     */
    virtual ~BielmeierForce()
    {
    }

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculate the force on each node in the vertex-based cell population based on the energy function in Bielmeier et al's paper.
     *
     * @param rCellPopulation reference to the cell population
     */
    virtual void AddForceContribution(AbstractCellPopulation<3>& rCellPopulation);

    /**
     * Overridden OutputForceParameters() method.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(BielmeierForce)

#endif /*BIELMEIERFORCE_HPP_*/
