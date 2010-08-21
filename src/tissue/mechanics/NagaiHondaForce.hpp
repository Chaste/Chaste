/*

Copyright (C) University of Oxford, 2005-2010

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
#ifndef NAGAIHONDAFORCE_HPP_
#define NAGAIHONDAFORCE_HPP_


#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "AbstractForce.hpp"
#include "VertexBasedTissue.hpp"

#include <iostream>

/**
 * A force class for use in vertex-based tissue simulations, based on a mechanical
 * model proposed by T. Nagai and H. Honda ("A dynamic cell model for the formation
 * of epithelial tissues", Philosophical Magazine Part B 81:699-719).
 * 
 * Each of the model parameter member variables are rescaled such that mDampingConstantNormal
 * takes the default value 1, whereas Nagai and Honda (who denote the parameter by
 * nu) take the value 0.01.
 */
template<unsigned DIM>
class NagaiHondaForce  : public AbstractForce<DIM>
{
friend class TestForcesNotForRelease;

private:

    /*
     * The following four parameters are used in vertex-based tissue simulations
     * based on the mechanical model proposed by T. Nagai and H. Honda ("A dynamic
     * cell model for the formation of epithelial tissues", Philosophical Magazine
     * Part B 81:699-719). They are rescaled such that mDampingConstantNormal takes
     * the default value 1, whereas Nagai and Honda (who denote the parameter by nu)
     * take the value 0.01.
     */

    /**
     * Cell deformation energy parameter. Has units of kg s^-2 (cell size at equilibrium rest length)^-1.
     */
    double mNagaiHondaDeformationEnergyParameter;

    /**
     * Cell membrane energy parameter. Has units of kg (cell size at equilibrium rest length) s^-2.
     */
    double mNagaiHondaMembraneSurfaceEnergyParameter;

    /**
     * Cell-cell adhesion energy parameter. Has has units of kg (cell size at equilibrium rest length)^2 s^-2.
     */
    double mNagaiHondaCellCellAdhesionEnergyParameter;

    /**
     * Cell-boundary adhesion energy parameter. Has units of kg (cell size at equilibrium rest length)^2 s^-2.
     */
    double mNagaiHondaCellBoundaryAdhesionEnergyParameter;

    /**
     * Non-dimensional target area of a mature (fully-grown) TissueCell.
     */
    double mMatureCellTargetArea;

    friend class boost::serialization::access;
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then '&' resolves to '<<'
        // If Archive is an input archive, then '&' resolves to '>>'
        archive & boost::serialization::base_object<AbstractForce<DIM> >(*this);
        archive & mNagaiHondaDeformationEnergyParameter;
        archive & mNagaiHondaMembraneSurfaceEnergyParameter;
        archive & mNagaiHondaCellCellAdhesionEnergyParameter;
        archive & mNagaiHondaCellBoundaryAdhesionEnergyParameter;
        archive & mMatureCellTargetArea;
    }

public:

    /**
     * Constructor.
     */
    NagaiHondaForce();

    /**
     * Destructor.
     */
    ~NagaiHondaForce();

    /**
     * Overridden AddForceContribution() method.
     *
     * Calculates the force on each node in the vertex-based tissue based on the
     * Nagai Honda model.
     *
     * @param rForces reference to vector of forces on nodes
     * @param rTissue reference to the tissue
     */
    void AddForceContribution(std::vector<c_vector<double, DIM> >& rForces, AbstractTissue<DIM>& rTissue);

    /**
     * Get the adhesion parameter for the edge between two given nodes.
     *
     * @param pNodeA one node
     * @param pNodeB the other node
     *
     * @return the adhesion parameter for this edge.
     */
    double GetAdhesionParameter(Node<DIM>* pNodeA, Node<DIM>* pNodeB);

    /**
     * @return mNagaiHondaDeformationEnergyParameter
     */
    double GetNagaiHondaDeformationEnergyParameter();

    /**
     * @return mNagaiHondaMembraneSurfaceEnergyParameter
     */
    double GetNagaiHondaMembraneSurfaceEnergyParameter();

    /**
     * @return mCellCellAdhesionEnergyParameter
     */
    double GetNagaiHondaCellCellAdhesionEnergyParameter();

    /**
     * @return mNagaiHondaCellBoundaryAdhesionEnergyParameter
     */
    double GetNagaiHondaCellBoundaryAdhesionEnergyParameter();

    /**
     * Set mNagaiHondaDeformationEnergyParameter.
     * 
     * @param nagaiHondaDeformationEnergyParameter the new value of mNagaiHondaDeformationEnergyParameter
     */
    void SetNagaiHondaDeformationEnergyParameter(double nagaiHondaDeformationEnergyParameter);

    /**
     * Set mNagaiHondaMembraneSurfaceEnergyParameter.
     * 
     * @param nagaiHondaMembraneSurfaceEnergyParameter the new value of mNagaiHondaMembraneSurfaceEnergyParameter
     */
    void SetNagaiHondaMembraneSurfaceEnergyParameter(double nagaiHondaMembraneSurfaceEnergyParameter);

    /**
     * Set mNagaiHondaCellCellAdhesionEnergyParameter.
     * 
     * @param nagaiHondaCellCellAdhesionEnergyEnergyParameter the new value of mNagaiHondaCellCellAdhesionEnergyParameter
     */
    void SetNagaiHondaCellCellAdhesionEnergyParameter(double nagaiHondaCellCellAdhesionEnergyEnergyParameter);

    /**
     * Set mNagaiHondaCellBoundaryAdhesionEnergyParameter.
     * 
     * @param nagaiHondaCellBoundaryAdhesionEnergyParameter the new value of mNagaiHondaCellBoundaryAdhesionEnergyParameter
     */
    void SetNagaiHondaCellBoundaryAdhesionEnergyParameter(double nagaiHondaCellBoundaryAdhesionEnergyParameter);

    /**
     * Get the target area of a given cell. This grows linearly from
     * 0.5*A to A during the G1 phase of the cell cycle, then remains
     * at A for the rest of the cell cycle, where A denotes the TissueConfig
     * member variable mMatureCellTargetArea.
     *
     * @param pCell the cell
     * @return the cell's target area
     */
    double GetTargetAreaOfCell(const TissueCellPtr pCell) const;

    /**
     * @return mMatureCellTargetArea
     */
    double GetMatureCellTargetArea() const;

    /**
     * Set mMatureCellTargetArea.
     * 
     * @param matureCellTargetArea the new value of mMatureCellTargetArea
     */
    void SetMatureCellTargetArea(double matureCellTargetArea);

    /**
     * Outputs force Parameters to file
     *
     * As this method is pure virtual, it must be overridden
     * in subclasses.
     *
     * @param rParamsFile the file stream to which the parameters are output
     */
    void OutputForceParameters(out_stream& rParamsFile);
};


#include "SerializationExportWrapper.hpp"

EXPORT_TEMPLATE_CLASS_SAME_DIMS(NagaiHondaForce)

#endif /*NAGAIHONDAFORCE_HPP_*/
