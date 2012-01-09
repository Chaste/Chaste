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

#ifndef CARDIACELECTROMECHPROBREGULARGEOM_HPP_
#define CARDIACELECTROMECHPROBREGULARGEOM_HPP_

#include "CardiacElectroMechanicsProblem.hpp"

/**
 *  Child class of CardiacElectroMechanicsProblem for setting up cardiac electromechanics
 *  problems on a square (currently just 2d). The user just has to specify the number
 *  of elements in each direction.
 *
 *  This class can only be used to set up highly-restricted problems (square geometry, fibres
 *  in X-direction, X=0 side fixed in space, no tractions, pressures or MEF) - use the more general class
 *  CardiacElectroMechanicsProblem if you want to do anything more complex.
 *
 *  Note: the X=0 surface is fixed in the deformation.
 */
template<unsigned DIM>
class CardiacElectroMechProbRegularGeom : public CardiacElectroMechanicsProblem<DIM>
{
public:

    /**
     * Constructor.
     *
     * @param compressibilityType Should be either INCOMPRESSIBLE or COMPRESSIBLE
     * @param width Width and height of the square
     * @param numMechanicsElementsEachDir Number of elements in each direction in the mechanics mesh
     * @param numElectricsElementsEachDir Number of elements in each direction in the electrics mesh
     * @param pCellFactory factory to use to create cells
     * @param contractionModel contraction model (see the enum "ContractionModel" for the options).
     * @param mechanicsSolveTimestep how often the mechanics is solved for (should be a multiple of electrics PDE timestep)
     * @param contractionModelOdeTimeStep Step size for contraction model (of active tension in cardiac cells) being used.
     * @param outputDirectory the output directory
     */
    CardiacElectroMechProbRegularGeom(CompressibilityType compressibilityType,
                                      double width,
                                      unsigned numMechanicsElementsEachDir,
                                      unsigned numElectricsElementsEachDir,
                                      AbstractCardiacCellFactory<DIM>* pCellFactory,
                                      ContractionModelName contractionModel,
                                      double mechanicsSolveTimestep,
                                      double contractionModelOdeTimeStep,
                                      std::string outputDirectory = "")
        : CardiacElectroMechanicsProblem<DIM>(compressibilityType,
                                              NULL, NULL,
                                              pCellFactory,
                                              NULL,
                                              outputDirectory)
    {
        assert(DIM==2); // the below assumes DIM==2

        assert(width > 0.0);
        assert(numMechanicsElementsEachDir > 0);
        assert(numElectricsElementsEachDir > 0);

        // create electrics mesh
        this->mpElectricsMesh = new TetrahedralMesh<DIM,DIM>();

        this->mpElectricsMesh->ConstructRegularSlabMesh(width/numElectricsElementsEachDir, width, width);

        // create mechanics mesh
        this->mpMechanicsMesh = new QuadraticMesh<DIM>(width/numMechanicsElementsEachDir, width, width);
        LOG(2, "Width of meshes is " << width);
        LOG(2, "Num nodes in electrical and mechanical meshes are: " << this->mpElectricsMesh->GetNumNodes() << ", " << this->mpMechanicsMesh->GetNumNodes() << "\n");

        this->mpProblemDefinition = new ElectroMechanicsProblemDefinition<DIM>(*(this->mpMechanicsMesh));

        // fix the nodes on x=0
        std::vector<unsigned> fixed_nodes;
        for(unsigned i=0; i<this->mpMechanicsMesh->GetNumNodes(); i++)
        {
            if( fabs(this->mpMechanicsMesh->GetNode(i)->rGetLocation()[0])<1e-6)
            {
                fixed_nodes.push_back(i);
            }
        }
        this->mpProblemDefinition->SetZeroDisplacementNodes(fixed_nodes);

        LOG(2, "Fixed the " << fixed_nodes.size() << " nodes on x=0");

        this->mpProblemDefinition->SetUseDefaultCardiacMaterialLaw(compressibilityType);
        this->mpProblemDefinition->SetContractionModel(contractionModel,contractionModelOdeTimeStep);
        this->mpProblemDefinition->SetMechanicsSolveTimestep(mechanicsSolveTimestep);
    }

    ~CardiacElectroMechProbRegularGeom()
    {
        delete this->mpElectricsMesh;
        delete this->mpMechanicsMesh;
        delete this->mpProblemDefinition;
    }
};

#endif /*CARDIACELECTROMECHPROBREGULARGEOM_HPP_*/
