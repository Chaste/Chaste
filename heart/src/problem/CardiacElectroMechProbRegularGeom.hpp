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
        : CardiacElectroMechanicsProblem<DIM,1>(compressibilityType,
                                              MONODOMAIN,
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

        // Fix the nodes on x=0
        std::vector<unsigned> fixed_nodes;
        for (unsigned i=0; i<this->mpMechanicsMesh->GetNumNodes(); i++)
        {
            if (fabs(this->mpMechanicsMesh->GetNode(i)->rGetLocation()[0])<1e-6)
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
