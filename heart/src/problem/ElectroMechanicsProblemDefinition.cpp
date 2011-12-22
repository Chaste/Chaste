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

#include "ElectroMechanicsProblemDefinition.hpp"

template<unsigned DIM>
ElectroMechanicsProblemDefinition<DIM>::ElectroMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh)
    : SolidMechanicsProblemDefinition<DIM>(rMesh),
      mContractionModelOdeTimeStep(-1.0),
      mMechanicsSolveTimestep(-1.0),
      mDeformationAffectsConductivity(false),
      mDeformationAffectsCellModels(false),
      mpDefaultMaterialLaw(NULL),
      mReadFibreSheetInformationFromFile(false)
{
}

template<unsigned DIM>
ElectroMechanicsProblemDefinition<DIM>::~ElectroMechanicsProblemDefinition()
{
    if(mpDefaultMaterialLaw)
    {
        delete mpDefaultMaterialLaw;
    }
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetContractionModel(ContractionModelName contractionModel, double timestep)
{
    assert(timestep > 0.0);
    mContractionModel = contractionModel;
    mContractionModelOdeTimeStep = timestep;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetUseDefaultCardiacMaterialLaw(CompressibilityType compressibilityType)
{
    if(mpDefaultMaterialLaw)
    {
        delete mpDefaultMaterialLaw;
    }

    if(compressibilityType == INCOMPRESSIBLE)
    {
        mpDefaultMaterialLaw = new NashHunterPoleZeroLaw<DIM>();
        this->SetMaterialLaw(INCOMPRESSIBLE, mpDefaultMaterialLaw);
    }
    else
    {
        mpDefaultMaterialLaw = new CompressibleExponentialLaw<DIM>();
        this->SetMaterialLaw(COMPRESSIBLE, mpDefaultMaterialLaw);
    }
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetDeformationAffectsElectrophysiology(bool deformationAffectsConductivity, bool deformationAffectsCellModels)
{
    mDeformationAffectsConductivity = deformationAffectsConductivity;
    mDeformationAffectsCellModels = deformationAffectsCellModels;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetMechanicsSolveTimestep(double timestep)
{
    assert(timestep > 0.0);
    mMechanicsSolveTimestep = timestep;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetVariableFibreSheetDirectionsFile(std::string fibreSheetDirectionsFile, bool definedPerQuadraturePoint)
{
    mReadFibreSheetInformationFromFile = true;
    mFibreSheetDirectionsFile = fibreSheetDirectionsFile;
    mFibreSheetDirectionsDefinedPerQuadraturePoint = definedPerQuadraturePoint;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::Validate()
{
    SolidMechanicsProblemDefinition<DIM>::Validate();

    if(mMechanicsSolveTimestep < 0.0)
    {
        EXCEPTION("Timestep for mechanics solve hasn't been set yet");
    }

    if(mContractionModelOdeTimeStep < 0.0)
    {
        EXCEPTION("Contraction model hasn't been set yet");
    }

    if(mDeformationAffectsConductivity && this->GetCompressibilityType()==COMPRESSIBLE)
    {
        // the conductivity depends on the deformation gradient and also scales in some way with
        // J=det(F), which is not equal to 1 in the compressible case. The F dependence
        // is implemented but the J dependence is not yet.
        EXCEPTION("Deformation affecting the conductivity is currently not implemented fully for compressible problems");
    }

    if(mDeformationAffectsCellModels && mReadFibreSheetInformationFromFile && mFibreSheetDirectionsDefinedPerQuadraturePoint)
    {
        // This combination is not allowed. For explanation see doxygen for SetDeformationAffectsElectrophysiology()
        std::stringstream message;
        message << "Deformation affecting cell models cannot be done when fibres-sheet information is defined for each quadrature point.";
        message << "Define fibre-sheet information for each element instead.";
        EXCEPTION(message.str());
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class ElectroMechanicsProblemDefinition<2>;
template class ElectroMechanicsProblemDefinition<3>;
