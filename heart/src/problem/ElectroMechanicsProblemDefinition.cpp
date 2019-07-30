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

#include "ElectroMechanicsProblemDefinition.hpp"
#include "LabelBasedContractionCellFactory.hpp"

template<unsigned DIM>
ElectroMechanicsProblemDefinition<DIM>::ElectroMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh)
    : SolidMechanicsProblemDefinition<DIM>(rMesh),
      mContractionModelOdeTimeStep(-1.0),
      mMechanicsSolveTimestep(-1.0),
      mDeformationAffectsConductivity(false),
      mDeformationAffectsCellModels(false),
      mpDefaultMaterialLaw(NULL),
      mReadFibreSheetInformationFromFile(false),
      mNumIncrementsForInitialDeformation(1),
      mApplyCrossFibreTension(false),
      mSheetTensionFraction(DOUBLE_UNSET),
      mSheetNormalTensionFraction(DOUBLE_UNSET),
      mpContractionCellFactory(NULL),
      mWeMadeCellFactory(false),
      mSolverType(IMPLICIT) // default solver is implicit
{
}

template<unsigned DIM>
ElectroMechanicsProblemDefinition<DIM>::~ElectroMechanicsProblemDefinition()
{
    if (mpDefaultMaterialLaw)
    {
        delete mpDefaultMaterialLaw;
    }

    if (mWeMadeCellFactory)
    {
        delete mpContractionCellFactory;
    }
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetContractionModel(ContractionModelName contractionModel, double timestep)
{
    assert(timestep > 0.0);
    SetContractionModelOdeTimestep(timestep);

    if (contractionModel == NASH2004 || contractionModel == CONSTANT)
    {
        // These models can use an Explicit solver, default is Implicit.
        SetSolverType(EXPLICIT);
    }

    // Make sure we aren't overwriting a problem that has been set up with a cell factory.
    assert(mpContractionCellFactory==NULL);

    AbstractContractionCellFactory<DIM>* p_factory = new LabelBasedContractionCellFactory<DIM>(contractionModel);
    mWeMadeCellFactory = true;
    SetContractionCellFactory(p_factory);
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetUseDefaultCardiacMaterialLaw(CompressibilityType compressibilityType)
{
    if (mpDefaultMaterialLaw)
    {
        delete mpDefaultMaterialLaw;
    }

    if (compressibilityType == INCOMPRESSIBLE)
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
void ElectroMechanicsProblemDefinition<DIM>::SetVariableFibreSheetDirectionsFile(const FileFinder& rFibreSheetDirectionsFile, bool definedPerQuadraturePoint)
{
    mReadFibreSheetInformationFromFile = true;
    mFibreSheetDirectionsFile = rFibreSheetDirectionsFile;
    mFibreSheetDirectionsDefinedPerQuadraturePoint = definedPerQuadraturePoint;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetApplyIsotropicCrossFibreTension(bool applyCrossFibreTension, double crossFibreTensionFraction)
{
    mApplyCrossFibreTension = applyCrossFibreTension;
    mSheetTensionFraction = crossFibreTensionFraction;
    mSheetNormalTensionFraction = crossFibreTensionFraction;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetApplyAnisotropicCrossFibreTension(bool applyCrossFibreTension, double sheetTensionFraction, double sheetNormalTensionFraction)
{
    if (DIM!=3)
    {
        EXCEPTION("You can only apply anisotropic cross fibre tensions in a 3D simulation.");
    }
    mApplyCrossFibreTension = applyCrossFibreTension;
    mSheetTensionFraction = sheetTensionFraction;
    mSheetNormalTensionFraction = sheetNormalTensionFraction;
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::SetContractionCellFactory(AbstractContractionCellFactory<DIM>* pCellFactory)
{
    // Make sure we aren't overwriting a problem that has been set up with a cell factory already.
    assert(mpContractionCellFactory == NULL);

    mpContractionCellFactory = pCellFactory;
    mpContractionCellFactory->SetMechanicsMesh(static_cast<QuadraticMesh<DIM>*>(&(this->mrMesh)));
}

template<unsigned DIM>
void ElectroMechanicsProblemDefinition<DIM>::Validate()
{
    SolidMechanicsProblemDefinition<DIM>::Validate();

    if (mMechanicsSolveTimestep < 0.0)
    {
        EXCEPTION("Timestep for mechanics solve hasn't been set yet");
    }

    if (mContractionModelOdeTimeStep < 0.0)
    {
        std::string message =  "Contraction model or contraction model ODE timestep have not been set. "
                               "Make sure SetContractionModel(), or SetContractionCellFactory() AND SetContractionModelOdeTimestep "
                               "are called. (Pass in a timestep even if contraction model is algebraic and won't use it). ";
        EXCEPTION(message);
    }

    if (mDeformationAffectsConductivity && this->GetCompressibilityType()==COMPRESSIBLE)
    {
        // the conductivity depends on the deformation gradient and also scales in some way with
        // J=det(F), which is not equal to 1 in the compressible case. The F dependence
        // is implemented but the J dependence is not yet.
        EXCEPTION("Deformation affecting the conductivity is currently not implemented fully for compressible problems");
    }

    if (mDeformationAffectsCellModels && mReadFibreSheetInformationFromFile && mFibreSheetDirectionsDefinedPerQuadraturePoint)
    {
        // This combination is not allowed. For explanation see doxygen for SetDeformationAffectsElectrophysiology()
        std::stringstream message;
        message << "Deformation affecting cell models cannot be done when fibres-sheet information is defined for each quadrature point.";
        message << "Define fibre-sheet information for each element instead.";
        EXCEPTION(message.str());
    }
}

// Explicit instantiation
template class ElectroMechanicsProblemDefinition<2>;
template class ElectroMechanicsProblemDefinition<3>;
