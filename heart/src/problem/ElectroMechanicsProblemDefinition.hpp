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

#ifndef ELECTROMECHANICSPROBLEMDEFINITION_HPP_
#define ELECTROMECHANICSPROBLEMDEFINITION_HPP_

#include "SolidMechanicsProblemDefinition.hpp"
#include "ContractionModelName.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "CompressibleExponentialLaw.hpp"



/**
 *  Subclass of SolidMechanicsProblemDefinition with some cardiac-electro-mechanics-specific
 *  methods.
 */
template<unsigned DIM>
class ElectroMechanicsProblemDefinition : public SolidMechanicsProblemDefinition<DIM>
{
private:
    /**
     *  The contraction model used (ContractionModelName is an enumeration containing all contraction
     *  models implemented.
     */
    ContractionModelName mContractionModel;

    /** Timestep to use when solving contraction models */
    double mContractionModelOdeTimeStep;

    /** How often a mechanics solve should be done */
    double mMechanicsSolveTimestep;

    /**
     *  Whether the deformation should affect the electrical physiological conductivity
     *  (or whether this effect is neglected)
     */
    bool mDeformationAffectsConductivity;

    /**
     *  Whether the deformation should affect the cardiac cell models, for example if there
     *  are stretch-activated channels in the cell model.
     */
    bool mDeformationAffectsCellModels;

    /**
     *  This member variable is used if SetDefaultCardiacMateriawLaw() is called.
     */
    AbstractMaterialLaw<DIM>* mpDefaultMaterialLaw;


    /** Whether to read fibre-sheet information from file */
    bool mReadFibreSheetInformationFromFile;

    /**
     * .ortho/.orthoquad file from which to read element-wise, or quadrature-point-wise
     * fibre-sheet-normal-directions
     */
    std::string mFibreSheetDirectionsFile;

    /**
     * Whether the mFibreSheetDirectionsFile file gives the fibre-sheet info for each element
     * or for each quadrature point
     */
    bool mFibreSheetDirectionsDefinedPerQuadraturePoint;

public:
    /**
     * Constructor
     * @param rMesh the mesh
     */
    ElectroMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh);

    /** Destructor */
    ~ElectroMechanicsProblemDefinition();

    /**
     * Set the contraction model to be used (throughout the tissue).
     *
     * Note the timestep should be set to a (typical) ODE time-step even if the contraction
     * model is not going to solve ODEs.
     *
     * @param contractionModel contraction model (from the enumeration ContractionModelName)
     * @param timestep timestep to be used in solving (ODE-based) contraction models.
     *
     */
    void SetContractionModel(ContractionModelName contractionModel, double timestep);

    /**
     * Use the default material law (NashHunter in the incompressible case, exponential in the
     *  compressible case), throughout the tissue.
     * @param compressibilityType Either INCOMPRESSIBLE or COMPRESSIBLE
     */
    void SetUseDefaultCardiacMaterialLaw(CompressibilityType compressibilityType);

    /**
     * Set if and how the deformation should affect the electro-physiology.
     *
     * @param deformationAffectsConductivity Whether the deformation should affect the electrical
     *   physiological conductivity (or whether this effect is neglected)
     * @param deformationAffectsCellModels Whether the deformation should affect the cardiac cell
     *   models, for example if there are stretch-activated channels in the cell model.
     *
     * Several important things to note:
     *   (i) this can't be called if fibre-sheet directions have been defined from file for each quadrature
     *  point (as opposed to each mechanics element) - this is because if the stretch is to be passed back to
     *  the electric mesh nodes, the fibre direction has to be defined at those nodes
     *   (ii) currently the set-up stage (computing mechanics mesh elements and weights for electrics mesh
     *  nodes) is inefficiently implemented - setup will be very slow for big meshes
     *   (iii) if deformationAffectsCellModels is true, the cell model ought to be one for which
     *   AbstractCardiacCell::SetStretch() has been implemented to do something (i.e. not generated automatically
     *   from CellML).
     *   (iv) deformationAffectsConductivity is not currently allowed in the compressible material law case
     *   as the effect of the determinant of the deformation gradient on the conductivity has not currently been
     *   implemented.
     *
     */
    void SetDeformationAffectsElectrophysiology(bool deformationAffectsConductivity, bool deformationAffectsCellModels);

    /**
     *  Set how often the mechanics is solved for.
     *  @param timestep timestep
     */
    void SetMechanicsSolveTimestep(double timestep);


    /**
     *  Set a variable fibre-sheet-normal direction (matrices), from file.
     *  If the second parameter is false, there should be one fibre-sheet definition for each element; otherwise
     *  there should be one fibre-sheet definition for each *quadrature point* in the mesh.
     *  In the first case, the file should be a standard .ortho file (ie each line has the fibre dir, sheet dir,
     *  normal dir for that element), in the second it should have .orthoquad as the format.
     *
     *  If this method is not called, the default fibre-sheet directions are used - ie fibres parallel to
     *  X-axis, sheets parallel to Y-axis.
     *
     *  @param fibreSheetDirectionsFile the file containing the fibre/sheet directions
     *  @param definedPerQuadPoint whether the fibre-sheet definitions are for each quadrature point in the mesh
     *   (if not, one for each element is assumed).
     */
    void SetVariableFibreSheetDirectionsFile(std::string fibreSheetDirectionsFile, bool definedPerQuadPoint);


    /**
     *  Get the contraction model
     */
    ContractionModelName GetContractionModel()
    {
        assert(mContractionModelOdeTimeStep>0.0); // if this fails SetContractionModel() probably hasn't been called - call Validate() to check
        return mContractionModel;
    }

    /**
     *  Get the contraction model timestep
     */
    double GetContractionModelOdeTimestep()
    {
        assert(mContractionModelOdeTimeStep>0.0); // if this fails SetContractionModel() probably hasn't been called - call Validate() to check
        return mContractionModelOdeTimeStep;
    }

    /**
     *  Get how often the mechanics is solved
     */
    double GetMechanicsSolveTimestep()
    {
        assert(mMechanicsSolveTimestep>0.0); // if this fails SetMechanicsSolveTimestep() probably hasn't been called - call Validate() to check
        return mMechanicsSolveTimestep;
    }

    /**
     *  Get whether the deformation affects the electrical physiological conductivity
     *  (or whether this effect is neglected).
     */
    bool GetDeformationAffectsConductivity()
    {
        return mDeformationAffectsConductivity;
    }

    /**
     *  Get whether the deformation affects the cardiac cell models, for example if
     *  there are stretch-activated channels in the cell model.
     */
    bool GetDeformationAffectsCellModels()
    {
        return mDeformationAffectsCellModels;
    }

    /**
     *  Whether the fibre-sheet info should be read from file (if not the defaults should be used).
     */
    bool ReadFibreSheetDirectionsFromFile()
    {
        return mReadFibreSheetInformationFromFile;
    }

    /**
     *  Get the fibre-sheet file (should only be called if ReadFibreSheetDirectionsFromFile() returns true).
     */
    std::string GetFibreSheetDirectionsFile()
    {
        assert(mReadFibreSheetInformationFromFile);
        assert(mFibreSheetDirectionsFile!="");
        return mFibreSheetDirectionsFile;
    }

    /**
     *  Get whether the fibre-sheet info is defined for each quadrature point in the mesh (if not,
     *  if it defined for each element in the mesh). (Should only be called if
     *  ReadFibreSheetDirectionsFromFile() returns true).
     */
    bool GetFibreSheetDirectionsDefinedPerQuadraturePoint()
    {
        assert(mReadFibreSheetInformationFromFile);
        return mFibreSheetDirectionsDefinedPerQuadraturePoint;
    }

    /**
     * Check all variables are set appropriately. Exceptions are thrown if any are not.
     * Derived classes can override but should call this version as well.
     */
    virtual void Validate();
};

#endif // ELECTROMECHANICSPROBLEMDEFINITION_HPP_
