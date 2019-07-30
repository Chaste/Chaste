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

#ifndef ELECTROMECHANICSPROBLEMDEFINITION_HPP_
#define ELECTROMECHANICSPROBLEMDEFINITION_HPP_

#include "SolidMechanicsProblemDefinition.hpp"
#include "ContractionModelName.hpp"
#include "NashHunterPoleZeroLaw.hpp"
#include "CompressibleExponentialLaw.hpp"
#include "FileFinder.hpp"
#include "AbstractContractionCellFactory.hpp"
#include "SolverType.hpp"

/**
 *  Subclass of SolidMechanicsProblemDefinition with some cardiac-electro-mechanics-specific
 *  methods.
 */
template<unsigned DIM>
class ElectroMechanicsProblemDefinition : public SolidMechanicsProblemDefinition<DIM>
{
private:
    friend class TestExplicitCardiacMechanicsSolver;
    friend class TestImplicitCardiacMechanicsSolver;
    friend class TestElectroMechanicsProblemDefinition;

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
    FileFinder mFibreSheetDirectionsFile;

    /**
     * Whether the mFibreSheetDirectionsFile file gives the fibre-sheet info for each element
     * or for each quadrature point
     */
    bool mFibreSheetDirectionsDefinedPerQuadraturePoint;

    /**
     *  The first deformation (to find the equilibrium state given the loading)
     *  may require the loading to be incremented, in order for the Solve() to
     *  converge - this stores the number of increments to be used (initialised to 1)
     */
    unsigned mNumIncrementsForInitialDeformation;

    /**
     * Whether to apply active tension in the cross-fibre directions (either iso- or aniso-tropic).
     */
    bool mApplyCrossFibreTension;

    /**
     * The fraction of the fibre tension to apply in the sheet direction.
     */
    double mSheetTensionFraction;

    /**
     * The fraction of the fibre tension to apply in the sheet-normal direction.
     */
    double mSheetNormalTensionFraction;

    /** The contraction cell factory, provides the contraction model to use on each element */
    AbstractContractionCellFactory<DIM>* mpContractionCellFactory;

    /** Whether we made the cell factory and therefore need to delete it! */
    bool mWeMadeCellFactory;

    /** Whether to use an explicit or implicit solver */
    SolverType mSolverType;

public:
    /**
     * Constructor
     * @param rMesh the mesh
     */
    ElectroMechanicsProblemDefinition(QuadraticMesh<DIM>& rMesh);

    /** Destructor */
    virtual ~ElectroMechanicsProblemDefinition();

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
     * We provide an interface to set heterogeneous contraction cells via a
     * cell factory, which allows different models (or parameters) to be set
     * per element.
     *
     * Note you also need to call #SetContractionModelOdeTimestep() before solving.
     *
     * @param pCellFactory  The contraction cell factory to be used.
     */
    void SetContractionCellFactory(AbstractContractionCellFactory<DIM>* pCellFactory);

    /**
     * Use the default material law (NashHunter in the incompressible case, exponential in the
     *  compressible case), throughout the tissue.
     * @param compressibilityType Either INCOMPRESSIBLE or COMPRESSIBLE
     */
    void SetUseDefaultCardiacMaterialLaw(CompressibilityType compressibilityType);

    /**
     * Decide which type of solver to attempt to use. Default value is
     * IMPLICIT if this method is not called.
     *
     * @param solver The solver type to use (EXPLICIT or IMPLICIT).
     */
    void SetSolverType(SolverType solver)
    {
        mSolverType = solver;
    }

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
     *  @param rFibreSheetDirectionsFile the file containing the fibre/sheet directions
     *  @param definedPerQuadPoint whether the fibre-sheet definitions are for each quadrature point in the mesh
     *   (if not, one for each element is assumed).
     */
    void SetVariableFibreSheetDirectionsFile(const FileFinder& rFibreSheetDirectionsFile, bool definedPerQuadPoint);


    /**
     * Set if active tension should be applied in the cross-fibre directions. This can be done in 2D or 3D.
     *
     * By default active tension is only applied in the fibre direction. This method allows the user to specify
     * that a proportion of the active tension should also be applied in the cross-fibre directions (both
     * in the sheet and sheet-normal directions equally).
     * The fraction of the active tension that is applied in the cross-fibre direction should be specified.
     *
     * @param  applyCrossFibreTension active tension is applied in the cross-fibre directions if set to true.
     * @param  crossFibreTensionFraction The fraction of the active tension to apply in the cross-fibre directions.
     */
    void SetApplyIsotropicCrossFibreTension(bool applyCrossFibreTension, double crossFibreTensionFraction);

    /**
     * Set if active tension should be applied in the cross-fibre directions. This only makes sense in 3D.
     *
     * By default active tension is only applied in the fibre direction. This method allows the user to specify
     * that a proportion of the active tension should also be applied in the cross-fibre directions.
     * The fraction of the active tension that is applied in the two cross-fibre direction should be specified.
     *
     * @param applyCrossFibreTension  active tension is applied in the cross-fibre directions if set to true.
     * @param sheetTensionFraction  The fraction of the active tension to apply in the sheet direction.
     * @param sheetNormalTensionFraction The fraction of the active tension to apply in the sheet-normal direction.
     */
    void SetApplyAnisotropicCrossFibreTension(bool applyCrossFibreTension, double sheetTensionFraction, double sheetNormalTensionFraction);

    /**
     * @return A pointer to the contraction cell factory that should be applied to the mechanics mesh.
     */
    AbstractContractionCellFactory<DIM>* GetContractionCellFactory()
    {
        assert(mpContractionCellFactory);
        return mpContractionCellFactory;
    }

    /**
     * Set the timestep to be used when solving contraction models.
     * Note that you will need to use the finest timestep that any of
     * the models requires, it is not set on a per-model basis.
     *
     * @param timestep  Timestep to use.
     */
    void SetContractionModelOdeTimestep(double timestep)
    {
        mContractionModelOdeTimeStep = timestep;
    }

    /**
     *  @return the contraction model time step.
     */
    double GetContractionModelOdeTimestep()
    {
        assert(mContractionModelOdeTimeStep>0.0); // see exception message in Validate()
        return mContractionModelOdeTimeStep;
    }

    /**
     *  @return how often the mechanics PDEs are solved.
     */
    double GetMechanicsSolveTimestep()
    {
        assert(mMechanicsSolveTimestep>0.0); // if this fails SetMechanicsSolveTimestep() probably hasn't been called - call Validate() to check
        return mMechanicsSolveTimestep;
    }

    /**
     *  @return whether the deformation affects the electrical physiological conductivity
     *  (or whether this effect is neglected).
     */
    bool GetDeformationAffectsConductivity()
    {
        return mDeformationAffectsConductivity;
    }

    /**
     *  @return whether the deformation affects the cardiac cell models, for example if
     *  there are stretch-activated channels in the cell model.
     */
    bool GetDeformationAffectsCellModels()
    {
        return mDeformationAffectsCellModels;
    }

    /**
     * @return Whether to use an implicit or explicit mechanics solver.
     */
    SolverType GetSolverType()
    {
        return mSolverType;
    }

    /**
     *  @return whether the fibre-sheet info should be read from file (if not the defaults should be used).
     */
    bool ReadFibreSheetDirectionsFromFile()
    {
        return mReadFibreSheetInformationFromFile;
    }

    /**
     *  @return the fibre-sheet file (should only be called if ReadFibreSheetDirectionsFromFile() returns true).
     */
    FileFinder GetFibreSheetDirectionsFile()
    {
        assert(mReadFibreSheetInformationFromFile);
        assert(mFibreSheetDirectionsFile.GetAbsolutePath()!="");
        return mFibreSheetDirectionsFile;
    }

    /**
     *  @return whether the fibre-sheet info is defined for each quadrature point in the mesh (if not,
     *  if it defined for each element in the mesh). (Should only be called if
     *  ReadFibreSheetDirectionsFromFile() returns true).
     */
    bool GetFibreSheetDirectionsDefinedPerQuadraturePoint()
    {
        assert(mReadFibreSheetInformationFromFile);
        return mFibreSheetDirectionsDefinedPerQuadraturePoint;
    }

    /**
     *  The first deformation (to find the equilibrium state given the loading)
     *  may require the loading to be incremented, in order for the Solve() to
     *  converge. Set the number of increments to be used.
     *  @param numIncrements number of increments
     */
    void SetNumIncrementsForInitialDeformation(unsigned numIncrements)
    {
        if (numIncrements==0)
        {
            EXCEPTION("Number of increments for initial deformation must be 1 or more");
        }
        mNumIncrementsForInitialDeformation = numIncrements;
    }

    /**
     *  @return the number of increments to be used in the initial deformation
     *  (see SetNumIncrementsForInitialDeformation() for more details).
     */
    unsigned GetNumIncrementsForInitialDeformation()
    {
        return mNumIncrementsForInitialDeformation;
    }

    /**
     * @return true if cross-fibre tension is applied.
     */
    double GetApplyCrossFibreTension()
    {
        return mApplyCrossFibreTension;
    }

    /**
     * @return the value of the cross-fibre sheet tension fraction.
     */
    double GetSheetTensionFraction()
    {
        return mSheetTensionFraction;
    }

    /**
     * @return the value of the cross-fibre sheet-normal tension fraction.
     */
    double GetSheetNormalTensionFraction()
    {
        return mSheetNormalTensionFraction;
    }


    /**
     * Check all variables are set appropriately. Exceptions are thrown if any are not.
     * Derived classes can override but should call this version as well.
     */
    virtual void Validate();
};

#endif // ELECTROMECHANICSPROBLEMDEFINITION_HPP_
