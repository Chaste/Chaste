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


#ifndef ABSTRACTCARDIACMECHANICSSOLVERINTERFACE_HPP_
#define ABSTRACTCARDIACMECHANICSSOLVERINTERFACE_HPP_

#include "AbstractNonlinearElasticitySolver.hpp"


/**
 * This class declares all the main public methods in AbstractCardiacMechanicsSolver.
 *
 * The methods declared here are all pure virtual, and implemented in AbstractCardiacMechanicsSolver
 * or it's children. The point of this class is that it is nothing templated over elasticity solver
 * (unlike AbstractCardiacMechanicsSolver), so is more convenient to use, in particular the cardiac
 * electro-mechanics problem class can own a pointer to one of these (ie, a
 * AbstractCardiacMechanicsSolverInterface<DIM>*), and then the problem does not also require
 * templating over solver.
 *
 */
template<unsigned DIM>
class AbstractCardiacMechanicsSolverInterface
{
public:
    /** Constructor, does nothing */
    AbstractCardiacMechanicsSolverInterface()
    {
    }

    /** Destructor, does nothing */
    virtual ~AbstractCardiacMechanicsSolverInterface()
    {
    }

    /** Get the total number of quad points in the mesh. Pure, implemented in concrete solver */
    virtual unsigned GetTotalNumQuadPoints()=0;

    /** Get the quadrature rule used in the elements. */
    virtual GaussianQuadratureRule<DIM>* GetQuadratureRule()=0;


    /**
     *  Set a constant fibre-sheet-normal direction (a matrix) to something other than the default (fibres in X-direction,
     *  sheet in the XY plane)
     *  @param rFibreSheetMatrix The fibre-sheet-normal matrix (fibre dir the first column, normal-to-fibre-in sheet in second
     *  column, sheet-normal in third column).
     */
    virtual void SetConstantFibreSheetDirections(const c_matrix<double,DIM,DIM>& rFibreSheetMatrix)=0;

    /**
     *  Set a variable fibre-sheet-normal direction (matrices), from file.
     *  If the second parameter is false, there should be one fibre-sheet definition for each element; otherwise
     *  there should be one fibre-sheet definition for each *quadrature point* in the mesh.
     *  In the first case, the file should be a .ortho file (ie each line has the fibre dir, sheet dir, normal dir
     *  for that element), in the second it should have .orthoquad as the format.
     *
     *  @param orthoFile the file containing the fibre/sheet directions
     *  @param definedPerQuadraturePoint whether the fibre-sheet definitions are for each quadrature point in the mesh
     *   (if not, one for each element is assumed).
     */
    virtual void SetVariableFibreSheetDirections(std::string orthoFile, bool definedPerQuadraturePoint)=0;



    /**
     *  Set the intracellular Calcium concentrations and voltages at each quad point. Pure.
     *
     *  Implicit solvers (for contraction models which are functions of stretch (and maybe
     *  stretch rate) would integrate the contraction model with this Ca/V/t using the current
     *  stretch (ie inside AssembleOnElement, ie inside GetActiveTensionAndTensionDerivs).
     *  Explicit solvers (for contraction models which are NOT functions of stretch can immediately
     *  integrate the contraction models to get the active tension.
     *
     *  @param rCalciumConcentrations Reference to a vector of intracellular calcium concentrations at each quadrature point
     *  @param rVoltages Reference to a vector of voltages at each quadrature point
     */
    virtual void SetCalciumAndVoltage(std::vector<double>& rCalciumConcentrations,
                                      std::vector<double>& rVoltages)=0;

    /**
     *  Solve for the deformation, integrating the contraction model ODEs.
     *
     *  @param time the current time
     *  @param nextTime the next time
     *  @param odeTimestep the ODE timestep
     */
    virtual void Solve(double time, double nextTime, double odeTimestep)=0;



    /**
     *  Compute the deformation gradient, and stretch in the fibre direction, for each element in the mesh.
     *  Note: using quadratic interpolation for position, the deformation gradients and stretches
     *  actually vary linearly in each element. However, for computational efficiency reasons, when computing
     *  deformation gradients and stretches to pass back to the electrophysiology solver, we just assume
     *  they are constant in each element (ie ignoring the quadratic correction to the displacement). This means
     *  that  the (const) deformation gradient and stretch for each element can be computed in advance and
     *  stored, and we don't have to worry about interpolation onto the precise location of the cell-model (electrics-mesh)
     *  node, just which element it is in, and ditto the electric mesh element centroid.
     *
     *  To compute this (elementwise-)constant F (and from it the constant stretch), we just have to compute
     *  F using the deformed positions at the vertices only, with linear bases, rather than all the
     *  nodes and quadratic bases.
     *
     *  @param rDeformationGradients A reference of a std::vector in which the deformation gradient
     *  in each element will be returned. Must be allocated prior to being passed in.
     *  @param rStretches A reference of a std::vector in which the stretch in each element will be returned.
     *  Must be allocated prior to being passed in.
     */
    virtual void ComputeDeformationGradientAndStretchInEachElement(std::vector<c_matrix<double,DIM,DIM> >& rDeformationGradients,
                                                                   std::vector<double>& rStretches)=0;
};

#endif /* ABSTRACTCARDIACMECHANICSSOLVERINTERFACE_HPP_ */
