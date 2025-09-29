/*

Copyright (c) 2005-2025, University of Oxford.
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


#ifndef CARDIACELECTROMECHANICSVTKHANDLER_HPP_
#define CARDIACELECTROMECHANICSVTKHANDLER_HPP_

#ifdef CHASTE_VTK

#include "AbstractCardiacProblem.hpp"
#include "AbstractNonlinearElasticitySolver.hpp"
#include "FineCoarseMeshPair.hpp"
#include "ReplicatableVector.hpp"
#include "ElectroMechanicsProblemDefinition.hpp"
#include "VtkDeformedMeshWriter.hpp"
#include "VoltageInterpolaterOntoMechanicsMesh.hpp"
#include "VtkNonlinearElasticitySolutionWriter.hpp"

/**
 * A convenience class to handle the VTK output of cardiac
 * electromechanics simulations. It makes use of various
 * other classes to compute mechanics qunatities, 
 * interpolate electrics quantities, and actually generate
 * files. It outputs voltage, displacements, and deformation gradient in VTU
 * format to be shown in the deformed mechanics mesh.
 */
template<unsigned DIM, unsigned ELEC_PROB_DIM=1>
class CardiacElectroMechanicsVtkHandler
{

private: 
    /** Cache for the mechanics solver (deformed solution is taken from this object)*/
    AbstractNonlinearElasticitySolver<DIM>& mrMechanicsSolver;
    /** Cache for the electrics problem used in the EM problem*/
    AbstractCardiacProblem<DIM,DIM,ELEC_PROB_DIM>& mrElectricsProblem;
       
    /**
     * A cache for the interpolated voltages from electrics to mechanics mesh but node-wise.
     * Memory is allocated within the constructor. Filled in when WriteSolution() is called.
     */
    std::vector<double> mInterpolatedVoltagesNodeWise;
    
    /** Vector of displacements to be written to file */
    std::vector<c_vector<double,DIM> > mDisplacements;

    /** vector to store strains to be printed out. */
    std::vector<c_matrix<double,DIM,DIM> > mStrains;

    /** Pointer to the mechanics solution writer class. Used to calculate some mechanics quantities */
    VtkNonlinearElasticitySolutionWriter<DIM>* mpVtkElastictyWriter;

    /** Poiunter to the actual mesh VTK writer. Initialized upon construction */
    VtkDeformedMeshWriter<DIM>* mpVtkWriter;

    /** Used to interpolate electrics solution onto mechanics mesh for VTK output*/
    VoltageInterpolaterOntoMechanicsMesh<DIM>* mpInterpolater;

    /** 
     * Local cache of the mechanics mesh. 
     * Object created within the constructor as a local copy of the
     * mesh onto which apply deformation by changing node coordinates.
     */
    QuadraticMesh<DIM>* mpVtkOutputMesh; 

public: 

    /**
     * Constructor. It creates all the necessary objects based on the input
     * @param rMechanicsSolver the mechanics solver used by the EM problem. The deformation
              that this class uses to write out a VRK deformed mesh is taken from this solver.
     * @param rQuadMesh The mesh used by the mechanics problem. A copy is made here.
     * @param rElectricsMesh The mesh used by the electrcis problem. 
     * @param rElectricsProblem The actual electrics problem (used for initial conditions)
     * @param rOutputDir The output directory where to write. 
     *        Files will be written (relative to CHASTE_TEST_OUTPUT) to rOutputDir/vtk/
     *        It will be cleared upon construction.
     */
    CardiacElectroMechanicsVtkHandler(AbstractNonlinearElasticitySolver<DIM>& rMechanicsSolver,
                                      QuadraticMesh<DIM>& rQuadMesh,
                                      TetrahedralMesh<DIM,DIM>& rElectricsMesh,
                                      AbstractCardiacProblem<DIM,DIM,ELEC_PROB_DIM>& rElectricsProblem,
                                      const std::string& rOutputDir);
    /**
     * Destructor
     */
     ~CardiacElectroMechanicsVtkHandler();
     
     /**
      * Write the solution to file. It will write to rOutputDir/vtk/deformed_mechanics_mesh_X.vtu where X=counter.
      * it writes Voltage ("V"), displacements ("displacements"),
      * and deformation gradient F ("deformation_gradient_F") using the deformed mechanics mesh.
      *
      * @param counter Used to name the file 
      * @param rElectricsSolution The solution of the electrcis problem. It will be interpolated
      *        onto the mechanics mesh within this method to outoput voltage V. 
      */
     void WriteSolution(unsigned counter, ReplicatableVector& rElectricsSolution);
};

#endif //CHATE_VTK

#endif //CARDIACELECTROMECHANICSVTKHANDLER_HPP_
