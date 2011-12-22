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




#ifndef VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_
#define VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_


#include <vector>
#include <string>
#include "UblasIncludes.hpp"
#include "TetrahedralMesh.hpp"
#include "QuadraticMesh.hpp"
#include "FineCoarseMeshPair.hpp"
#include "ReplicatableVector.hpp"
#include "HeartConfig.hpp"
#include "Hdf5ToCmguiConverter.hpp"
#include "Hdf5DataReader.hpp"
#include "PetscTools.hpp"
#include "Hdf5DataWriter.hpp"


/**
 *  Very simple one-method class which can be used to convert the voltage from an electrics
 *  (or electromechanics) simulation onto a coarser mechanics mesh, by interpolation. The
 *  class outputs a HDF5 file corresponding to nodes on the mechanics mesh, and converts it to
 *  CMGUI output.
 */
template<unsigned DIM>
class VoltageInterpolaterOntoMechanicsMesh
{
public:
    /**
     *  Constructor, also the main method of the class
     *
     *  @param rElectricsMesh The electrics mesh
     *  @param rMechanicsMesh The mechanics mesh
     *  @param directory Directory the voltage file is in
     *  @param inputFileNamePrefix Filename (without ".h5") of the electrics solution HDF5 file
     */
    VoltageInterpolaterOntoMechanicsMesh(TetrahedralMesh<DIM,DIM>& rElectricsMesh,
                                         QuadraticMesh<DIM>& rMechanicsMesh,
                                         std::string directory,
                                         std::string inputFileNamePrefix);
};

#endif /*VOLTAGEINTERPOLATERONTOMECHANICSMESH_HPP_*/
