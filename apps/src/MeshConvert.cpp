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

/*
 * Note: Do not put any VTK-specific functionality in this file, as we
 * don't ever test it with VTK support turned off!
 */

// Most of the work is done by this class.  It must be included first.
//#include "CardiacSimulation.hpp"

#include <string>
#include <libgen.h>

#include "ExecutableSupport.hpp"
#include "Exception.hpp"
#include "GenericMeshReader.hpp"
#include "DistributedTetrahedralMesh.hpp"
#include "TrianglesMeshWriter.hpp"

int main(int argc, char *argv[])
{
    ExecutableSupport::StandardStartup(&argc, &argv);

    int exit_code = ExecutableSupport::EXIT_OK;

    try
    {
        if (argc<2)
        {
            ExecutableSupport::PrintError("Usage: MeshConvert mesh_3d_file_base_name", true);
            exit_code = ExecutableSupport::EXIT_BAD_ARGUMENTS;
        }
        else
        {
            std::string base=basename(argv[1]);
#ifdef CHASTE_VTK
            ExecutableSupport::Print("Note: for VTK reading, give the full file path (including '.vtu' extension)");
#else
            ExecutableSupport::Print("Note: VTK reading is not supported");
#endif

            ExecutableSupport::Print("Opening "+base+" mesh file(s).");

            std::auto_ptr<AbstractMeshReader<3,3> > p_mesh_reader = GenericMeshReader<3,3>(argv[1]);
            //We have to make a mesh so that we can get the node connectivity list back
            DistributedTetrahedralMesh<3,3> mesh;
            mesh.ConstructFromMeshReader(*p_mesh_reader);

            //Find a dot
            std::string base_for_output=base;
            size_t pos = base.find('.');
            if (pos != std::string::npos)
            {
                //If dot found, then make the string smaller
                base_for_output.resize(pos);
            }
            base_for_output=base_for_output+"_bin";
            TrianglesMeshWriter<3,3> mesh_writer("", base_for_output);
            ExecutableSupport::Print("Writing  "+base_for_output+".node etc. mesh file in "+mesh_writer.GetOutputDirectory());
            mesh_writer.SetWriteFilesAsBinary();
            mesh_writer.WriteFilesUsingMesh(mesh);
            ExecutableSupport::Print("Done.");
        }
    }
    catch (const Exception& e)
    {
        ExecutableSupport::PrintError(e.GetMessage());
        exit_code = ExecutableSupport::EXIT_ERROR;
    }

    ExecutableSupport::FinalizePetsc();
    return exit_code;
}
