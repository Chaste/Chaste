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
#include "FileFinder.hpp"
#include "FibreConverter.hpp"

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

            std::shared_ptr<AbstractMeshReader<3,3> > p_mesh_reader = GenericMeshReader<3,3>(argv[1]);
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
            base_for_output = base_for_output + "_bin";
            TrianglesMeshWriter<3,3> mesh_writer("", base_for_output);
            ExecutableSupport::Print("Writing  " + base_for_output + ".node etc. mesh file in " + mesh_writer.GetOutputDirectory());
            mesh_writer.SetWriteFilesAsBinary();
            mesh_writer.WriteFilesUsingMesh(mesh);
            // Convert fibres if present
            FibreConverter fibre_converter;
            FileFinder mesh_file(argv[1], RelativeTo::AbsoluteOrCwd);
            fibre_converter.Convert(mesh_file, "");
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
