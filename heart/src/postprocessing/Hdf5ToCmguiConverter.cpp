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

#include "Hdf5ToCmguiConverter.hpp"
#include "CmguiMeshWriter.hpp"
#include "UblasCustomFunctions.hpp"
#include "HeartConfig.hpp"
#include "PetscTools.hpp"
#include "Exception.hpp"
#include "ReplicatableVector.hpp"
#include "DistributedVector.hpp"
#include "DistributedVectorFactory.hpp"
#include "Version.hpp"
#include "GenericMeshReader.hpp"

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Hdf5ToCmguiConverter<ELEMENT_DIM,SPACE_DIM>::Write(std::string type)
{
    out_stream p_file = out_stream(NULL);

    unsigned num_nodes = this->mpReader->GetNumberOfRows();
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();

    DistributedVectorFactory factory(num_nodes);

    Vec data = factory.CreateVec();//for V
    Vec data_phie = factory.CreateVec();//for phi_e
    Vec data_second_cell = factory.CreateVec();//for the V of the second cell, used in extended bidomain problems.

    for (unsigned time_step=0; time_step<num_timesteps; time_step++)
    {
        // Create the file for this time step
        std::stringstream time_step_string;

        // unsigned to string
        time_step_string << time_step;
        if (PetscTools::AmMaster())
        {
            p_file = this->mpOutputFileHandler->OpenOutputFile(this->mFileBaseName + "_" + time_step_string.str() + ".exnode");

            // Check how many digits are to be output in the solution (0 goes to default value of digits)
            unsigned int num_digits = HeartConfig::Instance()->GetVisualizerOutputPrecision();
            if (num_digits != 0)
            {
               p_file->precision(num_digits);
            }
        }

        std::vector<ReplicatableVector*> all_data;
        unsigned num_vars = this->mpReader->GetVariableNames().size();
        for (unsigned var=0; var<num_vars; var++)
        {
            // Read the data for this time step
            this->mpReader->GetVariableOverNodes(data, this->mpReader->GetVariableNames()[var], time_step);
            ReplicatableVector* p_repl_data = new ReplicatableVector(data);
            assert(p_repl_data->GetSize()==num_nodes);
            all_data.push_back(p_repl_data);
        }

        if (PetscTools::AmMaster())
        {
            // Write provenance info
            std::string comment = "! " + ChasteBuildInfo::GetProvenanceString();
            *p_file << comment;
            // The header first
            *p_file << "Group name: " << this->mFileBaseName << "\n";
            *p_file << "#Fields=" << num_vars << "\n";
            for (unsigned var=0; var<num_vars; var++)
            {
                *p_file << " " << var+1 << ") " <<this->mpReader->GetVariableNames()[var]<< " , field, rectangular cartesian, #Components=1" << "\n" << "x.  Value index=1, #Derivatives=0, #Versions=1"<<"\n";
                if (var != num_vars-1)
                {
                    *p_file << "\n";
                }
            }

            // Write the data
            for (unsigned i=0; i<num_nodes; i++)
            {
                // cmgui counts nodes from 1
                *p_file << "Node: "<< i+1 << "\n";
                for (unsigned var=0; var<num_vars; var++)
                {
                    *p_file  << (*(all_data[var]))[i] << "\n";
                }
            }
        }

        for (unsigned var=0; var<num_vars; var++)
        {
           delete all_data[var];
        }
    }
    PetscTools::Destroy(data);
    PetscTools::Destroy(data_phie);
    PetscTools::Destroy(data_second_cell);

    if (PetscTools::AmMaster())
    {
        p_file->close();
    }
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
Hdf5ToCmguiConverter<ELEMENT_DIM,SPACE_DIM>::Hdf5ToCmguiConverter(std::string inputDirectory,
                                                                  std::string fileBaseName,
                                                                  AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                                                  bool hasBath)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(inputDirectory, fileBaseName, pMesh, "cmgui_output")
{
    // Write the node data out
    Write("");

    // Write mesh in a suitable form for cmgui
    std::string output_directory =  HeartConfig::Instance()->GetOutputDirectory() + "/cmgui_output";

    CmguiMeshWriter<ELEMENT_DIM,SPACE_DIM> cmgui_mesh_writer(output_directory, HeartConfig::Instance()->GetOutputFilenamePrefix(), false);

    // Used to inform the mesh of the data names
    std::vector<std::string> field_names = this->mpReader->GetVariableNames();
    cmgui_mesh_writer.SetAdditionalFieldNames(field_names);
    if (hasBath)
    {
        std::vector<std::string> names;
        names.push_back("tissue");
        names.push_back("bath");
        cmgui_mesh_writer.SetRegionNames(names);
    }

    // Normally the in-memory mesh is converted:
    if (HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering() == false)
    {
        cmgui_mesh_writer.WriteFilesUsingMesh(*(this->mpMesh), false);
    }
    else
    {
        // In this case we expect the mesh to have been read in from file
        ///\todo What if the mesh has been scaled, translated or rotated?
        // Note that the next line will throw if the mesh has not been read from file
        std::string original_file=this->mpMesh->GetMeshFileBaseName();
        std::auto_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_original_mesh_reader
            = GenericMeshReader<ELEMENT_DIM, SPACE_DIM>(original_file);
        cmgui_mesh_writer.WriteFilesUsingMeshReader(*p_original_mesh_reader);
    }

    WriteCmguiScript();
    PetscTools::Barrier("Hdf5ToCmguiConverter");
}

template <unsigned ELEMENT_DIM, unsigned SPACE_DIM>
void Hdf5ToCmguiConverter<ELEMENT_DIM,SPACE_DIM>::WriteCmguiScript()
{
    unsigned num_timesteps = this->mpReader->GetUnlimitedDimensionValues().size();
    assert(this->mpReader->GetVariableNames().size() > 0); // seg fault guard
    std::string variable_name = this->mpReader->GetVariableNames()[0];

    if (PetscTools::AmMaster())
    {
        out_stream p_script_file = this->mpOutputFileHandler->OpenOutputFile("script.com");

        // Write provenance info, note the # instead of ! because this is - essentially - a PERL script that Cmgui interprets
        std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();
        *p_script_file << comment;

        *p_script_file << "# Read the mesh \n"
                       << "gfx read node "<<HeartConfig::Instance()->GetOutputFilenamePrefix()<<".exnode \n"
                       << "gfx read elem "<<HeartConfig::Instance()->GetOutputFilenamePrefix()<<".exelem generate_faces_and_lines \n" // note the mesh file name is taken from HeartConfig...
                       << "# Create a window \n"
                       << "gfx cre win 1 \n"
                       << "# Modify the scene (obtained by gfx list g_element XXXX commands) to visualize first var on lines and nodes \n"
                       << "gfx modify g_element "<< HeartConfig::Instance()->GetOutputFilenamePrefix()<<" general clear circle_discretization 6 default_coordinate coordinates element_discretization \"4*4*4\" native_discretization none; \n"
                       << "gfx modify g_element "<< HeartConfig::Instance()->GetOutputFilenamePrefix()<<" lines select_on material default data "<<variable_name<<" spectrum default selected_material default_selected; \n"
                       << "gfx modify g_element "<< HeartConfig::Instance()->GetOutputFilenamePrefix()<<" node_points glyph point general size \"1*1*1\" centre 0,0,0 font default select_on material default data "<<variable_name<<" spectrum default selected_material default_selected; \n"
                       << "# Load the data \n"
                       << "for ($i=0; $i<" << num_timesteps << "; $i++) { \n"
                       << "    gfx read node " << this->mFileBaseName << "_$i.exnode time $i\n" // ...while the data file from mFileBaseName...
                       << "}\n";
        p_script_file->close();
    }
}

/////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////

template class Hdf5ToCmguiConverter<1,1>;
template class Hdf5ToCmguiConverter<1,2>;
template class Hdf5ToCmguiConverter<2,2>;
template class Hdf5ToCmguiConverter<1,3>;
template class Hdf5ToCmguiConverter<2,3>;
template class Hdf5ToCmguiConverter<3,3>;
