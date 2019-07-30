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
            if (this->mPrecision != 0)
            {
               p_file->precision(this->mPrecision);
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
Hdf5ToCmguiConverter<ELEMENT_DIM,SPACE_DIM>::Hdf5ToCmguiConverter(const FileFinder& rInputDirectory,
                                                                  const std::string& rFileBaseName,
                                                                  AbstractTetrahedralMesh<ELEMENT_DIM,SPACE_DIM>* pMesh,
                                                                  bool hasBath,
                                                                  unsigned precision)
    : AbstractHdf5Converter<ELEMENT_DIM,SPACE_DIM>(rInputDirectory, rFileBaseName, pMesh, "cmgui_output", precision)
{
    ///\todo #1660 at present this converter is hardcoded to work with "Data" using the below statement
    while (this->mDatasetNames[this->mOpenDatasetIndex] != "Data")
    {
        bool next_open = this->MoveOntoNextDataset();
        UNUSED_OPT(next_open);
        assert(next_open);
    }

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
    if (HeartConfig::Instance()->GetOutputUsingOriginalNodeOrdering() == false || !this->mpMesh->IsMeshOnDisk())
    {
        cmgui_mesh_writer.WriteFilesUsingMesh(*(this->mpMesh), false);
    }
    else
    {
        // In this case we expect the mesh to have been read in from file
        ///\todo What if the mesh has been scaled, translated or rotated?
        // Note that the next line will throw if the mesh has not been read from file
        std::string original_file=this->mpMesh->GetMeshFileBaseName();
        std::shared_ptr<AbstractMeshReader<ELEMENT_DIM, SPACE_DIM> > p_original_mesh_reader
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
        out_stream p_script_file = this->mpOutputFileHandler->OpenOutputFile("LoadSolutions.com");

        // Write provenance info, note the # instead of ! because this is - essentially - a PERL script that Cmgui interprets
        std::string comment = "# " + ChasteBuildInfo::GetProvenanceString();
        *p_script_file << comment;

        *p_script_file << "# Read the mesh \n"
                       << "gfx read node "<<HeartConfig::Instance()->GetOutputFilenamePrefix()<<".exnode \n"
                       << "gfx read elem "<<HeartConfig::Instance()->GetOutputFilenamePrefix()<<".exelem \n" // note the mesh file name is taken from HeartConfig...
                       << "gfx define faces egroup "<<HeartConfig::Instance()->GetOutputFilenamePrefix()<<"\n"
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

// Explicit instantiation
template class Hdf5ToCmguiConverter<1,1>;
template class Hdf5ToCmguiConverter<1,2>;
template class Hdf5ToCmguiConverter<2,2>;
template class Hdf5ToCmguiConverter<1,3>;
template class Hdf5ToCmguiConverter<2,3>;
template class Hdf5ToCmguiConverter<3,3>;
