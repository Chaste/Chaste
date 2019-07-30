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

#ifndef CMGUIWRITER_HPP_
#define CMGUIWRITER_HPP_

#include "AbstractTetrahedralMeshWriter.hpp"
#include "OutputFileHandler.hpp"

/**
 * Header for node base file in 3D (.exnode)
 */
static const char CmguiNodeFileHeader3D[] = " #Fields=1\n\
 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  Value index= 1, #Derivatives= 0\n\
   y.  Value index= 2, #Derivatives= 0\n\
   z.  Value index= 3, #Derivatives= 0\n";

/**
 * Header for node base file in 2D (.exnode)
 */
static const char CmguiNodeFileHeader2D[] = " #Fields=1\n\
 1) coordinates, coordinate, rectangular cartesian, #Components=2\n\
   x.  Value index= 1, #Derivatives= 0\n\
   y.  Value index= 2, #Derivatives= 0\n";


/**
 * Header for node base file in 1D (.exnode)
 */
static const char CmguiNodeFileHeader1D[] = " #Fields=1\n\
 1) coordinates, coordinate, rectangular cartesian, #Components=1\n\
   x.  Value index= 1, #Derivatives= 0\n";

/**
 * Header for element base file in 3D (.exelem)
 */
static const char CmguiElementFileHeader3D[] = "Shape.  Dimension=3, simplex(2;3)*simplex*simplex\n\
 #Scale factor sets= 0\n\
 #Nodes= 4\n";

/**
 * Header for element base file in 3D (.exelem) for quadratic visualisation
 */
static const char CmguiElementFileHeader3DQuadratic[] = "Shape.  Dimension=3, simplex(2;3)*simplex*simplex\n\
 #Scale factor sets= 0\n\
 #Nodes= 10\n";


/**
 * Header for element base file in 2D (.exelem)
 */
static const char CmguiElementFileHeader2D[] = "Shape.  Dimension=2, simplex(2)*simplex\n\
 #Scale factor sets= 0\n\
 #Nodes= 3\n";

/**
 * Header for element base file in 2D (.exelem) for quadratic visualisation
 */
static const char CmguiElementFileHeader2DQuadratic[] = "Shape.  Dimension=2, simplex(2)*simplex\n\
 #Scale factor sets= 0\n\
 #Nodes= 6\n";

/**
 * Header for element base file in 1D (.exelem)
 */
static const char CmguiElementFileHeader1D[] = "Shape.  Dimension=1, line\n\
 #Scale factor sets= 0\n\
 #Nodes= 2\n";

/**
 * Header for element base file in 1D (.exelem)
 */
static const char CmguiElementFileHeader1DQuadratic[] = "Shape.  Dimension=1, line\n\
 #Scale factor sets= 0\n\
 #Nodes= 3\n";


/**
 * Header for element base file in 3D (.exelem), this comes after the definition of the number of fields
 */
static const char CmguiCoordinatesFileHeader3D[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
   y.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
   z.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n";


/**
 * Header for element base file in 3D (.exelem) (quadratic), this comes after the definition of the number of fields
 */
static const char CmguiCoordinatesFileHeader3DQuadratic[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=3\n\
   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.\n\
     #Nodes= 10\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n\
      7.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   7\n\
      8.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   8\n\
      9.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   9\n\
      10.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   10\n\
   y.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.\n\
     #Nodes= 10\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n\
      7.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   7\n\
      8.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   8\n\
      9.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   9\n\
      10.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   10\n\
   z.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.\n\
     #Nodes= 10\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n\
      7.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   7\n\
      8.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   8\n\
      9.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   9\n\
      10.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   10\n";


/**
 * Header for element base file in 2D (.exelem), this comes after the definition of the number of fields
 */
static const char CmguiCoordinatesFileHeader2D[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=2\n\
   x.  l.simplex(2)*l.simplex, no modify, standard node based.\n\
     #Nodes= 3\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
   y.  l.simplex(2)*l.simplex, no modify, standard node based.\n\
     #Nodes= 3\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n";


/**
 * Header for element base file in 2D (.exelem) (quadratic version), this comes after
 * the definition of the number of fields
 */
static const char CmguiCoordinatesFileHeader2DQuadratic[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=2\n\
   x.  q.simplex(2)*q.simplex, no modify, standard node based.\n\
     #Nodes= 6\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n\
   y.  q.simplex(2)*q.simplex, no modify, standard node based.\n\
     #Nodes= 6\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n";

/**
 * Header for element base file in 1D (.exelem), this comes after the definition of the number of fields
 * Note that in 1D the simplex doesn't seem to work, we use Lagrange instead
 */
static const char CmguiCoordinatesFileHeader1D[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=1\n\
   x.  l.Lagrange, no modify, standard node based.\n\
     #Nodes= 2\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n";


/**
 * Header for element base file in 1D (.exelem) (quadratic visualisation), this comes after the definition of the number of fields
 * Note that in 1D the simplex doesn't seem to work, we use Lagrange instead
 */
static const char CmguiCoordinatesFileHeader1DQuadratic[] = " 1) coordinates, coordinate, rectangular cartesian, #Components=1\n\
   x.  q.Lagrange, no modify, standard node based.\n\
     #Nodes= 3\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n";

/**
 * Header for additional fields in the element base file in 3D (.exelem),
 * Here we assume all additional fields will be interpolated by cmgui in the same way
 */
static const char CmguiAdditionalFieldHeader3D[] = " field, rectangular cartesian, #Components=1\n\
   x.  l.simplex(2;3)*l.simplex*l.simplex, no modify, standard node based.\n\
     #Nodes= 4\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n";


/**
 * Header for additional fields in the element base file in 3D (.exelem) (quadratic
 * visualisation). Here we assume all additional fields will be interpolated by cmgui
 * in the same way
 */
static const char CmguiAdditionalFieldHeader3DQuadratic[] = " field, rectangular cartesian, #Components=1\n\
   x.  q.simplex(2;3)*q.simplex*q.simplex, no modify, standard node based.\n\
     #Nodes= 10\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n\
      7.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   7\n\
      8.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   8\n\
      9.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   9\n\
      10.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   10\n";

/**
 * Header for additional fields in the element base file in 2D (.exelem),
 * Here we assume all additional fields will be interpolated by cmgui in the same way
 */
static const char CmguiAdditionalFieldHeader2D[] = " field, rectangular cartesian, #Components=1\n\
   x.  l.simplex(2)*l.simplex, no modify, standard node based.\n\
     #Nodes= 3\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n";

/**
 * Header for additional fields in the element base file in 2D (.exelem)  (quadratic version),
 * Here we assume all additional fields will be interpolated by cmgui in the same way
 */
static const char CmguiAdditionalFieldHeader2DQuadratic[] = " field, rectangular cartesian, #Components=1\n\
   x.  q.simplex(2)*q.simplex, no modify, standard node based.\n\
     #Nodes= 6\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n\
      4.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   4\n\
      5.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   5\n\
      6.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   6\n";


/**
 * Header for additional fields in the element base file in 1D (.exelem),
 * Here we assume all additional fields will be interpolated by cmgui in the same way
 */
static const char CmguiAdditionalFieldHeader1D[] = " field, rectangular cartesian, #Components=1\n\
   x.  l.Lagrange, no modify, standard node based.\n\
     #Nodes= 2\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n";

/**
 * Header for additional fields in the element base file in 1D (.exelem) (quadratic visualisation),
 * Here we assume all additional fields will be interpolated by cmgui in the same way
 */
static const char CmguiAdditionalFieldHeader1DQuadratic[] = " field, rectangular cartesian, #Components=1\n\
   x.  q.Lagrange, no modify, standard node based.\n\
     #Nodes= 2\n\
      1.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   1\n\
      2.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   2\n\
      3.  #Values=1\n\
       Value indices:     1\n\
       Scale factor indices:   3\n";


/**
 *  CmguiMeshWriter
 *
 *  Writes a mesh in Cmgui (the visualisation frontend of CMISS) format. Creates an exnode
 *  file and a exelem file. Note that the lines and faces are not written in the exelem
 *  file, so to load the data in Cmgui, you must use 'generate_faces_and_lines', i.e.
 *
 *  gfx read node base_file
 *  gfx read elem base_file generate_faces_and_lines
 *  gfx cr win
 */
template<unsigned ELEMENT_DIM, unsigned SPACE_DIM>
class CmguiMeshWriter : public AbstractTetrahedralMeshWriter<ELEMENT_DIM,SPACE_DIM>
{
protected:

    /**
     * For storage of names of additional fields.
     */
    std::vector<std::string> mAdditionalFieldNames;


    /**
     * Stored the names of different regions that we want Cmgui to be able to visualize separately.
     * There will be one .exelem file per region.
     */
    std::vector<std::string> mRegionNames;// {tissue, bath, whatever}

    /**
     * The group name to give in the output files. Defaults to the same as the
     * base name. The CmguiDeformedSolutionsWriter prepends a counter to the base name
     * (eg "solution_8.exnodes" instead of just "solution.exnodes"), but we would like
     * the group name to stay as "solution", hence this separate variable
     */
    std::string mGroupName;

    /**
     *  String which is set to either CmguiElementFileHeader2D or CmguiElementFileHeader2DQuadratic
     *  as appropriate
     */
    std::string mElementFileHeader;

    /**
     *  String which is set to either CmguiCoordinatesFileHeader2D or
     *  CmguiCoordinatesFileHeader2DQuadratic as appropriate
     */

    std::string mCoordinatesFileHeader;

    /**
     *  String which is set to either CmguiAdditionalFieldHeader2D or
     *  CmguiAdditionalFieldHeader2DQuadratic as appropriate
     */
    std::string mAdditionalFieldHeader;

    /**
     *  Number of nodes per element, eg, in 2D, 3 for linear visualisation
     *  and 6 for quadratic visualisation
     */
    unsigned mNumNodesPerElement;

    /**
     *  Ordering of the elements nodes, from Chaste convention to CMGUI convention
     */
    std::vector<unsigned> mReordering;

    /**
     * @return the mode to use when opening files.
     *
     * @param append  whether to append to the file, or overwrite it
     */
    std::ios_base::openmode GetOpenMode(bool append);

    /**
     * Open the file node information is written to.
     *
     * @return file handler
     * @param append whether to clear (append=false) or to append (append=true). False by default.
     */
    out_stream OpenNodeFile(bool append = false);

    /**
     * Helper method that open all the elements files (one per region in this case)
     * @return vector of file handlers
     *
     * @param append whether to clear (append=false) or to append (append=true). False by default. It applies to all files.
     */
    std::vector<boost::shared_ptr<std::ofstream> > OpenElementFiles(bool append = false);

    /**
     *  Write the header part of a node file, depending on the dimension. Short helper method,
     *  also called in CmguiDeformedSolutionsWriter. (Note, without the & below this method
     *  seg faults).
     *  @param rpNodeFile reference to the out_stream used for the node file
     */
    void WriteNodeFileHeader(out_stream& rpNodeFile);

    /**
     * Write the headers of each element file (as many as number of regions defined).
     *
     * @param rElemFiles vector of pointers to file streams
     */
    void WriteElementsFileHeader(std::vector<boost::shared_ptr<std::ofstream> >& rElemFiles);


    /**
     * Create output files and add headers.
     */
    void CreateFilesWithHeaders();

    /**
     * Append local mesh data to output files.
     */
    void AppendLocalDataToFiles();

    /**
     * Append footers to output files.
     */
    void WriteFilesFooter();

public:

    /**
     * Constructor.
     *
     * @param rDirectory  the directory in which to write the mesh to file
     * @param rBaseName  the base name of the files in which to write the mesh data
     * @param cleanDirectory  whether to clean the directory (defaults to true)
     */
    CmguiMeshWriter(const std::string& rDirectory,
                    const std::string& rBaseName,
                    bool cleanDirectory=true);

    /**
     * Write mesh data to files.
     */
    void WriteFiles();

    /**
     * Set any additional field that we want cmgui to visualize (interpolated over) elements and surfaces
     *
     * @param rFieldNames is a reference to a vector of string containing the names of each additional field name
     */
    void SetAdditionalFieldNames(std::vector<std::string>& rFieldNames);

    /**
     * Set the region names to be used when generating multiple element files
     *
     * @param rRegionNames is a reference to a vector of string containing the names of each region defined in the mesh
     */
    void SetRegionNames(std::vector<std::string>& rRegionNames);

    /**
     * Destructor.
     */
    virtual ~CmguiMeshWriter()
    {}

    // A method called CompareCmguiFiles() has been removed, please use FileComparison class instead.
};

#endif /*CMGUIWRITER_HPP_*/
