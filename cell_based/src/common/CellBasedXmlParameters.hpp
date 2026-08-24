/*

Copyright (c) 2005-2026, University of Oxford.
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

#ifndef CELLBASEDXMLPARAMETERS_HPP_
#define CELLBASEDXMLPARAMETERS_HPP_

#include <string>

/**
 * @file CellBasedXmlParameters.hpp
 *
 * Macros for writing XML parameter output in Output*Parameters() methods.
 *
 * All member variables that appear in XML output must follow the 'm' prefix
 * convention (e.g. mMyParam).  The macros below derive the XML tag name from
 * the variable name by stripping the leading 'm', so
 *
 *     CHASTE_PARAM(rParamsFile, level, mMyParam)
 *
 * emits  \<MyParam\>value\</MyParam\>  at the requested tab-indentation level.
 *
 * The indentation level is the number of tab characters to prepend:
 *  - level 2  simulation / population parameters
 *  - level 3  cell-cycle model / force / modifier / update-rule / numerical-method /
 *             boundary-condition parameters
 *
 * Declare the level once per Output*Parameters() function:
 * @code
 *     const unsigned level = 3;
 *     CHASTE_PARAM(rParamsFile, level, mFoo);
 *     CHASTE_PARAM(rParamsFile, level, mBar);
 * @endcode
 */

/**
 * Write a scalar member variable as an XML element.
 *
 * Strips the leading 'm' from the variable name to form the XML tag, e.g.
 * CHASTE_PARAM(rParamsFile, 3, mMyParam) emits
 *     `\t\t\t`\<MyParam\>value\</MyParam\>
 *
 * @param stream   the out_stream to write to (e.g. rParamsFile)
 * @param level    number of leading tabs (unsigned integer expression)
 * @param mVar     the member variable (must start with 'm')
 */
#define CHASTE_PARAM(stream, level, mVar) \
    *(stream) << std::string((level), '\t') \
              << "<" << (&(#mVar)[1]) << ">" << (mVar) << "</" << (&(#mVar)[1]) << ">\n"

/**
 * Write a c_vector member variable as a comma-separated XML element.
 *
 * Strips the leading 'm' from the variable name to form the XML tag, e.g.
 * CHASTE_PARAM_CVECTOR(rParamsFile, 3, mMyVec) emits
 *     `\t\t\t`\<MyVec\>v0,v1,...,vN\</MyVec\>
 *
 * @param stream   the out_stream to write to
 * @param level    number of leading tabs (unsigned integer expression)
 * @param mVar     the c_vector member variable (must start with 'm')
 */
#define CHASTE_PARAM_CVECTOR(stream, level, mVar) \
    do { \
        *(stream) << std::string((level), '\t') << "<" << (&(#mVar)[1]) << ">"; \
        for (unsigned chaste_cvec_i_ = 0; chaste_cvec_i_ + 1u < (unsigned)(mVar).size(); ++chaste_cvec_i_) \
        { \
            *(stream) << (mVar)[chaste_cvec_i_] << ","; \
        } \
        *(stream) << (mVar)[(mVar).size() - 1u] << "</" << (&(#mVar)[1]) << ">\n"; \
    } while (false)

/**
 * Write an arbitrary expression as an XML element with an explicit tag name.
 *
 * Use this for parameters that are not stored as a direct member variable,
 * e.g. values obtained from a helper object or singleton.
 *
 * CHASTE_PARAM_EXPR(rParamsFile, 3, MyTag, SomeObject::Instance()->GetValue())
 * emits   `\t\t\t`\<MyTag\>value\</MyTag\>
 *
 * @param stream   the out_stream to write to
 * @param level    number of leading tabs (unsigned integer expression)
 * @param tag      the XML tag name (written as an identifier, not a string)
 * @param expr     the expression whose value is written
 */
#define CHASTE_PARAM_EXPR(stream, level, tag, expr) \
    *(stream) << std::string((level), '\t') \
              << "<" << #tag << ">" << (expr) << "</" << #tag << ">\n"

#endif // CELLBASEDXMLPARAMETERS_HPP_
