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

#include "HeartFileFinder.hpp"
#include "Exception.hpp"


HeartFileFinder::HeartFileFinder(const cp::path_type& rPath)
{
    std::string relative_path(rPath);
    if (rPath.relative_to() == cp::relative_to_type::this_file)
    {
        SetPath(rPath, HeartConfig::Instance()->GetParametersFilePath());
    }
    else
    {
        RelativeTo::Value relative_to;
        switch (rPath.relative_to())
        {
            case cp::relative_to_type::cwd:
                relative_to = RelativeTo::CWD;
                break;
            case cp::relative_to_type::chaste_test_output:
                relative_to = RelativeTo::ChasteTestOutput;
                break;
            case cp::relative_to_type::chaste_source_root:
                relative_to = RelativeTo::ChasteSourceRoot;
                break;
            case cp::relative_to_type::chaste_build_root:
                relative_to = RelativeTo::ChasteBuildRoot;
                break;
            case cp::relative_to_type::absolute:
                relative_to = RelativeTo::Absolute;
                break;
            default:
                NEVER_REACHED;
                break;
        }
        SetPath(relative_path, relative_to);
    }
}
