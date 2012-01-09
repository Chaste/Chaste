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
