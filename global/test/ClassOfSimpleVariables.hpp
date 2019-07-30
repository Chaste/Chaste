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

#ifndef CLASSOFSIMPLEVARIABLES_HPP_
#define CLASSOFSIMPLEVARIABLES_HPP_

#include <vector>
#include <string>
#include <boost/serialization/vector.hpp>

class ClassOfSimpleVariables
{
private:
    friend class boost::serialization::access;
    friend class TestObjectCommunicator;

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & mNumber;
        archive & mString;
        archive & mVectorOfDoubles; // include <boost/serialization/vector.hpp> for this
        archive & mVectorOfBools;
    }

    int mNumber;
    std::string mString;
    std::vector<double> mVectorOfDoubles;
    std::vector<bool> mVectorOfBools;

public:

    ClassOfSimpleVariables()
    {
        // Do nothing. Used when loading into a pointer.
    }
    ClassOfSimpleVariables(int initial,
                           std::string string,
                           std::vector<double> doubles,
                           std::vector<bool> bools)
        : mString(string),
          mVectorOfDoubles(doubles),
          mVectorOfBools(bools)
    {
        mNumber = initial;
    }

    int GetNumber() const
    {
        return mNumber;
    }

    std::string GetString()
    {
        return mString;
    }

    std::vector<double>& GetVectorOfDoubles()
    {
        return mVectorOfDoubles;
    }

    std::vector<bool>& GetVectorOfBools()
    {
        return mVectorOfBools;
    }
};

#endif /*CLASSOFSIMPLEVARIABLES_HPP_*/
