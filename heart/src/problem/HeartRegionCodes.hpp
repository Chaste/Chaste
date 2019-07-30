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
#ifndef HEARTREGIONCODES_HPP_
#define HEARTREGIONCODES_HPP_

#include <climits>

/** Type for region codes */
typedef unsigned HeartRegionType;

/**
 * Codes that can be used to annotate regions of a cardiac mesh.
 *
 * See Node::GetRegion, Node::SetRegion, AbstractElement::GetUnsignedAttribute, AbstractElement::SetAttribute.
 *
 * Note: these constants are set explicitly to be of type unsigned, so as to match
 * the above methods.  Hence why we use a class instead of an enum - you can't
 * (until C++0x) specify the underlying type of an enum.
 */
class HeartRegionCode
{

public:
    /** @return a valid tissue identifier */
    static HeartRegionType GetValidTissueId();

    /** @return a valid bath identifier */
    static HeartRegionType GetValidBathId();

    /**
     *  @return For a given region identifier, determines whether it is a tissue identifier
     *
     *  @param regionId region identifier
     */
    static bool IsRegionTissue(HeartRegionType regionId);

    /**
     *  @return For a given region identifier, determines whether it is a bath identifier
     *
     *  @param regionId region identifier
     */
    static bool IsRegionBath(HeartRegionType regionId);

private:
    /** No instances of this class should be created. */
    HeartRegionCode();
    /** No instances of this class should be created. */
    HeartRegionCode(const HeartRegionCode&);
    /** No instances of this class should be created.
     * @return reference by language convention*/
    HeartRegionCode& operator=(const HeartRegionCode&);
};

#endif /*HEARTREGIONCODES_HPP_*/
