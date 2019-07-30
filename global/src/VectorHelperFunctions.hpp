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

#ifndef VECTORHELPERFUNCTIONS_HPP_
#define VECTORHELPERFUNCTIONS_HPP_

/**
 * @file
 *
 * A selection of helper functions to be able to access std::vector<double>
 * and CVODE's N_Vector types using the same interface.  These are used by
 * AbstractParameterisedSystem and some tests.
 */

#include <cassert>
#include <vector>

#ifdef CHASTE_CVODE
// CVODE headers
#include <nvector/nvector_serial.h>

#endif

/**
 * @return The vector component
 * Helper function to get a vector component.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector to access
 * @param index  the index of the component to get
 */
template <typename VECTOR>
inline double GetVectorComponent(const VECTOR& rVec, unsigned index);

/**
 * Helper function to set a vector component.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector to modify
 * @param index  the index of the component to set
 * @param value  the new value
 */
template <typename VECTOR>
inline void SetVectorComponent(VECTOR& rVec, unsigned index, double value);

/**
 * Helper function to determine a vector's size.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector
 * @return  its size
 */
template <typename VECTOR>
inline unsigned GetVectorSize(const VECTOR& rVec);

/**
 * Helper function to initialise a vector to be empty/unset.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector
 */
template <typename VECTOR>
inline void InitialiseEmptyVector(VECTOR& rVec);

/**
 * If the given vector is empty, set it to have a particular size, allocating memory if necessary.
 * It is a no-op if the vector is non-empty.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the empty vector
 * @param size  the size to create it as
 */
template <typename VECTOR>
inline void CreateVectorIfEmpty(VECTOR& rVec, unsigned size);

/**
 * Helper function to create a new empty/unset vector, useful for defining a default parameter value.
 *
 * @return a new vector
 */
template <typename VECTOR>
inline VECTOR CreateEmptyVector();

/**
 * Helper function to test whether a vector is empty/unset.
 *
 * @param rVec  the vector
 * @return true if empty
 */
template <typename VECTOR>
inline bool IsEmptyVector(VECTOR& rVec);

/**
 * Helper function to delete a vector.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector
 */
template <typename VECTOR>
inline void DeleteVector(VECTOR& rVec);

/**
 * Helper function to create a fresh copy of a vector.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rVec  the vector to copy
 * @return A copy
 */
template <typename VECTOR>
inline VECTOR CopyVector(VECTOR& rVec);

/**
 * A helper function to copy a VECTOR into a std::vector<double>.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rSrc  source vector
 * @param rDest  destination vector; will be resized and filled
 */
template <typename VECTOR>
inline void CopyToStdVector(const VECTOR& rSrc, std::vector<double>& rDest);

/**
 * A helper function to copy a std::vector<double> into a VECTOR.
 *
 * This isn't a member so that we can specialise it without having to
 * specialise the whole class.
 *
 * @param rSrc  source vector
 * @param rDest  destination vector; must exist and be the correct size
 */
template <typename VECTOR>
inline void CopyFromStdVector(const std::vector<double>& rSrc, VECTOR& rDest);

// Specialisations for std::vector<double>

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @param index
 * @return vector component
 */
template <>
inline double GetVectorComponent(const std::vector<double>& rVec, unsigned index)
{
    assert(index < rVec.size());
    return rVec[index];
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @param index
 * @param value
 */
template <>
inline void SetVectorComponent(std::vector<double>& rVec, unsigned index, double value)
{
    assert(index < rVec.size());
    rVec[index] = value;
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @return size
 */
template <>
inline unsigned GetVectorSize(const std::vector<double>& rVec)
{
    return rVec.size();
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 */
template <>
inline void InitialiseEmptyVector(std::vector<double>& rVec)
{
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @param size
 */
template <>
inline void CreateVectorIfEmpty(std::vector<double>& rVec, unsigned size)
{
    if (rVec.empty())
    {
        rVec.resize(size);
    }
}

/**
 * Specialisation for std::vector<double>.
 * @return empty vector
 */
template <>
inline std::vector<double> CreateEmptyVector()
{
    return std::vector<double>();
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @return true if empty
 */
template <>
inline bool IsEmptyVector(std::vector<double>& rVec)
{
    return rVec.empty();
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 */
template <>
inline void DeleteVector(std::vector<double>& rVec)
{
}

/**
 * Specialisation for std::vector<double>.
 * @param rVec
 * @return copy
 */
template <>
inline std::vector<double> CopyVector(std::vector<double>& rVec)
{
    return rVec;
}

/**
 * Specialisation for std::vector<double>.
 * @param rSrc
 * @param rDest
 */
template <>
inline void CopyToStdVector(const std::vector<double>& rSrc, std::vector<double>& rDest)
{
    rDest = rSrc;
}

/**
 * Specialisation for std::vector<double>.
 * @param rSrc
 * @param rDest
 */
template <>
inline void CopyFromStdVector(const std::vector<double>& rSrc, std::vector<double>& rDest)
{
    rDest = rSrc;
}

// Specialisations for N_Vector

#ifdef CHASTE_CVODE

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @param index
 * @return the vector component
 */
template <>
inline double GetVectorComponent(const N_Vector& rVec, unsigned index)
{
    assert(rVec != nullptr);
    return NV_Ith_S(rVec, index);
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @param index
 * @param value
 */
template <>
inline void SetVectorComponent(N_Vector& rVec, unsigned index, double value)
{
    assert(rVec != nullptr);
    NV_Ith_S(rVec, index) = value;
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @return size
 */
template <>
inline unsigned GetVectorSize(const N_Vector& rVec)
{
    assert(rVec != nullptr);
    return NV_LENGTH_S(rVec);
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 */
template <>
inline void InitialiseEmptyVector(N_Vector& rVec)
{
    rVec = nullptr;
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @param size
 */
template <>
inline void CreateVectorIfEmpty(N_Vector& rVec, unsigned size)
{
    if (rVec == nullptr)
    {
        rVec = N_VNew_Serial(size);
    }
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @return empty vector
 */
template <>
inline N_Vector CreateEmptyVector()
{
    return nullptr;
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @return true if empty
 */
template <>
inline bool IsEmptyVector(N_Vector& rVec)
{
    return rVec == nullptr;
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 */
template <>
inline void DeleteVector(N_Vector& rVec)
{
    if (rVec)
    {
        rVec->ops->nvdestroy(rVec);
        rVec = nullptr;
    }
}

/**
 * Specialisation for CVODE's N_Vector type.
 * @param rVec
 * @return copy
 */
template <>
inline N_Vector CopyVector(N_Vector& rVec)
{
    N_Vector copy = nullptr;
    if (rVec)
    {
        copy = N_VClone(rVec);
        unsigned size = NV_LENGTH_S(rVec);
        for (unsigned i = 0; i < size; i++)
        {
            NV_Ith_S(copy, i) = NV_Ith_S(rVec, i);
        }
    }
    return copy;
}

/**
 * A helper function to copy an N_Vector into a std::vector<double>.
 *
 * @param rSrc  source vector
 * @param rDest  destination vector; will be resized and filled
 */
template <>
inline void CopyToStdVector(const N_Vector& rSrc, std::vector<double>& rDest)
{
    // Check for no-op
    realtype* p_src = NV_DATA_S(rSrc);
    if (!rDest.empty() && p_src == &(rDest[0]))
        return;
    // Set dest size
    long size = NV_LENGTH_S(rSrc);
    rDest.resize(size);
    // Copy data
    for (long i = 0; i < size; i++)
    {
        rDest[i] = p_src[i];
    }
}

/**
 * A helper function to copy a std::vector<double> into an N_Vector.
 *
 * @param rSrc  source vector
 * @param rDest  destination vector; must exist and be the correct size
 */
template <>
inline void CopyFromStdVector(const std::vector<double>& rSrc, N_Vector& rDest)
{
    // Check for no-op
    realtype* p_dest = NV_DATA_S(rDest);
    if (p_dest == &(rSrc[0])) return;

    // Check dest size
    long size = NV_LENGTH_S(rDest);
    assert(size == (long)rSrc.size());

    // Copy data
    for (long i = 0; i < size; i++)
    {
        p_dest[i] = rSrc[i];
    }
}

/**
 * Make a standard vector from an N_Vector
 *
 * @param v  A CVODE N_Vector to copy the entries from
 * @return a std::vector of the entries of the N_Vector.
 */
inline std::vector<double> MakeStdVec(N_Vector v)
{
    std::vector<double> sv;
    CopyToStdVector(v, sv);
    return sv;
}

/**
 * Make an N_Vector from a standard vector
 *
 * @param rSrc a std::vector to copy the entries from
 * @return an N_Vector of the entries from the std::vector.
 */
inline N_Vector MakeNVector(const std::vector<double>& rSrc)
{
    N_Vector nv = nullptr;
    CreateVectorIfEmpty(nv, rSrc.size());
    CopyFromStdVector(rSrc, nv);
    return nv;
}

#endif // CHASTE_CVODE

// End of helper functions

#endif /*VECTORHELPERFUNCTIONS_HPP_*/
