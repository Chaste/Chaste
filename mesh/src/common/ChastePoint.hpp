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


#ifndef _CHASTEPOINT_HPP_
#define _CHASTEPOINT_HPP_

#include "ChasteSerialization.hpp"
#include "UblasVectorInclude.hpp"
#include <vector>


/**
 * A ChastePoint class, templated over spatial dimension.
 */
template<unsigned DIM>
class ChastePoint
{
    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the member variables.
     *
     * @param archive
     * @param version
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        //archive & mLocation; //earlier versions of boost are unable to do this. See #1709
    }

private:

    /** The location of the Point. */
    c_vector<double, DIM> mLocation;

public:

    /**
     * Create a Point object.
     * There are 3 optional arguments, which can be used to specify the values
     * of the first 3 dimensions, if present.
     *
     * Point now uses a ublas vector to store its location. The
     * rGetLocation method returns a reference to this vector.
     * Use of this method together with ublas operations
     * is the perfered way to use this class.
     *
     * @param v1  the point's x-coordinate (defaults to 0)
     * @param v2  the point's y-coordinate (defaults to 0)
     * @param v3  the point's z-coordinate (defaults to 0)
     */
    ChastePoint(double v1=0, double v2=0, double v3=0);

    /**
     * Create a Point object.
     * This constructor takes a vector giving the coordinates of the point.
     * The length of the vector must be at least the dimension of the point.
     *
     * @param coords  a std::vector storing the point's coordinates
     */
    ChastePoint(std::vector<double> coords);

    /**
     * Alternative constructor which takes in a c_vector.
     *
     * @param location  a c_vector storing the point's coordinates
     */
    ChastePoint(c_vector<double, DIM> location);

    /**
     * Get the location of the Point.
     */
    c_vector<double, DIM>& rGetLocation();

    /**
     * Access the vector mLocation.
     *
     * @param i the index of the vector to return
     */
    double operator[] (unsigned i) const;

    /**
     * Get a co-ordinate, returning the default value if the co-ordinate doesn't exist.
     * @param i  the co-ordinate to get
     * @param def  the default value
     */
    double GetWithDefault(unsigned i, double def=0.0) const;

    /**
     * Set one of the coordinates of the Point.
     *
     * @param i the index of the coordinate
     * @param value the value of the coordinate
     */
    void SetCoordinate(unsigned i, double value);

    /**
     * Checks whether one chaste point is the same as the one constructed
     *
     * @param rPoint the point to be checked
     * @return true if the are the same
     */
    bool IsSamePoint(const ChastePoint<DIM>& rPoint) const;
};


#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
EXPORT_TEMPLATE_CLASS_SAME_DIMS(ChastePoint)

namespace boost
{
namespace serialization
{

template<class Archive, unsigned SPACE_DIM>
inline void save_construct_data(
    Archive & ar, const ChastePoint<SPACE_DIM> * t, const BOOST_PFTO unsigned int file_version)
{
    for (unsigned i = 0; i < SPACE_DIM; i ++)
    {
        //we archive coordinates of mLocation one by one
        //this is because earlier version of boost (<1.40, I think) cannot archive c_vectors
        double coord = t->GetWithDefault(i);
        ar & coord;
    }
}

/**
 * Allow us to not need a default constructor, by specifying how Boost should
 * instantiate a ChastePoint instance (using existing constructor)
 */
template<class Archive,unsigned SPACE_DIM>
inline void load_construct_data(
    Archive & ar, ChastePoint<SPACE_DIM> * t, const unsigned int file_version)
{
    std::vector<double> coords;
    coords.resize(SPACE_DIM);
    for (unsigned i=0; i<SPACE_DIM; i++)
    {
        double coordinate;
        ar & coordinate;//resume coordinates one by one
        coords[i] = coordinate;
    }
    //use constructor with standard vectors to re-build the object
     ::new(t)ChastePoint<SPACE_DIM>(coords);
}
}
} // namespace ...


/**
 * A zero-dimensional ChastePoint class.
 * We need to specialise the entire class to avoid nonsense methods
 * that don't make sense for 0d.
 */
template<>
class ChastePoint<0>
{
public:
    /**
     * Create a zero-dimensional Point object.
     * There are 3 optional arguments, which should not be used.
     *
     * @param v1  the point's x-coordinate (defaults to 0)
     * @param v2  the point's y-coordinate (defaults to 0)
     * @param v3  the point's z-coordinate (defaults to 0)
     */
    ChastePoint(double v1=0, double v2=0, double v3=0);

    /**
     * Access the vector mLocation.  Actually raises an exception, since
     * a 0d point has no location.
     *
     * @param i the index of the vector to return
     */
    double operator[] (unsigned i) const;
};

#endif //_CHASTEPOINT_HPP_
