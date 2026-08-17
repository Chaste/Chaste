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

#ifndef _NODEATTRIBUTES_HPP_
#define _NODEATTRIBUTES_HPP_

#include "UblasVectorInclude.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/vector.hpp>
#include <boost/serialization/set.hpp>
#include <bitset>
#include <iterator>

/**
 * A container for attributes associated with the Node class.
 */
#define numNodes 4096

template <std::size_t N>
class MyBitset
{
public:
  using value_type      = unsigned;
  using size_type       = std::size_t;
  using difference_type = std::ptrdiff_t;
  using reference       = unsigned;
  using const_reference = unsigned;

  class iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type        = unsigned;
    using difference_type   = std::ptrdiff_t;
    using reference         = unsigned;
    using pointer           = void;

    iterator(MyBitset* c, size_type pos): container_(c), pos_(pos)   {    }

    reference operator*() const          { return pos_; }
    iterator& operator++()          {  nextBit();      return *this;      }
    iterator operator++(int)        {  iterator tmp = *this;   nextBit();      return tmp;      }
    friend bool operator==(const iterator& a, const iterator& b)      {    return a.pos_ == b.pos_;      }
    friend bool operator!=(const iterator& a, const iterator& b)      {    return !(a == b);      }
  private:
    void nextBit( )
    {
      while( ++pos_ < N && !container_->bits_.test( pos_ ) );
    }
    MyBitset* container_;
    size_type pos_;
  };
  class const_iterator
  {
  public:
    using iterator_category = std::forward_iterator_tag;
    using value_type        = unsigned;
    using difference_type   = std::ptrdiff_t;
    using reference         = unsigned;
    using pointer           = void;

    const_iterator(const MyBitset* c, size_type pos): container_(c), pos_(pos)      {      }
    reference operator*() const      { return pos_;        }
    const_iterator& operator++()     { nextBit(); return *this;  }
    const_iterator operator++(int)   { auto tmp = *this; nextBit(); return tmp;      }
    friend bool operator==(const const_iterator& a,
               const const_iterator& b)     { return a.pos_ == b.pos_;   }
    friend bool operator!=(const const_iterator& a,
               const const_iterator& b)     { return !(a == b);      }
  private:
    void nextBit( )
    {
      while( ++pos_ < N && !container_->bits_.test( pos_ ) );
    }
    const MyBitset* container_;
    size_type pos_;
  };
  private:
  iterator nextBit( size_type pos = 0 )
  {
    if( bits_.none() ) return iterator( this, N );
    for( size_type i = pos; i < numNodes; ++i )
      if ( bits_.test( i ))
    return iterator( this, i );
    return iterator( this, N );
  }
private:
  std::bitset<N> bits_;
  friend class boost::serialization::access;

  /**
   * not sure where this is use of if it is needed, but without it compile fails
   */
  template< class Archive >
  void serialize( Archive& ar, const unsigned int version )
  {
    std::vector<unsigned> temp;
    for( int x = 0; x < 1024; ++x )
      {
      if ( bits_.test(x) )
    temp.push_back( x );
      }
    ar & temp;
  }
public:
  bool test( size_t x ) const { return bits_.test( x ); }
  iterator begin()
  {
    return nextBit();
  }
  size_t         count() { return bits_.count(); }
  iterator               end()           { return iterator(this, N);      }
  const_iterator         begin() const   { return const_iterator(this, 0);    }
  const_iterator         end() const     { return const_iterator(this, N);    }

  const_iterator         cbegin() const  { return begin();    }
  const_iterator         cend() const    { return end();    }
  void                   insert( unsigned idx )  {    bits_.set( idx );  }
  void             clear()        { bits_.reset(); }
  bool             empty() const     { return bits_.none(); }
  unsigned erase( unsigned idx )
  {
    unsigned ret = bits_.test( idx );
    bits_.reset( idx );
    return ret;
  }
  auto count() const  {    return bits_.count();  }
  unsigned count(unsigned idx) const   {    return bits_.test( idx ) ? 1 : 0;  }
  iterator find( unsigned idx )
  {
    if ( bits_.test( idx ) )
      return iterator( this, idx );
    else
      return iterator( this, N );
  }
  const_iterator find( const unsigned idx ) const
  {
    if ( bits_.test( idx ) )
      return const_iterator( this, idx );
    else
      return const_iterator( this, N );
  }
};

template<unsigned SPACE_DIM>
class NodeAttributes
{
private:

    /** Arbitrary attributes that a user gives meaning to */
    std::vector<double> mAttributes;

    /** The ID of the region of mesh in which the Node lies */
  unsigned mRegion{ 0u };

    /** For mutable nodes in OffLatticeSimulations, a container for the force accumulated on this node. */
  c_vector<double, SPACE_DIM> mAppliedForce{ zero_vector<double>(SPACE_DIM) };

    /** The radius associated with the Node */
  double mRadius{ 0.0 };

    /** Vector of indices corresponding to neighbouring nodes. */
  MyBitset<numNodes>  mNeighbourIndices;

  /** A bool indicating whether the neighbours of this node have been calculated yet. */
  bool mNeighboursSetUp{ false };

    /** Whether the node represents a particle or not: Used for NodeBasedCellPopulationWithParticles */
  bool mIsParticle{ false };

    /** Needed for serialization. */
    friend class boost::serialization::access;

    /**
     * Archive the member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & mAttributes;
        archive & mRegion;
        archive & mRadius;
        archive & mNeighbourIndices;
        archive & mNeighboursSetUp;
        archive & mIsParticle;

        for (unsigned d = 0; d < SPACE_DIM; d++)
        {
            archive & mAppliedForce[d];
        }
    }

public:

    /**
     * Defaults all variables.
     */
  NodeAttributes() = default;

    /**
     * @return mAttributes
     */
    std::vector<double>& rGetAttributes();

    /**
     * Push an attribute back onto mAttributes
     *
     * @param attribute the value of the attribute.
     */
    void AddAttribute(double attribute);

    /**
     * Get the region ID
     *
     * @return mRegion
     */
    unsigned GetRegion();

    /**
     * Set the region ID
     *
     * @param region the value to to assign to mRegion.
     */
    void SetRegion(unsigned region);

    /**
     * Get the current value of the applied force on the node.
     *
     * @return mAppliedForce
     */
    c_vector<double, SPACE_DIM>& rGetAppliedForce();

    /**
     * Add a contribution to the force vector
     *
     * @param rForceContribution the contribution to add to mAppliedForce
     */
    void AddAppliedForceContribution(const c_vector<double, SPACE_DIM>& rForceContribution);

    /**
     * Set mAppliedForce to a zero vector.
     */
    void ClearAppliedForce();

    /**
     * Add a neighbour to this node's vector of neighbouring node  indices.
     *
     * @param index of the node to add.
     */
    void AddNeighbour(unsigned index);

    /**
     * Clear this node's vector of neighbour indices.
     */
    void ClearNeighbours();

    /**
     * Remove duplicates from the vector of node neighbour indices.
     */
  void RemoveDuplicateNeighbours(){}

    /**
     * Check whether the node neighbours collection is empty.
     *
     * @return whether this node has any neighbours.
     */
    bool NeighboursIsEmpty();

    /**
     * Sets a flag to indicate that the neighbours of this node have/have not been updated.
     *
     * @param flag whether the neighbours are set up or not
     */
    void SetNeighboursSetUp(bool flag);

    /**
     * @return a flag to indicate that the neighbours of this node have/have not been updated.
     */
    bool GetNeighboursSetUp();

    /**
     * @return this node's vector of neighbour indices.
     */
  std::vector<unsigned> rGetNeighbours() const;
  /**
   *  @return begin iterator to neighbours
   */
  MyBitset<numNodes>::iterator getNeighboursBegin() { return mNeighbourIndices.begin(); }

  /**
   * @return end iterator to neighbours
   */
  MyBitset<numNodes>::iterator getNeighboursEnd() { return mNeighbourIndices.end(); }

  /**
   * @return number of neighbours
   */
  size_t getNeighboursCount() { return mNeighbourIndices.count(); }

    /**
     * Get whether this node is a particle, or not.
     *
     * @return mIsParticle
     */
    bool IsParticle();

    /**
     * Set the flag mIsParticle.
     * @param isParticle whether this node is particle or not.
     */
    void SetIsParticle(bool isParticle);

    /**
     * Return the radius associated with the Node
     *
     * @return mRadius
     */
    double GetRadius();

    /**
     * Set the value of the radius.
     *
     * @param radius the value to assign to mRadius. Must be >= 0.0
     */
    void SetRadius(double radius);
};

#endif //_NODEATTRIBUTES_HPP_
