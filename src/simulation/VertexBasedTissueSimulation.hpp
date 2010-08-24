#ifndef VERTEXBASEDTISSUESIMULATION_HPP_
#define VERTEXBASEDTISSUESIMULATION_HPP_

#include <map>
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>

#include "TissueSimulation.hpp"
#include "VertexBasedTissue.hpp"
#include "OutputFileHandler.hpp"

/**
 * A Vertex-based tissue simulation object.
 */
// template<unsigned DIM>
class VertexBasedTissueSimulation : public TissueSimulation<2>
{
    // Allow tests to access private members, in order to test computation of
    // private functions eg. DoCellBirth
    friend class TestVertexBasedTissueSimulation;

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Serialize the object and its member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */

    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        // If Archive is an output archive, then & resolves to <<
        // If Archive is an input archive, then & resolves to >>
        archive & boost::serialization::base_object<TissueSimulation<2> >(*this);
     }

        /**
     *  Overridden PostSolve() method.
     */
    void PostSolve();

public:

    VertexBasedTissueSimulation(AbstractTissue<2>& rTissue,
                      std::vector<AbstractForce<2>*> forceCollection,
                      bool deleteTissueAndForceCollection=false,
                      bool initialiseCells=true);

    /**
     * Destructor.
     */
    ~VertexBasedTissueSimulation();

};

// Declare identifier for the serializer
#include "SerializationExportWrapper.hpp"
CHASTE_CLASS_EXPORT(VertexBasedTissueSimulation)

namespace boost
{
namespace serialization
{
/**
 * Serialize information required to construct a VertexBasedTissueSimulation.
 */
template<class Archive>//, unsigned DIM>
inline void save_construct_data(
    Archive & ar, const VertexBasedTissueSimulation * t, const BOOST_PFTO unsigned int file_version)
{
    // Save data required to construct instance
    const AbstractTissue<2> * p_tissue = &(t->rGetTissue());
    ar & p_tissue;
    const std::vector<AbstractForce<2>*> force_collection = t->rGetForceCollection();
    ar & force_collection;
}

/**
 * De-serialize constructor parameters and initialise a VertexBasedTissueSimulation.
 */
template<class Archive>//, unsigned DIM>
inline void load_construct_data(
    Archive & ar, VertexBasedTissueSimulation * t, const unsigned int file_version)
{
    // Retrieve data from archive required to construct new instance
    AbstractTissue<2>* p_tissue;
    ar & p_tissue;
    std::vector<AbstractForce<2>*> force_collection;
    ar & force_collection;

    // Invoke inplace constructor to initialise instance
    ::new(t)VertexBasedTissueSimulation(*p_tissue, force_collection, true, false);
}
}
} // namespace

#endif /*VERTEXBASEDTISSUESIMULATION_HPP_*/
