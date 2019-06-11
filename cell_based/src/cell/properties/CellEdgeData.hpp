
#ifndef CELLEDGEDATA_HPP_
#define CELLEDGEDATA_HPP_


#include <boost/shared_ptr.hpp>
#include <map>
#include <string>
#include <vector>

#include "AbstractCellProperty.hpp"
#include "ChasteSerialization.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/map.hpp>
#include "Exception.hpp"


class CellEdgeData  : public AbstractCellProperty{
private:

    /**
     * The cell data.
     */
    std::map<std::string, std::vector<double>> mCellEdgeData;

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
        archive & boost::serialization::base_object<AbstractCellProperty>(*this);
        archive & mCellEdgeData;
    }

public:

    /**
     * We need the empty virtual destructor in this class to ensure Boost
     * serialization works correctly with static libraries.
     */
    virtual ~CellEdgeData();

    /**
     * This assigns the cell data.
     *
     * @param rVariableName the name of the data to be set.
     * @param data the value to set it to.
     */
    void SetItem(const std::string& rVariableName, std::vector<double> data);

    /**
     * @return data.
     *
     * @param rVariableName the index of the data required.
     * throws if rVariableName has not been stored
     */
    std::vector<double> GetItem(const std::string& rVariableName) const;

    double GetItemAtIndex(const std::string& rVariableName, const unsigned int index);

    /**
     * @return number of data items
     */
    unsigned GetNumItems() const;

    /**
     * @return all keys.
     *
     * According to STL these are sorted in lexicographical/alphabetic order (so that the ordering here is predictable).
     */
    std::vector<std::string> GetKeys() const;
};

#include "SerializationExportWrapper.hpp"
// Declare identifier for the serializer
CHASTE_CLASS_EXPORT(CellEdgeData)


#endif //CELLEDGEDATA_HPP_
