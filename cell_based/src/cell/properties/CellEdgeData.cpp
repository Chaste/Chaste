//
// Created by twin on 05/03/19.
//

#include "CellEdgeData.hpp"


CellEdgeData::~CellEdgeData() {

}

void CellEdgeData::SetItem(const std::string &rVariableName, std::vector<double> data) {
    this->mCellEdgeData[rVariableName] = data;
}

std::vector<double> CellEdgeData::GetItem(const std::string &rVariableName) const {
    /*
     * Note that mCellData[rVariableName] is not const. If rVariableName is not
     * a key, then mCellData[rVariableName] will create a new item in the map
     * and increase the size by one.  Using a const_iterator ensures that the
     * map remains const.
     */
    auto it = mCellEdgeData.find(rVariableName);
    if (it == mCellEdgeData.end())
    {
        EXCEPTION("The item " << rVariableName << " is not stored");
    }
    return(it->second);
}

double CellEdgeData::GetItemAtIndex(const std::string &rVariableName, const unsigned int index) {
    /*
     * Note that mCellData[rVariableName] is not const. If rVariableName is not
     * a key, then mCellData[rVariableName] will create a new item in the map
     * and increase the size by one.  Using a const_iterator ensures that the
     * map remains const.
     */
    auto it = mCellEdgeData.find(rVariableName);
    if (it == mCellEdgeData.end())
    {
        EXCEPTION("The item " << rVariableName << " is not stored");
    }
    if(it->second.size() <= index)
    {
        EXCEPTION("The item " << rVariableName << " does not have index " << index);
    }

    return(it->second[index]);
}

unsigned CellEdgeData::GetNumItems() const {
    return mCellEdgeData.size();
}

std::vector<std::string> CellEdgeData::GetKeys() const {
    return std::vector<std::string>();
}
