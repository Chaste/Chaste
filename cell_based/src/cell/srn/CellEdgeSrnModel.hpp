//
// Created by twin on 22/01/19.
//

#ifndef CHASTE_CELLEGESRNMODEL_H
#define CHASTE_CELLEGESRNMODEL_H

#include <vector>
#include "ChasteSerialization.hpp"
#include "ClassIsAbstract.hpp"
#include <boost/serialization/base_object.hpp>
#include <boost/serialization/vector.hpp>
#include "AbstractSrnModel.hpp"
#include "CellCycleModelOdeHandler.hpp"
#include "SimulationTime.hpp"


class CellEdgeSrnModel : public AbstractSrnModel {

private:

    /** Needed for serialization. */
    friend class boost::serialization::access;
    /**
     * Archive the SRN model and member variables.
     *
     * @param archive the archive
     * @param version the current version of this class
     */
    template<class Archive>
    void serialize(Archive & archive, const unsigned int version)
    {
        archive & boost::serialization::base_object<AbstractSrnModel>(*this);
        archive & boost::serialization::base_object<CellCycleModelOdeHandler>(*this);
//        archive & mInitialConditions;
//        archive & mStateSize;
    }

    std::vector<boost::shared_ptr<AbstractSrnModel>> mEdgeSrnModels;

protected:

    using AbstractSrnModel::Initialise;


    void Initialise() override {

    }

    /*  */
    CellEdgeSrnModel(const CellEdgeSrnModel &rModel) : AbstractSrnModel(rModel) {

        //Makes a copy of all SRN models inside the system
        for(auto srnModel: rModel.mEdgeSrnModels){
            this->mEdgeSrnModels.push_back(boost::shared_ptr<AbstractSrnModel>(srnModel->CreateSrnModel()));
        }
    }

public:

    CellEdgeSrnModel() {}

    ~CellEdgeSrnModel(){}

    void SimulateToCurrentTime() override {

        for(auto srnModel: mEdgeSrnModels){
            srnModel->SimulateToCurrentTime();
        }

    }

    AbstractSrnModel* CreateSrnModel() override {
        return new CellEdgeSrnModel(*this);
    }

    void AddEdgeSrn(std::vector<boost::shared_ptr<AbstractSrnModel>> edgeSrn){
        mEdgeSrnModels = edgeSrn;
    }

    void AddEdgeSrn(boost::shared_ptr<AbstractSrnModel> edgeSrn){
        edgeSrn->SetEdgeLocalIndex(mEdgeSrnModels.size());
        mEdgeSrnModels.push_back(edgeSrn);

    }

    unsigned GetNumEdgeSrn()
    {
        return mEdgeSrnModels.size();
    }

    boost::shared_ptr<AbstractSrnModel> GetEdgeSrn(unsigned index)
    {
        assert(index < mEdgeSrnModels.size());
        return mEdgeSrnModels[index];
    }


};


#endif //CHASTE_CELLEGESRNMODEL_H
