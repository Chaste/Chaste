#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOnLatticeCellPopulation.hpp"

#include "AbstractOnLatticeCellPopulation3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOnLatticeCellPopulation<3 > AbstractOnLatticeCellPopulation3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_3_gt__gt__gt_const;

class AbstractOnLatticeCellPopulation3_Overrides : public AbstractOnLatticeCellPopulation3{
    public:
    using AbstractOnLatticeCellPopulation3::AbstractOnLatticeCellPopulation;
    void UpdateCellLocations(double dt) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOnLatticeCellPopulation3,
            UpdateCellLocations,
                    dt);
    }
    void SetNode(unsigned int index, ::ChastePoint<3> & rNewLocation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation3,
            SetNode,
                    index,
        rNewLocation);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            AbstractOnLatticeCellPopulation3,
            GetNeighbouringNodeIndices,
                    index);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE(
            double,
            AbstractOnLatticeCellPopulation3,
            GetDefaultTimeStep,
            );
    }
    void AddUpdateRule(::boost::shared_ptr<AbstractUpdateRule<3>> pUpdateRule) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOnLatticeCellPopulation3,
            AddUpdateRule,
                    pUpdateRule);
    }
    void RemoveAllUpdateRules() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation3,
            RemoveAllUpdateRules,
            );
    }
    ::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const GetUpdateRuleCollection() const  override {
        PYBIND11_OVERRIDE(
            _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_3_gt__gt__gt_const,
            AbstractOnLatticeCellPopulation3,
            GetUpdateRuleCollection,
            );
    }

};
void register_AbstractOnLatticeCellPopulation3_class(py::module &m){
py::class_<AbstractOnLatticeCellPopulation3 , AbstractOnLatticeCellPopulation3_Overrides , boost::shared_ptr<AbstractOnLatticeCellPopulation3 >  , AbstractCellPopulation<3>  >(m, "AbstractOnLatticeCellPopulation3")
        .def(
            "UpdateCellLocations",
            (void(AbstractOnLatticeCellPopulation3::*)(double)) &AbstractOnLatticeCellPopulation3::UpdateCellLocations,
            " " , py::arg("dt") )
        .def(
            "GetUpdateNodesInRandomOrder",
            (bool(AbstractOnLatticeCellPopulation3::*)()) &AbstractOnLatticeCellPopulation3::GetUpdateNodesInRandomOrder,
            " "  )
        .def(
            "SetUpdateNodesInRandomOrder",
            (void(AbstractOnLatticeCellPopulation3::*)(bool)) &AbstractOnLatticeCellPopulation3::SetUpdateNodesInRandomOrder,
            " " , py::arg("updateNodesInRandomOrder") )
        .def(
            "SetIterateRandomlyOverUpdateRuleCollection",
            (void(AbstractOnLatticeCellPopulation3::*)(bool)) &AbstractOnLatticeCellPopulation3::SetIterateRandomlyOverUpdateRuleCollection,
            " " , py::arg("iterateRandomly") )
        .def(
            "GetIterateRandomlyOverUpdateRuleCollection",
            (bool(AbstractOnLatticeCellPopulation3::*)()) &AbstractOnLatticeCellPopulation3::GetIterateRandomlyOverUpdateRuleCollection,
            " "  )
        .def(
            "SetNode",
            (void(AbstractOnLatticeCellPopulation3::*)(unsigned int, ::ChastePoint<3> &)) &AbstractOnLatticeCellPopulation3::SetNode,
            " " , py::arg("index"), py::arg("rNewLocation") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(AbstractOnLatticeCellPopulation3::*)(unsigned int)) &AbstractOnLatticeCellPopulation3::GetNeighbouringNodeIndices,
            " " , py::arg("index") )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractOnLatticeCellPopulation3::*)(::out_stream &)) &AbstractOnLatticeCellPopulation3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetDefaultTimeStep",
            (double(AbstractOnLatticeCellPopulation3::*)()) &AbstractOnLatticeCellPopulation3::GetDefaultTimeStep,
            " "  )
        .def(
            "AddUpdateRule",
            (void(AbstractOnLatticeCellPopulation3::*)(::boost::shared_ptr<AbstractUpdateRule<3>>)) &AbstractOnLatticeCellPopulation3::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "RemoveAllUpdateRules",
            (void(AbstractOnLatticeCellPopulation3::*)()) &AbstractOnLatticeCellPopulation3::RemoveAllUpdateRules,
            " "  )
        .def(
            "GetUpdateRuleCollection",
            (::std::vector<boost::shared_ptr<AbstractUpdateRule<3>>> const(AbstractOnLatticeCellPopulation3::*)() const ) &AbstractOnLatticeCellPopulation3::GetUpdateRuleCollection,
            " "  )
    ;
}
