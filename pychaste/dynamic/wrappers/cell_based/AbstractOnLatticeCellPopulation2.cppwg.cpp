#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOnLatticeCellPopulation.hpp"

#include "AbstractOnLatticeCellPopulation2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOnLatticeCellPopulation<2 > AbstractOnLatticeCellPopulation2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::vector<boost::shared_ptr<AbstractUpdateRule<2>>> const _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_2_gt__gt__gt_const;

class AbstractOnLatticeCellPopulation2_Overrides : public AbstractOnLatticeCellPopulation2{
    public:
    using AbstractOnLatticeCellPopulation2::AbstractOnLatticeCellPopulation;
    void UpdateCellLocations(double dt) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOnLatticeCellPopulation2,
            UpdateCellLocations,
                    dt);
    }
    void SetNode(unsigned int index, ::ChastePoint<2> & rNewLocation) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation2,
            SetNode,
                    index,
        rNewLocation);
    }
    ::std::set<unsigned int> GetNeighbouringNodeIndices(unsigned int index) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            AbstractOnLatticeCellPopulation2,
            GetNeighbouringNodeIndices,
                    index);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation2,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE(
            double,
            AbstractOnLatticeCellPopulation2,
            GetDefaultTimeStep,
            );
    }
    void AddUpdateRule(::boost::shared_ptr<AbstractUpdateRule<2>> pUpdateRule) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOnLatticeCellPopulation2,
            AddUpdateRule,
                    pUpdateRule);
    }
    void RemoveAllUpdateRules() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOnLatticeCellPopulation2,
            RemoveAllUpdateRules,
            );
    }
    ::std::vector<boost::shared_ptr<AbstractUpdateRule<2>>> const GetUpdateRuleCollection() const  override {
        PYBIND11_OVERRIDE(
            _std_vector_lt_boost_shared_ptr_lt_AbstractUpdateRule_lt_2_gt__gt__gt_const,
            AbstractOnLatticeCellPopulation2,
            GetUpdateRuleCollection,
            );
    }

};
void register_AbstractOnLatticeCellPopulation2_class(py::module &m){
py::class_<AbstractOnLatticeCellPopulation2 , AbstractOnLatticeCellPopulation2_Overrides , boost::shared_ptr<AbstractOnLatticeCellPopulation2 >  , AbstractCellPopulation<2>  >(m, "AbstractOnLatticeCellPopulation2")
        .def(
            "UpdateCellLocations",
            (void(AbstractOnLatticeCellPopulation2::*)(double)) &AbstractOnLatticeCellPopulation2::UpdateCellLocations,
            " " , py::arg("dt") )
        .def(
            "GetUpdateNodesInRandomOrder",
            (bool(AbstractOnLatticeCellPopulation2::*)()) &AbstractOnLatticeCellPopulation2::GetUpdateNodesInRandomOrder,
            " "  )
        .def(
            "SetUpdateNodesInRandomOrder",
            (void(AbstractOnLatticeCellPopulation2::*)(bool)) &AbstractOnLatticeCellPopulation2::SetUpdateNodesInRandomOrder,
            " " , py::arg("updateNodesInRandomOrder") )
        .def(
            "SetIterateRandomlyOverUpdateRuleCollection",
            (void(AbstractOnLatticeCellPopulation2::*)(bool)) &AbstractOnLatticeCellPopulation2::SetIterateRandomlyOverUpdateRuleCollection,
            " " , py::arg("iterateRandomly") )
        .def(
            "GetIterateRandomlyOverUpdateRuleCollection",
            (bool(AbstractOnLatticeCellPopulation2::*)()) &AbstractOnLatticeCellPopulation2::GetIterateRandomlyOverUpdateRuleCollection,
            " "  )
        .def(
            "SetNode",
            (void(AbstractOnLatticeCellPopulation2::*)(unsigned int, ::ChastePoint<2> &)) &AbstractOnLatticeCellPopulation2::SetNode,
            " " , py::arg("index"), py::arg("rNewLocation") )
        .def(
            "GetNeighbouringNodeIndices",
            (::std::set<unsigned int>(AbstractOnLatticeCellPopulation2::*)(unsigned int)) &AbstractOnLatticeCellPopulation2::GetNeighbouringNodeIndices,
            " " , py::arg("index") )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractOnLatticeCellPopulation2::*)(::out_stream &)) &AbstractOnLatticeCellPopulation2::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetDefaultTimeStep",
            (double(AbstractOnLatticeCellPopulation2::*)()) &AbstractOnLatticeCellPopulation2::GetDefaultTimeStep,
            " "  )
        .def(
            "AddUpdateRule",
            (void(AbstractOnLatticeCellPopulation2::*)(::boost::shared_ptr<AbstractUpdateRule<2>>)) &AbstractOnLatticeCellPopulation2::AddUpdateRule,
            " " , py::arg("pUpdateRule") )
        .def(
            "RemoveAllUpdateRules",
            (void(AbstractOnLatticeCellPopulation2::*)()) &AbstractOnLatticeCellPopulation2::RemoveAllUpdateRules,
            " "  )
        .def(
            "GetUpdateRuleCollection",
            (::std::vector<boost::shared_ptr<AbstractUpdateRule<2>>> const(AbstractOnLatticeCellPopulation2::*)() const ) &AbstractOnLatticeCellPopulation2::GetUpdateRuleCollection,
            " "  )
    ;
}
