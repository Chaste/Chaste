#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"

#include "AbstractOffLatticeCellPopulation2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOffLatticeCellPopulation<2,2 > AbstractOffLatticeCellPopulation2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractOffLatticeCellPopulation2_2_Overrides : public AbstractOffLatticeCellPopulation2_2{
    public:
    using AbstractOffLatticeCellPopulation2_2::AbstractOffLatticeCellPopulation;
    unsigned int AddNode(::Node<2> * pNewNode) override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractOffLatticeCellPopulation2_2,
            AddNode,
                    pNewNode);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<2> & rNewLocation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation2_2,
            SetNode,
                    nodeIndex,
        rNewLocation);
    }
    void UpdateNodeLocations(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation2_2,
            UpdateNodeLocations,
                    dt);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 2> & rDisplacement, double dt) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation2_2,
            CheckForStepSizeException,
                    nodeIndex,
        rDisplacement,
        dt);
    }
    double GetDampingConstant(unsigned int nodeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractOffLatticeCellPopulation2_2,
            GetDampingConstant,
                    nodeIndex);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation2_2,
            OutputCellPopulationParameters,
                    rParamsFile);
    }

};
void register_AbstractOffLatticeCellPopulation2_2_class(py::module &m){
py::class_<AbstractOffLatticeCellPopulation2_2 , AbstractOffLatticeCellPopulation2_2_Overrides , boost::shared_ptr<AbstractOffLatticeCellPopulation2_2 >  , AbstractCellPopulation<2>  >(m, "AbstractOffLatticeCellPopulation2_2")
        .def(
            "AddNode",
            (unsigned int(AbstractOffLatticeCellPopulation2_2::*)(::Node<2> *)) &AbstractOffLatticeCellPopulation2_2::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "SetNode",
            (void(AbstractOffLatticeCellPopulation2_2::*)(unsigned int, ::ChastePoint<2> &)) &AbstractOffLatticeCellPopulation2_2::SetNode,
            " " , py::arg("nodeIndex"), py::arg("rNewLocation") )
        .def(
            "UpdateNodeLocations",
            (void(AbstractOffLatticeCellPopulation2_2::*)(double)) &AbstractOffLatticeCellPopulation2_2::UpdateNodeLocations,
            " " , py::arg("dt") )
        .def(
            "CheckForStepSizeException",
            (void(AbstractOffLatticeCellPopulation2_2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2> &, double)) &AbstractOffLatticeCellPopulation2_2::CheckForStepSizeException,
            " " , py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt") )
        .def(
            "GetDampingConstant",
            (double(AbstractOffLatticeCellPopulation2_2::*)(unsigned int)) &AbstractOffLatticeCellPopulation2_2::GetDampingConstant,
            " " , py::arg("nodeIndex") )
        .def(
            "SetDampingConstantNormal",
            (void(AbstractOffLatticeCellPopulation2_2::*)(double)) &AbstractOffLatticeCellPopulation2_2::SetDampingConstantNormal,
            " " , py::arg("dampingConstantNormal") )
        .def(
            "SetDampingConstantMutant",
            (void(AbstractOffLatticeCellPopulation2_2::*)(double)) &AbstractOffLatticeCellPopulation2_2::SetDampingConstantMutant,
            " " , py::arg("dampingConstantMutant") )
        .def(
            "SetAbsoluteMovementThreshold",
            (void(AbstractOffLatticeCellPopulation2_2::*)(double)) &AbstractOffLatticeCellPopulation2_2::SetAbsoluteMovementThreshold,
            " " , py::arg("absoluteMovementThreshold") )
        .def(
            "GetAbsoluteMovementThreshold",
            (double(AbstractOffLatticeCellPopulation2_2::*)()) &AbstractOffLatticeCellPopulation2_2::GetAbsoluteMovementThreshold,
            " "  )
        .def(
            "GetDampingConstantNormal",
            (double(AbstractOffLatticeCellPopulation2_2::*)()) &AbstractOffLatticeCellPopulation2_2::GetDampingConstantNormal,
            " "  )
        .def(
            "GetDampingConstantMutant",
            (double(AbstractOffLatticeCellPopulation2_2::*)()) &AbstractOffLatticeCellPopulation2_2::GetDampingConstantMutant,
            " "  )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractOffLatticeCellPopulation2_2::*)(::out_stream &)) &AbstractOffLatticeCellPopulation2_2::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
    ;
}
