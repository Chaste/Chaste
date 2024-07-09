#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractOffLatticeCellPopulation.hpp"

#include "AbstractOffLatticeCellPopulation3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractOffLatticeCellPopulation<3,3 > AbstractOffLatticeCellPopulation3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;

class AbstractOffLatticeCellPopulation3_3_Overrides : public AbstractOffLatticeCellPopulation3_3{
    public:
    using AbstractOffLatticeCellPopulation3_3::AbstractOffLatticeCellPopulation;
    unsigned int AddNode(::Node<3> * pNewNode) override {
        PYBIND11_OVERRIDE_PURE(
            unsignedint,
            AbstractOffLatticeCellPopulation3_3,
            AddNode,
                    pNewNode);
    }
    void SetNode(unsigned int nodeIndex, ::ChastePoint<3> & rNewLocation) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation3_3,
            SetNode,
                    nodeIndex,
        rNewLocation);
    }
    void UpdateNodeLocations(double dt) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation3_3,
            UpdateNodeLocations,
                    dt);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 3> & rDisplacement, double dt) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractOffLatticeCellPopulation3_3,
            CheckForStepSizeException,
                    nodeIndex,
        rDisplacement,
        dt);
    }
    double GetDampingConstant(unsigned int nodeIndex) override {
        PYBIND11_OVERRIDE_PURE(
            double,
            AbstractOffLatticeCellPopulation3_3,
            GetDampingConstant,
                    nodeIndex);
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractOffLatticeCellPopulation3_3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }

};
void register_AbstractOffLatticeCellPopulation3_3_class(py::module &m){
py::class_<AbstractOffLatticeCellPopulation3_3 , AbstractOffLatticeCellPopulation3_3_Overrides , boost::shared_ptr<AbstractOffLatticeCellPopulation3_3 >  , AbstractCellPopulation<3>  >(m, "AbstractOffLatticeCellPopulation3_3")
        .def(
            "AddNode",
            (unsigned int(AbstractOffLatticeCellPopulation3_3::*)(::Node<3> *)) &AbstractOffLatticeCellPopulation3_3::AddNode,
            " " , py::arg("pNewNode") )
        .def(
            "SetNode",
            (void(AbstractOffLatticeCellPopulation3_3::*)(unsigned int, ::ChastePoint<3> &)) &AbstractOffLatticeCellPopulation3_3::SetNode,
            " " , py::arg("nodeIndex"), py::arg("rNewLocation") )
        .def(
            "UpdateNodeLocations",
            (void(AbstractOffLatticeCellPopulation3_3::*)(double)) &AbstractOffLatticeCellPopulation3_3::UpdateNodeLocations,
            " " , py::arg("dt") )
        .def(
            "CheckForStepSizeException",
            (void(AbstractOffLatticeCellPopulation3_3::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 3> &, double)) &AbstractOffLatticeCellPopulation3_3::CheckForStepSizeException,
            " " , py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt") )
        .def(
            "GetDampingConstant",
            (double(AbstractOffLatticeCellPopulation3_3::*)(unsigned int)) &AbstractOffLatticeCellPopulation3_3::GetDampingConstant,
            " " , py::arg("nodeIndex") )
        .def(
            "SetDampingConstantNormal",
            (void(AbstractOffLatticeCellPopulation3_3::*)(double)) &AbstractOffLatticeCellPopulation3_3::SetDampingConstantNormal,
            " " , py::arg("dampingConstantNormal") )
        .def(
            "SetDampingConstantMutant",
            (void(AbstractOffLatticeCellPopulation3_3::*)(double)) &AbstractOffLatticeCellPopulation3_3::SetDampingConstantMutant,
            " " , py::arg("dampingConstantMutant") )
        .def(
            "SetAbsoluteMovementThreshold",
            (void(AbstractOffLatticeCellPopulation3_3::*)(double)) &AbstractOffLatticeCellPopulation3_3::SetAbsoluteMovementThreshold,
            " " , py::arg("absoluteMovementThreshold") )
        .def(
            "GetAbsoluteMovementThreshold",
            (double(AbstractOffLatticeCellPopulation3_3::*)()) &AbstractOffLatticeCellPopulation3_3::GetAbsoluteMovementThreshold,
            " "  )
        .def(
            "GetDampingConstantNormal",
            (double(AbstractOffLatticeCellPopulation3_3::*)()) &AbstractOffLatticeCellPopulation3_3::GetDampingConstantNormal,
            " "  )
        .def(
            "GetDampingConstantMutant",
            (double(AbstractOffLatticeCellPopulation3_3::*)()) &AbstractOffLatticeCellPopulation3_3::GetDampingConstantMutant,
            " "  )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractOffLatticeCellPopulation3_3::*)(::out_stream &)) &AbstractOffLatticeCellPopulation3_3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
    ;
}
