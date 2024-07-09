#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include "MeshBasedCellPopulation.hpp"
#include "PottsBasedCellPopulation.hpp"
#include "CaBasedCellPopulation.hpp"
#include "ImmersedBoundaryCellPopulation.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "VertexBasedCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "NodeLocationWriter.hpp"

#include "NodeLocationWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef NodeLocationWriter<3,3 > NodeLocationWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NodeLocationWriter3_3_Overrides : public NodeLocationWriter3_3{
    public:
    using NodeLocationWriter3_3::NodeLocationWriter;
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeLocationWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeLocationWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeLocationWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeLocationWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeLocationWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeLocationWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_NodeLocationWriter3_3_class(py::module &m){
py::class_<NodeLocationWriter3_3 , NodeLocationWriter3_3_Overrides , boost::shared_ptr<NodeLocationWriter3_3 >  , AbstractCellPopulationWriter<3, 3>  >(m, "NodeLocationWriter3_3")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(NodeLocationWriter3_3::*)(::AbstractCellPopulation<3> *)) &NodeLocationWriter3_3::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeLocationWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &NodeLocationWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeLocationWriter3_3::*)(::CaBasedCellPopulation<3> *)) &NodeLocationWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeLocationWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &NodeLocationWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeLocationWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &NodeLocationWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeLocationWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &NodeLocationWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeLocationWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &NodeLocationWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
