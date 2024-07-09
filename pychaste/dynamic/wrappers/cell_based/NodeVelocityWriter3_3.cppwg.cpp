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
#include "NodeVelocityWriter.hpp"

#include "NodeVelocityWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef NodeVelocityWriter<3,3 > NodeVelocityWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NodeVelocityWriter3_3_Overrides : public NodeVelocityWriter3_3{
    public:
    using NodeVelocityWriter3_3::NodeVelocityWriter;
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_NodeVelocityWriter3_3_class(py::module &m){
py::class_<NodeVelocityWriter3_3 , NodeVelocityWriter3_3_Overrides , boost::shared_ptr<NodeVelocityWriter3_3 >  , AbstractCellPopulationWriter<3, 3>  >(m, "NodeVelocityWriter3_3")
        .def(py::init< >())
        .def(
            "Visit",
            (void(NodeVelocityWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &NodeVelocityWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter3_3::*)(::CaBasedCellPopulation<3> *)) &NodeVelocityWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &NodeVelocityWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &NodeVelocityWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &NodeVelocityWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &NodeVelocityWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
