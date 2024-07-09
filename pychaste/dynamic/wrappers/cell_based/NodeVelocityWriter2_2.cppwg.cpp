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

#include "NodeVelocityWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef NodeVelocityWriter<2,2 > NodeVelocityWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class NodeVelocityWriter2_2_Overrides : public NodeVelocityWriter2_2{
    public:
    using NodeVelocityWriter2_2::NodeVelocityWriter;
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            NodeVelocityWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_NodeVelocityWriter2_2_class(py::module &m){
py::class_<NodeVelocityWriter2_2 , NodeVelocityWriter2_2_Overrides , boost::shared_ptr<NodeVelocityWriter2_2 >  , AbstractCellPopulationWriter<2, 2>  >(m, "NodeVelocityWriter2_2")
        .def(py::init< >())
        .def(
            "Visit",
            (void(NodeVelocityWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &NodeVelocityWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter2_2::*)(::CaBasedCellPopulation<2> *)) &NodeVelocityWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &NodeVelocityWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &NodeVelocityWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &NodeVelocityWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(NodeVelocityWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &NodeVelocityWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
