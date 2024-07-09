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
#include "BoundaryNodeWriter.hpp"

#include "BoundaryNodeWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef BoundaryNodeWriter<3,3 > BoundaryNodeWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class BoundaryNodeWriter3_3_Overrides : public BoundaryNodeWriter3_3{
    public:
    using BoundaryNodeWriter3_3::BoundaryNodeWriter;
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BoundaryNodeWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BoundaryNodeWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BoundaryNodeWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BoundaryNodeWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BoundaryNodeWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            BoundaryNodeWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_BoundaryNodeWriter3_3_class(py::module &m){
py::class_<BoundaryNodeWriter3_3 , BoundaryNodeWriter3_3_Overrides , boost::shared_ptr<BoundaryNodeWriter3_3 >  , AbstractCellPopulationWriter<3, 3>  >(m, "BoundaryNodeWriter3_3")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(BoundaryNodeWriter3_3::*)(::AbstractCellPopulation<3> *)) &BoundaryNodeWriter3_3::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(BoundaryNodeWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &BoundaryNodeWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(BoundaryNodeWriter3_3::*)(::CaBasedCellPopulation<3> *)) &BoundaryNodeWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(BoundaryNodeWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &BoundaryNodeWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(BoundaryNodeWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &BoundaryNodeWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(BoundaryNodeWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &BoundaryNodeWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(BoundaryNodeWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &BoundaryNodeWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
