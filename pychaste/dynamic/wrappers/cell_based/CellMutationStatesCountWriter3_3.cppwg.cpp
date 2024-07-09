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
#include "CellMutationStatesCountWriter.hpp"

#include "CellMutationStatesCountWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellMutationStatesCountWriter<3,3 > CellMutationStatesCountWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellMutationStatesCountWriter3_3_Overrides : public CellMutationStatesCountWriter3_3{
    public:
    using CellMutationStatesCountWriter3_3::CellMutationStatesCountWriter;
    void WriteHeader(::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            WriteHeader,
                    pCellPopulation);
    }
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_CellMutationStatesCountWriter3_3_class(py::module &m){
py::class_<CellMutationStatesCountWriter3_3 , CellMutationStatesCountWriter3_3_Overrides , boost::shared_ptr<CellMutationStatesCountWriter3_3 >  , AbstractCellPopulationCountWriter<3, 3>  >(m, "CellMutationStatesCountWriter3_3")
        .def(py::init< >())
        .def(
            "WriteHeader",
            (void(CellMutationStatesCountWriter3_3::*)(::AbstractCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::WriteHeader,
            " " , py::arg("pCellPopulation") )
        .def(
            "VisitAnyPopulation",
            (void(CellMutationStatesCountWriter3_3::*)(::AbstractCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter3_3::*)(::CaBasedCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &CellMutationStatesCountWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
