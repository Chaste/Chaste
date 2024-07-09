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

#include "CellMutationStatesCountWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellMutationStatesCountWriter<2,2 > CellMutationStatesCountWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellMutationStatesCountWriter2_2_Overrides : public CellMutationStatesCountWriter2_2{
    public:
    using CellMutationStatesCountWriter2_2::CellMutationStatesCountWriter;
    void WriteHeader(::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            WriteHeader,
                    pCellPopulation);
    }
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesCountWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_CellMutationStatesCountWriter2_2_class(py::module &m){
py::class_<CellMutationStatesCountWriter2_2 , CellMutationStatesCountWriter2_2_Overrides , boost::shared_ptr<CellMutationStatesCountWriter2_2 >  , AbstractCellPopulationCountWriter<2, 2>  >(m, "CellMutationStatesCountWriter2_2")
        .def(py::init< >())
        .def(
            "WriteHeader",
            (void(CellMutationStatesCountWriter2_2::*)(::AbstractCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::WriteHeader,
            " " , py::arg("pCellPopulation") )
        .def(
            "VisitAnyPopulation",
            (void(CellMutationStatesCountWriter2_2::*)(::AbstractCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter2_2::*)(::CaBasedCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellMutationStatesCountWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &CellMutationStatesCountWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
