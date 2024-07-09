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
#include "CellPopulationAdjacencyMatrixWriter.hpp"

#include "CellPopulationAdjacencyMatrixWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellPopulationAdjacencyMatrixWriter<2,2 > CellPopulationAdjacencyMatrixWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellPopulationAdjacencyMatrixWriter2_2_Overrides : public CellPopulationAdjacencyMatrixWriter2_2{
    public:
    using CellPopulationAdjacencyMatrixWriter2_2::CellPopulationAdjacencyMatrixWriter;
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationAdjacencyMatrixWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationAdjacencyMatrixWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationAdjacencyMatrixWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationAdjacencyMatrixWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationAdjacencyMatrixWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationAdjacencyMatrixWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_CellPopulationAdjacencyMatrixWriter2_2_class(py::module &m){
py::class_<CellPopulationAdjacencyMatrixWriter2_2 , CellPopulationAdjacencyMatrixWriter2_2_Overrides , boost::shared_ptr<CellPopulationAdjacencyMatrixWriter2_2 >  , AbstractCellPopulationWriter<2, 2>  >(m, "CellPopulationAdjacencyMatrixWriter2_2")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::AbstractCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::CaBasedCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationAdjacencyMatrixWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &CellPopulationAdjacencyMatrixWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
