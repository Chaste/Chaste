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
#include "CellPopulationElementWriter.hpp"

#include "CellPopulationElementWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellPopulationElementWriter<2,2 > CellPopulationElementWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellPopulationElementWriter2_2_Overrides : public CellPopulationElementWriter2_2{
    public:
    using CellPopulationElementWriter2_2::CellPopulationElementWriter;
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_CellPopulationElementWriter2_2_class(py::module &m){
py::class_<CellPopulationElementWriter2_2 , CellPopulationElementWriter2_2_Overrides , boost::shared_ptr<CellPopulationElementWriter2_2 >  , AbstractCellPopulationWriter<2, 2>  >(m, "CellPopulationElementWriter2_2")
        .def(py::init< >())
        .def(
            "Visit",
            (void(CellPopulationElementWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &CellPopulationElementWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter2_2::*)(::CaBasedCellPopulation<2> *)) &CellPopulationElementWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &CellPopulationElementWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &CellPopulationElementWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &CellPopulationElementWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &CellPopulationElementWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
