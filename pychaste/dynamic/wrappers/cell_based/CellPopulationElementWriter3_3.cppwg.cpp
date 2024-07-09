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

#include "CellPopulationElementWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellPopulationElementWriter<3,3 > CellPopulationElementWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellPopulationElementWriter3_3_Overrides : public CellPopulationElementWriter3_3{
    public:
    using CellPopulationElementWriter3_3::CellPopulationElementWriter;
    void Visit(::MeshBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter3_3,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellPopulationElementWriter3_3,
            Visit,
                    pCellPopulation);
    }

};
void register_CellPopulationElementWriter3_3_class(py::module &m){
py::class_<CellPopulationElementWriter3_3 , CellPopulationElementWriter3_3_Overrides , boost::shared_ptr<CellPopulationElementWriter3_3 >  , AbstractCellPopulationWriter<3, 3>  >(m, "CellPopulationElementWriter3_3")
        .def(py::init< >())
        .def(
            "Visit",
            (void(CellPopulationElementWriter3_3::*)(::MeshBasedCellPopulation<3> *)) &CellPopulationElementWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter3_3::*)(::CaBasedCellPopulation<3> *)) &CellPopulationElementWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter3_3::*)(::NodeBasedCellPopulation<3> *)) &CellPopulationElementWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter3_3::*)(::PottsBasedCellPopulation<3> *)) &CellPopulationElementWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter3_3::*)(::VertexBasedCellPopulation<3> *)) &CellPopulationElementWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellPopulationElementWriter3_3::*)(::ImmersedBoundaryCellPopulation<3> *)) &CellPopulationElementWriter3_3::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
