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
#include "CellRemovalLocationsWriter.hpp"

#include "CellRemovalLocationsWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellRemovalLocationsWriter<2,2 > CellRemovalLocationsWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellRemovalLocationsWriter2_2_Overrides : public CellRemovalLocationsWriter2_2{
    public:
    using CellRemovalLocationsWriter2_2::CellRemovalLocationsWriter;
    void Visit(::MeshBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRemovalLocationsWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::CaBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRemovalLocationsWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::NodeBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRemovalLocationsWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::PottsBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRemovalLocationsWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::VertexBasedCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRemovalLocationsWriter2_2,
            Visit,
                    pCellPopulation);
    }
    void Visit(::ImmersedBoundaryCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellRemovalLocationsWriter2_2,
            Visit,
                    pCellPopulation);
    }

};
void register_CellRemovalLocationsWriter2_2_class(py::module &m){
py::class_<CellRemovalLocationsWriter2_2 , CellRemovalLocationsWriter2_2_Overrides , boost::shared_ptr<CellRemovalLocationsWriter2_2 >  , AbstractCellPopulationEventWriter<2, 2>  >(m, "CellRemovalLocationsWriter2_2")
        .def(py::init< >())
        .def(
            "VisitAnyPopulation",
            (void(CellRemovalLocationsWriter2_2::*)(::AbstractCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::VisitAnyPopulation,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellRemovalLocationsWriter2_2::*)(::MeshBasedCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellRemovalLocationsWriter2_2::*)(::CaBasedCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellRemovalLocationsWriter2_2::*)(::NodeBasedCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellRemovalLocationsWriter2_2::*)(::PottsBasedCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellRemovalLocationsWriter2_2::*)(::VertexBasedCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
        .def(
            "Visit",
            (void(CellRemovalLocationsWriter2_2::*)(::ImmersedBoundaryCellPopulation<2> *)) &CellRemovalLocationsWriter2_2::Visit,
            " " , py::arg("pCellPopulation") )
    ;
}
