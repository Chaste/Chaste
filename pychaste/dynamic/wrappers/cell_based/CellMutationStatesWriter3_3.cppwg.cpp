#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellMutationStatesWriter.hpp"

#include "CellMutationStatesWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellMutationStatesWriter<3,3 > CellMutationStatesWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellMutationStatesWriter3_3_Overrides : public CellMutationStatesWriter3_3{
    public:
    using CellMutationStatesWriter3_3::CellMutationStatesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellMutationStatesWriter3_3,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellMutationStatesWriter3_3_class(py::module &m){
py::class_<CellMutationStatesWriter3_3 , CellMutationStatesWriter3_3_Overrides , boost::shared_ptr<CellMutationStatesWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellMutationStatesWriter3_3")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellMutationStatesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellMutationStatesWriter3_3::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellMutationStatesWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellMutationStatesWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
