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

#include "CellMutationStatesWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellMutationStatesWriter<2,2 > CellMutationStatesWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellMutationStatesWriter2_2_Overrides : public CellMutationStatesWriter2_2{
    public:
    using CellMutationStatesWriter2_2::CellMutationStatesWriter;
    double GetCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            double,
            CellMutationStatesWriter2_2,
            GetCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellMutationStatesWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellMutationStatesWriter2_2_class(py::module &m){
py::class_<CellMutationStatesWriter2_2 , CellMutationStatesWriter2_2_Overrides , boost::shared_ptr<CellMutationStatesWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellMutationStatesWriter2_2")
        .def(py::init< >())
        .def(
            "GetCellDataForVtkOutput",
            (double(CellMutationStatesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellMutationStatesWriter2_2::GetCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellMutationStatesWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellMutationStatesWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
