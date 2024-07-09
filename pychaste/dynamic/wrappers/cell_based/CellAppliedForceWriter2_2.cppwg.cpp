#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "AbstractCellPopulation.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellAppliedForceWriter.hpp"

#include "CellAppliedForceWriter2_2.cppwg.hpp"

namespace py = pybind11;
typedef CellAppliedForceWriter<2,2 > CellAppliedForceWriter2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class CellAppliedForceWriter2_2_Overrides : public CellAppliedForceWriter2_2{
    public:
    using CellAppliedForceWriter2_2::CellAppliedForceWriter;
    ::boost::numeric::ublas::c_vector<double, 2> GetVectorCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            CellAppliedForceWriter2_2,
            GetVectorCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<2> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellAppliedForceWriter2_2,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellAppliedForceWriter2_2_class(py::module &m){
py::class_<CellAppliedForceWriter2_2 , CellAppliedForceWriter2_2_Overrides , boost::shared_ptr<CellAppliedForceWriter2_2 >  , AbstractCellWriter<2, 2>  >(m, "CellAppliedForceWriter2_2")
        .def(py::init< >())
        .def(
            "GetVectorCellDataForVtkOutput",
            (::boost::numeric::ublas::c_vector<double, 2>(CellAppliedForceWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellAppliedForceWriter2_2::GetVectorCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellAppliedForceWriter2_2::*)(::CellPtr, ::AbstractCellPopulation<2> *)) &CellAppliedForceWriter2_2::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
