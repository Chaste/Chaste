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

#include "CellAppliedForceWriter3_3.cppwg.hpp"

namespace py = pybind11;
typedef CellAppliedForceWriter<3,3 > CellAppliedForceWriter3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class CellAppliedForceWriter3_3_Overrides : public CellAppliedForceWriter3_3{
    public:
    using CellAppliedForceWriter3_3::CellAppliedForceWriter;
    ::boost::numeric::ublas::c_vector<double, 3> GetVectorCellDataForVtkOutput(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            CellAppliedForceWriter3_3,
            GetVectorCellDataForVtkOutput,
                    pCell,
        pCellPopulation);
    }
    void VisitCell(::CellPtr pCell, ::AbstractCellPopulation<3> * pCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            CellAppliedForceWriter3_3,
            VisitCell,
                    pCell,
        pCellPopulation);
    }

};
void register_CellAppliedForceWriter3_3_class(py::module &m){
py::class_<CellAppliedForceWriter3_3 , CellAppliedForceWriter3_3_Overrides , boost::shared_ptr<CellAppliedForceWriter3_3 >  , AbstractCellWriter<3, 3>  >(m, "CellAppliedForceWriter3_3")
        .def(py::init< >())
        .def(
            "GetVectorCellDataForVtkOutput",
            (::boost::numeric::ublas::c_vector<double, 3>(CellAppliedForceWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellAppliedForceWriter3_3::GetVectorCellDataForVtkOutput,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
        .def(
            "VisitCell",
            (void(CellAppliedForceWriter3_3::*)(::CellPtr, ::AbstractCellPopulation<3> *)) &CellAppliedForceWriter3_3::VisitCell,
            " " , py::arg("pCell"), py::arg("pCellPopulation") )
    ;
}
