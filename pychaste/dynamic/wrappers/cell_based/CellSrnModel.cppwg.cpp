#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellSrnModel.hpp"

#include "CellSrnModel.cppwg.hpp"

namespace py = pybind11;
typedef CellSrnModel CellSrnModel;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::AbstractSrnModel * _AbstractSrnModelPtr;

class CellSrnModel_Overrides : public CellSrnModel{
    public:
    using CellSrnModel::CellSrnModel;
    void Initialise() override {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            Initialise,
            );
    }
    void ResetForDivision() override {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            ResetForDivision,
            );
    }
    void SimulateToCurrentTime() override {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            SimulateToCurrentTime,
            );
    }
    ::AbstractSrnModel * CreateSrnModel() override {
        PYBIND11_OVERRIDE(
            _AbstractSrnModelPtr,
            CellSrnModel,
            CreateSrnModel,
            );
    }
    void SetCell(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            void,
            CellSrnModel,
            SetCell,
                    pCell);
    }

};
void register_CellSrnModel_class(py::module &m){
py::class_<CellSrnModel , CellSrnModel_Overrides , boost::shared_ptr<CellSrnModel >  , AbstractSrnModel  >(m, "CellSrnModel")
        .def(py::init< >())
        .def(
            "begin",
            (::CellSrnModel::iterator(CellSrnModel::*)()) &CellSrnModel::begin,
            " "  )
        .def(
            "end",
            (::CellSrnModel::iterator(CellSrnModel::*)()) &CellSrnModel::end,
            " "  )
        .def(
            "begin",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const ) &CellSrnModel::begin,
            " "  )
        .def(
            "end",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const ) &CellSrnModel::end,
            " "  )
        .def(
            "cbegin",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const ) &CellSrnModel::cbegin,
            " "  )
        .def(
            "cend",
            (::CellSrnModel::const_iterator(CellSrnModel::*)() const ) &CellSrnModel::cend,
            " "  )
        .def(
            "Initialise",
            (void(CellSrnModel::*)()) &CellSrnModel::Initialise,
            " "  )
        .def(
            "ResetForDivision",
            (void(CellSrnModel::*)()) &CellSrnModel::ResetForDivision,
            " "  )
        .def(
            "SimulateToCurrentTime",
            (void(CellSrnModel::*)()) &CellSrnModel::SimulateToCurrentTime,
            " "  )
        .def(
            "CreateSrnModel",
            (::AbstractSrnModel *(CellSrnModel::*)()) &CellSrnModel::CreateSrnModel,
            " "  , py::return_value_policy::reference)
        .def(
            "AddEdgeSrn",
            (void(CellSrnModel::*)(::std::vector<boost::shared_ptr<AbstractSrnModel>>)) &CellSrnModel::AddEdgeSrn,
            " " , py::arg("edgeSrns") )
        .def(
            "AddEdgeSrnModel",
            (void(CellSrnModel::*)(::AbstractSrnModelPtr)) &CellSrnModel::AddEdgeSrnModel,
            " " , py::arg("pEdgeSrn") )
        .def(
            "GetNumEdgeSrn",
            (unsigned int(CellSrnModel::*)() const ) &CellSrnModel::GetNumEdgeSrn,
            " "  )
        .def(
            "GetEdgeSrn",
            (::AbstractSrnModelPtr(CellSrnModel::*)(unsigned int) const ) &CellSrnModel::GetEdgeSrn,
            " " , py::arg("index") )
        .def(
            "GetEdges",
            (::std::vector<boost::shared_ptr<AbstractSrnModel>> const &(CellSrnModel::*)() const ) &CellSrnModel::GetEdges,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "SetInteriorSrnModel",
            (void(CellSrnModel::*)(::AbstractSrnModelPtr)) &CellSrnModel::SetInteriorSrnModel,
            " " , py::arg("pInteriorSrn") )
        .def(
            "GetInteriorSrn",
            (::AbstractSrnModelPtr(CellSrnModel::*)() const ) &CellSrnModel::GetInteriorSrn,
            " "  )
        .def(
            "SetCell",
            (void(CellSrnModel::*)(::CellPtr)) &CellSrnModel::SetCell,
            " " , py::arg("pCell") )
    ;
}
