#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "ChemotacticForce.hpp"

#include "ChemotacticForce3.cppwg.hpp"

namespace py = pybind11;
typedef ChemotacticForce<3 > ChemotacticForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class ChemotacticForce3_Overrides : public ChemotacticForce3{
    public:
    using ChemotacticForce3::ChemotacticForce;
    void AddForceContribution(::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            void,
            ChemotacticForce3,
            AddForceContribution,
                    rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            ChemotacticForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_ChemotacticForce3_class(py::module &m){
py::class_<ChemotacticForce3 , ChemotacticForce3_Overrides , boost::shared_ptr<ChemotacticForce3 >  , AbstractForce<3>  >(m, "ChemotacticForce3")
        .def(py::init< >())
        .def(
            "AddForceContribution",
            (void(ChemotacticForce3::*)(::AbstractCellPopulation<3> &)) &ChemotacticForce3::AddForceContribution,
            " " , py::arg("rCellPopulation") )
        .def(
            "OutputForceParameters",
            (void(ChemotacticForce3::*)(::out_stream &)) &ChemotacticForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
