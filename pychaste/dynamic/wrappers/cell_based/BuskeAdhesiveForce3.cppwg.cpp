#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BuskeAdhesiveForce.hpp"

#include "BuskeAdhesiveForce3.cppwg.hpp"

namespace py = pybind11;
typedef BuskeAdhesiveForce<3 > BuskeAdhesiveForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class BuskeAdhesiveForce3_Overrides : public BuskeAdhesiveForce3{
    public:
    using BuskeAdhesiveForce3::BuskeAdhesiveForce;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            BuskeAdhesiveForce3,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeAdhesiveForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_BuskeAdhesiveForce3_class(py::module &m){
py::class_<BuskeAdhesiveForce3 , BuskeAdhesiveForce3_Overrides , boost::shared_ptr<BuskeAdhesiveForce3 >  , AbstractTwoBodyInteractionForce<3>  >(m, "BuskeAdhesiveForce3")
        .def(py::init< >())
        .def(
            "GetAdhesionEnergyParameter",
            (double(BuskeAdhesiveForce3::*)()) &BuskeAdhesiveForce3::GetAdhesionEnergyParameter,
            " "  )
        .def(
            "SetAdhesionEnergyParameter",
            (void(BuskeAdhesiveForce3::*)(double)) &BuskeAdhesiveForce3::SetAdhesionEnergyParameter,
            " " , py::arg("adhesionEnergyParameter") )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 3>(BuskeAdhesiveForce3::*)(unsigned int, unsigned int, ::AbstractCellPopulation<3> &)) &BuskeAdhesiveForce3::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "GetMagnitudeOfForce",
            (double(BuskeAdhesiveForce3::*)(double, double, double)) &BuskeAdhesiveForce3::GetMagnitudeOfForce,
            " " , py::arg("distanceBetweenNodes"), py::arg("radiusOfCellOne"), py::arg("radiusOfCellTwo") )
        .def(
            "OutputForceParameters",
            (void(BuskeAdhesiveForce3::*)(::out_stream &)) &BuskeAdhesiveForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
