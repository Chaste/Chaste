#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BuskeAdhesiveForce.hpp"

#include "BuskeAdhesiveForce2.cppwg.hpp"

namespace py = pybind11;
typedef BuskeAdhesiveForce<2 > BuskeAdhesiveForce2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;

class BuskeAdhesiveForce2_Overrides : public BuskeAdhesiveForce2{
    public:
    using BuskeAdhesiveForce2::BuskeAdhesiveForce;
    ::boost::numeric::ublas::c_vector<double, 2> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<2> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            BuskeAdhesiveForce2,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeAdhesiveForce2,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_BuskeAdhesiveForce2_class(py::module &m){
py::class_<BuskeAdhesiveForce2 , BuskeAdhesiveForce2_Overrides , boost::shared_ptr<BuskeAdhesiveForce2 >  , AbstractTwoBodyInteractionForce<2>  >(m, "BuskeAdhesiveForce2")
        .def(py::init< >())
        .def(
            "GetAdhesionEnergyParameter",
            (double(BuskeAdhesiveForce2::*)()) &BuskeAdhesiveForce2::GetAdhesionEnergyParameter,
            " "  )
        .def(
            "SetAdhesionEnergyParameter",
            (void(BuskeAdhesiveForce2::*)(double)) &BuskeAdhesiveForce2::SetAdhesionEnergyParameter,
            " " , py::arg("adhesionEnergyParameter") )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 2>(BuskeAdhesiveForce2::*)(unsigned int, unsigned int, ::AbstractCellPopulation<2> &)) &BuskeAdhesiveForce2::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "GetMagnitudeOfForce",
            (double(BuskeAdhesiveForce2::*)(double, double, double)) &BuskeAdhesiveForce2::GetMagnitudeOfForce,
            " " , py::arg("distanceBetweenNodes"), py::arg("radiusOfCellOne"), py::arg("radiusOfCellTwo") )
        .def(
            "OutputForceParameters",
            (void(BuskeAdhesiveForce2::*)(::out_stream &)) &BuskeAdhesiveForce2::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
