#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "BuskeElasticForce.hpp"

#include "BuskeElasticForce3.cppwg.hpp"

namespace py = pybind11;
typedef BuskeElasticForce<3 > BuskeElasticForce3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;

class BuskeElasticForce3_Overrides : public BuskeElasticForce3{
    public:
    using BuskeElasticForce3::BuskeElasticForce;
    ::boost::numeric::ublas::c_vector<double, 3> CalculateForceBetweenNodes(unsigned int nodeAGlobalIndex, unsigned int nodeBGlobalIndex, ::AbstractCellPopulation<3> & rCellPopulation) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            BuskeElasticForce3,
            CalculateForceBetweenNodes,
                    nodeAGlobalIndex,
        nodeBGlobalIndex,
        rCellPopulation);
    }
    void OutputForceParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            BuskeElasticForce3,
            OutputForceParameters,
                    rParamsFile);
    }

};
void register_BuskeElasticForce3_class(py::module &m){
py::class_<BuskeElasticForce3 , BuskeElasticForce3_Overrides , boost::shared_ptr<BuskeElasticForce3 >  , AbstractTwoBodyInteractionForce<3>  >(m, "BuskeElasticForce3")
        .def(py::init< >())
        .def(
            "GetDeformationEnergyParameter",
            (double(BuskeElasticForce3::*)()) &BuskeElasticForce3::GetDeformationEnergyParameter,
            " "  )
        .def(
            "SetDeformationEnergyParameter",
            (void(BuskeElasticForce3::*)(double)) &BuskeElasticForce3::SetDeformationEnergyParameter,
            " " , py::arg("deformationEnergyParameter") )
        .def(
            "CalculateForceBetweenNodes",
            (::boost::numeric::ublas::c_vector<double, 3>(BuskeElasticForce3::*)(unsigned int, unsigned int, ::AbstractCellPopulation<3> &)) &BuskeElasticForce3::CalculateForceBetweenNodes,
            " " , py::arg("nodeAGlobalIndex"), py::arg("nodeBGlobalIndex"), py::arg("rCellPopulation") )
        .def(
            "GetMagnitudeOfForce",
            (double(BuskeElasticForce3::*)(double, double, double)) &BuskeElasticForce3::GetMagnitudeOfForce,
            " " , py::arg("distanceBetweenNodes"), py::arg("radiusOfCellOne"), py::arg("radiusOfCellTwo") )
        .def(
            "OutputForceParameters",
            (void(BuskeElasticForce3::*)(::out_stream &)) &BuskeElasticForce3::OutputForceParameters,
            " " , py::arg("rParamsFile") )
    ;
}
