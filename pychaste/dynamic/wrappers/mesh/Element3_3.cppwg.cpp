#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Element.hpp"

#include "Element3_3.cppwg.hpp"

namespace py = pybind11;
typedef Element<3,3 > Element3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class Element3_3_Overrides : public Element3_3{
    public:
    using Element3_3::Element;
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE(
            void,
            Element3_3,
            RegisterWithNodes,
            );
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE(
            void,
            Element3_3,
            MarkAsDeleted,
            );
    }
    void UpdateNode(unsigned int const & rIndex, ::Node<3> * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            Element3_3,
            UpdateNode,
                    rIndex,
        pNode);
    }

};
void register_Element3_3_class(py::module &m){
py::class_<Element3_3 , Element3_3_Overrides , boost::shared_ptr<Element3_3 >   >(m, "Element3_3")
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const &, bool >(), py::arg("index"), py::arg("rNodes"), py::arg("registerWithNodes") = true)
        .def(py::init<::Element<3, 3> const &, unsigned int const >(), py::arg("rElement"), py::arg("index"))
        .def(
            "RegisterWithNodes",
            (void(Element3_3::*)()) &Element3_3::RegisterWithNodes,
            " "  )
        .def(
            "MarkAsDeleted",
            (void(Element3_3::*)()) &Element3_3::MarkAsDeleted,
            " "  )
        .def(
            "UpdateNode",
            (void(Element3_3::*)(unsigned int const &, ::Node<3> *)) &Element3_3::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ResetIndex",
            (void(Element3_3::*)(unsigned int)) &Element3_3::ResetIndex,
            " " , py::arg("index") )
        .def(
            "CalculateCircumsphere",
            (::boost::numeric::ublas::c_vector<double, 4>(Element3_3::*)(::boost::numeric::ublas::c_matrix<double, 3, 3> &, ::boost::numeric::ublas::c_matrix<double, 3, 3> &)) &Element3_3::CalculateCircumsphere,
            " " , py::arg("rJacobian"), py::arg("rInverseJacobian") )
        .def(
            "CalculateQuality",
            (double(Element3_3::*)()) &Element3_3::CalculateQuality,
            " "  )
        .def(
            "CalculateMinMaxEdgeLengths",
            (::boost::numeric::ublas::c_vector<double, 2>(Element3_3::*)()) &Element3_3::CalculateMinMaxEdgeLengths,
            " "  )
        .def(
            "CalculateInterpolationWeights",
            (::boost::numeric::ublas::c_vector<double, 4>(Element3_3::*)(::ChastePoint<3> const &)) &Element3_3::CalculateInterpolationWeights,
            " " , py::arg("rTestPoint") )
        .def(
            "CalculateInterpolationWeightsWithProjection",
            (::boost::numeric::ublas::c_vector<double, 4>(Element3_3::*)(::ChastePoint<3> const &)) &Element3_3::CalculateInterpolationWeightsWithProjection,
            " " , py::arg("rTestPoint") )
        .def(
            "CalculateXi",
            (::boost::numeric::ublas::c_vector<double, 3>(Element3_3::*)(::ChastePoint<3> const &)) &Element3_3::CalculateXi,
            " " , py::arg("rTestPoint") )
        .def(
            "IncludesPoint",
            (bool(Element3_3::*)(::ChastePoint<3> const &, bool)) &Element3_3::IncludesPoint,
            " " , py::arg("rTestPoint"), py::arg("strict") = false )
    ;
}
