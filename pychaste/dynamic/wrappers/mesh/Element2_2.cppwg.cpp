#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "Element.hpp"

#include "Element2_2.cppwg.hpp"

namespace py = pybind11;
typedef Element<2,2 > Element2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class Element2_2_Overrides : public Element2_2{
    public:
    using Element2_2::Element;
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE(
            void,
            Element2_2,
            RegisterWithNodes,
            );
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE(
            void,
            Element2_2,
            MarkAsDeleted,
            );
    }
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE(
            void,
            Element2_2,
            UpdateNode,
                    rIndex,
        pNode);
    }

};
void register_Element2_2_class(py::module &m){
py::class_<Element2_2 , Element2_2_Overrides , boost::shared_ptr<Element2_2 >   >(m, "Element2_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const &, bool >(), py::arg("index"), py::arg("rNodes"), py::arg("registerWithNodes") = true)
        .def(py::init<::Element<2, 2> const &, unsigned int const >(), py::arg("rElement"), py::arg("index"))
        .def(
            "RegisterWithNodes",
            (void(Element2_2::*)()) &Element2_2::RegisterWithNodes,
            " "  )
        .def(
            "MarkAsDeleted",
            (void(Element2_2::*)()) &Element2_2::MarkAsDeleted,
            " "  )
        .def(
            "UpdateNode",
            (void(Element2_2::*)(unsigned int const &, ::Node<2> *)) &Element2_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ResetIndex",
            (void(Element2_2::*)(unsigned int)) &Element2_2::ResetIndex,
            " " , py::arg("index") )
        .def(
            "CalculateCircumsphere",
            (::boost::numeric::ublas::c_vector<double, 3>(Element2_2::*)(::boost::numeric::ublas::c_matrix<double, 2, 2> &, ::boost::numeric::ublas::c_matrix<double, 2, 2> &)) &Element2_2::CalculateCircumsphere,
            " " , py::arg("rJacobian"), py::arg("rInverseJacobian") )
        .def(
            "CalculateQuality",
            (double(Element2_2::*)()) &Element2_2::CalculateQuality,
            " "  )
        .def(
            "CalculateMinMaxEdgeLengths",
            (::boost::numeric::ublas::c_vector<double, 2>(Element2_2::*)()) &Element2_2::CalculateMinMaxEdgeLengths,
            " "  )
        .def(
            "CalculateInterpolationWeights",
            (::boost::numeric::ublas::c_vector<double, 3>(Element2_2::*)(::ChastePoint<2> const &)) &Element2_2::CalculateInterpolationWeights,
            " " , py::arg("rTestPoint") )
        .def(
            "CalculateInterpolationWeightsWithProjection",
            (::boost::numeric::ublas::c_vector<double, 3>(Element2_2::*)(::ChastePoint<2> const &)) &Element2_2::CalculateInterpolationWeightsWithProjection,
            " " , py::arg("rTestPoint") )
        .def(
            "CalculateXi",
            (::boost::numeric::ublas::c_vector<double, 2>(Element2_2::*)(::ChastePoint<2> const &)) &Element2_2::CalculateXi,
            " " , py::arg("rTestPoint") )
        .def(
            "IncludesPoint",
            (bool(Element2_2::*)(::ChastePoint<2> const &, bool)) &Element2_2::IncludesPoint,
            " " , py::arg("rTestPoint"), py::arg("strict") = false )
    ;
}
