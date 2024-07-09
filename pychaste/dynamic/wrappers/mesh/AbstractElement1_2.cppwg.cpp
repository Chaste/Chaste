#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractElement.hpp"

#include "AbstractElement1_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractElement<1,2 > AbstractElement1_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractElement1_2_Overrides : public AbstractElement1_2{
    public:
    using AbstractElement1_2::AbstractElement;
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement1_2,
            UpdateNode,
                    rIndex,
        pNode);
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement1_2,
            MarkAsDeleted,
            );
    }
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement1_2,
            RegisterWithNodes,
            );
    }

};
void register_AbstractElement1_2_class(py::module &m){
py::class_<AbstractElement1_2 , AbstractElement1_2_Overrides , boost::shared_ptr<AbstractElement1_2 >   >(m, "AbstractElement1_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(py::init<unsigned int >(), py::arg("index") = ::INDEX_IS_NOT_USED)
        .def(
            "UpdateNode",
            (void(AbstractElement1_2::*)(unsigned int const &, ::Node<2> *)) &AbstractElement1_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ReplaceNode",
            (void(AbstractElement1_2::*)(::Node<2> *, ::Node<2> *)) &AbstractElement1_2::ReplaceNode,
            " " , py::arg("pOldNode"), py::arg("pNewNode") )
        .def(
            "MarkAsDeleted",
            (void(AbstractElement1_2::*)()) &AbstractElement1_2::MarkAsDeleted,
            " "  )
        .def(
            "RegisterWithNodes",
            (void(AbstractElement1_2::*)()) &AbstractElement1_2::RegisterWithNodes,
            " "  )
        .def(
            "GetNodeLocation",
            (double(AbstractElement1_2::*)(unsigned int, unsigned int) const ) &AbstractElement1_2::GetNodeLocation,
            " " , py::arg("localIndex"), py::arg("dimension") )
        .def(
            "GetNodeLocation",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractElement1_2::*)(unsigned int) const ) &AbstractElement1_2::GetNodeLocation,
            " " , py::arg("localIndex") )
        .def(
            "GetNodeGlobalIndex",
            (unsigned int(AbstractElement1_2::*)(unsigned int) const ) &AbstractElement1_2::GetNodeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetNode",
            (::Node<2> *(AbstractElement1_2::*)(unsigned int) const ) &AbstractElement1_2::GetNode,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(AbstractElement1_2::*)() const ) &AbstractElement1_2::GetNumNodes,
            " "  )
        .def(
            "AddNode",
            (void(AbstractElement1_2::*)(::Node<2> *)) &AbstractElement1_2::AddNode,
            " " , py::arg("pNode") )
        .def(
            "IsDeleted",
            (bool(AbstractElement1_2::*)() const ) &AbstractElement1_2::IsDeleted,
            " "  )
        .def(
            "GetIndex",
            (unsigned int(AbstractElement1_2::*)() const ) &AbstractElement1_2::GetIndex,
            " "  )
        .def(
            "SetIndex",
            (void(AbstractElement1_2::*)(unsigned int)) &AbstractElement1_2::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetOwnership",
            (bool(AbstractElement1_2::*)() const ) &AbstractElement1_2::GetOwnership,
            " "  )
        .def(
            "SetOwnership",
            (void(AbstractElement1_2::*)(bool)) &AbstractElement1_2::SetOwnership,
            " " , py::arg("ownership") )
        .def(
            "SetAttribute",
            (void(AbstractElement1_2::*)(double)) &AbstractElement1_2::SetAttribute,
            " " , py::arg("attribute") )
        .def(
            "GetAttribute",
            (double(AbstractElement1_2::*)()) &AbstractElement1_2::GetAttribute,
            " "  )
        .def(
            "GetUnsignedAttribute",
            (unsigned int(AbstractElement1_2::*)()) &AbstractElement1_2::GetUnsignedAttribute,
            " "  )
        .def(
            "AddElementAttribute",
            (void(AbstractElement1_2::*)(double)) &AbstractElement1_2::AddElementAttribute,
            " " , py::arg("attribute") )
        .def(
            "rGetElementAttributes",
            (::std::vector<double> &(AbstractElement1_2::*)()) &AbstractElement1_2::rGetElementAttributes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumElementAttributes",
            (unsigned int(AbstractElement1_2::*)()) &AbstractElement1_2::GetNumElementAttributes,
            " "  )
    ;
}
