#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractElement.hpp"

#include "AbstractElement3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractElement<3,3 > AbstractElement3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractElement3_3_Overrides : public AbstractElement3_3{
    public:
    using AbstractElement3_3::AbstractElement;
    void UpdateNode(unsigned int const & rIndex, ::Node<3> * pNode) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement3_3,
            UpdateNode,
                    rIndex,
        pNode);
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement3_3,
            MarkAsDeleted,
            );
    }
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement3_3,
            RegisterWithNodes,
            );
    }

};
void register_AbstractElement3_3_class(py::module &m){
py::class_<AbstractElement3_3 , AbstractElement3_3_Overrides , boost::shared_ptr<AbstractElement3_3 >   >(m, "AbstractElement3_3")
        .def(py::init<unsigned int, ::std::vector<Node<3> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(py::init<unsigned int >(), py::arg("index") = ::INDEX_IS_NOT_USED)
        .def(
            "UpdateNode",
            (void(AbstractElement3_3::*)(unsigned int const &, ::Node<3> *)) &AbstractElement3_3::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ReplaceNode",
            (void(AbstractElement3_3::*)(::Node<3> *, ::Node<3> *)) &AbstractElement3_3::ReplaceNode,
            " " , py::arg("pOldNode"), py::arg("pNewNode") )
        .def(
            "MarkAsDeleted",
            (void(AbstractElement3_3::*)()) &AbstractElement3_3::MarkAsDeleted,
            " "  )
        .def(
            "RegisterWithNodes",
            (void(AbstractElement3_3::*)()) &AbstractElement3_3::RegisterWithNodes,
            " "  )
        .def(
            "GetNodeLocation",
            (double(AbstractElement3_3::*)(unsigned int, unsigned int) const ) &AbstractElement3_3::GetNodeLocation,
            " " , py::arg("localIndex"), py::arg("dimension") )
        .def(
            "GetNodeLocation",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractElement3_3::*)(unsigned int) const ) &AbstractElement3_3::GetNodeLocation,
            " " , py::arg("localIndex") )
        .def(
            "GetNodeGlobalIndex",
            (unsigned int(AbstractElement3_3::*)(unsigned int) const ) &AbstractElement3_3::GetNodeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetNode",
            (::Node<3> *(AbstractElement3_3::*)(unsigned int) const ) &AbstractElement3_3::GetNode,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(AbstractElement3_3::*)() const ) &AbstractElement3_3::GetNumNodes,
            " "  )
        .def(
            "AddNode",
            (void(AbstractElement3_3::*)(::Node<3> *)) &AbstractElement3_3::AddNode,
            " " , py::arg("pNode") )
        .def(
            "IsDeleted",
            (bool(AbstractElement3_3::*)() const ) &AbstractElement3_3::IsDeleted,
            " "  )
        .def(
            "GetIndex",
            (unsigned int(AbstractElement3_3::*)() const ) &AbstractElement3_3::GetIndex,
            " "  )
        .def(
            "SetIndex",
            (void(AbstractElement3_3::*)(unsigned int)) &AbstractElement3_3::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetOwnership",
            (bool(AbstractElement3_3::*)() const ) &AbstractElement3_3::GetOwnership,
            " "  )
        .def(
            "SetOwnership",
            (void(AbstractElement3_3::*)(bool)) &AbstractElement3_3::SetOwnership,
            " " , py::arg("ownership") )
        .def(
            "SetAttribute",
            (void(AbstractElement3_3::*)(double)) &AbstractElement3_3::SetAttribute,
            " " , py::arg("attribute") )
        .def(
            "GetAttribute",
            (double(AbstractElement3_3::*)()) &AbstractElement3_3::GetAttribute,
            " "  )
        .def(
            "GetUnsignedAttribute",
            (unsigned int(AbstractElement3_3::*)()) &AbstractElement3_3::GetUnsignedAttribute,
            " "  )
        .def(
            "AddElementAttribute",
            (void(AbstractElement3_3::*)(double)) &AbstractElement3_3::AddElementAttribute,
            " " , py::arg("attribute") )
        .def(
            "rGetElementAttributes",
            (::std::vector<double> &(AbstractElement3_3::*)()) &AbstractElement3_3::rGetElementAttributes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumElementAttributes",
            (unsigned int(AbstractElement3_3::*)()) &AbstractElement3_3::GetNumElementAttributes,
            " "  )
    ;
}
