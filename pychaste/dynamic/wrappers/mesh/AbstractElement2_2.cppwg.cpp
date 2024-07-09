#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractElement.hpp"

#include "AbstractElement2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractElement<2,2 > AbstractElement2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class AbstractElement2_2_Overrides : public AbstractElement2_2{
    public:
    using AbstractElement2_2::AbstractElement;
    void UpdateNode(unsigned int const & rIndex, ::Node<2> * pNode) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement2_2,
            UpdateNode,
                    rIndex,
        pNode);
    }
    void MarkAsDeleted() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement2_2,
            MarkAsDeleted,
            );
    }
    void RegisterWithNodes() override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractElement2_2,
            RegisterWithNodes,
            );
    }

};
void register_AbstractElement2_2_class(py::module &m){
py::class_<AbstractElement2_2 , AbstractElement2_2_Overrides , boost::shared_ptr<AbstractElement2_2 >   >(m, "AbstractElement2_2")
        .def(py::init<unsigned int, ::std::vector<Node<2> *> const & >(), py::arg("index"), py::arg("rNodes"))
        .def(py::init<unsigned int >(), py::arg("index") = ::INDEX_IS_NOT_USED)
        .def(
            "UpdateNode",
            (void(AbstractElement2_2::*)(unsigned int const &, ::Node<2> *)) &AbstractElement2_2::UpdateNode,
            " " , py::arg("rIndex"), py::arg("pNode") )
        .def(
            "ReplaceNode",
            (void(AbstractElement2_2::*)(::Node<2> *, ::Node<2> *)) &AbstractElement2_2::ReplaceNode,
            " " , py::arg("pOldNode"), py::arg("pNewNode") )
        .def(
            "MarkAsDeleted",
            (void(AbstractElement2_2::*)()) &AbstractElement2_2::MarkAsDeleted,
            " "  )
        .def(
            "RegisterWithNodes",
            (void(AbstractElement2_2::*)()) &AbstractElement2_2::RegisterWithNodes,
            " "  )
        .def(
            "GetNodeLocation",
            (double(AbstractElement2_2::*)(unsigned int, unsigned int) const ) &AbstractElement2_2::GetNodeLocation,
            " " , py::arg("localIndex"), py::arg("dimension") )
        .def(
            "GetNodeLocation",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractElement2_2::*)(unsigned int) const ) &AbstractElement2_2::GetNodeLocation,
            " " , py::arg("localIndex") )
        .def(
            "GetNodeGlobalIndex",
            (unsigned int(AbstractElement2_2::*)(unsigned int) const ) &AbstractElement2_2::GetNodeGlobalIndex,
            " " , py::arg("localIndex") )
        .def(
            "GetNode",
            (::Node<2> *(AbstractElement2_2::*)(unsigned int) const ) &AbstractElement2_2::GetNode,
            " " , py::arg("localIndex") , py::return_value_policy::reference)
        .def(
            "GetNumNodes",
            (unsigned int(AbstractElement2_2::*)() const ) &AbstractElement2_2::GetNumNodes,
            " "  )
        .def(
            "AddNode",
            (void(AbstractElement2_2::*)(::Node<2> *)) &AbstractElement2_2::AddNode,
            " " , py::arg("pNode") )
        .def(
            "IsDeleted",
            (bool(AbstractElement2_2::*)() const ) &AbstractElement2_2::IsDeleted,
            " "  )
        .def(
            "GetIndex",
            (unsigned int(AbstractElement2_2::*)() const ) &AbstractElement2_2::GetIndex,
            " "  )
        .def(
            "SetIndex",
            (void(AbstractElement2_2::*)(unsigned int)) &AbstractElement2_2::SetIndex,
            " " , py::arg("index") )
        .def(
            "GetOwnership",
            (bool(AbstractElement2_2::*)() const ) &AbstractElement2_2::GetOwnership,
            " "  )
        .def(
            "SetOwnership",
            (void(AbstractElement2_2::*)(bool)) &AbstractElement2_2::SetOwnership,
            " " , py::arg("ownership") )
        .def(
            "SetAttribute",
            (void(AbstractElement2_2::*)(double)) &AbstractElement2_2::SetAttribute,
            " " , py::arg("attribute") )
        .def(
            "GetAttribute",
            (double(AbstractElement2_2::*)()) &AbstractElement2_2::GetAttribute,
            " "  )
        .def(
            "GetUnsignedAttribute",
            (unsigned int(AbstractElement2_2::*)()) &AbstractElement2_2::GetUnsignedAttribute,
            " "  )
        .def(
            "AddElementAttribute",
            (void(AbstractElement2_2::*)(double)) &AbstractElement2_2::AddElementAttribute,
            " " , py::arg("attribute") )
        .def(
            "rGetElementAttributes",
            (::std::vector<double> &(AbstractElement2_2::*)()) &AbstractElement2_2::rGetElementAttributes,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetNumElementAttributes",
            (unsigned int(AbstractElement2_2::*)()) &AbstractElement2_2::GetNumElementAttributes,
            " "  )
    ;
}
