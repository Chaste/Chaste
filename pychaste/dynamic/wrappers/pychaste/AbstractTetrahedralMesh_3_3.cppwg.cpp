/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

// This file is auto-generated. Manual changes will be overwritten. For changes
// to persist, update the configuration in pychaste/dynamic/config.yaml.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractTetrahedralMesh.hpp"

#include "AbstractTetrahedralMesh_3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractTetrahedralMesh<3, 3> AbstractTetrahedralMesh_3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef unsigned int unsignedint;
typedef unsigned int unsignedint;

class AbstractTetrahedralMesh_3_3_Overrides : public AbstractTetrahedralMesh_3_3
{
public:
    using AbstractTetrahedralMesh_3_3::AbstractTetrahedralMesh;
    unsigned int GetNumElements() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetNumElements,
            );
    }
    unsigned int GetNumLocalElements() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetNumLocalElements,
            );
    }
    unsigned int GetNumBoundaryElements() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetNumBoundaryElements,
            );
    }
    unsigned int GetNumLocalBoundaryElements() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetNumLocalBoundaryElements,
            );
    }
    unsigned int GetNumCableElements() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetNumCableElements,
            );
    }
    unsigned int GetNumVertices() const override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetNumVertices,
            );
    }
    unsigned int GetMaximumNodeIndex() override
    {
        PYBIND11_OVERRIDE(
            unsignedint,
            AbstractTetrahedralMesh_3_3,
            GetMaximumNodeIndex,
            );
    }
    void ConstructFromMeshReader(::AbstractMeshReader<3, 3> & rMeshReader) override
    {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractTetrahedralMesh_3_3,
            ConstructFromMeshReader,
            rMeshReader);
    }
    void GetInverseJacobianForElement(unsigned int elementIndex, ::boost::numeric::ublas::c_matrix<double, 3, 3> & rJacobian, double & rJacobianDeterminant, ::boost::numeric::ublas::c_matrix<double, 3, 3> & rInverseJacobian) const override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            GetInverseJacobianForElement,
            elementIndex,
            rJacobian,
            rJacobianDeterminant,
            rInverseJacobian);
    }
    void GetWeightedDirectionForBoundaryElement(unsigned int elementIndex, ::boost::numeric::ublas::c_vector<double, 3> & rWeightedDirection, double & rJacobianDeterminant) const override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            GetWeightedDirectionForBoundaryElement,
            elementIndex,
            rWeightedDirection,
            rJacobianDeterminant);
    }
    void ConstructLinearMesh(unsigned int width) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            ConstructLinearMesh,
            width);
    }
    void ConstructRectangularMesh(unsigned int width, unsigned int height, bool stagger) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            ConstructRectangularMesh,
            width,
            height,
            stagger);
    }
    void ConstructCuboid(unsigned int width, unsigned int height, unsigned int depth) override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            ConstructCuboid,
            width,
            height,
            depth);
    }
    bool CalculateDesignatedOwnershipOfBoundaryElement(unsigned int faceIndex) override
    {
        PYBIND11_OVERRIDE(
            bool,
            AbstractTetrahedralMesh_3_3,
            CalculateDesignatedOwnershipOfBoundaryElement,
            faceIndex);
    }
    bool CalculateDesignatedOwnershipOfElement(unsigned int elementIndex) override
    {
        PYBIND11_OVERRIDE(
            bool,
            AbstractTetrahedralMesh_3_3,
            CalculateDesignatedOwnershipOfElement,
            elementIndex);
    }
    void GetHaloNodeIndices(::std::vector<unsigned int> & rHaloIndices) const override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            GetHaloNodeIndices,
            rHaloIndices);
    }
    ::boost::numeric::ublas::c_vector<double, 2> CalculateMinMaxEdgeLengths() override
    {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractTetrahedralMesh_3_3,
            CalculateMinMaxEdgeLengths,
            );
    }
    void SetElementOwnerships() override
    {
        PYBIND11_OVERRIDE(
            void,
            AbstractTetrahedralMesh_3_3,
            SetElementOwnerships,
            );
    }
};

void register_AbstractTetrahedralMesh_3_3_class(py::module &m)
{
    py::class_<AbstractTetrahedralMesh_3_3, AbstractTetrahedralMesh_3_3_Overrides, boost::shared_ptr<AbstractTetrahedralMesh_3_3>, AbstractMesh<3, 3>>(m, "AbstractTetrahedralMesh_3_3")
        .def("GetElementIteratorBegin",
            (::AbstractTetrahedralMesh<3, 3>::ElementIterator(AbstractTetrahedralMesh_3_3::*)(bool)) &AbstractTetrahedralMesh_3_3::GetElementIteratorBegin,
            " ", py::arg("skipDeletedElements") = true)
        .def("GetElementIteratorEnd",
            (::AbstractTetrahedralMesh<3, 3>::ElementIterator(AbstractTetrahedralMesh_3_3::*)()) &AbstractTetrahedralMesh_3_3::GetElementIteratorEnd,
            " ")
        .def("GetNumElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumElements,
            " ")
        .def("GetNumLocalElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumLocalElements,
            " ")
        .def("GetNumBoundaryElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumBoundaryElements,
            " ")
        .def("GetNumLocalBoundaryElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumLocalBoundaryElements,
            " ")
        .def("GetNumAllElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumAllElements,
            " ")
        .def("GetNumAllBoundaryElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumAllBoundaryElements,
            " ")
        .def("GetNumCableElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumCableElements,
            " ")
        .def("GetNumVertices",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetNumVertices,
            " ")
        .def("GetMaximumNodeIndex",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)()) &AbstractTetrahedralMesh_3_3::GetMaximumNodeIndex,
            " ")
        .def("GetElement",
            (::Element<3, 3> *(AbstractTetrahedralMesh_3_3::*)(unsigned int) const) &AbstractTetrahedralMesh_3_3::GetElement,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("GetBoundaryElement",
            (::BoundaryElement<2, 3> *(AbstractTetrahedralMesh_3_3::*)(unsigned int) const) &AbstractTetrahedralMesh_3_3::GetBoundaryElement,
            " ", py::arg("index"), py::return_value_policy::reference)
        .def("ConstructFromMeshReader",
            (void(AbstractTetrahedralMesh_3_3::*)(::AbstractMeshReader<3, 3> &)) &AbstractTetrahedralMesh_3_3::ConstructFromMeshReader,
            " ", py::arg("rMeshReader"))
        .def("ConstructFromMesh",
            (void(AbstractTetrahedralMesh_3_3::*)(::AbstractTetrahedralMesh<3, 3> &)) &AbstractTetrahedralMesh_3_3::ConstructFromMesh,
            " ", py::arg("rOtherMesh"))
        .def("GetBoundaryElementIteratorBegin",
            (::AbstractTetrahedralMesh<3, 3>::BoundaryElementIterator(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetBoundaryElementIteratorBegin,
            " ")
        .def("GetBoundaryElementIteratorEnd",
            (::AbstractTetrahedralMesh<3, 3>::BoundaryElementIterator(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::GetBoundaryElementIteratorEnd,
            " ")
        .def("GetInverseJacobianForElement",
            (void(AbstractTetrahedralMesh_3_3::*)(unsigned int, ::boost::numeric::ublas::c_matrix<double, 3, 3> &, double &, ::boost::numeric::ublas::c_matrix<double, 3, 3> &) const) &AbstractTetrahedralMesh_3_3::GetInverseJacobianForElement,
            " ", py::arg("elementIndex"), py::arg("rJacobian"), py::arg("rJacobianDeterminant"), py::arg("rInverseJacobian"))
        .def("GetWeightedDirectionForBoundaryElement",
            (void(AbstractTetrahedralMesh_3_3::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 3> &, double &) const) &AbstractTetrahedralMesh_3_3::GetWeightedDirectionForBoundaryElement,
            " ", py::arg("elementIndex"), py::arg("rWeightedDirection"), py::arg("rJacobianDeterminant"))
        .def("CheckOutwardNormals",
            (void(AbstractTetrahedralMesh_3_3::*)()) &AbstractTetrahedralMesh_3_3::CheckOutwardNormals,
            " ")
        .def("ConstructLinearMesh",
            (void(AbstractTetrahedralMesh_3_3::*)(unsigned int)) &AbstractTetrahedralMesh_3_3::ConstructLinearMesh,
            " ", py::arg("width"))
        .def("ConstructRectangularMesh",
            (void(AbstractTetrahedralMesh_3_3::*)(unsigned int, unsigned int, bool)) &AbstractTetrahedralMesh_3_3::ConstructRectangularMesh,
            " ", py::arg("width"), py::arg("height"), py::arg("stagger") = true)
        .def("ConstructCuboid",
            (void(AbstractTetrahedralMesh_3_3::*)(unsigned int, unsigned int, unsigned int)) &AbstractTetrahedralMesh_3_3::ConstructCuboid,
            " ", py::arg("width"), py::arg("height"), py::arg("depth"))
        .def("ConstructRegularSlabMesh",
            (void(AbstractTetrahedralMesh_3_3::*)(double, double, double, double)) &AbstractTetrahedralMesh_3_3::ConstructRegularSlabMesh,
            " ", py::arg("spaceStep"), py::arg("width"), py::arg("height") = 0, py::arg("depth") = 0)
        .def("ConstructRegularSlabMeshWithDimensionSplit",
            (void(AbstractTetrahedralMesh_3_3::*)(unsigned int, double, double, double, double)) &AbstractTetrahedralMesh_3_3::ConstructRegularSlabMeshWithDimensionSplit,
            " ", py::arg("dimension"), py::arg("spaceStep"), py::arg("width"), py::arg("height") = 0, py::arg("depth") = 0)
        .def("CalculateDesignatedOwnershipOfBoundaryElement",
            (bool(AbstractTetrahedralMesh_3_3::*)(unsigned int)) &AbstractTetrahedralMesh_3_3::CalculateDesignatedOwnershipOfBoundaryElement,
            " ", py::arg("faceIndex"))
        .def("CalculateDesignatedOwnershipOfElement",
            (bool(AbstractTetrahedralMesh_3_3::*)(unsigned int)) &AbstractTetrahedralMesh_3_3::CalculateDesignatedOwnershipOfElement,
            " ", py::arg("elementIndex"))
        .def("CalculateMaximumNodeConnectivityPerProcess",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)() const) &AbstractTetrahedralMesh_3_3::CalculateMaximumNodeConnectivityPerProcess,
            " ")
        .def("GetHaloNodeIndices",
            (void(AbstractTetrahedralMesh_3_3::*)(::std::vector<unsigned int> &) const) &AbstractTetrahedralMesh_3_3::GetHaloNodeIndices,
            " ", py::arg("rHaloIndices"))
        .def("CalculateNodeExchange",
            (void(AbstractTetrahedralMesh_3_3::*)(::std::vector<std::vector<unsigned int>> &, ::std::vector<std::vector<unsigned int>> &)) &AbstractTetrahedralMesh_3_3::CalculateNodeExchange,
            " ", py::arg("rNodesToSendPerProcess"), py::arg("rNodesToReceivePerProcess"))
        .def("CalculateMinMaxEdgeLengths",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractTetrahedralMesh_3_3::*)()) &AbstractTetrahedralMesh_3_3::CalculateMinMaxEdgeLengths,
            " ")
        .def("GetContainingElementIndex",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)(::ChastePoint<3> const &, bool, ::std::set<unsigned int>, bool)) &AbstractTetrahedralMesh_3_3::GetContainingElementIndex,
            " ", py::arg("rTestPoint"), py::arg("strict") = false, py::arg("testElements") = std::set<unsigned int>(), py::arg("onlyTryWithTestElements") = false)
        .def("GetNearestElementIndexFromTestElements",
            (unsigned int(AbstractTetrahedralMesh_3_3::*)(::ChastePoint<3> const &, ::std::set<unsigned int>)) &AbstractTetrahedralMesh_3_3::GetNearestElementIndexFromTestElements,
            " ", py::arg("rTestPoint"), py::arg("testElements"))
    ;
}
