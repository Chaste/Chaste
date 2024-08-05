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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

#include "AbstractCentreBasedCellPopulation2_2.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCentreBasedCellPopulation<2,2 > AbstractCentreBasedCellPopulation2_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 2> _boost_numeric_ublas_c_vector_lt_double_2_gt_;
typedef ::CellPtr _CellPtr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::vector<std::pair<Node<2> *, Node<2> *>> & _std_vector_lt_std_pair_lt_Node_lt_2_gt_Ptr_Node_lt_2_gt_Ptr_gt__gt_Ref;

class AbstractCentreBasedCellPopulation2_2_Overrides : public AbstractCentreBasedCellPopulation2_2{
    public:
    using AbstractCentreBasedCellPopulation2_2::AbstractCentreBasedCellPopulation;
    ::boost::numeric::ublas::c_vector<double, 2> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_2_gt_,
            AbstractCentreBasedCellPopulation2_2,
            GetLocationOfCellCentre,
                    pCell);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCentreBasedCellPopulation2_2,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            AbstractCentreBasedCellPopulation2_2,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCentreBasedCellPopulation2_2,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            AbstractCentreBasedCellPopulation2_2,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 2> & rDisplacement, double dt) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedCellPopulation2_2,
            CheckForStepSizeException,
                    nodeIndex,
        rDisplacement,
        dt);
    }
    double GetDampingConstant(unsigned int nodeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCentreBasedCellPopulation2_2,
            GetDampingConstant,
                    nodeIndex);
    }
    bool IsGhostNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCentreBasedCellPopulation2_2,
            IsGhostNode,
                    index);
    }
    bool IsParticle(unsigned int index) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCentreBasedCellPopulation2_2,
            IsParticle,
                    index);
    }
    ::std::vector<std::pair<Node<2> *, Node<2> *>> & rGetNodePairs() override {
        PYBIND11_OVERRIDE_PURE(
            _std_vector_lt_std_pair_lt_Node_lt_2_gt_Ptr_Node_lt_2_gt_Ptr_gt__gt_Ref,
            AbstractCentreBasedCellPopulation2_2,
            rGetNodePairs,
            );
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedCellPopulation2_2,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCentreBasedCellPopulation2_2,
            GetDefaultTimeStep,
            );
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCentreBasedCellPopulation2_2,
            WriteVtkResultsToFile,
                    rDirectory);
    }
    void AcceptCellWritersAcrossPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedCellPopulation2_2,
            AcceptCellWritersAcrossPopulation,
            );
    }

};
void register_AbstractCentreBasedCellPopulation2_2_class(py::module &m){
py::class_<AbstractCentreBasedCellPopulation2_2 , AbstractCentreBasedCellPopulation2_2_Overrides , boost::shared_ptr<AbstractCentreBasedCellPopulation2_2 >  , AbstractOffLatticeCellPopulation<2>  >(m, "AbstractCentreBasedCellPopulation2_2")
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 2>(AbstractCentreBasedCellPopulation2_2::*)(::CellPtr)) &AbstractCentreBasedCellPopulation2_2::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "GetNodeCorrespondingToCell",
            (::Node<2> *(AbstractCentreBasedCellPopulation2_2::*)(::CellPtr)) &AbstractCentreBasedCellPopulation2_2::GetNodeCorrespondingToCell,
            " " , py::arg("pCell") , py::return_value_policy::reference)
        .def(
            "GetCellDataItemAtPdeNode",
            (double(AbstractCentreBasedCellPopulation2_2::*)(unsigned int, ::std::string &, bool, double)) &AbstractCentreBasedCellPopulation2_2::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def(
            "AddCell",
            (::CellPtr(AbstractCentreBasedCellPopulation2_2::*)(::CellPtr, ::CellPtr)) &AbstractCentreBasedCellPopulation2_2::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "CreateCellPair",
            (::std::pair<boost::shared_ptr<Cell>, boost::shared_ptr<Cell>>(AbstractCentreBasedCellPopulation2_2::*)(::CellPtr, ::CellPtr)) &AbstractCentreBasedCellPopulation2_2::CreateCellPair,
            " " , py::arg("pCell1"), py::arg("pCell2") )
        .def(
            "IsMarkedSpring",
            (bool(AbstractCentreBasedCellPopulation2_2::*)(::std::pair<boost::shared_ptr<Cell>, boost::shared_ptr<Cell>> const &)) &AbstractCentreBasedCellPopulation2_2::IsMarkedSpring,
            " " , py::arg("rCellPair") )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(AbstractCentreBasedCellPopulation2_2::*)(::CellPtr)) &AbstractCentreBasedCellPopulation2_2::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(AbstractCentreBasedCellPopulation2_2::*)(::CellPtr)) &AbstractCentreBasedCellPopulation2_2::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "CheckForStepSizeException",
            (void(AbstractCentreBasedCellPopulation2_2::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 2> &, double)) &AbstractCentreBasedCellPopulation2_2::CheckForStepSizeException,
            " " , py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt") )
        .def(
            "GetDampingConstant",
            (double(AbstractCentreBasedCellPopulation2_2::*)(unsigned int)) &AbstractCentreBasedCellPopulation2_2::GetDampingConstant,
            " " , py::arg("nodeIndex") )
        .def(
            "IsGhostNode",
            (bool(AbstractCentreBasedCellPopulation2_2::*)(unsigned int)) &AbstractCentreBasedCellPopulation2_2::IsGhostNode,
            " " , py::arg("index") )
        .def(
            "IsParticle",
            (bool(AbstractCentreBasedCellPopulation2_2::*)(unsigned int)) &AbstractCentreBasedCellPopulation2_2::IsParticle,
            " " , py::arg("index") )
        .def(
            "rGetNodePairs",
            (::std::vector<std::pair<Node<2> *, Node<2> *>> &(AbstractCentreBasedCellPopulation2_2::*)()) &AbstractCentreBasedCellPopulation2_2::rGetNodePairs,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetMeinekeDivisionSeparation",
            (double(AbstractCentreBasedCellPopulation2_2::*)()) &AbstractCentreBasedCellPopulation2_2::GetMeinekeDivisionSeparation,
            " "  )
        .def(
            "SetMeinekeDivisionSeparation",
            (void(AbstractCentreBasedCellPopulation2_2::*)(double)) &AbstractCentreBasedCellPopulation2_2::SetMeinekeDivisionSeparation,
            " " , py::arg("divisionSeparation") )
        .def(
            "GetCentreBasedDivisionRule",
            (::boost::shared_ptr<AbstractCentreBasedDivisionRule<2, 2>>(AbstractCentreBasedCellPopulation2_2::*)()) &AbstractCentreBasedCellPopulation2_2::GetCentreBasedDivisionRule,
            " "  )
        .def(
            "SetCentreBasedDivisionRule",
            (void(AbstractCentreBasedCellPopulation2_2::*)(::boost::shared_ptr<AbstractCentreBasedDivisionRule<2, 2>>)) &AbstractCentreBasedCellPopulation2_2::SetCentreBasedDivisionRule,
            " " , py::arg("pCentreBasedDivisionRule") )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractCentreBasedCellPopulation2_2::*)(::out_stream &)) &AbstractCentreBasedCellPopulation2_2::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetDefaultTimeStep",
            (double(AbstractCentreBasedCellPopulation2_2::*)()) &AbstractCentreBasedCellPopulation2_2::GetDefaultTimeStep,
            " "  )
    ;
}
