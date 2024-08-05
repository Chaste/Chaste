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

#include "AbstractCentreBasedCellPopulation3_3.cppwg.hpp"

namespace py = pybind11;
typedef AbstractCentreBasedCellPopulation<3,3 > AbstractCentreBasedCellPopulation3_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);
typedef ::boost::numeric::ublas::c_vector<double, 3> _boost_numeric_ublas_c_vector_lt_double_3_gt_;
typedef ::CellPtr _CellPtr;
typedef ::std::set<unsigned int> _std_set_lt_unsignedint_gt_;
typedef ::std::vector<std::pair<Node<3> *, Node<3> *>> & _std_vector_lt_std_pair_lt_Node_lt_3_gt_Ptr_Node_lt_3_gt_Ptr_gt__gt_Ref;

class AbstractCentreBasedCellPopulation3_3_Overrides : public AbstractCentreBasedCellPopulation3_3{
    public:
    using AbstractCentreBasedCellPopulation3_3::AbstractCentreBasedCellPopulation;
    ::boost::numeric::ublas::c_vector<double, 3> GetLocationOfCellCentre(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _boost_numeric_ublas_c_vector_lt_double_3_gt_,
            AbstractCentreBasedCellPopulation3_3,
            GetLocationOfCellCentre,
                    pCell);
    }
    double GetCellDataItemAtPdeNode(unsigned int pdeNodeIndex, ::std::string & rVariableName, bool dirichletBoundaryConditionApplies, double dirichletBoundaryValue) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCentreBasedCellPopulation3_3,
            GetCellDataItemAtPdeNode,
                    pdeNodeIndex,
        rVariableName,
        dirichletBoundaryConditionApplies,
        dirichletBoundaryValue);
    }
    ::CellPtr AddCell(::CellPtr pNewCell, ::CellPtr pParentCell) override {
        PYBIND11_OVERRIDE(
            _CellPtr,
            AbstractCentreBasedCellPopulation3_3,
            AddCell,
                    pNewCell,
        pParentCell);
    }
    bool IsCellAssociatedWithADeletedLocation(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCentreBasedCellPopulation3_3,
            IsCellAssociatedWithADeletedLocation,
                    pCell);
    }
    ::std::set<unsigned int> GetNeighbouringLocationIndices(::CellPtr pCell) override {
        PYBIND11_OVERRIDE(
            _std_set_lt_unsignedint_gt_,
            AbstractCentreBasedCellPopulation3_3,
            GetNeighbouringLocationIndices,
                    pCell);
    }
    void CheckForStepSizeException(unsigned int nodeIndex, ::boost::numeric::ublas::c_vector<double, 3> & rDisplacement, double dt) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedCellPopulation3_3,
            CheckForStepSizeException,
                    nodeIndex,
        rDisplacement,
        dt);
    }
    double GetDampingConstant(unsigned int nodeIndex) override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCentreBasedCellPopulation3_3,
            GetDampingConstant,
                    nodeIndex);
    }
    bool IsGhostNode(unsigned int index) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCentreBasedCellPopulation3_3,
            IsGhostNode,
                    index);
    }
    bool IsParticle(unsigned int index) override {
        PYBIND11_OVERRIDE(
            bool,
            AbstractCentreBasedCellPopulation3_3,
            IsParticle,
                    index);
    }
    ::std::vector<std::pair<Node<3> *, Node<3> *>> & rGetNodePairs() override {
        PYBIND11_OVERRIDE_PURE(
            _std_vector_lt_std_pair_lt_Node_lt_3_gt_Ptr_Node_lt_3_gt_Ptr_gt__gt_Ref,
            AbstractCentreBasedCellPopulation3_3,
            rGetNodePairs,
            );
    }
    void OutputCellPopulationParameters(::out_stream & rParamsFile) override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedCellPopulation3_3,
            OutputCellPopulationParameters,
                    rParamsFile);
    }
    double GetDefaultTimeStep() override {
        PYBIND11_OVERRIDE(
            double,
            AbstractCentreBasedCellPopulation3_3,
            GetDefaultTimeStep,
            );
    }
    void WriteVtkResultsToFile(::std::string const & rDirectory) override {
        PYBIND11_OVERRIDE_PURE(
            void,
            AbstractCentreBasedCellPopulation3_3,
            WriteVtkResultsToFile,
                    rDirectory);
    }
    void AcceptCellWritersAcrossPopulation() override {
        PYBIND11_OVERRIDE(
            void,
            AbstractCentreBasedCellPopulation3_3,
            AcceptCellWritersAcrossPopulation,
            );
    }

};
void register_AbstractCentreBasedCellPopulation3_3_class(py::module &m){
py::class_<AbstractCentreBasedCellPopulation3_3 , AbstractCentreBasedCellPopulation3_3_Overrides , boost::shared_ptr<AbstractCentreBasedCellPopulation3_3 >  , AbstractOffLatticeCellPopulation<3>  >(m, "AbstractCentreBasedCellPopulation3_3")
        .def(
            "GetLocationOfCellCentre",
            (::boost::numeric::ublas::c_vector<double, 3>(AbstractCentreBasedCellPopulation3_3::*)(::CellPtr)) &AbstractCentreBasedCellPopulation3_3::GetLocationOfCellCentre,
            " " , py::arg("pCell") )
        .def(
            "GetNodeCorrespondingToCell",
            (::Node<3> *(AbstractCentreBasedCellPopulation3_3::*)(::CellPtr)) &AbstractCentreBasedCellPopulation3_3::GetNodeCorrespondingToCell,
            " " , py::arg("pCell") , py::return_value_policy::reference)
        .def(
            "GetCellDataItemAtPdeNode",
            (double(AbstractCentreBasedCellPopulation3_3::*)(unsigned int, ::std::string &, bool, double)) &AbstractCentreBasedCellPopulation3_3::GetCellDataItemAtPdeNode,
            " " , py::arg("pdeNodeIndex"), py::arg("rVariableName"), py::arg("dirichletBoundaryConditionApplies") = false, py::arg("dirichletBoundaryValue") = 0. )
        .def(
            "AddCell",
            (::CellPtr(AbstractCentreBasedCellPopulation3_3::*)(::CellPtr, ::CellPtr)) &AbstractCentreBasedCellPopulation3_3::AddCell,
            " " , py::arg("pNewCell"), py::arg("pParentCell") = ::CellPtr( ) )
        .def(
            "CreateCellPair",
            (::std::pair<boost::shared_ptr<Cell>, boost::shared_ptr<Cell>>(AbstractCentreBasedCellPopulation3_3::*)(::CellPtr, ::CellPtr)) &AbstractCentreBasedCellPopulation3_3::CreateCellPair,
            " " , py::arg("pCell1"), py::arg("pCell2") )
        .def(
            "IsMarkedSpring",
            (bool(AbstractCentreBasedCellPopulation3_3::*)(::std::pair<boost::shared_ptr<Cell>, boost::shared_ptr<Cell>> const &)) &AbstractCentreBasedCellPopulation3_3::IsMarkedSpring,
            " " , py::arg("rCellPair") )
        .def(
            "IsCellAssociatedWithADeletedLocation",
            (bool(AbstractCentreBasedCellPopulation3_3::*)(::CellPtr)) &AbstractCentreBasedCellPopulation3_3::IsCellAssociatedWithADeletedLocation,
            " " , py::arg("pCell") )
        .def(
            "GetNeighbouringLocationIndices",
            (::std::set<unsigned int>(AbstractCentreBasedCellPopulation3_3::*)(::CellPtr)) &AbstractCentreBasedCellPopulation3_3::GetNeighbouringLocationIndices,
            " " , py::arg("pCell") )
        .def(
            "CheckForStepSizeException",
            (void(AbstractCentreBasedCellPopulation3_3::*)(unsigned int, ::boost::numeric::ublas::c_vector<double, 3> &, double)) &AbstractCentreBasedCellPopulation3_3::CheckForStepSizeException,
            " " , py::arg("nodeIndex"), py::arg("rDisplacement"), py::arg("dt") )
        .def(
            "GetDampingConstant",
            (double(AbstractCentreBasedCellPopulation3_3::*)(unsigned int)) &AbstractCentreBasedCellPopulation3_3::GetDampingConstant,
            " " , py::arg("nodeIndex") )
        .def(
            "IsGhostNode",
            (bool(AbstractCentreBasedCellPopulation3_3::*)(unsigned int)) &AbstractCentreBasedCellPopulation3_3::IsGhostNode,
            " " , py::arg("index") )
        .def(
            "IsParticle",
            (bool(AbstractCentreBasedCellPopulation3_3::*)(unsigned int)) &AbstractCentreBasedCellPopulation3_3::IsParticle,
            " " , py::arg("index") )
        .def(
            "rGetNodePairs",
            (::std::vector<std::pair<Node<3> *, Node<3> *>> &(AbstractCentreBasedCellPopulation3_3::*)()) &AbstractCentreBasedCellPopulation3_3::rGetNodePairs,
            " "  , py::return_value_policy::reference_internal)
        .def(
            "GetMeinekeDivisionSeparation",
            (double(AbstractCentreBasedCellPopulation3_3::*)()) &AbstractCentreBasedCellPopulation3_3::GetMeinekeDivisionSeparation,
            " "  )
        .def(
            "SetMeinekeDivisionSeparation",
            (void(AbstractCentreBasedCellPopulation3_3::*)(double)) &AbstractCentreBasedCellPopulation3_3::SetMeinekeDivisionSeparation,
            " " , py::arg("divisionSeparation") )
        .def(
            "GetCentreBasedDivisionRule",
            (::boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>>(AbstractCentreBasedCellPopulation3_3::*)()) &AbstractCentreBasedCellPopulation3_3::GetCentreBasedDivisionRule,
            " "  )
        .def(
            "SetCentreBasedDivisionRule",
            (void(AbstractCentreBasedCellPopulation3_3::*)(::boost::shared_ptr<AbstractCentreBasedDivisionRule<3, 3>>)) &AbstractCentreBasedCellPopulation3_3::SetCentreBasedDivisionRule,
            " " , py::arg("pCentreBasedDivisionRule") )
        .def(
            "OutputCellPopulationParameters",
            (void(AbstractCentreBasedCellPopulation3_3::*)(::out_stream &)) &AbstractCentreBasedCellPopulation3_3::OutputCellPopulationParameters,
            " " , py::arg("rParamsFile") )
        .def(
            "GetDefaultTimeStep",
            (double(AbstractCentreBasedCellPopulation3_3::*)()) &AbstractCentreBasedCellPopulation3_3::GetDefaultTimeStep,
            " "  )
    ;
}
