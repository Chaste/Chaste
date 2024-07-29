#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "NoCellCycleModel.hpp"
#include "UniformCellCycleModel.hpp"
#include "SimpleOxygenBasedCellCycleModel.hpp"
#include "UniformG1GenerationalCellCycleModel.hpp"
#include "BiasedBernoulliTrialCellCycleModel.hpp"
#include "LabelDependentBernoulliTrialCellCycleModel.hpp"
#include "AlwaysDivideCellCycleModel.hpp"
#include "ContactInhibitionCellCycleModel.hpp"
#include "StochasticOxygenBasedCellCycleModel.hpp"
#include "GammaG1CellCycleModel.hpp"
#include "ExponentialG1GenerationalCellCycleModel.hpp"
#include "TysonNovakCellCycleModel.hpp"
#include "Alarcon2004OxygenBasedCellCycleModel.hpp"
#include "FixedSequenceCellCycleModel.hpp"
#include "BernoulliTrialCellCycleModel.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "CellsGenerator.hpp"

#include "CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2.cppwg.hpp"

namespace py = pybind11;
typedef CellsGenerator<LabelDependentBernoulliTrialCellCycleModel,2 > CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2_Overrides : CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2{
    public:
    using CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2;
    std::vector<CellPtr> GenerateBasic(unsigned numCells, const std::vector<unsigned> locationIndices=std::vector<unsigned>(), boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {
        std::vector<CellPtr> cells;
        CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::GenerateBasic(cells, numCells, locationIndices, pCellProliferativeType);
        return cells;
    };
    std::vector<CellPtr> GenerateBasicRandom(unsigned numCells, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {
        std::vector<CellPtr> cells;
        CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::GenerateBasicRandom(cells, numCells, pCellProliferativeType);
        return cells;
    };
    std::vector<CellPtr> GenerateGivenLocationIndices(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {
        std::vector<CellPtr> cells;
        CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::GenerateGivenLocationIndices(cells, locationIndices, pCellProliferativeType);
        return cells;
    };
};


void register_CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2_class(py::module &m){
py::class_<CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2  , boost::shared_ptr<CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2 >   >(m, "CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2")
        .def(py::init< >())
        .def(
            "GenerateBasic",
            (void(CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::*)(::std::vector<boost::shared_ptr<Cell>> &, unsigned int, ::std::vector<unsigned int> const, ::boost::shared_ptr<AbstractCellProperty>)) &CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::GenerateBasic,
            " " , py::arg("rCells"), py::arg("numCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>() )
        .def(
            "GenerateBasicRandom",
            (void(CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::*)(::std::vector<boost::shared_ptr<Cell>> &, unsigned int, ::boost::shared_ptr<AbstractCellProperty>)) &CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::GenerateBasicRandom,
            " " , py::arg("rCells"), py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>() )
        .def(
            "GenerateGivenLocationIndices",
            (void(CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::*)(::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, ::boost::shared_ptr<AbstractCellProperty>)) &CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::GenerateGivenLocationIndices,
            " " , py::arg("rCells"), py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>() )

        .def("GenerateBasic",
            (std::vector<CellPtr>(CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::*)(unsigned int, const std::vector<unsigned>, boost::shared_ptr<AbstractCellProperty>)) 
            &CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2_Overrides::GenerateBasic, " " , py::arg("numCells"),  py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateBasicRandom",
            (std::vector<CellPtr>(CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::*)(unsigned int, boost::shared_ptr<AbstractCellProperty>)) 
            &CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2_Overrides::GenerateBasicRandom, " " , py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateGivenLocationIndices",
            (std::vector<CellPtr>(CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2::*)(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty>)) 
            &CellsGeneratorLabelDependentBernoulliTrialCellCycleModel_2_Overrides::GenerateGivenLocationIndices, " " , py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
    ;
}
