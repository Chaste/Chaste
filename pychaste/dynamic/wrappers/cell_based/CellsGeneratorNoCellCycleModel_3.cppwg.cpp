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

#include "CellsGeneratorNoCellCycleModel_3.cppwg.hpp"

namespace py = pybind11;
typedef CellsGenerator<NoCellCycleModel,3 > CellsGeneratorNoCellCycleModel_3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class CellsGeneratorNoCellCycleModel_3_Overrides : CellsGeneratorNoCellCycleModel_3{
    public:
    using CellsGeneratorNoCellCycleModel_3::CellsGeneratorNoCellCycleModel_3;
    std::vector<CellPtr> GenerateBasic(unsigned numCells, const std::vector<unsigned> locationIndices=std::vector<unsigned>(), boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {
        std::vector<CellPtr> cells;
        CellsGeneratorNoCellCycleModel_3::GenerateBasic(cells, numCells, locationIndices, pCellProliferativeType);
        return cells;
    };
    std::vector<CellPtr> GenerateBasicRandom(unsigned numCells, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {
        std::vector<CellPtr> cells;
        CellsGeneratorNoCellCycleModel_3::GenerateBasicRandom(cells, numCells, pCellProliferativeType);
        return cells;
    };
    std::vector<CellPtr> GenerateGivenLocationIndices(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {
        std::vector<CellPtr> cells;
        CellsGeneratorNoCellCycleModel_3::GenerateGivenLocationIndices(cells, locationIndices, pCellProliferativeType);
        return cells;
    };
};


void register_CellsGeneratorNoCellCycleModel_3_class(py::module &m){
py::class_<CellsGeneratorNoCellCycleModel_3  , boost::shared_ptr<CellsGeneratorNoCellCycleModel_3 >   >(m, "CellsGeneratorNoCellCycleModel_3")
        .def(py::init< >())
        .def(
            "GenerateBasic",
            (void(CellsGeneratorNoCellCycleModel_3::*)(::std::vector<boost::shared_ptr<Cell>> &, unsigned int, ::std::vector<unsigned int> const, ::boost::shared_ptr<AbstractCellProperty>)) &CellsGeneratorNoCellCycleModel_3::GenerateBasic,
            " " , py::arg("rCells"), py::arg("numCells"), py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>() )
        .def(
            "GenerateBasicRandom",
            (void(CellsGeneratorNoCellCycleModel_3::*)(::std::vector<boost::shared_ptr<Cell>> &, unsigned int, ::boost::shared_ptr<AbstractCellProperty>)) &CellsGeneratorNoCellCycleModel_3::GenerateBasicRandom,
            " " , py::arg("rCells"), py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>() )
        .def(
            "GenerateGivenLocationIndices",
            (void(CellsGeneratorNoCellCycleModel_3::*)(::std::vector<boost::shared_ptr<Cell>> &, ::std::vector<unsigned int> const, ::boost::shared_ptr<AbstractCellProperty>)) &CellsGeneratorNoCellCycleModel_3::GenerateGivenLocationIndices,
            " " , py::arg("rCells"), py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>() )

        .def("GenerateBasic",
            (std::vector<CellPtr>(CellsGeneratorNoCellCycleModel_3::*)(unsigned int, const std::vector<unsigned>, boost::shared_ptr<AbstractCellProperty>)) 
            &CellsGeneratorNoCellCycleModel_3_Overrides::GenerateBasic, " " , py::arg("numCells"),  py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateBasicRandom",
            (std::vector<CellPtr>(CellsGeneratorNoCellCycleModel_3::*)(unsigned int, boost::shared_ptr<AbstractCellProperty>)) 
            &CellsGeneratorNoCellCycleModel_3_Overrides::GenerateBasicRandom, " " , py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateGivenLocationIndices",
            (std::vector<CellPtr>(CellsGeneratorNoCellCycleModel_3::*)(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty>)) 
            &CellsGeneratorNoCellCycleModel_3_Overrides::GenerateGivenLocationIndices, " " , py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
    ;
}
