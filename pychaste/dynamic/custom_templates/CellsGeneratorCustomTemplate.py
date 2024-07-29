import cppwg.templates.custom


class CellsGeneratorCustomTemplate(cppwg.templates.custom.Custom):

    def __init__(self):
        pass

    def get_class_cpp_pre_code(self, class_name):
        """
        Creates custom pybind11 trampoline wrapper code for CellsGenerator.
        """

        code = f"""
class {class_name}_Overrides : {class_name}{{
    public:
    using {class_name}::{class_name};
    std::vector<CellPtr> GenerateBasic(unsigned numCells, const std::vector<unsigned> locationIndices=std::vector<unsigned>(), boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {{
        std::vector<CellPtr> cells;
        {class_name}::GenerateBasic(cells, numCells, locationIndices, pCellProliferativeType);
        return cells;
    }};
    std::vector<CellPtr> GenerateBasicRandom(unsigned numCells, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {{
        std::vector<CellPtr> cells;
        {class_name}::GenerateBasicRandom(cells, numCells, pCellProliferativeType);
        return cells;
    }};
    std::vector<CellPtr> GenerateGivenLocationIndices(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty> pCellProliferativeType=boost::shared_ptr<AbstractCellProperty>())
    {{
        std::vector<CellPtr> cells;
        {class_name}::GenerateGivenLocationIndices(cells, locationIndices, pCellProliferativeType);
        return cells;
    }};
}};

"""
        return code

    def get_class_cpp_def_code(self, class_name):
        """
        Creates custom pybind11 class wrapper code for CellsGenerator.
        """

        code = f"""
        .def("GenerateBasic",
            (std::vector<CellPtr>({class_name}::*)(unsigned int, const std::vector<unsigned>, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overrides::GenerateBasic, " " , py::arg("numCells"),  py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateBasicRandom",
            (std::vector<CellPtr>({class_name}::*)(unsigned int, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overrides::GenerateBasicRandom, " " , py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateGivenLocationIndices",
            (std::vector<CellPtr>({class_name}::*)(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overrides::GenerateGivenLocationIndices, " " , py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
"""
        return code
