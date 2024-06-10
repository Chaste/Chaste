
import cppwg.templates.custom


class CellsGenerator_custom(cppwg.templates.custom.Custom):

    def __init__(self):
        pass

    def get_class_cpp_pre_code(self, class_name):

        replacements = {"class_name": class_name}
        code = """
class {class_name}_Overloads : {class_name}{{
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
        return code.format(**replacements)

    def get_class_cpp_def_code(self, class_name):

        replacements = {"class_name": class_name}
        code = """        .def("GenerateBasic",
            (std::vector<CellPtr>({class_name}::*)(unsigned int, const std::vector<unsigned>, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overloads::GenerateBasic, " " , py::arg("numCells"),  py::arg("locationIndices") = std::vector<unsigned int>(), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateBasicRandom",
            (std::vector<CellPtr>({class_name}::*)(unsigned int, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overloads::GenerateBasicRandom, " " , py::arg("numCells"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
        .def("GenerateGivenLocationIndices",
            (std::vector<CellPtr>({class_name}::*)(const std::vector<unsigned> locationIndices, boost::shared_ptr<AbstractCellProperty>)) 
            &{class_name}_Overloads::GenerateGivenLocationIndices, " " , py::arg("locationIndices"), py::arg("pCellProliferativeType") = boost::shared_ptr<AbstractCellProperty>())
"""
        return code.format(**replacements)
