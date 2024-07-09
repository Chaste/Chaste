#include <pybind11/pybind11.h>
#include "VtkScene2.cppwg.hpp"
#include "VtkScene3.cppwg.hpp"
#include "AbstractPyChasteActorGenerator2.cppwg.hpp"
#include "AbstractPyChasteActorGenerator3.cppwg.hpp"
#include "CellPopulationPyChasteActorGenerator2.cppwg.hpp"
#include "CellPopulationPyChasteActorGenerator3.cppwg.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pychaste_visualization, m)
{
    register_VtkScene2_class(m);
    register_VtkScene3_class(m);
    register_AbstractPyChasteActorGenerator2_class(m);
    register_AbstractPyChasteActorGenerator3_class(m);
    register_CellPopulationPyChasteActorGenerator2_class(m);
    register_CellPopulationPyChasteActorGenerator3_class(m);
}
