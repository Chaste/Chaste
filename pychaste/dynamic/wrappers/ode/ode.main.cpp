#include <pybind11/pybind11.h>
#include "AbstractOdeSystemInformation.cppwg.hpp"
#include "AbstractOdeSystem.cppwg.hpp"
#include "DeltaNotchOdeSystem.cppwg.hpp"
#include "DeltaNotchEdgeOdeSystem.cppwg.hpp"
#include "DeltaNotchInteriorOdeSystem.cppwg.hpp"
#include "Alarcon2004OxygenBasedCellCycleOdeSystem.cppwg.hpp"
#include "Goldbeter1991OdeSystem.cppwg.hpp"
#include "TysonNovak2001OdeSystem.cppwg.hpp"
#include "AbstractPythonOdeSystemInformation.cppwg.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_pychaste_ode, m)
{
    register_AbstractOdeSystemInformation_class(m);
    register_AbstractOdeSystem_class(m);
    register_DeltaNotchOdeSystem_class(m);
    register_DeltaNotchEdgeOdeSystem_class(m);
    register_DeltaNotchInteriorOdeSystem_class(m);
    register_Alarcon2004OxygenBasedCellCycleOdeSystem_class(m);
    register_Goldbeter1991OdeSystem_class(m);
    register_TysonNovak2001OdeSystem_class(m);
    register_AbstractPythonOdeSystemInformation_class(m);
}
