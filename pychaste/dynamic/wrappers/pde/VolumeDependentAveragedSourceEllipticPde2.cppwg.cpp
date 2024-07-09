#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"

#include "VolumeDependentAveragedSourceEllipticPde2.cppwg.hpp"

namespace py = pybind11;
typedef VolumeDependentAveragedSourceEllipticPde<2 > VolumeDependentAveragedSourceEllipticPde2;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeDependentAveragedSourceEllipticPde2_Overrides : public VolumeDependentAveragedSourceEllipticPde2{
    public:
    using VolumeDependentAveragedSourceEllipticPde2::VolumeDependentAveragedSourceEllipticPde;
    void SetupSourceTerms(::TetrahedralMesh<2, 2> & rCoarseMesh, ::std::map<boost::shared_ptr<Cell>, unsigned int> * pCellPdeElementMap) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeDependentAveragedSourceEllipticPde2,
            SetupSourceTerms,
                    rCoarseMesh,
        pCellPdeElementMap);
    }

};
void register_VolumeDependentAveragedSourceEllipticPde2_class(py::module &m){
py::class_<VolumeDependentAveragedSourceEllipticPde2 , VolumeDependentAveragedSourceEllipticPde2_Overrides , boost::shared_ptr<VolumeDependentAveragedSourceEllipticPde2 >  , AveragedSourceEllipticPde<2>  >(m, "VolumeDependentAveragedSourceEllipticPde2")
        .def(py::init<::AbstractCellPopulation<2> &, double >(), py::arg("rCellPopulation"), py::arg("coefficient") = 0.)
        .def(
            "SetupSourceTerms",
            (void(VolumeDependentAveragedSourceEllipticPde2::*)(::TetrahedralMesh<2, 2> &, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &VolumeDependentAveragedSourceEllipticPde2::SetupSourceTerms,
            " " , py::arg("rCoarseMesh"), py::arg("pCellPdeElementMap") = nullptr )
    ;
}
