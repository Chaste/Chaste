#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <set>
#include <vector>
#include <string>
#include <map>
#include "SmartPointers.hpp"
#include "UblasIncludes.hpp"
#include "VolumeDependentAveragedSourceEllipticPde.hpp"

#include "VolumeDependentAveragedSourceEllipticPde3.cppwg.hpp"

namespace py = pybind11;
typedef VolumeDependentAveragedSourceEllipticPde<3 > VolumeDependentAveragedSourceEllipticPde3;
PYBIND11_DECLARE_HOLDER_TYPE(T, boost::shared_ptr<T>);

class VolumeDependentAveragedSourceEllipticPde3_Overrides : public VolumeDependentAveragedSourceEllipticPde3{
    public:
    using VolumeDependentAveragedSourceEllipticPde3::VolumeDependentAveragedSourceEllipticPde;
    void SetupSourceTerms(::TetrahedralMesh<3, 3> & rCoarseMesh, ::std::map<boost::shared_ptr<Cell>, unsigned int> * pCellPdeElementMap) override {
        PYBIND11_OVERRIDE(
            void,
            VolumeDependentAveragedSourceEllipticPde3,
            SetupSourceTerms,
                    rCoarseMesh,
        pCellPdeElementMap);
    }

};
void register_VolumeDependentAveragedSourceEllipticPde3_class(py::module &m){
py::class_<VolumeDependentAveragedSourceEllipticPde3 , VolumeDependentAveragedSourceEllipticPde3_Overrides , boost::shared_ptr<VolumeDependentAveragedSourceEllipticPde3 >  , AveragedSourceEllipticPde<3>  >(m, "VolumeDependentAveragedSourceEllipticPde3")
        .def(py::init<::AbstractCellPopulation<3> &, double >(), py::arg("rCellPopulation"), py::arg("coefficient") = 0.)
        .def(
            "SetupSourceTerms",
            (void(VolumeDependentAveragedSourceEllipticPde3::*)(::TetrahedralMesh<3, 3> &, ::std::map<boost::shared_ptr<Cell>, unsigned int> *)) &VolumeDependentAveragedSourceEllipticPde3::SetupSourceTerms,
            " " , py::arg("rCoarseMesh"), py::arg("pCellPdeElementMap") = nullptr )
    ;
}
