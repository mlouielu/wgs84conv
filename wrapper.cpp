#include <armadillo>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

#include "carma/carma.h"
#include "wgs84conv.hh"

PYBIND11_MODULE(wgs84conv, m) {
  m.def("ecef2lla", &ecef2lla, py::arg("vec3"));
  m.def("lla2ecef", &lla2ecef, py::arg("vec3"));
}
