#ifndef WGS84_WGS84CONV_HH
#define WGS84_WGS84CONV_HH

#include <armadillo>

#include "wgs84const.hh"

arma::vec3 lla2ecef(const arma::vec3 &lla);
arma::vec3 ecef2lla(const arma::vec3 &ecef);

#endif /* WGS84_WGS84CONV_HH */
