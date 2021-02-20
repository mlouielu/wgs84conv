wgs84conv
=========

wgs84conv is a accurate & fast Earth-Fixed Earth-Centered (ECEF) to Geodetic (lla)
convert module based on:

> Karl Osen. Accurate Conversion of Earth-Fixed Earth-Centered Coordinates to Geodetic Coordinates.[Research Report] Norwegian University of Science and Technology. 2017. hal-01704943v2

And with `pybind11` & `carma`, this module also support Python binding.


Prerequirements
---------------

* armadillo
* pybind11

Build
-----

* mkdir build && cd build
* cmake -G"Ninja" ..

API
---

### ecef2lla

* Input: ECEF position in meter
* Output: WGS84 LLA in radius, meter

### lla2ecef

* Input: WGS84 LLA in radius
* Output: ECEF position in meter


Compare
-------

Compare to pyproj convertion

```
$ time python gps_track_pybind11.py test.tlm
    0.42s user 0.04s system 99% cpu 0.459 total

$ time python gps_track.py test.tlm
    34.54s user 3.20s system 98% cpu 38.256 total
```
