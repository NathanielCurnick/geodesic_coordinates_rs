# Geodesic Coordinates 

Crate for calculating coordinates in various reference frames on or close to the earth

## Basic Usage

Convert a WGS84 latitude and longitude to ECEF frame 

```
let reference_point = WGS84Coord::new_from_degrees(50.0_f64, 10.0_f64, 0.0);

let ecef = ECEF::new_from_wgs84(&reference_point);
```

Convert ECEF to NED 

```
let lat = 50.0_f64;
let lon = 10.0_f64;
let alt = 0.0_f64;
let reference_point = WGS84Coord::new_from_degrees(lat, lon, alt);

let rotation = generate_ecef_to_ned_matrix(&reference_point);

let reference_point = ECEF::new_from_wgs84(&reference_point);

let pos = ECEF {x: 3004296, y: 1093474, z: 5500477};

let ned = NED::new_from_ecef_rot(&pos, &rotation, &reference_point);
```