use crate::{
    constants::{DEG_TO_RAD, EARTH_ECCENTRICITY_SQUARED, EARTH_MAJOR},
    types::{Degrees, Metres, Radians},
    utils::transpose_times_vec,
};

use super::{ecef::ECEF, ned::NED};

#[derive(Debug, Clone)]
pub struct WGS84Coord {
    lat: Radians,
    lon: Radians,
    alt: Metres,
}

// Structs for WGS84 system coordinates
impl WGS84Coord {
    // pub fn new_from_astrocoord(astro: &AstroCoord) -> WGS84Coord {}
    pub fn new_from_radians(latitude: Radians, longitude: Radians, altitude: Metres) -> WGS84Coord {
        // given a latitude in radians, a longitude in radians and an altitude in m gives a WGS84Coord
        return WGS84Coord {
            lat: latitude,
            lon: longitude,
            alt: altitude,
        };
    }
    pub fn new_from_degrees(latitude: Degrees, longitude: Degrees, altitude: Metres) -> WGS84Coord {
        return WGS84Coord {
            lat: latitude * DEG_TO_RAD,
            lon: longitude * DEG_TO_RAD,
            alt: altitude,
        };
    }

    pub fn new_from_ecef_struct(ecef: &ECEF) -> WGS84Coord {
        Self::new_from_ecef(ecef.x, ecef.y, ecef.z)
    }
    pub fn new_from_ecef(x: f64, y: f64, z: f64) -> WGS84Coord {
        // Given an x, y, z in the ECEF frame, produces a WGS84Coord
        let p = (x.powf(2_f64) + y.powf(2_f64)).sqrt();
        let lambda_lon = f64::atan2(y, x);
        let mut phi_lat = f64::atan2(z, p * (1_f64 - EARTH_ECCENTRICITY_SQUARED));
        let mut precision = 1_f64;
        let mut nutation: f64 = 0_f64;
        while precision > 1e-12 {
            nutation = EARTH_MAJOR
                / (1_f64 - (EARTH_ECCENTRICITY_SQUARED * phi_lat.sin().powf(2_f64))).sqrt();
            let new_phi_lat = f64::atan2(
                z + (EARTH_ECCENTRICITY_SQUARED * nutation * phi_lat.sin()),
                p,
            );
            precision = (phi_lat - new_phi_lat).abs();
            phi_lat = new_phi_lat;
        }
        let h = p / phi_lat.cos() - nutation;
        return WGS84Coord::new_from_radians(phi_lat, lambda_lon, h);
    }

    pub fn new_from_ned(
        ned: &NED,
        rotation_matrix: &Vec<f64>,
        reference_point: &ECEF,
    ) -> WGS84Coord {
        let ned_matrix = vec![ned.n, ned.e, ned.d];

        let temp_matrix = transpose_times_vec(rotation_matrix, &ned_matrix);

        let ecef = ECEF {
            x: temp_matrix[0] + reference_point.x,
            y: temp_matrix[1] + reference_point.y,
            z: temp_matrix[2] + reference_point.z,
        };

        return WGS84Coord::new_from_ecef_struct(&ecef);
    }

    pub fn get_lat_radians(&self) -> f64 {
        return self.lat;
    }

    pub fn get_lat_degrees(&self) -> f64 {
        return self.lat.to_degrees();
    }

    pub fn get_lon_radians(&self) -> f64 {
        return self.lon;
    }

    pub fn get_lon_degrees(&self) -> f64 {
        return self.lon.to_degrees();
    }

    pub fn get_altitude(&self) -> f64 {
        return self.alt;
    }

    pub fn set_altitude(&mut self, alt: f64) {
        self.alt = alt
    }
}
