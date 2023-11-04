use chrono::NaiveDateTime;

use crate::{
    constants::{EARTH_ECCENTRICITY_SQUARED, EARTH_MAJOR},
    utils::{get_polar_motion_matrix, old_maybe_broken_jday, transpose_times_vec},
};

use super::{
    ned::{NEDVel, NED},
    pef::{PEFVel, PEF},
    teme::{TEMEVel, TEME},
    wgs84::WGS84Coord,
};

#[derive(Debug, Clone)]
pub struct ECEF {
    // Wrapper for ECEF coordinates
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl ECEF {
    pub fn default() -> Self {
        Self {
            x: 0.0,
            y: 0.0,
            z: 0.0,
        }
    }
    pub fn new_from_wgs84(wgs84: &WGS84Coord) -> ECEF {
        // Given a reference to a WGS84Coord produces an ECEF coord
        let lat = wgs84.get_lat_radians();
        let lon = wgs84.get_lon_radians();
        let alt = wgs84.get_altitude();

        let lat_cos = lat.cos();
        let lat_sin = lat.sin();

        let nutation =
            EARTH_MAJOR / (1_f64 - (EARTH_ECCENTRICITY_SQUARED * lat_sin.powf(2_f64))).sqrt();
        let x = (nutation + alt) * lat_cos * lon.cos();
        let y = (nutation + alt) * lat_cos * lon.sin();
        let z = (((1_f64 - EARTH_ECCENTRICITY_SQUARED) * nutation) + alt) * lat_sin;

        return ECEF { x, y, z };
    }

    pub fn new_from_pef(pef: &PEF, utc_time: &NaiveDateTime) -> ECEF {
        // Given a reference to a PEF coordinate and a reference to a NaiveDateTime produces an ECEF
        // Since PEF and ECEF rotate with respect to one another, a time is necessary
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        let (jday, jfrac) = old_maybe_broken_jday(utc_time);
        let julian = jday + jfrac;
        let polar_motion_matrix = get_polar_motion_matrix(julian);

        // polar_motion_matrix^T * pef

        let x = polar_motion_matrix[0] * pef.x
            + polar_motion_matrix[3] * pef.y
            + polar_motion_matrix[6] * pef.z;
        let y = polar_motion_matrix[1] * pef.x
            + polar_motion_matrix[4] * pef.y
            + polar_motion_matrix[7] * pef.z;
        let z = polar_motion_matrix[2] * pef.x
            + polar_motion_matrix[5] * pef.y
            + polar_motion_matrix[8] * pef.z;

        return ECEF { x, y, z };
    }

    // ! Broken
    // pub fn new_from_enu(enu: &ENU, wgs84_ref: &WGS84Coord) -> ECEF {
    //     // wgs84_ref is a constant reference point - MUST NOT VARY OVER THE COURSE OF THE PROGRAM
    //     let e = enu.e;
    //     let n = enu.n;
    //     let u = enu.u;

    //     let lambda = wgs84_ref.lat;
    //     let phi = wgs84_ref.lon;

    //     let sin_lambda = lambda.sin();
    //     let cos_lambda = lambda.cos();

    //     let sin_phi = phi.sin();
    //     let cos_phi = phi.cos();

    //     let x = -sin_lambda * e + cos_lambda * n;
    //     let y = -cos_lambda * sin_phi * e
    //         + -sin_lambda * sin_phi * n
    //         + cos_phi * u;
    //     let z =
    //         cos_lambda * cos_phi * e + sin_lambda * cos_phi * n + sin_phi * u;

    //     return ECEF { x, y, z };
    // }

    // ! Broken!
    // pub fn new_from_ned(ned: &NED, wgs84_ref: &WGS84Coord) -> ECEF {
    //     // wgs84_ref is a constant reference point - MUST NOT VARY OVER THE COURSE OF THE PROGRAM
    //     let enu = ENU::new_from_ned(ned);
    //     let ecef = ECEF::new_from_enu(&enu, wgs84_ref);
    //     return ecef;
    // }

    pub fn new_from_ned_rot(ned: &NED, rotation_matrix: &Vec<f64>, reference_point: &ECEF) -> ECEF {
        // rotation_matrix^T * ned

        let ned_v = vec![ned.n, ned.e, ned.d];

        let tmp = transpose_times_vec(&rotation_matrix, &ned_v);

        let x = tmp[0] + reference_point.x;
        let y = tmp[1] + reference_point.y;
        let z = tmp[2] + reference_point.z;

        return ECEF { x, y, z };
    }
}

#[derive(Debug, Clone)]
pub struct ECEFVel {
    // Wrapper for ECEF velocities
    pub x_vel: f64,
    pub y_vel: f64,
    pub z_vel: f64,
}

impl ECEFVel {
    pub fn new_from_teme_vel(teme: &TEME, teme_vel: &TEMEVel, utc_time: &NaiveDateTime) -> ECEFVel {
        // Given a reference to a TEME coordinate, a reference to a TEME velocity and a reference to a NaiveDateTime produces an ECEF velocity
        // ECEF and TEME rotate with respect to one another so both a location and a velocity are needed
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        let (jday, frac) = old_maybe_broken_jday(utc_time);
        let polar_motion_matrix = get_polar_motion_matrix(jday + frac);
        let velocity_pef = PEFVel::new_from_teme_vel(teme, teme_vel, utc_time);
        let velocity_pef = vec![velocity_pef.x_vel, velocity_pef.y_vel, velocity_pef.z_vel];

        let v_ecef = transpose_times_vec(&polar_motion_matrix, &velocity_pef);

        return ECEFVel {
            x_vel: v_ecef[0],
            y_vel: v_ecef[1],
            z_vel: v_ecef[2],
        };
    }

    pub fn new_from_wgs84(
        point1: &WGS84Coord,
        point1_time: &NaiveDateTime,
        point2: &WGS84Coord,
        point2_time: &NaiveDateTime,
    ) -> ECEFVel {
        // Given a reference to two WGS84 points and a reference to a NaiveDateTime corresponding to both those locations produces an ECEF velocity
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        // Note this assumes no acceleration between the points
        let ecef_1 = ECEF::new_from_wgs84(point1);
        let ecef_2 = ECEF::new_from_wgs84(point2);
        let seconds_diff = (point2_time.timestamp() - point1_time.timestamp()) as f64;
        let x_vel = (ecef_2.x - ecef_1.x) / seconds_diff;
        let y_vel = (ecef_2.y - ecef_1.y) / seconds_diff;
        let z_vel = (ecef_2.z - ecef_1.z) / seconds_diff;

        return ECEFVel {
            x_vel,
            y_vel,
            z_vel,
        };
    }

    pub fn new_from_wgs84_timediff(
        point1: &WGS84Coord,
        point2: &WGS84Coord,
        seconds_diff: f64,
    ) -> ECEFVel {
        let ecef_1 = ECEF::new_from_wgs84(point1);
        let ecef_2 = ECEF::new_from_wgs84(point2);
        let x_vel = (ecef_2.x - ecef_1.x) / seconds_diff;
        let y_vel = (ecef_2.y - ecef_1.y) / seconds_diff;
        let z_vel = (ecef_2.z - ecef_1.z) / seconds_diff;

        return ECEFVel {
            x_vel,
            y_vel,
            z_vel,
        };
    }

    pub fn new_from_ecef(
        point1: &ECEF,
        point1_time: &NaiveDateTime,
        point2: &ECEF,
        point2_time: &NaiveDateTime,
    ) -> ECEFVel {
        // Given a reference to two ECEF points and a reference to a NaiveDateTime corresponding to both those locations produces an ECEF velocity
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        let seconds_diff = (point2_time.timestamp() - point1_time.timestamp()) as f64;
        let x_vel = (point2.x - point1.x) / seconds_diff;
        let y_vel = (point2.y - point1.y) / seconds_diff;
        let z_vel = (point2.z - point1.z) / seconds_diff;

        return ECEFVel {
            x_vel,
            y_vel,
            z_vel,
        };
    }

    pub fn new_from_ecef_timestamp(point1: &ECEF, point2: &ECEF, seconds_diff: f64) -> ECEFVel {
        // Given a reference to two ECEF points and a reference to a NaiveDateTime corresponding to both those locations produces an ECEF velocity
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        let x_vel = (point2.x - point1.x) / seconds_diff;
        let y_vel = (point2.y - point1.y) / seconds_diff;
        let z_vel = (point2.z - point1.z) / seconds_diff;

        return ECEFVel {
            x_vel,
            y_vel,
            z_vel,
        };
    }

    pub fn new_from_ned_rot(ned: &NEDVel, rotation_matrix: &Vec<f64>) -> ECEFVel {
        // rotation_matrix^T * ned_vel

        let ned_v = vec![ned.n_vel, ned.e_vel, ned.d_vel];

        let tmp = transpose_times_vec(rotation_matrix, &ned_v);

        return ECEFVel {
            x_vel: tmp[0],
            y_vel: tmp[1],
            z_vel: tmp[2],
        };
    }

    pub fn get_speed(&self) -> f64 {
        return (self.x_vel.powf(2_f64) + self.y_vel.powf(2_f64) + self.z_vel.powf(2_f64)).sqrt();
    }
}

pub fn generate_ecef_to_ned_matrix(initial_point: &WGS84Coord) -> Vec<f64> {
    let sin_lon = initial_point.get_lon_radians().sin();
    let cos_lon = initial_point.get_lon_radians().cos();

    let sin_lat = initial_point.get_lat_radians().sin();
    let cos_lat = initial_point.get_lat_radians().cos();

    let rotation = vec![
        -sin_lat * cos_lon,
        -sin_lat * sin_lon,
        cos_lat,
        -sin_lon,
        cos_lon,
        0_f64,
        -cos_lat * cos_lon,
        -cos_lat * sin_lon,
        -sin_lat,
    ];

    return rotation;
}

pub fn construct_ecef_to_ned_jacobian(ecef_to_ned: &Vec<f64>) -> Vec<f64> {
    let jacobian_vec = vec![
        vec![ecef_to_ned[0], ecef_to_ned[1], ecef_to_ned[2]],
        vec![0.0; 3],
        vec![ecef_to_ned[3], ecef_to_ned[4], ecef_to_ned[5]],
        vec![0.0; 3],
        vec![ecef_to_ned[6], ecef_to_ned[7], ecef_to_ned[8]],
        vec![0.0; 3],
        vec![0.0; 3],
        vec![ecef_to_ned[0], ecef_to_ned[1], ecef_to_ned[2]],
        vec![0.0; 3],
        vec![ecef_to_ned[3], ecef_to_ned[4], ecef_to_ned[5]],
        vec![0.0; 3],
        vec![ecef_to_ned[6], ecef_to_ned[7], ecef_to_ned[8]],
    ];

    let jacobian_matrix = jacobian_vec.concat();

    return jacobian_matrix;
}
