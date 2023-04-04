use chrono::NaiveDateTime;
use peroxide::prelude::{matrix, Matrix, Shape::Row};

use crate::{
    constants::{EARTH_ECCENTRICITY_SQUARED, EARTH_MAJOR},
    utils::{get_polar_motion_matrix, old_maybe_broken_jday},
};

use super::{
    ned::{NEDVel, NED},
    pef::{PEFVel, PEF},
    teme::{TEMEVel, TEME},
    wgs84::WGS84Coord,
};

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
        let pef_matrix = matrix(vec![pef.x, pef.y, pef.z], 3, 1, Row);
        let efec_matrix = polar_motion_matrix.t() * pef_matrix;
        return ECEF {
            x: efec_matrix[(0, 0)],
            y: efec_matrix[(1, 0)],
            z: efec_matrix[(2, 0)],
        };
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

    pub fn new_from_ned_rot(ned: &NED, rotation_matrix: &Matrix, reference_point: &ECEF) -> ECEF {
        let ned_matrix = matrix(vec![ned.n, ned.e, ned.d], 3, 1, Row);
        let ecef_matrix = rotation_matrix.t() * ned_matrix;

        let x = ecef_matrix[(0, 0)] + reference_point.x;
        let y = ecef_matrix[(1, 0)] + reference_point.y;
        let z = ecef_matrix[(2, 0)] + reference_point.z;

        return ECEF { x, y, z };
    }
}

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
        let velocity_pef_matrix = matrix(
            vec![velocity_pef.x_vel, velocity_pef.y_vel, velocity_pef.z_vel],
            3,
            1,
            Row,
        );
        let velocity_ecef = polar_motion_matrix.t() * velocity_pef_matrix;
        return ECEFVel {
            x_vel: velocity_ecef[(0, 0)],
            y_vel: velocity_ecef[(1, 0)],
            z_vel: velocity_ecef[(2, 0)],
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

    pub fn new_from_ned_rot(ned: &NEDVel, rotation_matrix: &Matrix) -> ECEFVel {
        let ned_matrix = matrix(vec![ned.n_vel, ned.e_vel, ned.d_vel], 3, 1, Row);
        let ecef_matrix = rotation_matrix.t() * ned_matrix;
        return ECEFVel {
            x_vel: ecef_matrix[(0, 0)],
            y_vel: ecef_matrix[(1, 0)],
            z_vel: ecef_matrix[(2, 0)],
        };
    }

    pub fn get_speed(&self) -> f64 {
        return (self.x_vel.powf(2_f64) + self.y_vel.powf(2_f64) + self.z_vel.powf(2_f64)).sqrt();
    }
}

pub fn generate_ecef_to_ned_matrix(initial_point: &WGS84Coord) -> Matrix {
    let sin_lon = initial_point.get_lon_radians().sin();
    let cos_lon = initial_point.get_lon_radians().cos();

    let sin_lat = initial_point.get_lat_radians().sin();
    let cos_lat = initial_point.get_lat_radians().cos();

    let rotation = matrix(
        vec![
            -sin_lat * cos_lon,
            -sin_lat * sin_lon,
            cos_lat,
            -sin_lon,
            cos_lon,
            0_f64,
            -cos_lat * cos_lon,
            -cos_lat * sin_lon,
            -sin_lat,
        ],
        3,
        3,
        Row,
    );
    return rotation;
}

pub fn construct_ecef_to_ned_jacobian(ecef_to_ned: &Matrix) -> Matrix {
    assert_eq!(ecef_to_ned.row, 3);
    assert_eq!(ecef_to_ned.col, 3);

    // let rot = ecef_to_ned.data.clone();
    let jacobian_vec = vec![
        ecef_to_ned.row(0),
        vec![0.0; 3],
        ecef_to_ned.row(1),
        vec![0.0; 3],
        ecef_to_ned.row(2),
        vec![0.0; 3],
        vec![0.0; 3],
        ecef_to_ned.row(0),
        vec![0.0; 3],
        ecef_to_ned.row(1),
        vec![0.0; 3],
        ecef_to_ned.row(2),
    ];

    let jacobian_matrix = matrix(jacobian_vec.concat(), 6, 6, Row);

    return jacobian_matrix;
}
