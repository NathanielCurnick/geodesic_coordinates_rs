use chrono::NaiveDateTime;

use crate::utils::{
    get_pef_tod_matrix, get_polar_motion_matrix, invert_matrix, julian_to_gmst, matrix_times_vec,
    old_maybe_broken_jday, transpose, transpose_times_vec,
};

use super::{
    ecef::{ECEFVel, ECEF},
    teme::{TEMEVel, TEME},
};

#[derive(Debug, Clone)]
pub struct PEF {
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl PEF {
    // fn new_from_ecef(ecef: &ECEF) -> PEF {}
    pub fn new_from_teme(teme: &TEME, utc_time: &NaiveDateTime) -> PEF {
        let teme_matrix = vec![teme.x, teme.y, teme.z];
        let (jday, jfrac) = old_maybe_broken_jday(utc_time);
        let julian = jday + jfrac;
        let gmst = julian_to_gmst(julian);
        let pef_tod_matrix = get_pef_tod_matrix(gmst);

        let pef = transpose_times_vec(&pef_tod_matrix, &teme_matrix);
        return PEF {
            x: (pef[0]),
            y: (pef[1]),
            z: (pef[2]),
        };
    }

    pub fn new_from_ecef(ecef: &ECEF, utc_time: &NaiveDateTime) -> PEF {
        let ecef_matrix = vec![ecef.x, ecef.y, ecef.z];
        let (jday, jfrac) = old_maybe_broken_jday(utc_time);
        let julian = jday + jfrac;
        let polar_motion_matrix = get_polar_motion_matrix(julian);

        let trans = transpose(&polar_motion_matrix);
        let inv = invert_matrix(&trans);

        let pef = matrix_times_vec(&inv, &ecef_matrix);

        return PEF {
            x: (pef[0]),
            y: (pef[1]),
            z: (pef[2]),
        };
    }
}

#[derive(Debug, Clone)]
pub struct PEFVel {
    // Wrapper for the PEF coordinate system
    pub x_vel: f64,
    pub y_vel: f64,
    pub z_vel: f64,
}

impl PEFVel {
    pub fn new_from_teme_vel(teme: &TEME, teme_vel: &TEMEVel, utc_time: &NaiveDateTime) -> PEFVel {
        // Given a reference to a TEME coordinate, a reference to a TEME velocity and a reference to a NaiveDateTime produces a PEF velocity
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        let pef = PEF::new_from_teme(teme, utc_time);
        let teme_vel_matrix = vec![teme_vel.x_vel, teme_vel.y_vel, teme_vel.z_vel];

        let (jday, frac) = old_maybe_broken_jday(utc_time);
        let gmst = julian_to_gmst(jday + frac);
        let pef_tod_matrix = get_pef_tod_matrix(gmst);
        let omega_earth = 7.29211514670698e-05_f64 * (1.0_f64 - 0.002_f64 / 86400.0_f64);

        let velocity_pef_temp = transpose_times_vec(&pef_tod_matrix, &teme_vel_matrix);

        return PEFVel {
            x_vel: velocity_pef_temp[0] - omega_earth * pef.y,
            y_vel: velocity_pef_temp[1] - omega_earth * pef.x,
            z_vel: velocity_pef_temp[2],
        };
    }

    pub fn new_from_ecef_vel(ecef_vel: &ECEFVel, utc_time: &NaiveDateTime) -> PEFVel {
        let (jday, frac) = old_maybe_broken_jday(utc_time);
        let polar_motion_matrix = get_polar_motion_matrix(jday + frac);
        let velocity_ecef_matrix = vec![ecef_vel.x_vel, ecef_vel.y_vel, ecef_vel.z_vel];

        let trans = transpose(&polar_motion_matrix);
        let inv = invert_matrix(&trans);

        let velocity_pef_matrix = matrix_times_vec(&inv, &velocity_ecef_matrix);

        return PEFVel {
            x_vel: velocity_pef_matrix[0],
            y_vel: velocity_pef_matrix[1],
            z_vel: velocity_pef_matrix[2],
        };
    }
}
