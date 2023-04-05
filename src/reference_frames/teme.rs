use chrono::NaiveDateTime;
use peroxide::prelude::{matrix, Shape::Row, SimplerLinearAlgebra};

use crate::utils::{get_pef_tod_matrix, julian_to_gmst, old_maybe_broken_jday};

use super::pef::{PEFVel, PEF};

#[derive(Debug, Clone)]
pub struct TEME {
    // Wrapper for TEME coordinates
    // This is the coordinate system used by TLEs
    pub x: f64,
    pub y: f64,
    pub z: f64,
}

impl TEME {
    pub fn new_from_pef(pef: &PEF, utc_time: &NaiveDateTime) -> TEME {
        // Given a reference to a PEF and a reference to a NaiveDateTime produces a TEME location
        // Be careful of the NaiveDateTime time zones - I am assuming this is using UNIX seconds
        // In general, you are best of making all NaiveDateTimes from UNIX timestamps
        let pef_matrix = matrix(vec![pef.x, pef.y, pef.z], 3, 1, Row);
        let (jday, jfrac) = old_maybe_broken_jday(utc_time);
        let julian = jday + jfrac;
        let gmst = julian_to_gmst(julian);
        let pef_tod_matrix = get_pef_tod_matrix(gmst);
        let teme = pef_tod_matrix.transpose().inv() * pef_matrix;
        return TEME {
            x: (teme[(0, 0)]),
            y: (teme[(1, 0)]),
            z: (teme[(2, 0)]),
        };
    }
}

#[derive(Debug, Clone)]
pub struct TEMEVel {
    // Wrapper for TEME velocities
    pub x_vel: f64,
    pub y_vel: f64,
    pub z_vel: f64,
}
impl TEMEVel {
    pub fn new_from_pef_vel(pef: &PEF, pef_vel: &PEFVel, utc_time: &NaiveDateTime) -> TEMEVel {
        let pef_vel_matrix = matrix(vec![pef_vel.x_vel, pef_vel.y_vel, pef_vel.z_vel], 3, 1, Row);
        let (jday, frac) = old_maybe_broken_jday(utc_time);
        let gmst = julian_to_gmst(jday + frac);
        let pef_tod_matrix = get_pef_tod_matrix(gmst);
        let omega_earth = 7.29211514670698e-05_f64 * (1.0_f64 - 0.002_f64 / 86400.0_f64);
        let velocity_teme_temp = pef_tod_matrix.t().inv() * pef_vel_matrix;
        let velocity_teme = velocity_teme_temp
            + matrix(
                vec![omega_earth * pef.y, omega_earth * pef.x, 0_f64],
                3,
                1,
                Row,
            );
        return TEMEVel {
            x_vel: velocity_teme[(0, 0)],
            y_vel: velocity_teme[(1, 0)],
            z_vel: velocity_teme[(2, 0)],
        };
    }
}
