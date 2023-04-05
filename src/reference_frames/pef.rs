use chrono::NaiveDateTime;
use peroxide::prelude::{matrix, Shape::Row, SimplerLinearAlgebra};

use crate::utils::{
    get_pef_tod_matrix, get_polar_motion_matrix, julian_to_gmst, old_maybe_broken_jday,
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
        let teme_matrix = matrix(vec![teme.x, teme.y, teme.z], 3, 1, Row);
        let (jday, jfrac) = old_maybe_broken_jday(utc_time);
        let julian = jday + jfrac;
        let gmst = julian_to_gmst(julian);
        let pef_tod_matrix = get_pef_tod_matrix(gmst);
        let pef = pef_tod_matrix.transpose() * teme_matrix;
        // ! These matrices are not compatable, teme needs to be a column matrix - everything will be broken like this
        return PEF {
            x: (pef[(0, 0)]),
            y: (pef[(1, 0)]),
            z: (pef[(2, 0)]),
        };
        // ! accessing the pef matrix isn't going well
        // ! Looks like the index is out of trange
        // ! NOTE - (COLUMN, ROW)
    }

    pub fn new_from_ecef(ecef: &ECEF, utc_time: &NaiveDateTime) -> PEF {
        let ecef_matrix = matrix(vec![ecef.x, ecef.y, ecef.z], 3, 1, Row);
        let (jday, jfrac) = old_maybe_broken_jday(utc_time);
        let julian = jday + jfrac;
        let polar_motion_matrix = get_polar_motion_matrix(julian);
        println!(
            "polar motion m :{:?} | ecef matrix {:?}",
            polar_motion_matrix, ecef_matrix
        );
        let pef_matrix = polar_motion_matrix.t().inv() * ecef_matrix;

        return PEF {
            x: (pef_matrix[(0, 0)]),
            y: (pef_matrix[(1, 0)]),
            z: (pef_matrix[(2, 0)]),
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
        let teme_vel_matrix = matrix(
            vec![teme_vel.x_vel, teme_vel.y_vel, teme_vel.z_vel],
            3,
            1,
            Row,
        );
        let (jday, frac) = old_maybe_broken_jday(utc_time);
        let gmst = julian_to_gmst(jday + frac);
        let pef_tod_matrix = get_pef_tod_matrix(gmst);
        let omega_earth = 7.29211514670698e-05_f64 * (1.0_f64 - 0.002_f64 / 86400.0_f64);
        let velocity_pef_temp = pef_tod_matrix.t() * teme_vel_matrix;
        // println!("{:?}", velocity_pef_temp);
        let velocity_pef = velocity_pef_temp
            - matrix(
                vec![omega_earth * pef.y, omega_earth * pef.x, 0_f64],
                3,
                1,
                Row,
            );
        // println!("{:?}", velocity_pef);
        return PEFVel {
            x_vel: velocity_pef[(0, 0)],
            y_vel: velocity_pef[(1, 0)],
            z_vel: velocity_pef[(2, 0)],
        };
    }

    pub fn new_from_ecef_vel(ecef_vel: &ECEFVel, utc_time: &NaiveDateTime) -> PEFVel {
        let (jday, frac) = old_maybe_broken_jday(utc_time);
        let polar_motion_matrix = get_polar_motion_matrix(jday + frac);
        let velocity_ecef_matrix = matrix(
            vec![ecef_vel.x_vel, ecef_vel.y_vel, ecef_vel.z_vel],
            3,
            1,
            Row,
        );
        let velocity_pef_matrix = polar_motion_matrix.t().inv() * velocity_ecef_matrix;
        return PEFVel {
            x_vel: velocity_pef_matrix[(0, 0)],
            y_vel: velocity_pef_matrix[(1, 0)],
            z_vel: velocity_pef_matrix[(2, 0)],
        };
    }
}
