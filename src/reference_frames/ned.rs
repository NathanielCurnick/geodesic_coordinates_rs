use crate::utils::matrix_times_vec;

use super::{
    ecef::{ECEFVel, ECEF},
    enu::ENU,
};

#[derive(Debug, Clone)]
pub struct NED {
    // Wrapper for NED coordinates
    pub n: f64,
    pub e: f64,
    pub d: f64,
}

impl NED {
    pub fn default() -> NED {
        return NED {
            n: 0_f64,
            e: 0_f64,
            d: 0_f64,
        };
    }
    pub fn new_from_raw(n: f64, e: f64, d: f64) -> NED {
        return NED { n, e, d };
    }
    pub fn new_from_vec(ned_vec: Vec<f64>) -> Self {
        return NED::new_from_raw(ned_vec[0], ned_vec[1], ned_vec[2]);
    }

    pub fn new_from_enu(enu: &ENU) -> NED {
        // Given a reference to a ENU returns a NED
        // Z axis is flipped; x and y swap
        return NED {
            n: enu.e,
            e: enu.n,
            d: -enu.u,
        };
    }

    pub fn new_from_ecef_rot(
        ecef: &ECEF,
        rotation_matrix: &Vec<f64>,
        reference_point: &ECEF,
    ) -> NED {
        let x = ecef.x - reference_point.x;
        let y = ecef.y - reference_point.y;
        let z = ecef.z - reference_point.z;
        let ecef_vec = vec![x, y, z];

        let ned_vec = matrix_times_vec(rotation_matrix, &ecef_vec);

        return NED {
            n: ned_vec[0],
            e: ned_vec[1],
            d: ned_vec[2],
        };
    }
    pub fn distance_between(&self, other: &NED) -> f64 {
        return ((self.n - other.n).powf(2_f64)
            + (self.e - other.e).powf(2_f64)
            + (self.d - other.d).powf(2_f64))
        .sqrt();
    }

    pub fn heading_between_2d(&self, other: &NED) -> f64 {
        return (other.n - self.n).atan2(other.e - self.e);
    }

    pub fn to_vec(&self) -> Vec<f64> {
        return vec![self.n, self.e, self.d];
    }

    pub fn distance_between_points(&self, other_point: &NED) -> f64 {
        return ((self.n - other_point.n).powf(2.)
            + (self.e - other_point.e).powf(2.)
            + (self.d - other_point.d).powf(2.))
        .sqrt();
    }

    // pub fn new_from_ecef_ref(ecef: &ECEF, ref_point: &WGS84Coord) -> NED {} // TODO
}

#[derive(Debug, Clone)]
pub struct NEDVel {
    pub n_vel: f64,
    pub e_vel: f64,
    pub d_vel: f64,
}

impl NEDVel {
    /// Speed from a to b assuming no acceleration
    pub fn linear_speed(a: &NED, b: &NED, dt: f64) -> NEDVel {
        return NEDVel {
            n_vel: (b.n - a.n) / dt,
            e_vel: (b.e - a.e) / dt,
            d_vel: (b.d - a.d) / dt,
        };
    }

    pub fn new_from_ecef_rot(ecef: &ECEFVel, rotation_matrix: &Vec<f64>) -> NEDVel {
        let ecef_vec = vec![ecef.x_vel, ecef.y_vel, ecef.z_vel];

        let ned_vec = matrix_times_vec(rotation_matrix, &ecef_vec);

        return NEDVel {
            n_vel: ned_vec[0],
            e_vel: ned_vec[1],
            d_vel: ned_vec[2],
        };
    }
}

#[derive(Debug, Clone)]
pub struct NEDAccel {
    pub n_accel: f64,
    pub e_accel: f64,
    pub d_accel: f64,
}

impl NEDAccel {
    pub fn magnitude(&self) -> f64 {
        return (self.n_accel.powf(2_f64) + self.e_accel.powf(2_f64) + self.d_accel.powf(2_f64))
            .sqrt();
    }
}
