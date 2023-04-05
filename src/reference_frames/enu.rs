use super::{ecef::ECEF, ned::NED, wgs84::WGS84Coord};

#[derive(Debug, Clone)]
pub struct ENU {
    pub e: f64,
    pub n: f64,
    pub u: f64,
}

impl ENU {
    pub fn new_from_raw(e: f64, n: f64, u: f64) -> ENU {
        return ENU { e, n, u };
    }

    pub fn new_from_ecef(ecef: &ECEF) -> ENU {
        let x = ecef.x;
        let y = ecef.y;
        let z = ecef.z;
        let wgs84 = WGS84Coord::new_from_ecef_struct(ecef);
        let lambda = wgs84.get_lat_radians();
        let phi = wgs84.get_lon_radians();

        let sin_lambda = lambda.sin();
        let cos_lambda = lambda.cos();

        let sin_phi = phi.sin();
        let cos_phi = phi.cos();

        let e = -sin_lambda * x + cos_lambda * y;
        let n = -cos_lambda * sin_phi * x + -sin_lambda * sin_phi * y + cos_phi * z;
        let u = cos_lambda * cos_phi * x + sin_lambda * cos_phi * y + sin_phi * z;

        return ENU { e, n, u };
    }

    pub fn new_from_ned(ned: &NED) -> ENU {
        return ENU {
            e: ned.e,
            n: ned.n,
            u: -ned.d,
        };
    }
}
