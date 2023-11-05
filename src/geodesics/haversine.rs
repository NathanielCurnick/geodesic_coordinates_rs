use std::f64::consts::PI;

use crate::{
    constants::EARTH_MEAN_RADIUS,
    types::{DistBearing, LocBearing, Metres, Radians},
};

use wasm_bindgen::prelude::*;

#[wasm_bindgen]
pub fn haversine_distance_and_bearing(
    lat1: Radians,
    lon1: Radians,
    lat2: Radians,
    lon2: Radians,
) -> DistBearing {
    // Calc distance
    let delta_lon = lon2 - lon1;
    let lat1_cos = lat1.cos();
    let lat2_cos = lat2.cos();

    let a = ((lat2 - lat1) / 2.).sin().powf(2.)
        + (lat1_cos * lat2_cos) * (delta_lon / 2.).sin().powf(2.);
    let c = 2. * a.sqrt().atan2(1. - a);

    let distance = EARTH_MEAN_RADIUS * c;

    // Calc bearing
    let bearing = calc_initial_bearing(lat1_cos, lat2_cos, lat1.sin(), lat2.sin(), delta_lon);

    return DistBearing { distance, bearing };
}

#[wasm_bindgen]
pub fn haversine_location_and_bearing(
    lat1: Radians,
    lon1: Radians,
    bearing: Radians,
    distance: Metres,
) -> LocBearing {
    let lat1_cos = lat1.cos();
    let lat1_sin = lat1.sin();
    let dr_cos = (distance / EARTH_MEAN_RADIUS).cos();
    let dr_sin = (distance / EARTH_MEAN_RADIUS).sin();

    let lat2 = ((lat1_sin * dr_cos) + (lat1_cos * dr_sin * bearing.cos())).asin();
    let lon2 = lon1 + (bearing.sin() * dr_sin * lat1_cos).atan2(dr_cos - lat1_sin * lat2.sin());

    let mut bearing = calc_initial_bearing(lat2.cos(), lat1_cos, lat2.sin(), lat1_sin, lon1 - lon2);

    bearing = (bearing - PI) % (2.0 * PI);

    return LocBearing {
        lat: lat2,
        lon: lon2,
        bearing,
    };
}

fn calc_initial_bearing(
    lat1_cos: Radians,
    lat2_cos: Radians,
    lat1_sin: Radians,
    lat2_sin: Radians,
    delta_lon: Radians,
) -> f64 {
    let y = delta_lon.sin() * lat2_cos;
    let x = (lat1_cos * lat2_sin) - (lat1_sin * lat2_cos * delta_lon.cos());

    let theta = y.atan2(x);

    return (theta + (2. * PI)) % (2. * PI);
}
