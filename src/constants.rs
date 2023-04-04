use std::f64::consts::PI;

use crate::types::{Kilometres, Metres, Minutes};

pub const DEG_TO_RAD: f64 = PI / 180.0;
pub const EARTH_FLATTENING: f64 = 1.0 / 298.257223563;
pub const EARTH_FLATTENING_SQUARED: f64 = EARTH_FLATTENING * EARTH_FLATTENING;
pub const EARTH_MAJOR: Metres = 6378137.0;
pub const EARTH_MINOR: Metres = EARTH_MAJOR - (EARTH_FLATTENING * EARTH_MAJOR);
pub const EARTH_MEAN_RADIUS: Metres = 6371e3;
pub const EARTH_THIRD_FLATTENING: f64 = EARTH_FLATTENING / (2. - EARTH_FLATTENING);
pub const EARTH_ECCENTRICITY_SQUARED: f64 = EARTH_FLATTENING * (2. - EARTH_FLATTENING);
pub const EARTH_SECOND_ECCENTRICITY_SQUARED: f64 =
    EARTH_ECCENTRICITY_SQUARED / (1. - EARTH_ECCENTRICITY_SQUARED);

pub const KEPLER_CALC_ACCURACY: f64 = DEG_TO_RAD / 3600.; // Arcseconds

pub const SUN_SEMI_MAJOR_AXIS: Kilometres = 149598845.0;

pub const SUN_RADIUS: Kilometres = 695008.0;
pub const EARTH_RADIUS: Kilometres = 6378.16; // Maybe make this more accurate?

pub const MINUTES_PER_DAY: Minutes = 24. * 60.;

pub const GM: f64 = 398600.0; // Kilometers^3/seconds^2 .
pub const SIDEREAL_SOLAR: f64 = 1.0027379093;

pub const LTLIM: f64 = 14.5;

pub const GELIM: f64 = 15.5;

pub const DEFAULT_PRESSURE: f64 = 1010.;

pub const DEFAULT_TEMP: f64 = 15.0;
pub const GRAVITY: f64 = 9.81_f64;
