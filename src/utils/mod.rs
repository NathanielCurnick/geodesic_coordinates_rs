use std::{
    f64::consts::PI,
    ops::{Div, Neg},
};

use chrono::{Datelike, NaiveDateTime, Timelike};

pub fn get_pef_tod_matrix(gmst: f64) -> Vec<f64> {
    // Given a gmst time, returns a PEF TOD matrix
    // TOD (True of Date) is ANOTHER reference system - https://mycoordinates.org/tracking-satellite-footprints-on-earth%E2%80%99s-surface/
    let pef_tod = vec![
        gmst.cos(),
        gmst.sin().neg(),
        0_f64,
        gmst.sin(),
        gmst.cos(),
        0_f64,
        0_f64,
        0_f64,
        1_f64,
    ];

    return pef_tod;
}

pub fn get_polar_motion_matrix(julian: f64) -> Vec<f64> {
    // Given a julian day and fraction (in one number) produces a polar motion matrix
    // The Earf wobbles on its path - https://www.iers.org/IERS/EN/Science/EarthRotation/PolarMotion.html
    let mjd = julian - 2400000.5;
    let a = 2.0 * PI * (mjd - 57226_f64);
    // let b = a / 365.25_f64; // TODO
    let c = a / 435_f64;

    let xp = (0.1033_f64
        + 0.0494_f64 * a.cos()
        + 0.0482_f64 * a.sin()
        + 0.0297_f64 * c.cos()
        + 0.0307_f64 * c.sin())
        * 4.84813681e-6_f64;

    let yp = (0.3498_f64 + 0.0441_f64 * a.cos() - 0.0393_f64 * a.sin() + 0.0307_f64 * c.cos()
        - 0.0297_f64 * c.sin())
        * 4.84813681e-6_f64;

    let ma = vec![
        xp.cos(),
        0_f64,
        xp.sin().neg(),
        xp.sin() * yp.sin(),
        yp.cos(),
        xp.cos() * yp.sin(),
        xp.sin() * yp.cos(),
        yp.sin().neg(),
        xp.cos() * yp.cos(),
    ];

    return ma;
}

pub fn old_maybe_broken_jday(utc_time: &NaiveDateTime) -> (f64, f64) {
    let mut year = utc_time.year() as f64;
    let mut month = utc_time.month() as f64;
    let day = utc_time.day() as f64;
    let hour = utc_time.hour() as f64;
    let minute = utc_time.minute() as f64;
    let second = utc_time.second() as f64;

    let term_1 = 367_f64 * year;
    // TODO: Why divide by 1.0?
    let term_2 = (7_f64 * ((month + 9_f64).div(12_f64).floor()) * 0.25)
        .div(1.0_f64)
        .floor();
    let term_3 = (275_f64 * (month / 9_f64)).div(1.0_f64).floor();
    let term_4 = day + 1721013.5;

    let jd = term_1 - term_2 + term_3 + term_4;

    let fr = (second + minute * 60_f64 + hour * 3600_f64) / 86400_f64;

    // Julian day and Julian day fraction
    // e.g. jday 100 and jfrac 0.5 would be 12 hours into the 100th Julian day

    return (jd, fr);
}

pub fn jday(utc_time: &NaiveDateTime) -> f64 {
    let mut year = utc_time.year() as f64;
    let mut month = utc_time.month() as f64;
    let mut day = utc_time.day() as f64;
    let hour = utc_time.hour() as f64;
    let minute = utc_time.minute() as f64;
    let second = utc_time.second() as f64;

    // Add decimals to day of month. Divide by hours in day, minutes in day, seconds in day.
    day += (hour / 24.0) + (minute / 1440.0) + (second / 86400.0);

    if month <= 2.0 {
        year -= 1.0;
        month += 12.0;
    }

    let a = (year / 100.0).floor();

    let b = 2.0 - a + (a / 4.0).floor();

    let jd =
        (365.25 * (year + 4716.0)).floor() + (30.6001 * (month + 1.0)).floor() + day + b - 1524.5;

    return jd;
}

pub fn julian_to_gmst(jday: f64) -> f64 {
    let tut1 = (jday - 2451545.0_f64) / 36525.0_f64;
    let mut gmst = -6.2e-6_f64 * tut1.powf(3.0)
        + 0.093104_f64 * tut1.powf(2_f64)
        + (876600.0_f64 * 3600.0_f64 + 8640184.812866_f64) * tut1
        + 67310.54841_f64;

    if gmst < 0_f64 {
        gmst += 2_f64 * PI;
    }

    return gmst;
}

pub fn transpose_times_vec(m: &Vec<f64>, v: &Vec<f64>) -> Vec<f64> {
    // Takes a 3X3 matrix and multiplies its transpose by a 3 element vector

    let x = m[0] * v[0] + m[3] * v[1] + m[6] * v[2];
    let y = m[1] * v[0] + m[4] * v[1] + m[7] * v[2];
    let z = m[2] * v[0] + m[5] * v[1] + m[8] * v[2];

    return vec![x, y, z];
}

pub fn matrix_times_vec(m: &Vec<f64>, v: &Vec<f64>) -> Vec<f64> {
    let x = m[0] * v[0] + m[1] * v[1] + m[2] * v[2];
    let y = m[3] * v[0] + m[4] * v[1] + m[5] * v[2];
    let z = m[6] * v[0] + m[7] * v[1] + m[8] * v[2];

    return vec![x, y, z];
}

pub fn invert_matrix(m: &Vec<f64>) -> Vec<f64> {
    // invert 3x3 matrix
    // TODO: check length 9

    let a = m[0];
    let b = m[1];
    let c = m[2];
    let d = m[3];
    let e = m[4];
    let f = m[5];
    let g = m[6];
    let h = m[7];
    let i = m[8];

    let determinant = a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);

    // TODO: check determinant exists

    let inv_det = 1.0 / determinant;

    let mut matrix = Vec::with_capacity(9);

    matrix[0] = (e * i - f * h) * inv_det;
    matrix[1] = (c * h - b * i) * inv_det;
    matrix[2] = (b * f - c * e) * inv_det;
    matrix[3] = (f * g - d * i) * inv_det;
    matrix[4] = (a * i - c * g) * inv_det;
    matrix[5] = (c * d - a * f) * inv_det;
    matrix[6] = (d * h - e * g) * inv_det;
    matrix[7] = (b * g - a * h) * inv_det;
    matrix[8] = (a * e - b * d) * inv_det;

    return matrix;
}

pub fn transpose(m: &Vec<f64>) -> Vec<f64> {
    // Transpose 3x3 matrix
    return vec![m[0], m[3], m[6], m[1], m[4], m[7], m[2], m[5], m[8]];
}
