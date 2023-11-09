use crate::{
    constants::{
        EARTH_ECCENTRICITY_SQUARED, EARTH_FLATTENING, EARTH_MAJOR, EARTH_MINOR,
        EARTH_SECOND_ECCENTRICITY_SQUARED, EARTH_THIRD_FLATTENING,
    },
    types::{DistBearing, LocBearing, Metres, Radians},
};

use wasm_bindgen::prelude::*;

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn karney_location_and_bearing(
    lat1: Radians,
    lon1: Radians,
    bearing: Radians,
    distance: Metres,
) -> LocBearing {
    return location_and_bearing(lat1, lon1, bearing, distance);
}

pub fn location_and_bearing(
    lat1: Radians,
    lon1: Radians,
    bearing: Radians,
    distance: Metres,
) -> LocBearing {
    let beta: f64 = ((1. - EARTH_FLATTENING) * lat1.tan()).atan();

    let sin_alpha = bearing.sin();
    let cos_alpha = bearing.cos();
    let sin_beta = beta.sin();
    let cos_beta = beta.cos();

    let alpha_0 = (sin_alpha * beta.cos())
        .atan2((cos_alpha.powf(2.) + (beta.sin().powf(2.) * sin_alpha.powf(2.))).powf(0.5));

    let sin_alpha_0 = alpha_0.sin();
    let cos_alpha_0 = alpha_0.cos();

    let sigma = sin_beta.atan2(cos_beta * cos_alpha);

    let omega = (sin_alpha_0 * sigma.sin()).atan2(sigma.cos());

    let k_squared = EARTH_SECOND_ECCENTRICITY_SQUARED * cos_alpha_0.powf(2.);

    let one_plus_k_squared_square_rooted = (1. + k_squared).sqrt();

    let epsilon = (one_plus_k_squared_square_rooted - 1.) / (one_plus_k_squared_square_rooted + 1.);

    let epsilon_pow2 = epsilon.powf(2.);
    let epsilon_pow3 = epsilon.powf(3.);
    let epsilon_pow4 = epsilon.powf(4.);
    let epsilon_pow5 = epsilon.powf(5.);
    let epsilon_pow6 = epsilon.powf(6.);
    let epsilon_pow7 = epsilon.powf(7.);
    let epsilon_pow8 = epsilon.powf(8.);
    let epsilon_pow9 = epsilon.powf(9.);
    let epsilon_pow10 = epsilon.powf(10.);

    let epsilon_pow_tup = (
        epsilon,
        epsilon_pow2,
        epsilon_pow3,
        epsilon_pow4,
        epsilon_pow5,
        epsilon_pow6,
        epsilon_pow7,
        epsilon_pow8,
        epsilon_pow9,
        epsilon_pow10,
    );

    let a_1 = (1.
        + (1. / 4. * epsilon_pow2)
        + (1. / 64. * epsilon_pow4)
        + (1. / 256. * epsilon_pow6)
        + (25. / 16384. * epsilon_pow8)
        + (49. / 65536. * epsilon_pow10))
        / (1. - epsilon);

    let s_one = i_one_fourier_series(epsilon_pow_tup, sigma, a_1) * EARTH_MINOR;

    let s_two = s_one + distance;

    let tau_two = s_two / (EARTH_MINOR * a_1);

    let sigma_two = sigma_two_fourier_series(epsilon_pow_tup, tau_two);

    let alpha_two = sin_alpha_0.atan2(cos_alpha_0 * sigma_two.cos());

    let sin_sigma_two = sigma_two.sin();
    let cos_sigma_two = sigma_two.cos();

    let beta_two = (cos_alpha_0 * sin_sigma_two)
        .atan2(((cos_alpha_0.powf(2.) * cos_sigma_two.powf(2.)) + sin_alpha_0.powf(2.)).sqrt());

    let omega_two = (sin_alpha_0 * sin_sigma_two).atan2(cos_sigma_two);

    let a_3 = 1.
        - (1. / 2. - 1. / 2. * EARTH_THIRD_FLATTENING) * epsilon
        - (1. / 4. + 1. / 8. * EARTH_THIRD_FLATTENING - 3. / 8. * EARTH_THIRD_FLATTENING.powf(2.))
            * epsilon_pow2
        - (1. / 16.
            + 3. / 16. * EARTH_THIRD_FLATTENING
            + 1. / 16. * EARTH_THIRD_FLATTENING.powf(2.)
            - 5. / 16. * EARTH_THIRD_FLATTENING.powf(3.))
            * epsilon_pow3
        - (3. / 64.
            + 1. / 32. * EARTH_THIRD_FLATTENING
            + 5. / 32. * EARTH_THIRD_FLATTENING.powf(2.)
            + 5. / 128. * EARTH_THIRD_FLATTENING.powf(3.)
            - 35. / 128. * EARTH_THIRD_FLATTENING.powf(4.))
            * epsilon_pow4
        - (3. / 128.
            + 5. / 128. * EARTH_THIRD_FLATTENING
            + 5. / 256. * EARTH_THIRD_FLATTENING.powf(2.)
            + 35. / 256. * EARTH_THIRD_FLATTENING.powf(3.)
            + 7. / 256. * EARTH_THIRD_FLATTENING.powf(4.))
            * epsilon_pow5
        - (5. / 256.
            + 15. / 1024. * EARTH_THIRD_FLATTENING
            + 35. / 1024. * EARTH_THIRD_FLATTENING.powf(2.)
            + 7. / 512. * EARTH_THIRD_FLATTENING.powf(3.))
            * epsilon_pow6
        - (25. / 2048.
            + 35. / 2048. * EARTH_THIRD_FLATTENING
            + 21. / 2048. * EARTH_THIRD_FLATTENING.powf(2.))
            * epsilon_pow7
        - (175. / 16384. + 35. / 4096. * EARTH_THIRD_FLATTENING) * epsilon_pow8
        - 245. / 32768. * epsilon_pow9;

    let lambda_one = omega
        - (EARTH_FLATTENING * sin_alpha_0 * i_three_fourier_series(epsilon_pow_tup, sigma, a_3));
    let lambda_two = omega_two
        - (EARTH_FLATTENING
            * sin_alpha_0
            * i_three_fourier_series(epsilon_pow_tup, sigma_two, a_3));

    let lon = lambda_two - lambda_one;
    let lat = (beta_two.tan() / (1. - EARTH_FLATTENING)).atan();

    return LocBearing {
        lat,
        lon: lon + lon1,
        bearing: alpha_two,
    };
}

#[cfg(target_arch = "wasm32")]
#[wasm_bindgen]
pub fn karney_distance_and_bearing(
    lat1: Radians,
    lon1: Radians,
    lat2: Radians,
    lon2: Radians,
) -> DistBearing {
    return distance_and_bearing(lat1, lon1, lat2, lon2);
}

pub fn distance_and_bearing(
    lat1: Radians,
    lon1: Radians,
    lat2: Radians,
    lon2: Radians,
) -> DistBearing {
    let lambda_one_two = lon2 - lon1;
    let beta_one = ((1. - EARTH_FLATTENING) * lat1.tan()).atan();
    let beta_two = ((1. - EARTH_FLATTENING) * lat2.tan()).atan();

    let sin_beta_one = beta_one.sin();
    let sin_beta_two = beta_two.sin();
    let cos_beta_one = beta_one.cos();
    let cos_beta_two = beta_two.cos();

    let omega =
        (1. - EARTH_ECCENTRICITY_SQUARED * ((cos_beta_one + cos_beta_two) / 2.).powf(2.)).sqrt();

    let omega_one_two = lambda_one_two / omega;

    let cos_omega_one_two = omega_one_two.cos();
    let sin_omega_one_two = omega_one_two.sin();

    let abs_z_one = (((cos_beta_one * sin_beta_two)
        - (sin_beta_one * cos_beta_two * cos_omega_one_two))
        .powf(2.)
        + (cos_beta_two * sin_omega_one_two).powf(2.))
    .sqrt();

    let sigma_one_two = abs_z_one
        .atan2((sin_beta_one * sin_beta_two) + (cos_beta_one * cos_beta_two * cos_omega_one_two));

    let alpha_one = (cos_beta_two * sin_omega_one_two)
        .atan2(cos_beta_one * sin_beta_two - sin_beta_one * cos_beta_two * cos_omega_one_two);

    let distance = EARTH_MAJOR * omega * sigma_one_two;

    return DistBearing {
        distance,
        bearing: alpha_one,
    };
}

fn i_one_fourier_series(
    epsilon_pow_tup: (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
    sigma: f64,
    a_1: f64,
) -> f64 {
    // The Cs are coeffs of the fourier transform --> https://geographiclib.sourceforge.io/html/geodesic.html
    // and Karney 2011.
    let b_1: f64 = ((- 1./2. * epsilon_pow_tup.0   // C11
                    + 3./16. * epsilon_pow_tup.2
                    - 1./32. * epsilon_pow_tup.4
                    + 19./2048. * epsilon_pow_tup.6
                    - 3./4096. * epsilon_pow_tup.8) * (2.* sigma).sin())
                + ((- 1./16. * epsilon_pow_tup.1 // C12
                    + 1./32. * epsilon_pow_tup.3
                    - 9./2048. * epsilon_pow_tup.5
                    + 7./4096. * epsilon_pow_tup.7
                    + 1./65536. * epsilon_pow_tup.9) * (2. * sigma * 2.).sin())
                + ((- 1./48. * epsilon_pow_tup.2// C13
                    + 3./256. * epsilon_pow_tup.4
                    - 3./2048. * epsilon_pow_tup.6
                    + 17./24576. * epsilon_pow_tup.8) * (3. * sigma * 2.).sin())
                + ((- 5./512. * epsilon_pow_tup.3 // C14
                    + 3./512. * epsilon_pow_tup.5
                    - 11./16384. * epsilon_pow_tup.7
                    + 3./8192. * epsilon_pow_tup.9) * (4. * sigma * 2.).sin())
                + ((- 7./1280. * epsilon_pow_tup.4 //C15
                    + 7./2048. * epsilon_pow_tup.6
                    - 3./8192. * epsilon_pow_tup.8) * (5. * sigma * 2.).sin())
                + ((- 7./2048. * epsilon_pow_tup.5 // C16
                    + 9./4096. * epsilon_pow_tup.7
                    - 117./524288. * epsilon_pow_tup.9) * (6. * sigma * 2.).sin())
                + ((- 33./14336. * epsilon_pow_tup.6// C17
                    + 99./65536. * epsilon_pow_tup.8) * (7. * sigma * 2.).sin())
                + ((- 429./262144. * epsilon_pow_tup.7// C18
                    + 143./131072. * epsilon_pow_tup.9) * (8. * sigma * 2.).sin())
                + (( - 715./589824. * epsilon_pow_tup.8) * (9. * sigma * 2.).sin()) // C19
                + (( - 2431./2621440. * epsilon_pow_tup.9) * (10. * sigma * 2.).sin()); // C110
    return a_1 * (sigma + b_1);
}

fn sigma_two_fourier_series(
    epsilon_pow_tup: (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
    tau: f64,
) -> f64 {
    let d_1: f64 = ((1./2. * epsilon_pow_tup.0 //C'11
                    - 9./32. * epsilon_pow_tup.2
        + 205. / 1536. * epsilon_pow_tup.4
        - 4879. / 73728. * epsilon_pow_tup.6
        + 9039. / 327680. * epsilon_pow_tup.8)
        * (2. * tau).sin())
        + ((5./16. * epsilon_pow_tup.1 //C'12
                    - 37./96. * epsilon_pow_tup.3
            + 1335. / 4096. * epsilon_pow_tup.5
            - 86171. / 368640. * epsilon_pow_tup.7
            + 4119073. / 28311552. * epsilon_pow_tup.9)
            * (2. * 2. * tau).sin())
        + ((29./96. * epsilon_pow_tup.2 // C'13
                    - 75./128. * epsilon_pow_tup.4
            + 2901. / 4096. * epsilon_pow_tup.6
            - 443327. / 655360. * epsilon_pow_tup.8)
            * (2. * 3. * tau).sin())
        + ((539./1536. * epsilon_pow_tup.3 // C'14
                    - 2391./2560. * epsilon_pow_tup.5
            + 1082857. / 737280. * epsilon_pow_tup.7
            - 2722891. / 1548288. * epsilon_pow_tup.9)
            * (2. * 4. * tau).sin())
        + ((3467./7680. * epsilon_pow_tup.4 // C'15
                    - 28223./18432. * epsilon_pow_tup.6
            + 1361343. / 458752. * epsilon_pow_tup.8)
            * (2. * 5. * tau).sin())
        + ((38081./61440. * epsilon_pow_tup.5 // C'16
                    - 733437./286720. * epsilon_pow_tup.7
            + 10820079. / 1835008. * epsilon_pow_tup.9)
            * (2. * 6. * tau).sin())
        + ((459485./516096. * epsilon_pow_tup.6 //C'17
                    - 709743./163840. * epsilon_pow_tup.8)
            * (2. * 7. * tau).sin())
        + ((109167851./82575360. * epsilon_pow_tup.7 //C'18
                    - 550835669./74317824. * epsilon_pow_tup.9)
            * (2. * 8. * tau).sin())
        + ((83141299. / 41287680. * epsilon_pow_tup.8
            + 9303339907. / 2972712960. * epsilon_pow_tup.9)
            * (2. * 9. * tau).sin());

    return tau + d_1;
}

fn i_three_fourier_series(
    epsilon_pow_tup: (f64, f64, f64, f64, f64, f64, f64, f64, f64, f64),
    sigma: f64,
    a_3: f64,
) -> f64 {
    let e_1: f64 = (((1./4. - 1./4.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.0 // C31
                    + (1./8. - 1./8.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.1
                    + (3./64. + 3./64.*EARTH_THIRD_FLATTENING - 1./64.*EARTH_THIRD_FLATTENING.powf(2.) - 5./64.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.2
                    + (5./128. + 1./64.*EARTH_THIRD_FLATTENING + 1./64.*EARTH_THIRD_FLATTENING.powf(2.) - 1./64.*EARTH_THIRD_FLATTENING.powf(3.) - 7./128.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.3
                    + (3./128. + 11./512.*EARTH_THIRD_FLATTENING + 3./512.*EARTH_THIRD_FLATTENING.powf(2.) + 1./256.*EARTH_THIRD_FLATTENING.powf(3.) - 7./512.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.4
                    + (21./1024. + 5./512.*EARTH_THIRD_FLATTENING + 13./1024.*EARTH_THIRD_FLATTENING.powf(2.) + 1./512.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.5
                    + (243./16384. + 189./16384.*EARTH_THIRD_FLATTENING + 83./16384.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6
                    + (435./32768. + 109./16384.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (345./32768. * epsilon_pow_tup.8)) * (2. * sigma).sin())
                + (((1./16. - 3./32.*EARTH_THIRD_FLATTENING + 1./32.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.1 // C32
                    + (3./64. - 1./32.*EARTH_THIRD_FLATTENING - 3./64.*EARTH_THIRD_FLATTENING.powf(2.) + 1./32.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.2
                    + (3./128. + 1./128.*EARTH_THIRD_FLATTENING - 9./256.*EARTH_THIRD_FLATTENING.powf(2.) - 3./128.*EARTH_THIRD_FLATTENING.powf(3.) + 7./256.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.3
                    + (5./256. + 1./256.*EARTH_THIRD_FLATTENING - 1./128.*EARTH_THIRD_FLATTENING.powf(2.) - 7./256.*EARTH_THIRD_FLATTENING.powf(3.) - 3./256.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.4
                    + (27./2048. + 69./8192.*EARTH_THIRD_FLATTENING - 39./8192.*EARTH_THIRD_FLATTENING.powf(2.) - 47./4096.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.5
                    + (187./16384. + 39./8192.*EARTH_THIRD_FLATTENING + 31./16384.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6
                    + (287./32768. + 47./8192.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (255./32768. * epsilon_pow_tup.8)) * (2. * 2. * sigma).sin())
                + (((5./192. - 3./64.*EARTH_THIRD_FLATTENING + 5./192.*EARTH_THIRD_FLATTENING.powf(2.) - 1./192.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.2 // C33
                    + (3./128. - 5./192.*EARTH_THIRD_FLATTENING - 1./64.*EARTH_THIRD_FLATTENING.powf(2.) + 5./192.*EARTH_THIRD_FLATTENING.powf(3.) - 1./128.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.3
                    + (7./512. - 1./384.*EARTH_THIRD_FLATTENING - 77./3072.*EARTH_THIRD_FLATTENING.powf(2.) + 5./3072.*EARTH_THIRD_FLATTENING.powf(3.) + 65./3072.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.4
                    + (3./256. - 1./1024.*EARTH_THIRD_FLATTENING - 71./6144.*EARTH_THIRD_FLATTENING.powf(2.) - 47./3072.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.5
                    + (139./16384. + 143./49152.*EARTH_THIRD_FLATTENING - 383./49152.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6
                    + (243./32768. + 95./49152.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (581./98304. * epsilon_pow_tup.8)) * (2. * 3. * sigma).sin())
                + (((7./512. - 7./256.*EARTH_THIRD_FLATTENING + 5./256.*EARTH_THIRD_FLATTENING.powf(2.) - 7./1024.*EARTH_THIRD_FLATTENING.powf(3.) + 1./1024.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.3// C34
                    + (7./512. - 5./256.*EARTH_THIRD_FLATTENING - 7./2048.*EARTH_THIRD_FLATTENING.powf(2.)+ 9./512.*EARTH_THIRD_FLATTENING.powf(3.) - 21./2048.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.4
                    + (9./1024. - 43./8192.*EARTH_THIRD_FLATTENING - 129./8192.*EARTH_THIRD_FLATTENING.powf(2.) + 39./4096.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.5
                    + (127./16384. - 23./8192.*EARTH_THIRD_FLATTENING - 165./16384.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6
                    + (193./32768. + 3./8192.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (171./32768. * epsilon_pow_tup.8)) * (2. * 4. * sigma).sin()) // C35
                + (((21./2560. - 9./512.*EARTH_THIRD_FLATTENING + 15./1024.*EARTH_THIRD_FLATTENING.powf(2.) - 7./1024.*EARTH_THIRD_FLATTENING.powf(3.) + 9./5120.*EARTH_THIRD_FLATTENING.powf(4.)) * epsilon_pow_tup.4
                    + (9./1024. - 15./1024.*EARTH_THIRD_FLATTENING + 3./2048.*EARTH_THIRD_FLATTENING.powf(2.) + 57./5120.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.5
                    + (99./16384. - 91./16384.*EARTH_THIRD_FLATTENING - 781./81920.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6
                    + (179./32768. - 55./16384.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (141./32768. * epsilon_pow_tup.8)) * (2. * 5. * sigma).sin())
                + (((11./2048. - 99./8192.*EARTH_THIRD_FLATTENING + 275./24576.*EARTH_THIRD_FLATTENING.powf(2.) - 77./12288.*EARTH_THIRD_FLATTENING.powf(3.)) * epsilon_pow_tup.5 // C36
                    + (99./16384. - 275./24576.*EARTH_THIRD_FLATTENING + 55./16384.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6
                    + (143./32768. - 253./49152.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (33./8192. * epsilon_pow_tup.8)) * (2. * 6. * sigma).sin())
                + (((429./114688. - 143./16384.*EARTH_THIRD_FLATTENING + 143./16384.*EARTH_THIRD_FLATTENING.powf(2.)) * epsilon_pow_tup.6 // C37
                    + (143./32768. - 143./16384.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7
                    + (429./131072. * epsilon_pow_tup.8)) * (2. * 7. * sigma).sin())
                + (((715./262144. - 429./65536.*EARTH_THIRD_FLATTENING) * epsilon_pow_tup.7 + (429./131072. * epsilon_pow_tup.8)) * (2. * 8. * sigma).sin()) // C38
                + ((2431./1179648. * epsilon_pow_tup.8) * (2. * 9. * sigma).sin()); //C39

    return a_3 * (sigma + e_1);
}
