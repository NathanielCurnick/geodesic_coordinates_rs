use crate::{
    constants::{EARTH_FLATTENING, EARTH_MAJOR, EARTH_MINOR},
    types::{DistBearing, LocBearing, Metres, Radians},
};

pub fn location_and_bearing(
    lat1: Radians,
    lon1: Radians,
    bearing: Radians,
    distance: Metres,
) -> LocBearing {
    let sin_alpha1 = bearing.sin();
    let cos_alpha1 = bearing.cos();

    let tan_u1 = (1. - EARTH_FLATTENING) * lat1.tan();
    let cos_u1 = 1. / (1. + (tan_u1 * tan_u1)).sqrt();
    let sin_u1 = tan_u1 * cos_u1;

    let sigma_1 = tan_u1.atan2(cos_alpha1);
    let sin_alpha = cos_u1 * sin_alpha1;
    let cos_s_q_alpha = 1.0 - (sin_alpha * sin_alpha);

    let u_s_q = cos_s_q_alpha * ((EARTH_MAJOR * EARTH_MAJOR) - (EARTH_MINOR * EARTH_MINOR))
        / (EARTH_MINOR * EARTH_MINOR);

    let a_vin =
        1. + (u_s_q / 16384.0) * (4096.0 + (u_s_q * (-768.0 + u_s_q * (320.0 - 175.0 * u_s_q))));
    let b_vin = (u_s_q / 1024.0) * (256.0 + u_s_q * (-128.0 + u_s_q * (74.0 - 47.0 * u_s_q)));

    let sigma = distance / (a_vin * EARTH_MINOR);

    let mut cos2_sigma_m = ((2. * sigma_1) + sigma).cos();
    let mut sin_sigma = sigma.sin();
    let mut cos_sigma = sigma.cos();
    let mut cos_2_sigma_m_squared = cos2_sigma_m * cos2_sigma_m;

    let mut delta_sigma = b_vin
        * sin_sigma
        * (cos2_sigma_m
            + (b_vin / 4.0)
                * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m_squared)
                    - (b_vin / 6.0)
                        * cos2_sigma_m
                        * (-3.0 + 4.0 * sin_sigma.powf(2.))
                        * (-3.0 + 4.0 * cos_2_sigma_m_squared)));

    let mut sigma_prime = sigma;
    let mut sigma = (distance / (EARTH_MINOR * a_vin)) + delta_sigma;

    while (sigma - sigma_prime).abs() > 1e-12 {
        cos2_sigma_m = ((2. * sigma_1) + sigma).cos();
        sin_sigma = sigma.sin();
        cos_sigma = sigma.cos();
        cos_2_sigma_m_squared = cos2_sigma_m * cos2_sigma_m;

        delta_sigma = b_vin
            * sin_sigma
            * (cos2_sigma_m
                + (b_vin / 4.0)
                    * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m_squared)
                        - (b_vin / 6.0)
                            * cos2_sigma_m
                            * (-3.0 + 4.0 * sin_sigma.powf(2.))
                            * (-3.0 + 4.0 * cos_2_sigma_m_squared)));

        sigma_prime = sigma;
        sigma = (distance / (EARTH_MINOR * a_vin)) + delta_sigma;
    }
    let x = sin_u1 * sin_sigma - cos_u1 * cos_sigma * cos_alpha1;
    let lat2 = (sin_u1 * cos_sigma + cos_u1 * sin_sigma * cos_alpha1)
        .atan2((1. - EARTH_FLATTENING) * (sin_alpha.powf(2.) + x.powf(2.)).sqrt());

    let lambda =
        (sin_sigma * sin_alpha1).atan2(cos_u1 * cos_sigma - sin_u1 * sin_sigma * cos_alpha1);

    let c = EARTH_FLATTENING
        / (16. * cos_s_q_alpha * (4.0 + EARTH_FLATTENING * (4.0 - 3.0 * cos_s_q_alpha)));

    let l = lambda
        - ((1.0 - c)
            * EARTH_FLATTENING
            * sin_alpha
            * (sigma
                + c * sin_sigma
                    * (cos2_sigma_m + c * cos_sigma * (-1.0 + 2.0 * cos2_sigma_m.powf(2.)))));

    return LocBearing {
        lat: lat2,
        lon: lon1 + l,
        bearing: sin_alpha.atan2(-x),
    };
}

pub fn distance_and_bearing(
    lat1: Radians,
    lon1: Radians,
    lat2: Radians,
    lon2: Radians,
) -> DistBearing {
    let l = lon2 - lon1;

    let tan_u1 = (1. - EARTH_FLATTENING) * lat1.tan();
    let cos_u1 = 1. / (1. + (tan_u1 * tan_u1)).sqrt();
    let sin_u1 = tan_u1 * cos_u1;

    let tan_u2 = (1. - EARTH_FLATTENING) * lat2.tan();
    let cos_u2 = 1. / (1. + (tan_u2 * tan_u2)).sqrt();
    let sin_u2 = tan_u2 * cos_u2;

    let mut lambda = l;

    let mut sin_lambda = lambda.sin();
    let mut cos_lambda = lambda.cos();
    let mut sin_s_q_sigma =
        (cos_u2 * sin_lambda).powf(2.) + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda).powf(2.);
    let mut sin_sigma = sin_s_q_sigma.sqrt();
    let mut cos_sigma = (sin_u1 * sin_u2) + (cos_u1 * cos_u2 * cos_lambda);

    let mut sigma = sin_sigma.atan2(cos_sigma);
    let mut sin_alpha = (cos_u1 * cos_u2 * sin_lambda) / sin_sigma;
    let mut cos_s_q_alpha = 1. - (sin_alpha * sin_alpha);
    let mut cos_2_sigma_m = cos_sigma - 2. * sin_u1 * sin_u2 / cos_s_q_alpha;
    let mut c = (EARTH_FLATTENING / 16.)
        * (cos_s_q_alpha * (4. + EARTH_FLATTENING * (4. - 3. * cos_s_q_alpha)));
    let mut lambda_prime = lambda;

    lambda = l
        + (1. - c)
            * EARTH_FLATTENING
            * sin_alpha
            * (sigma
                + c * sin_sigma
                    * (cos_2_sigma_m
                        + c * cos_sigma * (-1. + 2. * (cos_2_sigma_m * cos_2_sigma_m))));

    while (lambda - lambda_prime).abs() > 1e-12 {
        sin_lambda = lambda.sin();
        cos_lambda = lambda.cos();
        sin_s_q_sigma = (cos_u2 * sin_lambda).powf(2.)
            + (cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda).powf(2.);
        sin_sigma = sin_s_q_sigma.sqrt();
        cos_sigma = sin_u1 * sin_u2 + cos_u1 * cos_u2 * cos_lambda;

        sigma = sin_sigma.atan2(cos_sigma);
        sin_alpha = (cos_u1 * cos_u2 * sin_lambda) / sin_sigma;
        cos_s_q_alpha = 1. - (sin_alpha * sin_alpha);
        cos_2_sigma_m = cos_sigma - ((2. * sin_u1 * sin_u2) / cos_s_q_alpha);
        c = (EARTH_FLATTENING / 16.)
            * (cos_s_q_alpha * (4. + EARTH_FLATTENING * (4. - (3. * cos_s_q_alpha))));
        lambda_prime = lambda;

        lambda = l
            + (1. - c)
                * EARTH_FLATTENING
                * sin_alpha
                * (sigma
                    + c * sin_sigma
                        * (cos_2_sigma_m
                            + c * cos_sigma * (-1. + 2. * (cos_2_sigma_m * cos_2_sigma_m))));
    }

    let u_s_q = cos_s_q_alpha
        * (((EARTH_MAJOR * EARTH_MAJOR) - (EARTH_MINOR * EARTH_MINOR))
            / (EARTH_MINOR * EARTH_MINOR));
    let a = 1. + (u_s_q / 16384.0) * (4096.0 + u_s_q * (-768.0 + u_s_q * (320.0 - 175.0 * u_s_q)));
    let b = (u_s_q / 1024.0) * (256.0 + u_s_q * (-128.0 + u_s_q * (74.0 - 47.0 * u_s_q)));
    let delta_sigma = b
        * sin_sigma
        * (cos_2_sigma_m
            + (b / 4.0)
                * (cos_sigma * (-1.0 + 2.0 * cos_2_sigma_m.powf(2.))
                    - (b / 6.0)
                        * cos_2_sigma_m
                        * (-3.0 + 4.0 * sin_sigma.powf(2.))
                        * (-3.0 + 4.0 * cos_2_sigma_m.powf(2.))));

    let distance = EARTH_MINOR * a * (sigma - delta_sigma);

    let alpha1 = (cos_u2 * sin_lambda).atan2(cos_u1 * sin_u2 - sin_u1 * cos_u2 * cos_lambda);

    return DistBearing {
        distance,
        bearing: alpha1,
    };
}
