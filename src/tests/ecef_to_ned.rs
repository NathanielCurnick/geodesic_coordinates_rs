use approx::assert_relative_eq;
use peroxide::prelude::matrix;
use peroxide::prelude::Shape::Row;

use crate::reference_frames::ecef::{
    construct_ecef_to_ned_jacobian, generate_ecef_to_ned_matrix, ECEFVel, ECEF,
};
use crate::reference_frames::ned::{NEDVel, NED};
use crate::reference_frames::wgs84::WGS84Coord;

#[test]
fn test_ecef_to_ned_jacobian() {
    let basic_matrix = matrix(vec![1.0; 9], 3, 3, Row);

    let jacobian = construct_ecef_to_ned_jacobian(&basic_matrix);

    println!("{}", jacobian);

    let correct = vec![
        1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0,
        0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0,
    ];

    let correct_matrix = matrix(correct.clone(), 6, 6, Row);

    println!("{}", correct_matrix);

    for (a, b) in correct.iter().zip(jacobian.data.iter()) {
        assert_eq!(a, b);
    }

    // =============================================

    let basic_matrix = matrix(vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 9.0], 3, 3, Row);

    let jacobian = construct_ecef_to_ned_jacobian(&basic_matrix);

    println!("{}", jacobian);

    let correct = vec![
        1.0, 2.0, 3.0, 0.0, 0.0, 0.0, //
        4.0, 5.0, 6.0, 0.0, 0.0, 0.0, //
        7.0, 8.0, 9.0, 0.0, 0.0, 0.0, //
        0.0, 0.0, 0.0, 1.0, 2.0, 3.0, //
        0.0, 0.0, 0.0, 4.0, 5.0, 6.0, //
        0.0, 0.0, 0.0, 7.0, 8.0, 9.0, //
    ];

    for (a, b) in correct.iter().zip(jacobian.data.iter()) {
        assert_eq!(a, b);
    }
}

#[test]
fn test_ecef_from_ned_rot() {
    let lat = 50.0_f64;
    let lon = 10.0_f64;
    let alt = 0.0_f64;
    let reference_point = WGS84Coord::new_from_radians(lat.to_radians(), lon.to_radians(), alt);

    let ecef = ECEF::new_from_wgs84(&reference_point);

    let ned = NED {
        n: 0.0,
        e: 0.0,
        d: 0.0,
    };

    let rotation = generate_ecef_to_ned_matrix(&reference_point);

    let final_ecef = ECEF::new_from_ned_rot(&ned, &rotation, &ecef);

    assert_relative_eq!(ecef.x, final_ecef.x, epsilon = 1e-4);
    assert_relative_eq!(ecef.y, final_ecef.y, epsilon = 1e-4);
    assert_relative_eq!(ecef.z, final_ecef.z, epsilon = 1e-4);
}

#[test]
fn test_ecef_to_ned_vel() {
    let lat = 50.0_f64;
    let lon = 10.0_f64;
    let alt = 0.0_f64;
    let reference_point = WGS84Coord::new_from_radians(lat.to_radians(), lon.to_radians(), alt);

    let ecef = ECEF::new_from_wgs84(&reference_point);

    // You can get lla to ecef test cases from
    // https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm

    assert_relative_eq!(ecef.x, 4045456.0_f64, epsilon = 1.0_f64);
    assert_relative_eq!(ecef.y, 713323.0_f64, epsilon = 1.0);
    assert_relative_eq!(ecef.z, 4862789_f64, epsilon = 1.0);

    let rotation = generate_ecef_to_ned_matrix(&reference_point);

    let start_ned = NED {
        n: 0.0,
        e: 0.0,
        d: 0.0,
    };

    let end_ned = NED {
        n: 100.0,
        e: 100.0,
        d: 0.0,
    };

    let dt = 10_f64;

    let ned_vel = NEDVel::linear_speed(&start_ned, &end_ned, dt);

    let ecef_vel = ECEFVel::new_from_ned_rot(&ned_vel, &rotation);

    let final_ecef_point = ECEF {
        x: ecef.x + ecef_vel.x_vel * dt,
        y: ecef.y + ecef_vel.y_vel * dt,
        z: ecef.z + ecef_vel.z_vel * dt,
    };

    let final_ecef_to_ned = NED::new_from_ecef_rot(&final_ecef_point, &rotation, &ecef);

    assert_relative_eq!(final_ecef_to_ned.n, end_ned.n, epsilon = 1e-6);
    assert_relative_eq!(final_ecef_to_ned.e, end_ned.e, epsilon = 1e-6);
    assert_relative_eq!(final_ecef_to_ned.d, end_ned.d, epsilon = 1e-6);
}

#[test]
fn test_ecef_to_ned() {
    let lat = 50.0_f64;
    let lon = 10.0_f64;
    let alt = 0.0_f64;
    let reference_point = WGS84Coord::new_from_radians(lat.to_radians(), lon.to_radians(), alt);

    let ecef = ECEF::new_from_wgs84(&reference_point);

    // You can get lla to ecef test cases from
    // https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm

    assert_relative_eq!(ecef.x, 4045456.0_f64, epsilon = 1.0_f64);
    assert_relative_eq!(ecef.y, 713323.0_f64, epsilon = 1.0);
    assert_relative_eq!(ecef.z, 4862789_f64, epsilon = 1.0);

    let rotation = generate_ecef_to_ned_matrix(&reference_point);

    let last = WGS84Coord::new_from_ned(
        &NED {
            n: 0.0,
            e: 0.0,
            d: 0.0,
        },
        &rotation,
        &ecef,
    );

    assert_relative_eq!(last.get_lat_radians(), lat.to_radians(), epsilon = 1e-8);
    assert_relative_eq!(last.get_lon_radians(), lon.to_radians(), epsilon = 1e-8);
    assert_relative_eq!(last.get_altitude(), alt, epsilon = 1e-8);
}

#[test]
fn test_lla_to_ecef() {
    let reference_point =
        WGS84Coord::new_from_radians(50.0_f64.to_radians(), 10.0_f64.to_radians(), 0.0);

    let ecef = ECEF::new_from_wgs84(&reference_point);

    // You can get lla to ecef test cases from
    // https://www.oc.nps.edu/oc2902w/coord/llhxyz.htm

    assert_relative_eq!(ecef.x, 4045456.0_f64, epsilon = 1.0);
    assert_relative_eq!(ecef.y, 713323.0_f64, epsilon = 1.0);
    assert_relative_eq!(ecef.z, 4862789_f64, epsilon = 1.0);
}
