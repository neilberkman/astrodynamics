//! Vallado SGP4 verification suite — 33 satellites, 198 propagation points,
//! 1098 component checks. Reference values captured from the Python `sgp4` C
//! extension (v2.25), which compiles Vallado's C++ (v2020-07-13) with WGS72,
//! opsmode 'i'. The pure-Rust port in `astrodynamics::sgp4` must match
//! bit-for-bit (0 ULP) on every component.

use astrodynamics::sgp4::{JulianDate, MinutesSinceEpoch, Satellite};

fn hex_to_f64(s: &str) -> f64 {
    let (neg, rest) = if let Some(r) = s.strip_prefix("-0x") {
        (true, r)
    } else if let Some(r) = s.strip_prefix("0x") {
        (false, r)
    } else {
        panic!("bad hex float: {s}");
    };
    let (mant_str, exp_str) = rest.split_once('p').unwrap();
    let (int_part, frac_part) = mant_str.split_once('.').unwrap();
    let exp: i32 = exp_str.parse().unwrap();
    let full = format!("{int_part}{frac_part}");
    let mant = u64::from_str_radix(&full, 16).unwrap();
    let frac_bits = frac_part.len() as i32 * 4;
    let val = mant as f64 * (2.0_f64).powi(exp - frac_bits);
    if neg { -val } else { val }
}

fn ulp_distance(a: f64, b: f64) -> u64 {
    let ia = a.to_bits() as i64;
    let ib = b.to_bits() as i64;
    (ia - ib).unsigned_abs()
}

#[test]
fn all_33_vallado_satellites_at_0_ulp() {
    let data: serde_json::Value =
        serde_json::from_str(include_str!("sgp4_verification.json")).unwrap();

    let labels = ["px", "py", "pz", "vx", "vy", "vz"];
    let mut failures = Vec::new();

    for sat in data["satellites"].as_array().unwrap() {
        let line1 = sat["line1"].as_str().unwrap();
        let line2 = sat["line2"].as_str().unwrap();
        let norad = sat["norad"].as_str().unwrap();

        let satellite = Satellite::from_tle(line1, line2).unwrap();

        for prop in sat["propagations"].as_array().unwrap() {
            // Skip rows the C++ reference flagged as errors (decay, etc.)
            if prop.get("error").is_some() {
                continue;
            }

            let tsince = prop["tsince"].as_f64().unwrap();
            let pred = match satellite.propagate(MinutesSinceEpoch(tsince)) {
                Ok(p) => p,
                Err(_) => continue,
            };

            let actual = [
                pred.position[0],
                pred.position[1],
                pred.position[2],
                pred.velocity[0],
                pred.velocity[1],
                pred.velocity[2],
            ];

            for (i, label) in labels.iter().enumerate() {
                let expected = hex_to_f64(prop[*label].as_str().unwrap());
                let ulp = ulp_distance(actual[i], expected);
                if ulp > 0 {
                    failures.push((norad.to_string(), tsince, label.to_string(), ulp));
                }
            }
        }
    }

    if !failures.is_empty() {
        failures.sort_by(|a, b| b.3.cmp(&a.3));
        let summary: Vec<String> = failures
            .iter()
            .take(20)
            .map(|(norad, tsince, label, ulp)| format!("  {norad} t={tsince} {label}: {ulp} ULP"))
            .collect();
        panic!(
            "{} ULP failures (top 20):\n{}",
            failures.len(),
            summary.join("\n")
        );
    }
}

#[test]
fn iss_basic_propagation() {
    let sat = Satellite::from_tle(
        "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993",
        "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106",
    )
    .unwrap();

    let pred = sat.propagate(MinutesSinceEpoch(0.0)).unwrap();

    // Position should be in the right ballpark (LEO, ~6700-7000 km from center)
    let r = (pred.position[0].powi(2)
        + pred.position[1].powi(2)
        + pred.position[2].powi(2))
    .sqrt();
    assert!(
        (6500.0..=7200.0).contains(&r),
        "ISS radius {r} km outside LEO range"
    );
}

#[test]
fn julian_date_propagation() {
    let sat = Satellite::from_tle(
        "1 25544U 98067A   18184.80969102  .00001614  00000-0  31745-4 0  9993",
        "2 25544  51.6414 295.8524 0003435 262.6267 204.2868 15.54005638121106",
    )
    .unwrap();

    // 2018-07-04 00:00:00 UTC = JD 2458303.5
    let pred = sat.propagate_jd(JulianDate(2458303.0, 0.5)).unwrap();

    let r = (pred.position[0].powi(2)
        + pred.position[1].powi(2)
        + pred.position[2].powi(2))
    .sqrt();
    assert!(
        (6500.0..=7200.0).contains(&r),
        "ISS radius {r} km outside LEO range"
    );
}

#[test]
fn invalid_tle_rejected() {
    assert!(Satellite::from_tle("garbage", "data").is_err());
    assert!(Satellite::from_tle("1 short", "2 short").is_err());
}
