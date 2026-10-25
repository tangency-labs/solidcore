//! Non-C++ regression coverage for the runner's actual `.sms` reader.

pub use difftest::HarnessError;

#[path = "../src/units/nlib_bezier_sms.rs"]
mod sms;

fn assert_rejected(text: &str, message: &str) {
    // A panic (including capacity overflow) fails the test, rather than counting
    // as rejection. No allocation should be driven by these header-only inputs.
    let result = sms::parse_sms(text, "malformed.sms");
    assert!(result.is_err(), "accepted {text}");
    if let Err(error) = result {
        let error = error.to_string();
        assert!(error.contains("malformed.sms"), "{error}");
        assert!(error.contains(message), "{error}");
    }
}

#[test]
fn saturated_header_count_is_rejected_without_panicking() {
    assert_rejected("1e100 0 0 0 0", "bad n count");
}

#[test]
fn exclusive_usize_limit_is_rejected() {
    let upper = 2.0_f64.powi(usize::BITS as i32);
    assert_rejected(&format!("{upper:.0} 0 0 0 0"), "bad n count");
}

#[test]
fn nonfinite_counts_are_rejected() {
    for value in ["NaN", "inf", "-inf"] {
        for header in [
            format!("{value} 0 0 0 0"),
            format!("0 {value} 0 0 0"),
            format!("0 0 {value} 0 0"),
            format!("0 0 0 {value} 0"),
        ] {
            assert_rejected(&header, "count");
        }
    }
}

#[test]
fn negative_and_fractional_counts_are_rejected() {
    for value in ["-1", "0.5"] {
        for header in [
            format!("{value} 0 0 0 0"),
            format!("0 {value} 0 0 0"),
            format!("0 0 {value} 0 0"),
            format!("0 0 0 {value} 0"),
        ] {
            assert_rejected(&header, "count");
        }
    }
}

#[test]
fn oversized_knot_dimensions_are_rejected() {
    let large = 2.0_f64.powi(usize::BITS as i32 - 1);
    assert_rejected(
        &format!("{large:.0} 0 {large:.0} 0 0"),
        "u knot count overflow",
    );
    assert_rejected(
        &format!("0 {large:.0} 0 {large:.0} 0"),
        "v knot count overflow",
    );
}

#[test]
fn point_token_arithmetic_overflow_is_rejected() {
    let large = 2.0_f64.powi(usize::BITS as i32 - 1);
    assert_rejected(
        &format!("{large:.0} {large:.0} 0 0 0"),
        "required token count overflow",
    );
}

#[test]
fn large_representable_dimensions_with_no_data_are_rejected() {
    assert_rejected("1000000000 0 0 0 0", "truncated");
    assert_rejected("0 1000000000 0 0 1", "truncated");
    assert_rejected("0 0 1000000000 0 0", "truncated");
    assert_rejected("0 0 0 1000000000 0", "truncated");
}

#[test]
fn truncated_headers_points_and_knots_are_rejected() {
    for text in [
        "",
        "0 0 0 0",
        "0 0 0 0 0 1 2",
        "0 0 0 0 1 1 2 3",
        "0 0 0 0 0 1 2 3 0",
        "0 0 0 0 0 1 2 3 0 1 0",
    ] {
        assert_rejected(text, "truncated");
    }
}

#[test]
fn minimal_nonrational_surface_decodes_exactly() -> Result<(), HarnessError> {
    let surface = sms::parse_sms("0e0 0.0 0 0 0 -0.0 2.5 -3e-2 -1 1 2 4", "valid.sms")?;
    assert_eq!(
        surface.net,
        vec![vec![[-0.0, 2.5, -0.03, solidcore::nurbs::math::NOW]]]
    );
    assert_eq!(surface.net[0][0][0].to_bits(), (-0.0_f64).to_bits());
    assert_eq!((surface.p, surface.q), (0, 0));
    assert_eq!(surface.u, vec![-1.0, 1.0]);
    assert_eq!(surface.v, vec![2.0, 4.0]);
    assert_eq!(surface.used, None);
    Ok(())
}

#[test]
fn minimal_rational_surface_decodes_exactly() -> Result<(), HarnessError> {
    // Weighted coordinates are stored verbatim, not divided by w.
    let surface = sms::parse_sms("0 0 0 0 1 2 -4 6 2 0 1 -2 3", "valid.sms")?;
    assert_eq!(surface.net, vec![vec![[2.0, -4.0, 6.0, 2.0]]]);
    assert_eq!((surface.p, surface.q), (0, 0));
    assert_eq!(surface.u, vec![0.0, 1.0]);
    assert_eq!(surface.v, vec![-2.0, 3.0]);
    assert_eq!(surface.used, None);
    Ok(())
}

#[test]
fn nonsquare_surface_preserves_row_order_degrees_and_knots() -> Result<(), HarnessError> {
    let surface = sms::parse_sms(
        "1 2 1 2 0 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 \
         0 0 1 1 -1 -1 -1 2 2 2 trailing",
        "valid.sms",
    )?;
    let w = solidcore::nurbs::math::NOW;
    assert_eq!(
        surface.net,
        vec![
            vec![[0.0, 1.0, 2.0, w], [3.0, 4.0, 5.0, w], [6.0, 7.0, 8.0, w]],
            vec![
                [9.0, 10.0, 11.0, w],
                [12.0, 13.0, 14.0, w],
                [15.0, 16.0, 17.0, w]
            ],
        ]
    );
    assert_eq!((surface.p, surface.q), (1, 2));
    assert_eq!(surface.u, vec![0.0, 0.0, 1.0, 1.0]);
    assert_eq!(surface.v, vec![-1.0, -1.0, -1.0, 2.0, 2.0, 2.0]);
    assert_eq!(surface.used, None);
    Ok(())
}

#[test]
fn file_entry_point_uses_the_guarded_parser() {
    // This Rust source is an existing non-SMS file: the file entry point must
    // report a contextual decoding error too, without requiring temp files.
    let path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/nlib_bezier_sms.rs");
    let result = sms::read_sms(path);
    assert!(result.is_err());
    if let Err(error) = result {
        assert!(error.to_string().contains(path));
    }
}
