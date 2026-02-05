-- ============================================================================
-- Caliper Metrology Engine: Unit Tests for GD&T Engine
-- Module: tests/futhark/test_gdt.fut
--
-- Tests ASME Y14.5 / ISO 1101 Flatness, Perpendicularity, and True Position
-- with Maximum Material Condition (MMC) bonus tolerance evaluation.
-- ============================================================================

import "../../src/futhark/gdt"
import "../../src/futhark/math"

def test_flatness_evaluation: bool =
  -- Planar grid with 0.03 mm peak-to-valley sinusoidal ripple (tool chatter)
  let n = 25i64
  let pts =
    tabulate n (\i ->
      let x = f64.i64 (i % 5) * 10.0
      let y = f64.i64 (i / 5) * 10.0
      let ripple = 0.015 * f64.sin (x * 0.5) -- -0.015 to +0.015 -> 0.030 mm total
      in (x, y, ripple)
    )
  let res_pass = check_flatness pts 0.050 -- 50 micron tolerance -> PASS
  let res_fail = check_flatness pts 0.010 -- 10 micron tolerance -> FAIL
  let ok_pass = res_pass.passed && res_pass.measured_value <= 0.035
  let ok_fail = !res_fail.passed
  in ok_pass && ok_fail

def test_perpendicularity: bool =
  -- Perfect 90 degree surface normal relative to datum +Z normal
  let datum_n: vec3 = (0.0, 0.0, 1.0)
  let feat_n_ortho: vec3 = (1.0, 0.0, 0.0)
  let length_mm = 100.0
  let res_ortho = check_perpendicularity feat_n_ortho datum_n length_mm 0.050
  -- Tilted surface by 1 degree
  let feat_n_tilt: vec3 = v3_normalize (1.0, 0.0, f64.sin (1.0 * f64.pi / 180.0))
  let res_tilt = check_perpendicularity feat_n_tilt datum_n length_mm 0.050
  let ok_ortho = res_ortho.passed && res_ortho.measured_value < 1e-6
  let ok_tilt = !res_tilt.passed -- 1 deg tilt across 100 mm is ~1.74 mm >> 0.05 mm
  in ok_ortho && ok_tilt

def test_true_position_mmc: bool =
  -- Nominal hole location (50.0, 50.0), nominal diameter 10.00 mm
  let nom_xy = (50.0, 50.0)
  -- Actual hole drilled at (50.02, 50.01) -> radial offset = sqrt(0.02^2 + 0.01^2) = 0.02236
  -- Diametral error = 2 * 0.02236 = 0.04472 mm
  let act_xy = (50.02, 50.01)
  let act_dia = 10.12 -- Internal hole produced larger than MMC limit
  let mmc_dia = 10.00
  let lmc_dia = 10.25
  let base_tol = 0.040 -- Base position tolerance 40 microns (error is 44.7 microns)

  -- Check under RFS (no bonus): 44.7 > 40 -> FAIL
  let res_rfs = check_true_position act_xy nom_xy act_dia mmc_dia lmc_dia 0 base_tol
  -- Check under MMC: bonus = 10.12 - 10.00 = 0.120 mm -> total tol = 0.160 mm -> PASS!
  let res_mmc = check_true_position act_xy nom_xy act_dia mmc_dia lmc_dia 1 base_tol

  let ok_rfs = !res_rfs.passed && res_rfs.bonus_tolerance == 0.0
  let ok_mmc = res_mmc.passed && f64.abs (res_mmc.bonus_tolerance - 0.12) < 1e-6
  in ok_rfs && ok_mmc

entry main: bool =
  test_flatness_evaluation && test_perpendicularity && test_true_position_mmc
