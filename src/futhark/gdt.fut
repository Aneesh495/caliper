-- ============================================================================
-- Caliper Metrology Engine: ASME Y14.5 / ISO 1101 GD&T Tolerance Engine
-- Module: gdt.fut
--
-- Formal evaluation of Form, Orientation, Location, and Runout tolerances,
-- Maximum Material Condition (MMC) / Least Material Condition (LMC) bonus
-- tolerance calculations, and automated Pass/Fail quality compliance certification.
-- ============================================================================

import "math"
import "fit"

-- ----------------------------------------------------------------------------
-- GD&T Evaluation Result Records
-- ----------------------------------------------------------------------------

type gdt_check_result =
  { measured_value: f64
  , specified_tolerance: f64
  , bonus_tolerance: f64
  , total_allowable_tolerance: f64
  , passed: bool
  , deviation_margin: f64
  }

-- ----------------------------------------------------------------------------
-- Form Tolerances (No Datum Reference Required)
-- ----------------------------------------------------------------------------

-- Flatness: Evaluates peak-to-valley separation between two parallel bounding planes
def check_flatness [n] (pts: [n]vec3) (tolerance_mm: f64): gdt_check_result =
  let fit = fit_plane pts
  let measured = fit.flatness_error
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- Cylindricity: Evaluates radial separation between two coaxial bounding cylinders
def check_cylindricity [n] (pts: [n]vec3) (tolerance_mm: f64): gdt_check_result =
  let fit = fit_cylinder pts
  let measured = fit.cylindricity_error
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- Circularity / Roundness: Radial variation of cross-sectional slice
def check_circularity [n] (pts: [n]vec3) (tolerance_mm: f64): gdt_check_result =
  let fit = fit_cylinder pts
  let measured = fit.cylindricity_error
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- ----------------------------------------------------------------------------
-- Orientation Tolerances (Relative to Datum Reference Plane or Axis)
-- ----------------------------------------------------------------------------

-- Perpendicularity: Surface normal must be 90 degrees relative to datum plane normal
-- Measured error across feature length L: L * |dot(n_feature, n_datum)|
def check_perpendicularity
    (feature_normal: vec3)
    (datum_normal: vec3)
    (feature_length_mm: f64)
    (tolerance_mm: f64): gdt_check_result =
  let cos_theta = f64.abs (v3_dot (v3_normalize feature_normal) (v3_normalize datum_normal))
  -- For perpendicular surfaces, cos(theta) should be 0.0 (theta = 90 deg)
  let measured = cos_theta * feature_length_mm
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- Parallelism: Surface normal must be 0 degrees relative to datum plane normal
-- Measured error across feature length L: L * ||n_feature x n_datum||
def check_parallelism
    (feature_normal: vec3)
    (datum_normal: vec3)
    (feature_length_mm: f64)
    (tolerance_mm: f64): gdt_check_result =
  let sin_theta = v3_norm (v3_cross (v3_normalize feature_normal) (v3_normalize datum_normal))
  let measured = sin_theta * feature_length_mm
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- Angularity: Surface normal at nominal design angle theta_nominal relative to datum
def check_angularity
    (feature_normal: vec3)
    (datum_normal: vec3)
    (theta_nominal_rad: f64)
    (feature_length_mm: f64)
    (tolerance_mm: f64): gdt_check_result =
  let dot_val = v3_dot (v3_normalize feature_normal) (v3_normalize datum_normal)
  let actual_angle = f64.acos (f64.max (-1.0) (f64.min 1.0 dot_val))
  let angle_diff = f64.abs (actual_angle - theta_nominal_rad)
  let measured = f64.sin angle_diff * feature_length_mm
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- ----------------------------------------------------------------------------
-- Location Tolerances: True Position & MMC / LMC Material Modifiers
-- ----------------------------------------------------------------------------

-- Material Condition Modifier: 0 = RFS (Regardless of Feature Size), 1 = MMC, 2 = LMC
-- For an internal cylindrical hole:
-- MMC diameter = min diameter limit (smallest hole = most material in part)
-- LMC diameter = max diameter limit (largest hole = least material in part)
def check_true_position
    (actual_center: (f64, f64))
    (nominal_center: (f64, f64))
    (actual_diameter: f64)
    (mmc_diameter: f64)
    (lmc_diameter: f64)
    (material_condition: i32)
    (base_tolerance_mm: f64): gdt_check_result =
  let dx = actual_center.0 - nominal_center.0
  let dy = actual_center.1 - nominal_center.1
  -- Diametral true position error: 2 * sqrt(dx^2 + dy^2)
  let measured = 2.0 * f64.sqrt (dx*dx + dy*dy)
  let bonus =
    if material_condition == 1 then
      -- MMC: bonus is excess departure from MMC limit
      f64.max 0.0 (actual_diameter - mmc_diameter)
    else if material_condition == 2 then
      -- LMC: bonus is departure from LMC limit
      f64.max 0.0 (lmc_diameter - actual_diameter)
    else
      0.0 -- RFS
  let total_tol = base_tolerance_mm + bonus
  let passed = measured <= total_tol
  in { measured_value = measured
     , specified_tolerance = base_tolerance_mm
     , bonus_tolerance = bonus
     , total_allowable_tolerance = total_tol
     , passed = passed
     , deviation_margin = total_tol - measured
     }

-- Coaxiality / Concentricity: Axis offset relative to Datum Axis
def check_coaxiality
    (feature_axis_pt: (f64, f64))
    (datum_axis_pt: (f64, f64))
    (tolerance_mm: f64): gdt_check_result =
  let dx = feature_axis_pt.0 - datum_axis_pt.0
  let dy = feature_axis_pt.1 - datum_axis_pt.1
  let measured = 2.0 * f64.sqrt (dx*dx + dy*dy)
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- Circular Runout: Cumulative radial displacement under 360 degree rotation
def check_runout [n]
    (pts: [n]vec3)
    (datum_axis_pt: vec3)
    (datum_axis_dir: vec3)
    (tolerance_mm: f64): gdt_check_result =
  let axis = v3_normalize datum_axis_dir
  let radial_dists =
    map (\p ->
      let diff = v3_sub p datum_axis_pt
      in v3_norm (v3_cross diff axis)
    ) pts
  let max_r = reduce f64.max (-1e30) radial_dists
  let min_r = reduce f64.min (1e30) radial_dists
  let measured = max_r - min_r
  let passed = measured <= tolerance_mm
  in { measured_value = measured
     , specified_tolerance = tolerance_mm
     , bonus_tolerance = 0.0
     , total_allowable_tolerance = tolerance_mm
     , passed = passed
     , deviation_margin = tolerance_mm - measured
     }

-- ----------------------------------------------------------------------------
-- Exported Entry Points
-- ----------------------------------------------------------------------------

entry fut_check_flatness [n]
    (xs: [n]f64, ys: [n]f64, zs: [n]f64)
    (tolerance_mm: f64): (f64, f64, bool, f64) =
  let pts = tabulate n (\i -> (xs[i], ys[i], zs[i]))
  let res = check_flatness pts tolerance_mm
  in (res.measured_value, res.total_allowable_tolerance, res.passed, res.deviation_margin)

entry fut_check_true_position
    (act_x: f64, act_y: f64, nom_x: f64, nom_y: f64)
    (act_dia: f64, mmc_dia: f64, lmc_dia: f64, condition: i32, tol: f64): (f64, f64, f64, bool) =
  let res = check_true_position (act_x, act_y) (nom_x, nom_y) act_dia mmc_dia lmc_dia condition tol
  in (res.measured_value, res.bonus_tolerance, res.total_allowable_tolerance, res.passed)
