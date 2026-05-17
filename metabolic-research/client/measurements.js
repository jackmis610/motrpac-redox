/* YOUR biomarker measurements — this is the one file you edit.
 *
 * window.MEASUREMENTS is a list of panels. Each panel is one point in time:
 *   { date: "YYYY-MM-DD", source: "where it came from", values: { ... } }
 * `values` maps biomarker id -> your result, in the units the dashboard shows
 * for that biomarker (open the dashboard; each row states its expected units).
 *
 * Submit as many or as few biomarkers as you have. Anything outside the
 * 135-biomarker vocabulary, or without a reference range yet, is still stored
 * and listed as "tracked" — it just is not scored.
 *
 * The dashboard is useful from a single panel; a second panel unlocks the
 * trajectory view. The two panels below are EXAMPLE data — replace them.
 */
window.MEASUREMENTS = [
  {
    date: "2025-11-01",
    source: "example — replace with your own",
    values: {
      apob: 110, ldl_c: 130, hdl_c: 48, triglycerides: 140, non_hdl_c: 152,
      hba1c: 5.6, fasting_glucose: 95, fasting_insulin: 11,
      hscrp: 1.8, vo2max: 42, grip_strength: 44,
      resting_bp: 128, waist_circumference: 92, vitamin_d_25oh: 28
    }
  },
  {
    date: "2026-02-15",
    source: "example — replace with your own",
    values: {
      apob: 88, ldl_c: 100, hdl_c: 52, triglycerides: 95, non_hdl_c: 116,
      hba1c: 5.4, fasting_glucose: 90, fasting_insulin: 8,
      hscrp: 1.1, vo2max: 46, grip_strength: 46,
      resting_bp: 122, waist_circumference: 89, vitamin_d_25oh: 39
    }
  }
];
