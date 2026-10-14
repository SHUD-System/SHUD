#ifndef MD_CONSISTENCY_DIAG_HPP
#define MD_CONSISTENCY_DIAG_HPP
/* Consistency diagnostics of the coupled solver loop (issue #14): an
 * acceptance tool for the fixes under #13 that does not depend on the model's
 * own output files. Default off. Driven by one-line hooks in SHUD().
 *
 * Two gates, strict `=1` as in MD_osc_diag.hpp, read once at construction:
 *   SHUD_DIAG_CONSISTENCY=1  read-only; every model output stays identical.
 *   SHUD_DIAG_RHS_REPEAT=1   calls f() twice before each CVode call; this
 *                            changes the results (see the limits below).
 * With both gates off every hook returns at once and no file is written.
 *
 * Output <outpath>/diag_consistency.txt, `key=value`, written at run end.
 * Totals over the run. Always:
 *   macro_steps          number of solver macro steps
 * SHUD_DIAG_CONSISTENCY:
 *   land_elements        elements with iLake <= 0
 *   prcp_m3              precipitation (rain and snow) on the land elements
 *   netprcp_m3           net precipitation passed on to the land surface
 *   canopy_evap_m3       evaporation from canopy interception, as ET() set it
 *   d_canopy_m3          change of canopy interception storage
 *   d_snow_m3            change of snow storage
 *   residual_m3          prcp - netprcp - canopy_evap - d_canopy - d_snow
 *   residual_frac_prcp   residual_m3 / prcp_m3 (0 when there is none)
 *   cvode_overshoot_steps    macro steps after which CVODE's internal time
 *                            is past the output time tnext
 *   cvode_max_overshoot_min  largest such excess, minutes
 *   landland_edges             edges shared by two land elements
 *   landland_sub_residual_m3   sum of |Q_ij + Q_ji| * dt, subsurface flux
 *   landland_sub_gross_m3      sum of |Q_ij| * dt, subsurface flux
 *   landland_surf_residual_m3  the same for the surface flux
 *   landland_surf_gross_m3
 *   lakes                      number of lakes
 *   landlake_sub_residual_m3   sum over lakes of |land side - lake side| * dt
 *   landlake_sub_gross_m3      sum over lakes of |lake side| * dt
 *   landlake_surf_residual_m3  the same for the surface flux
 *   landlake_surf_gross_m3
 *
 * Interface fluxes: QeleSubAt(i, j) and QeleSurfAt(i, j) are m3/min, positive
 * out of element i through edge j (rhs_apply subtracts their sum from the
 * storage of i). QLakeSub and QLakeSurf are m3/min, positive into the lake
 * (rhs_apply adds them to the lake storage). So on a lake edge positive is
 * land to lake on both sides. The fluxes are read once per macro step as
 * CVODE's last f() call left them and weighted by the step length; nothing
 * is recomputed. Both fluxes are antisymmetric by construction except for
 * the frozen-soil factor fu_Sub, which is 1 unless `cryosphere` is on, and
 * for the `depression` threshold of the surface flux, which is the element's
 * own (MD_ElementFlux.cpp:86-88) and today has the same default value in
 * every element (Element.hpp:89). With `cryosphere` off a residual must be
 * at rounding level and verifies the pairing of the two sides; with it on
 * the residual measures #20.
 *
 * SHUD_DIAG_RHS_REPEAT:
 *   rhs_repeat_calls      macro steps in which f() was called twice at the
 *                         same (t, y)
 *   rhs_repeat_differing  those in which any component of DY differed
 *   rhs_repeat_max_abs    largest |DY second - DY first|, in the unit of DY
 *                         (m/min for elements, rivers and lakes)
 *
 * Known limits:
 *   - The repeat check adds two f() calls per macro step. While f() depends
 *     on its call history (#16), a run with SHUD_DIAG_RHS_REPEAT=1 does not
 *     reproduce the default trajectory, and nfcall.txt counts the extra
 *     calls. With both gates on, the canopy and snow keys are unaffected
 *     (they depend on the forcing only); the cvode_* and interface keys and
 *     all model outputs then describe the perturbed run.
 *   - Coupled mode only; SHUD_uncouple (-g) is not instrumented (#25).
 *   - In dummy mode (no CVode call) only the ledger runs.
 *   - The lake keys are not covered by `make regress` (ccw has no lake).
 *   - Lake elements are left out of the ledger: they pass raw precipitation
 *     to the lake.
 *
 * Removal (when #13 closes): delete
 *   - this header;
 *   - the six lines in src/Model/shud.cpp that include it or name cdiag;
 *   - in tools/ci/runtime_checks.py: the constants DIAG_FILE, LEDGER_KEYS,
 *     COUNT_KEYS and INTERFACES with their comment; the functions read_keys,
 *     compare_ledger, check_closure and check_diagnostics; the call of
 *     check_diagnostics, so that gate_regress ends with `return failures`
 *     again; the parameter extra_env of run_case and the line that uses it;
 *     in gate_coverage, the `--exclude` of this header and the comment above
 *     the gcovr call, which fits on one line again;
 *   - the `ledger` dict and the "ledger tolerance" and "ledger closure" rows
 *     in tools/ci/test_guardrails.py;
 *   - regression.ledger_rel, with its comment, in constraints.yaml;
 *   - tests/reference/ccw/diag_consistency.ref. */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <cvode/cvode.h>
#include "f.hpp"
#include "Model_Data.hpp"
#include "MD_adjacency.hpp" /* lake_bank_edge_by_lake: (element, edge) pairs on each lake */
#include "MD_osc_diag.hpp" /* osc_diag_env_is_1() */

class ConsistencyDiag {
public:
    /* Constructed after the initial condition is loaded, so the snapshot
     * holds the canopy and snow stores before the first ET() call. */
    ConsistencyDiag(Model_Data *md, const char *out_dir)
        : MD(md), outpath(out_dir),
          consistency_on(osc_diag_env_is_1("SHUD_DIAG_CONSISTENCY")),
          repeat_on(osc_diag_env_is_1("SHUD_DIAG_RHS_REPEAT")) {
        if (!consistency_on) return;
        canopy_before.assign(MD->yEleIS, MD->yEleIS + MD->NumEle);
        snow_before.assign(MD->yEleSnow, MD->yEleSnow + MD->NumEle);
        for (int i = 0; i < MD->NumEle; i++) {
            if (is_land(i)) land_elements++;
            for (int j = 0; j < 3; j++) {
                if (is_landland_edge(i, j)) landland_edges++;
            }
        }
    }

    /* Canopy and snow water ledger of one macro step. Called right after
     * ET(t, tnext): ET() has just moved both stores over dt = tnext - t and
     * set qEleNetPrep and qEleE_IC for the step. No f() has run yet in this
     * step; f() overwrites qEleE_IC. Only ET() writes the two stores inside
     * the loop, so the values kept here are the ones before the next ET().
     * Rates are m/min and stores m; times the element area they are m3.
     * Also counts the macro step (for both gates) and keeps dt for
     * after_cvode(). */
    void after_et(double t, double tnext) {
        if (!consistency_on && !repeat_on) return;
        macro_steps++;
        if (!consistency_on) return;
        dt = tnext - t;
        for (int i = 0; i < MD->NumEle; i++) {
            if (!is_land(i)) continue;
            const double area = MD->Ele[i].area;
            prcp        += MD->qElePrep[i] * dt * area;
            netprcp     += MD->qEleNetPrep[i] * dt * area;
            canopy_evap += MD->qEleE_IC[i] * dt * area;
            d_canopy    += (MD->yEleIS[i] - canopy_before[i]) * area;
            d_snow      += (MD->yEleSnow[i] - snow_before[i]) * area;
            canopy_before[i] = MD->yEleIS[i];
            snow_before[i] = MD->yEleSnow[i];
        }
    }

    /* Called right before CVode(): evaluates the right-hand side twice at
     * the same (t, y) and compares the two DY. A right-hand side that is a
     * function of (t, y) only gives the same DY both times. `dy` is the
     * scratch vector of SHUD(), which nothing else reads. */
    void rhs_repeat(double t, N_Vector y, N_Vector dy) {
        if (!repeat_on) return;
        const double *DY = N_VGetArrayPointer(dy);
        f(t, y, dy, MD);
        first_dy.assign(DY, DY + MD->NumY);
        f(t, y, dy, MD);
        double max_abs = 0.;
        for (int i = 0; i < MD->NumY; i++) {
            if (fabs(DY[i] - first_dy[i]) > max_abs) max_abs = fabs(DY[i] - first_dy[i]);
        }
        repeat_calls++;
        if (max_abs > 0.) repeat_differing++;
        if (max_abs > repeat_max_abs) repeat_max_abs = max_abs;
    }

    /* Called right after CVode() returned for the output time tnext. In
     * CV_NORMAL mode CVODE steps past tnext and interpolates back, so the
     * model arrays written by its last f() call belong to a later time. */
    void after_cvode(void *mem, double tnext) {
        if (!consistency_on) return;
        realtype tcur = tnext;
        CVodeGetCurrentTime(mem, &tcur);
        if (tcur > tnext) {
            overshoot_steps++;
            if (tcur - tnext > max_overshoot) max_overshoot = tcur - tnext;
        }
        land_land();
        land_lake();
    }

    /* Write diag_consistency.txt. Called once, after the loop. */
    void finish() const {
        if (!consistency_on && !repeat_on) return;
        char path[MAXLEN];
        snprintf(path, sizeof(path), "%s/diag_consistency.txt", outpath);
        FILE *fp = fopen(path, "w");
        if (fp == NULL) {
            fprintf(stderr, "[DIAG_CONSISTENCY] WARNING: cannot open %s\n", path);
            return;
        }
        fprintf(fp, "macro_steps=%ld\n", macro_steps);
        if (consistency_on) write_consistency(fp);
        if (repeat_on) {
            fprintf(fp, "rhs_repeat_calls=%ld\n", repeat_calls);
            fprintf(fp, "rhs_repeat_differing=%ld\n", repeat_differing);
            put(fp, "rhs_repeat_max_abs", repeat_max_abs);
        }
        fclose(fp);
    }

private:
    void write_consistency(FILE *fp) const {
        const double residual = prcp - netprcp - canopy_evap - d_canopy - d_snow;
        fprintf(fp, "land_elements=%d\n", land_elements);
        put(fp, "prcp_m3", prcp);
        put(fp, "netprcp_m3", netprcp);
        put(fp, "canopy_evap_m3", canopy_evap);
        put(fp, "d_canopy_m3", d_canopy);
        put(fp, "d_snow_m3", d_snow);
        put(fp, "residual_m3", residual);
        put(fp, "residual_frac_prcp", prcp > 0. ? residual / prcp : 0.);
        fprintf(fp, "cvode_overshoot_steps=%ld\n", overshoot_steps);
        put(fp, "cvode_max_overshoot_min", max_overshoot);
        fprintf(fp, "landland_edges=%d\n", landland_edges);
        put(fp, "landland_sub_residual_m3", landland_sub_residual);
        put(fp, "landland_sub_gross_m3", landland_sub_gross);
        put(fp, "landland_surf_residual_m3", landland_surf_residual);
        put(fp, "landland_surf_gross_m3", landland_surf_gross);
        fprintf(fp, "lakes=%d\n", MD->NumLake);
        put(fp, "landlake_sub_residual_m3", landlake_sub_residual);
        put(fp, "landlake_sub_gross_m3", landlake_sub_gross);
        put(fp, "landlake_surf_residual_m3", landlake_surf_residual);
        put(fp, "landlake_surf_gross_m3", landlake_surf_gross);
    }

    Model_Data *MD;
    const char *outpath;
    const bool consistency_on;
    const bool repeat_on;

    long macro_steps = 0;
    double dt = 0.; /* min; set by after_et(), read by after_cvode() in the same macro step */

    int land_elements = 0;
    std::vector<double> canopy_before, snow_before; /* stores before ET(), m */
    double prcp = 0., netprcp = 0., canopy_evap = 0., d_canopy = 0., d_snow = 0.; /* m3 */

    long overshoot_steps = 0;
    double max_overshoot = 0.; /* min */

    int landland_edges = 0;
    double landland_sub_residual = 0., landland_sub_gross = 0.;   /* m3 */
    double landland_surf_residual = 0., landland_surf_gross = 0.; /* m3 */
    double landlake_sub_residual = 0., landlake_sub_gross = 0.;   /* m3 */
    double landlake_surf_residual = 0., landlake_surf_gross = 0.; /* m3 */

    std::vector<double> first_dy; /* DY of the first of the two f() calls */
    long repeat_calls = 0, repeat_differing = 0;
    double repeat_max_abs = 0.;

    bool is_land(int i) const { return MD->Ele[i].iLake <= 0; }

    /* Edge j of element i is shared with a land element of higher index, so
     * each interior edge is taken once. nabr and nabrToMe are 1-based;
     * nabrToMe is -1 when the neighbour does not list element i. */
    bool is_landland_edge(int i, int j) const {
        const _Element &e = MD->Ele[i];
        const int k = e.nabr[j] - 1;
        return is_land(i) && e.lakenabr[j] <= 0 && k > i && is_land(k) && e.nabrToMe[j] > 0;
    }

    /* Antisymmetry of the land-land fluxes: Q_ij of element i through its
     * edge j against Q_ji of the neighbour through the same edge. */
    void land_land() {
        for (int i = 0; i < MD->NumEle; i++) {
            for (int j = 0; j < 3; j++) {
                if (!is_landland_edge(i, j)) continue;
                const int k = MD->Ele[i].nabr[j] - 1;
                const int kj = MD->Ele[i].nabrToMe[j] - 1;
                landland_sub_residual += fabs(MD->QeleSubAt(i, j) + MD->QeleSubAt(k, kj)) * dt;
                landland_sub_gross += fabs(MD->QeleSubAt(i, j)) * dt;
                landland_surf_residual += fabs(MD->QeleSurfAt(i, j) + MD->QeleSurfAt(k, kj)) * dt;
                landland_surf_gross += fabs(MD->QeleSurfAt(i, j)) * dt;
            }
        }
    }

    /* Per lake: what the land elements book on their lake edges (the flux
     * their storage loses) against what the lake receives. */
    void land_lake() {
        for (int l = 0; l < MD->NumLake; l++) {
            double land_sub = 0., land_surf = 0.;
            for (const std::pair<int, int> &edge : lake_bank_edge_by_lake[l]) {
                land_sub += MD->QeleSubAt(edge.first, edge.second);
                land_surf += MD->QeleSurfAt(edge.first, edge.second);
            }
            landlake_sub_residual += fabs(land_sub - MD->QLakeSub[l]) * dt;
            landlake_sub_gross += fabs(MD->QLakeSub[l]) * dt;
            landlake_surf_residual += fabs(land_surf - MD->QLakeSurf[l]) * dt;
            landlake_surf_gross += fabs(MD->QLakeSurf[l]) * dt;
        }
    }

    static void put(FILE *fp, const char *key, double value) {
        fprintf(fp, "%s=%.12e\n", key, value);
    }
};

#endif /* MD_CONSISTENCY_DIAG_HPP */
