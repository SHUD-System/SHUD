#ifndef MD_CONSISTENCY_DIAG_HPP
#define MD_CONSISTENCY_DIAG_HPP
/* Consistency diagnostics of the coupled solver loop (issue #14): an
 * acceptance tool for the fixes under #13 that does not depend on the model's
 * own output files. Default off. Driven by one-line hooks in SHUD().
 *
 * Gate, strict `=1` as in MD_osc_diag.hpp, read once at construction:
 *   SHUD_DIAG_CONSISTENCY=1  read-only; every model output stays identical.
 * When the gate is off every hook returns at once and no file is written.
 *
 * Output <outpath>/diag_consistency.txt, `key=value`, written at run end.
 * Totals over the run:
 *   macro_steps          number of solver macro steps
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
 *
 * Interface fluxes: QeleSubAt(i, j) and QeleSurfAt(i, j) are m3/min, positive
 * out of element i through edge j (rhs_apply subtracts their sum from the
 * storage of i). They are read once per macro step as CVODE's last f() call
 * left them and weighted by the step length; nothing is recomputed. Both
 * fluxes are antisymmetric by construction except for the frozen-soil factor
 * fu_Sub, which is 1 unless `cryosphere` is on, and for the `depression`
 * threshold of the surface flux, which is the element's own
 * (MD_ElementFlux.cpp:86-88) and today has the same default value in every
 * element (Element.hpp:89). With `cryosphere` off a residual must be at
 * rounding level and verifies the pairing of the two sides; with it on the
 * residual measures #20.
 *
 * Known limits:
 *   - Coupled mode only; SHUD_uncouple (-g) is not instrumented (#25).
 *   - In dummy mode (no CVode call) only the ledger runs.
 *   - Lake elements are left out of the ledger: they pass raw precipitation
 *     to the lake.
 *
 * Removal (when #13 closes): delete this header and the lines in
 * src/Model/shud.cpp that include it or name cdiag. */

#include <stdio.h>
#include <math.h>
#include <vector>
#include <cvode/cvode.h>
#include "Model_Data.hpp"
#include "MD_osc_diag.hpp" /* osc_diag_env_is_1() */

class ConsistencyDiag {
public:
    /* Constructed after the initial condition is loaded, so the snapshot
     * holds the canopy and snow stores before the first ET() call. */
    ConsistencyDiag(Model_Data *md, const char *out_dir)
        : MD(md), outpath(out_dir),
          consistency_on(osc_diag_env_is_1("SHUD_DIAG_CONSISTENCY")) {
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
     * Also counts the macro step and keeps dt for after_cvode(). */
    void after_et(double t, double tnext) {
        if (!consistency_on) return;
        macro_steps++;
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
    }

    /* Write diag_consistency.txt. Called once, after the loop. */
    void finish() const {
        if (!consistency_on) return;
        char path[MAXLEN];
        snprintf(path, sizeof(path), "%s/diag_consistency.txt", outpath);
        FILE *fp = fopen(path, "w");
        if (fp == NULL) {
            fprintf(stderr, "[DIAG_CONSISTENCY] WARNING: cannot open %s\n", path);
            return;
        }
        fprintf(fp, "macro_steps=%ld\n", macro_steps);
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
        fclose(fp);
    }

private:
    Model_Data *MD;
    const char *outpath;
    const bool consistency_on;

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

    static void put(FILE *fp, const char *key, double value) {
        fprintf(fp, "%s=%.12e\n", key, value);
    }
};

#endif /* MD_CONSISTENCY_DIAG_HPP */
