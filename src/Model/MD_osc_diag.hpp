#ifndef MD_OSC_DIAG_HPP
#define MD_OSC_DIAG_HPP
/* =====================================================================
 * P11-osc PR-D1 (#434) — env-gated CVODE stepping diagnostics.
 *
 * Two independent, DEFAULT-OFF diagnostic emitters, both driven from the
 * shud.cpp coupled driver loop at the accepted CVode-return boundary
 * (one SolverStep interval per iteration):
 *
 *   SHUD_DIAG_DT=1  -> per-interval CVODE counter-delta trace
 *                      (<outpath>/diag_dt_trace.csv)
 *   SHUD_DIAG_OSC=1 -> per-entity state-delta sign-flip counters
 *                      (<outpath>/diag_osc_flips.csv + _daily.csv)
 *
 * BITWISE NEUTRALITY: the whole diagnostic path is behind a STRICT `=1`
 * env compare (strcmp == 0). Presence-only or `=0` MUST NOT enable — this
 * is deliberately stricter than the SHUD_DUMP_CV_Y presence-only precedent
 * (shud.cpp), whose block *placement* we reuse but whose predicate we do
 * not. When disabled the constructor sets both flags false and every hook
 * is a no-op, so the default hot path takes zero extra CVodeGet* calls and
 * zero I/O -> standard outputs stay byte-identical (keliya B0 SHA gate).
 *
 * SOURCE AUDIT (osc-flip-counters spec, scenario "instrumentation location
 * and source audit"): flip detection reads states ONLY as slices of the
 * accepted CVODE state vector via N_VGetArrayPointer(udata), layout
 * [sf n1 | us n1 | gw n1 | riv n2 | lake n3] per functions.hpp:83-90 with
 * n1=NumEle, n2=NumRiv, n3=NumLake. It NEVER touches the uYsf/uYus/uYgw
 * RHS scratch globals (P1e PR-B0 hazard: in the coupled driver those hold
 * the last internal f() trial evaluation at t_internal != tout). No code
 * here lives on the RHS f() call path; f.cpp / MD_rhs_core.cpp are
 * untouched. The lake slice is out of scope (no lake case in the matrix).
 * ===================================================================== */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <map>
#include <vector>
#include <string>
#include <cvode/cvode.h>
/* Macros.hpp unconditionally pulls in nvector/nvector_serial.h and documents
 * that N_VGetArrayPointer dispatches via the ops table (works for the
 * Serial and OpenMP backends alike), plus MAXLEN. */
#include "Macros.hpp"

/* Strict `=1` env predicate — the single gate semantics for both emitters.
 * getenv present AND exact string "1". NULL, "", "0", "true", "1 " etc. all
 * return false. */
inline bool osc_diag_env_is_1(const char *name) {
    const char *v = getenv(name);
    return (v != NULL) && (strcmp(v, "1") == 0);
}

class OscDiag {
public:
    /* Read the two strict `=1` gates once at construction. */
    OscDiag()
        : dt_on(osc_diag_env_is_1("SHUD_DIAG_DT")),
          osc_on(osc_diag_env_is_1("SHUD_DIAG_OSC")),
          dt_fp(NULL),
          prev_nst(0), prev_nfe(0), prev_ncfn(0), prev_netf(0),
          n_ele(0), n_riv(0), have_prev_state(false) {}

    ~OscDiag() {
        if (dt_fp != NULL) { fclose(dt_fp); dt_fp = NULL; }
    }

    bool any_on() const { return dt_on || osc_on; }

    /* Open files + write headers + allocate OSC buffers + capture the
     * initial (t = StartTime) state snapshot. Called ONCE, before the
     * driver's output-point loop. `solverstep_min` is the raw CS.SolverStep
     * (minutes); `project_name` is fout->projectname (the ./shud arg). */
    void begin(void *mem, N_Vector udata,
               const char *project_name, double solverstep_min,
               int num_ele, int num_riv, const char *outpath) {
        (void)mem;
        if (dt_on) {
            char path[MAXLEN];
            snprintf(path, sizeof(path), "%s/diag_dt_trace.csv", outpath);
            dt_fp = fopen(path, "w");
            if (dt_fp != NULL) {
                /* Header row carries project_name + solverstep_min so the
                 * analyzer derives interval-mean dt without reading cfg.para
                 * (dt-step-trace spec / design R5). */
                fprintf(dt_fp,
                        "# project_name=%s solverstep_min=%.10g\n",
                        project_name, solverstep_min);
                fprintf(dt_fp,
                        "t_next,delta_nst,delta_nfe,delta_ncfn,delta_netf,"
                        "h_last,h_cur\n");
            } else {
                fprintf(stderr,
                        "[DIAG_DT] WARNING: cannot open %s for write; "
                        "trace disabled for this run.\n", path);
            }
            /* Seed the cumulative-counter baseline so the first interval's
             * deltas are measured from run start (all counters start at 0
             * pre-integration, but read them to be robust). */
            read_counters(mem, &prev_nst, &prev_nfe, &prev_ncfn, &prev_netf);
        }
        if (osc_on) {
            n_ele = num_ele;
            n_riv = num_riv;
            osc_project = project_name ? project_name : "";
            osc_outpath = outpath ? outpath : "";
            osc_solverstep_min = solverstep_min;
            /* per-entity cumulative flip counts, one vector per family */
            flips_sf.assign((size_t)n_ele, 0);
            flips_us.assign((size_t)n_ele, 0);
            flips_gw.assign((size_t)n_ele, 0);
            flips_stage.assign((size_t)n_riv, 0);
            /* alternation state: sign of the PREVIOUS supra-epsilon delta.
             * 0 = none seen yet (sign-holding start). */
            sign_sf.assign((size_t)n_ele, 0);
            sign_us.assign((size_t)n_ele, 0);
            sign_gw.assign((size_t)n_ele, 0);
            sign_stage.assign((size_t)n_riv, 0);
            /* previous accepted-boundary state snapshot */
            const int ny = 3 * n_ele + n_riv; /* sf+us+gw+riv (lake excl.) */
            prev_state.assign((size_t)(ny > 0 ? ny : 0), 0.0);
            snapshot_state(udata);
            have_prev_state = true;
        }
    }

    /* Called once per SolverStep interval, AFTER CVode(...) returns for that
     * interval (accepted-boundary). t_next is the model time (minutes) at
     * this boundary. */
    void record(void *mem, N_Vector udata, double t_next) {
        if (dt_on && dt_fp != NULL) {
            long int nst, nfe, ncfn, netf;
            read_counters(mem, &nst, &nfe, &ncfn, &netf);
            realtype h_last = 0.0, h_cur = 0.0;
            CVodeGetLastStep(mem, &h_last);
            CVodeGetCurrentStep(mem, &h_cur);
            /* cumulative counters are monotonically non-decreasing, so
             * deltas are >= 0 by construction. */
            long int d_nst  = nst  - prev_nst;
            long int d_nfe  = nfe  - prev_nfe;
            long int d_ncfn = ncfn - prev_ncfn;
            long int d_netf = netf - prev_netf;
            fprintf(dt_fp,
                    "%.6f,%ld,%ld,%ld,%ld,%.17g,%.17g\n",
                    t_next, d_nst, d_nfe, d_ncfn, d_netf,
                    (double)h_last, (double)h_cur);
            prev_nst = nst; prev_nfe = nfe;
            prev_ncfn = ncfn; prev_netf = netf;
        }
        if (osc_on && have_prev_state) {
            update_flips(udata, t_next);
        }
    }

    /* Dump the two flip CSVs at run end. Called ONCE, after the loop. */
    void finish() {
        if (!osc_on) return;
        dump_per_entity();
        dump_daily();
    }

private:
    /* ---- configuration flags (strict `=1`) ---- */
    const bool dt_on;
    const bool osc_on;

    /* ---- DT trace state ---- */
    FILE *dt_fp;
    long int prev_nst, prev_nfe, prev_ncfn, prev_netf;

    /* per-day aggregate totals keyed by day_index = floor(t_next/1440) */
    struct DayAgg {
        long sf, us, gw, stage;
        DayAgg() : sf(0), us(0), gw(0), stage(0) {}
    };

    /* ---- OSC flip state ---- */
    int n_ele, n_riv;
    bool have_prev_state;
    double osc_solverstep_min;
    std::string osc_project;
    std::string osc_outpath;
    std::vector<double> prev_state;          /* [sf|us|gw|riv], len 3*n_ele+n_riv */
    std::vector<long> flips_sf, flips_us, flips_gw, flips_stage;
    std::vector<signed char> sign_sf, sign_us, sign_gw, sign_stage;
    std::map<long, DayAgg> daily;

    static const double EPSILON_M; /* 1e-6 m storage-height dead-floor */

    static void read_counters(void *mem, long int *nst, long int *nfe,
                              long int *ncfn, long int *netf) {
        CVodeGetNumSteps(mem, nst);
        CVodeGetNumRhsEvals(mem, nfe);
        CVodeGetNumNonlinSolvConvFails(mem, ncfn);
        CVodeGetNumErrTestFails(mem, netf);
    }

    void snapshot_state(N_Vector udata) {
        const double *y = N_VGetArrayPointer(udata);
        const int len = 3 * n_ele + n_riv;
        for (int i = 0; i < len; i++) prev_state[(size_t)i] = y[i];
    }

    /* +1 / -1 / 0(sub-epsilon) classification of an interval delta. */
    static signed char delta_sign(double d) {
        if (d >  EPSILON_M) return  1;
        if (d < -EPSILON_M) return -1;
        return 0; /* within dead-floor -> sign-holding */
    }

    /* Update one family's flip counter + alternation state for one entity.
     * sub-epsilon deltas (s==0) HOLD: neither increment nor reset prev_sign.
     * A flip is counted when the new supra-epsilon sign is opposite the
     * stored previous supra-epsilon sign. */
    static void step_family(double cur, double prev, signed char &prev_sign,
                            long &flip_count, long &day_bucket) {
        signed char s = delta_sign(cur - prev);
        if (s == 0) return; /* sign-holding: ignore sub-epsilon interval */
        if (prev_sign != 0 && s != prev_sign) {
            flip_count += 1;
            day_bucket += 1;
        }
        prev_sign = s;
    }

    void update_flips(N_Vector udata, double t_next) {
        const double *y = N_VGetArrayPointer(udata);
        long day = (long)floor(t_next / 1440.0);
        DayAgg &agg = daily[day];
        const int off_sf = 0;
        const int off_us = n_ele;
        const int off_gw = 2 * n_ele;
        const int off_rv = 3 * n_ele;
        for (int i = 0; i < n_ele; i++) {
            step_family(y[off_sf + i], prev_state[(size_t)(off_sf + i)],
                        sign_sf[(size_t)i], flips_sf[(size_t)i], agg.sf);
            step_family(y[off_us + i], prev_state[(size_t)(off_us + i)],
                        sign_us[(size_t)i], flips_us[(size_t)i], agg.us);
            step_family(y[off_gw + i], prev_state[(size_t)(off_gw + i)],
                        sign_gw[(size_t)i], flips_gw[(size_t)i], agg.gw);
        }
        for (int j = 0; j < n_riv; j++) {
            step_family(y[off_rv + j], prev_state[(size_t)(off_rv + j)],
                        sign_stage[(size_t)j], flips_stage[(size_t)j],
                        agg.stage);
        }
        /* advance the snapshot to this accepted boundary */
        const int len = 3 * n_ele + n_riv;
        for (int i = 0; i < len; i++) prev_state[(size_t)i] = y[i];
    }

    void write_flip_header(FILE *fp) const {
        fprintf(fp, "# project_name=%s solverstep_min=%.10g epsilon_m=%.10g\n",
                osc_project.c_str(), osc_solverstep_min, EPSILON_M);
    }

    void dump_per_entity() const {
        char path[MAXLEN];
        snprintf(path, sizeof(path), "%s/diag_osc_flips.csv",
                 osc_outpath.c_str());
        FILE *fp = fopen(path, "w");
        if (fp == NULL) {
            fprintf(stderr, "[DIAG_OSC] WARNING: cannot open %s\n", path);
            return;
        }
        write_flip_header(fp);
        fprintf(fp, "entity_type,entity_id,flips_sf,flips_us,flips_gw,"
                    "flips_stage\n");
        /* ele rows: element families, flips_stage=0; entity_id 1-based. */
        for (int i = 0; i < n_ele; i++) {
            fprintf(fp, "ele,%d,%ld,%ld,%ld,0\n",
                    i + 1, flips_sf[(size_t)i], flips_us[(size_t)i],
                    flips_gw[(size_t)i]);
        }
        /* riv rows: stage family only, element families=0; separate index
         * space, entity_id 1-based within rivers. */
        for (int j = 0; j < n_riv; j++) {
            fprintf(fp, "riv,%d,0,0,0,%ld\n", j + 1, flips_stage[(size_t)j]);
        }
        fclose(fp);
    }

    void dump_daily() const {
        char path[MAXLEN];
        snprintf(path, sizeof(path), "%s/diag_osc_flips_daily.csv",
                 osc_outpath.c_str());
        FILE *fp = fopen(path, "w");
        if (fp == NULL) {
            fprintf(stderr, "[DIAG_OSC] WARNING: cannot open %s\n", path);
            return;
        }
        write_flip_header(fp);
        fprintf(fp, "day_index,flips_sf,flips_us,flips_gw,flips_stage,"
                    "flips_total\n");
        /* std::map iterates day_index ascending. */
        for (std::map<long, DayAgg>::const_iterator it = daily.begin();
             it != daily.end(); ++it) {
            const DayAgg &a = it->second;
            long total = a.sf + a.us + a.gw + a.stage;
            fprintf(fp, "%ld,%ld,%ld,%ld,%ld,%ld\n",
                    it->first, a.sf, a.us, a.gw, a.stage, total);
        }
        fclose(fp);
    }
};

const double OscDiag::EPSILON_M = 1e-6;

#endif /* MD_OSC_DIAG_HPP */
