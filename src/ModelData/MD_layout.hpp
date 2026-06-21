#ifndef MD_LAYOUT_HPP
#define MD_LAYOUT_HPP
/* MD_layout.hpp — S5d.1 (#178) ElementHotData SoA container.
 *
 * Purpose: extract _Element multi-inheritance fat-AoS hot fields (master
 * plan §4.22.1) into a separate contiguous SoA container so RHS hot path
 * loads stay cache-friendly. Field set is THE SINGLE SOURCE OF TRUTH at
 * docs/s5d_hot_fields.yaml — any drift between yaml and this header MUST
 * be caught by the CI grep gate (tools/check_manifest/check_hot_fields.py)
 * and the DEBUG assert path in Model_Data::initialize_hot().
 *
 * Field selection method (audit):
 *   The roster is the set of distinct field names accessed via
 *   `Ele[<expr>].<field>` in the three RHS hot-path TUs
 *     SHUD/src/ModelData/MD_ElementFlux.cpp
 *     SHUD/src/ModelData/MD_f.cpp
 *     SHUD/src/ModelData/MD_ET.cpp
 *   determined by the grep
 *     grep -nE 'Ele\[[^]]+\]\.' src/ModelData/MD_*.cpp \
 *       | grep -oE 'Ele\[[^]]+\]\.[A-Za-z_0-9]+(\[[^]]+\])?' \
 *       | sort -u
 *   minus the four MEMBER METHODS
 *     Flux_Infiltration, Flux_Recharge, updateElement, updateLakeElement
 *   (kept as AoS-resident methods; they write _Element members in-place
 *    and the SoA is synced from AoS via sync_hot_dynamic() AFTER each
 *    invocation, so RHS subsequent reads of u_qi/u_qex/u_effKH/u_satn
 *    see the just-updated values from the SoA). The roster reflects the
 *    actual hot-path footprint, not the master-plan §4.22.1 estimate
 *    (which over-counts by listing fields that grep does not hit).
 *
 * Contract:
 *   - Static SoA fields (geometry / topology / BC / Soil-Layer-derived
 *     scalars) are POPULATED ONCE by Model_Data::initialize_hot() from
 *     _Element AoS values during Model_Data::initialize().
 *   - Dynamic SoA fields (u_qi, u_qex, u_effKH, u_satn) are RE-POPULATED
 *     after every _Element method invocation in the RHS hot path via
 *     Model_Data::sync_hot_dynamic(i) — a per-element inline helper
 *     invoked immediately after Ele[i].updateElement(), updateLakeElement(),
 *     Flux_Infiltration() or Flux_Recharge(). This is the price of the
 *     "double-track" design (D2): _Element AoS remains source-of-truth
 *     for the writers (methods) and the SoA is a cache-friendly read
 *     mirror for the consumers (geometry-driven flux loops).
 *   - _Element AoS stays intact for init / IO / calibration. RHS hot path
 *     (MD_ElementFlux.cpp / MD_f.cpp / MD_ET.cpp) reads ONLY this SoA.
 *     Non-hot paths read _Element.
 *   - Per-Element layout uses flat int/double arrays. NumEle-sized arrays
 *     hold scalars; arrays with 3 slots are flattened as
 *     `<name>_flat[NumEle * 3]` accessed at index (3*i + j) per index
 *     convention; matches existing PassValue deterministic_gather idiom
 *     for hot-loop predictability.
 *   - Bitwise contract: ElementHotData populated values SHALL match
 *     _Element source values bit-for-bit (no rounding, no cast loss).
 *     DEBUG builds assert via initialize_hot() and sync_hot_dynamic().
 */
#include "Macros.hpp"

struct ElementHotData {
    /* All members allocated by Model_Data::malloc_EleRiv() under
     * MD_layout-managed contiguous block (S5d.3 will switch to parallel
     * first-touch). Free()d by FreeData() symmetrically. */

    /* Geometry — from class Triangle */
    int    *nabr_flat;        /* nabr[NumEle][3] flat */
    int    *lakenabr_flat;    /* lakenabr[NumEle][3] flat */
    double *edge_flat;        /* edge[NumEle][3] flat */
    double *area;             /* per-element area */
    double *z_bottom;         /* per-element aquifer bottom elevation */
    double *z_surf;           /* per-element surface elevation */

    /* Topology indices — from class AttriuteIndex (RHS-touched only) */
    int    *iSoil;            /* index into Soil[] for soil-stress lookup */
    int    *iLC;              /* index into Landcover[] for LAI tsd lookup */
    int    *iMF;              /* meltFactor tsd column index */
    int    *iForc;            /* forcing tsd index */
    int    *iLake;            /* lake-cell membership (0 = land) */
    int    *iBC;              /* boundary-condition tag (sign-encoded) */
    int    *iSS;              /* source/sink tag (sign-encoded) */

    /* Direct _Element fields */
    double *Dist2Nabor_flat;  /* Dist2Nabor[NumEle][3] flat */
    double *Dist2Edge_flat;   /* Dist2Edge[NumEle][3] flat */
    double *avgRough_flat;    /* avgRough[NumEle][3] flat */
    double *FixPressure;      /* per-element atmospheric pressure */
    double *WetlandLevel;     /* Aquiferdepth - infD */
    double *RootReachLevel;   /* Aquiferdepth - RzD */
    double *depression;       /* per-element depression storage */
    double *QBC;              /* per-element BC flux */
    double *QSS;              /* per-element source/sink flux */
    double *windH;            /* per-element wind-measurement height */
    double *u_qi;             /* infiltration; written by Flux_Infiltration */
    double *u_qex;            /* exfiltration; written by Flux_Infiltration */
    double *u_effKH;          /* horizontal effective K; written by updateElement */
    double *u_satn;           /* saturation ratio; written by updateElement */

    /* From Soil_Layer parent (RHS-touched only) */
    double *Sy;               /* specific yield; DY scale in f_applyDY */

    /* From Landcover parent (RHS-touched only) */
    double *VegFrac;          /* vegetation fraction */
    double *Albedo;           /* surface albedo for net radiation */
    double *Rough;            /* surface roughness for boundary flux */
    double *ImpAF;            /* impervious area fraction */
};

#endif /* MD_LAYOUT_HPP */
