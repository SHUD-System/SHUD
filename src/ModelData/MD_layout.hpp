#ifndef MD_LAYOUT_HPP
#define MD_LAYOUT_HPP
/* MD_layout.hpp — ElementHotData structure-of-arrays (SoA) container.
 *
 * Purpose: copy the hot fields of the multi-inheritance, fat
 * array-of-structures (AoS) _Element into a separate contiguous SoA
 * container so RHS hot path loads stay cache-friendly. This header is
 * the definitive list of hot fields; malloc_EleRiv(), initialize_hot()
 * and FreeData() must be kept in step with it. DEBUG builds check the
 * copied values in Model_Data::initialize_hot().
 *
 * Field selection method:
 *   The roster is the set of distinct field names accessed via
 *   `Ele[<expr>].<field>` in the three RHS hot-path TUs
 *     src/ModelData/MD_ElementFlux.cpp
 *     src/ModelData/MD_f.cpp
 *     src/ModelData/MD_ET.cpp
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
 *    actual hot-path footprint: fields the grep does not hit are not
 *    mirrored.
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
 *     "double-track" design: _Element AoS remains source-of-truth
 *     for the writers (methods) and the SoA is a cache-friendly read
 *     mirror for the consumers (geometry-driven flux loops).
 *   - _Element AoS stays intact for init / IO / calibration. RHS hot path
 *     (MD_ElementFlux.cpp / MD_f.cpp / MD_ET.cpp) reads ONLY this SoA.
 *     Non-hot paths read _Element.
 *   - Per-Element layout uses flat int/double arrays. NumEle-sized arrays
 *     hold scalars; arrays with 3 slots are flattened as
 *     `<name>_flat[NumEle * 3]` accessed at index (3*i + j) per index
 *     convention; matches the existing rhs_deterministic_gather() row-major
 *     idiom for hot-loop predictability.
 *   - Bitwise contract: ElementHotData populated values SHALL match
 *     _Element source values bit-for-bit (no rounding, no cast loss).
 *     DEBUG builds assert via initialize_hot() and sync_hot_dynamic().
 */
#include "Macros.hpp"

struct ElementHotData {
    /* All members are allocated by Model_Data::malloc_EleRiv(), one
     * contiguous array each (with an optional parallel first-touch
     * there), and freed symmetrically by FreeData().
     *
     * Every pointer field defaults to nullptr so that
     * Model_Data::FreeData()'s unconditional `delete[] hot.*` chain is
     * a defined no-op when malloc_EleRiv() never ran or aborted
     * partway (e.g. out of memory in mid-allocation). Without these
     * defaults, the pointer slots would hold indeterminate bytes and
     * `delete[]` on those triggers `free(): invalid pointer` heap
     * corruption (glibc) or SEGV (other allocators). The defaults are
     * overwritten by the `hot.* = new ...` assignments in
     * Model_Data::malloc_EleRiv(), so they do not affect model
     * output. */

    /* Geometry — from class Triangle */
    int    *nabr_flat = nullptr;        /* nabr[NumEle][3] flat */
    int    *lakenabr_flat = nullptr;    /* lakenabr[NumEle][3] flat */
    double *edge_flat = nullptr;        /* edge[NumEle][3] flat */
    double *area = nullptr;             /* per-element area */
    double *z_bottom = nullptr;         /* per-element aquifer bottom elevation */
    double *z_surf = nullptr;           /* per-element surface elevation */

    /* Topology indices — from class AttriuteIndex (RHS-touched only) */
    int    *iSoil = nullptr;            /* index into Soil[] for soil-stress lookup */
    int    *iLC = nullptr;              /* index into Landcover[] for LAI tsd lookup */
    int    *iMF = nullptr;              /* meltFactor tsd column index */
    int    *iForc = nullptr;            /* forcing tsd index */
    int    *iLake = nullptr;            /* lake-cell membership (0 = land) */
    int    *iBC = nullptr;              /* boundary-condition tag (sign-encoded) */
    int    *iSS = nullptr;              /* source/sink tag (sign-encoded) */

    /* Direct _Element fields */
    double *Dist2Nabor_flat = nullptr;  /* Dist2Nabor[NumEle][3] flat */
    double *Dist2Edge_flat = nullptr;   /* Dist2Edge[NumEle][3] flat */
    double *avgRough_flat = nullptr;    /* avgRough[NumEle][3] flat */
    double *FixPressure = nullptr;      /* per-element atmospheric pressure */
    double *WetlandLevel = nullptr;     /* Aquiferdepth - infD */
    double *RootReachLevel = nullptr;   /* Aquiferdepth - RzD */
    double *depression = nullptr;       /* per-element depression storage */
    double *QBC = nullptr;              /* per-element BC flux */
    double *QSS = nullptr;              /* per-element source/sink flux */
    double *windH = nullptr;            /* per-element wind-measurement height */
    double *u_qi = nullptr;             /* infiltration; written by Flux_Infiltration */
    double *u_qex = nullptr;            /* exfiltration; written by Flux_Infiltration */
    double *u_effKH = nullptr;          /* horizontal effective K; written by updateElement */
    double *u_satn = nullptr;           /* saturation ratio; written by updateElement */

    /* From Soil_Layer parent (RHS-touched only) */
    double *Sy = nullptr;               /* specific yield; DY scale in f_applyDY */

    /* From Landcover parent (RHS-touched only) */
    double *VegFrac = nullptr;          /* vegetation fraction */
    double *Albedo = nullptr;           /* surface albedo for net radiation */
    double *Rough = nullptr;            /* surface roughness for boundary flux */
    double *ImpAF = nullptr;            /* impervious area fraction */
};

#endif /* MD_LAYOUT_HPP */
