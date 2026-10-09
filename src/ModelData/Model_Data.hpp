//  Model_Data.hpp
//  Created by Lele Shu (lele.shu@gmail.com) on 2018.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#ifndef Model_Data_hpp
#define Model_Data_hpp

#include <stdio.h>
#include "TimeSeriesData.hpp"
#include "ModelConfigure.hpp"
#include "IO.hpp"
#include "River.hpp"
#include "Element.hpp"
#include "Model_Control.hpp"
#include "TabularData.hpp"
#include "FloodAlert.hpp"
#include "Lake.hpp"
#include "is_sm_et.hpp"
#include "Flux_RiverElement.hpp"
#include "Macros.hpp"
#include "AccTemperature.hpp"
#include "MD_layout.hpp" /* ElementHotData SoA */
using namespace std;

/* ExecPolicy
 *
 * Selects the execution backend for `Model_Data::rhs_core(Y, DY, t,
 * policy)`. Defined here (rather than in MD_rhs_core.hpp) because
 * `Model_Data::rhs_core`'s declaration takes it by value and
 * MD_rhs_core.hpp `#include`s this file (circular dependency if it
 * lived only there). MD_rhs_core.hpp re-exposes it via its own
 * `#include "Model_Data.hpp"`.
 *
 * Design:
 *   - plain `enum class` + `switch(policy)` inside `rhs_core`
 *   - NO template specialization, NO virtual dispatch (compile-time-
 *     known policy keeps binary size flat + dispatch predictable)
 *   - Serial = plain call chain; StrictOMP = the same chain inside one
 *     OpenMP parallel region (only compiled with
 *     SHUD_ENABLE_OPENMP_RHS); ProductionOMP and any policy that is
 *     not compiled in call `std::abort()` (NOT `assert(false)` —
 *     `-DNDEBUG` strips assert to a no-op and would let execution
 *     fall through silently)
 */
enum class ExecPolicy { Serial, StrictOMP, ProductionOMP };

class Model_Data {        /* Model_data definition */
public:
    FileIn  *pf_in;
    FileOut *pf_out;
    double t0;  /* Time tag before current f loop */
    double t1;  /* Time tag before current iteration */
    double tnow;       /* Current time tag*/
    double dt;     /* DT = tnow - t1 */
    //int UnsatMode;        /* Unsat Mode */
    //int SurfMode;        /* Surface Overland Flow Mode */
    //int RivMode;        /* River Routing Mode */
//    unsigned long nFCall1;
//    unsigned long nFCall2;
    char file_debug[MAXLEN];
    unsigned long nFCall = 0;
    
    unsigned long nFCall1 = 0;
    unsigned long nFCall2 = 0;
    unsigned long nFCall3 = 0;
    unsigned long nFCall4 = 0;
    unsigned long nFCall5 = 0;
    
    double tic;
    int NumForc = 0;
    int NumEle = 0;            /* Number of Elements */
    int NumNode = 0;        /* Number of Nodes */
//    int NumY1;
//    int NumY2;
    int NumY = 0;
    int NumSSEle = 0;       /* Number of Souce/Sink for Elements */
    int NumBCEle1 = 0;      /* Number of Boundary Condition for Elements */
    int NumBCEle2 = 0;      /* Number of Boundary Condition for Elements */
    int NumBCRiv1 = 0;      /* Number of Boundary Condition for Rivers */
    int NumBCRiv2 = 0;      /* Number of Boundary Condition for Rivers */
    int NumBCLake1 = 0;      /* Number of Boundary Condition for Lakes */
    int NumBCLake2 = 0;      /* Number of Boundary Condition for Lakes */
    int NumRiv = 0;            /* Number of Rivere Reaches */
    
    
    int NumSoil = 0;        /* Number of Soils */
    int NumGeol = 0;        /* Number of Geologies */
    int NumLC = 0;            /* Number of Land Cover Index Data */
    int NumMeltF = 0;        /* Number of Melt Factor Time series */
    int NumRivType = 0;        /* Number of River Shape */
    int NumRivNode = 0;
    /* io_riv / io_lake are allocated conditionally in MD_readin.cpp
     * read_cfgout() (NumRiv>0 / NumLake>0). Initialize to nullptr so
     * FreeData()'s unconditional `delete[]` is a defined no-op when
     * the project has no river or no lake (delete[] on an
     * uninitialized pointer is undefined behavior). io_ele is always
     * allocated, but initialize for symmetry. */
    int *io_ele = nullptr;
    int *io_riv = nullptr;
    int *io_lake = nullptr; /* Wether Export the data of these elements */
    
    /* The heap pointer members below default to nullptr so that
     * Model_Data::FreeData()'s unconditional `delete[]` / `delete`
     * chain is a defined no-op when allocation never ran (e.g. out of
     * memory part-way through malloc_EleRiv; ctor short-circuit on
     * exception). Neither constructor sets them: the default ctor
     * `Model_Data()` is empty and `Model_Data(FileIn*, FileOut*)` sets
     * only pf_in/pf_out. They are allocated in
     * malloc_EleRiv()/malloc_Y()/MD_Lake.cpp/MD_readin.cpp, where the
     * assignment overrides the default, so the defaults do not affect
     * model output. */
    _TimeSeriesData *tsd_weather = nullptr;
    _TimeSeriesData tsd_LAI;
//    _TimeSeriesData tsd_RL;
    _TimeSeriesData tsd_MF;
    _TimeSeriesData tsd_eleSS; /* Element Source/Sink Term [L3/T] */
    _TimeSeriesData tsd_eyBC; /* Element Y BC */
    _TimeSeriesData tsd_eqBC; /* Element Q BC */
    _TimeSeriesData tsd_ryBC;
    _TimeSeriesData tsd_rqBC;
    _TimeSeriesData tsd_lyBC;
    _TimeSeriesData tsd_lqBC;
    int ieBC1 = 0;
    int ieBC2 = 0;
    int irBC1 = 0;
    int irBC2 = 0;
    int ilBC1 = 0;
    int ilBC2 = 0;
    int ieSS = 0;
    
    globalCal gc;
    Control_Data CS;
    
    _Element *Ele = nullptr;        /* Store Element Information */
    /* ElementHotData SoA hot field container. See
     * MD_layout.hpp. Populated by
     * initialize_hot() called from initialize() AFTER _Element AoS data
     * is fully loaded. Dynamic fields (u_qi, u_qex, u_effKH, u_satn) are
     * resynced by sync_hot_dynamic(i) after each Ele[i].updateElement /
     * updateLakeElement / Flux_Infiltration / Flux_Recharge call. RHS hot
     * path (MD_ElementFlux/MD_f/MD_ET) reads this; init/IO/calib path
     * reads _Element. */
    ElementHotData hot;
    _Node *Node = nullptr;        /* Store Node Information */
    //element_IC * Ele_IC;    /* Store Element Initial Condtion */
    Soil_Layer *Soil = nullptr;        /* Store Soil Information */
    Geol_Layer *Geol = nullptr;            /* Store Geology Information */
    Landcover *LandC = nullptr;        /* Store Land Cover Information */

    _River *Riv = nullptr;        /* Store River Reach Information */
    river_para *Riv_Type = nullptr;    /* Store River Shape Information */
    _Node *rivNode = nullptr;
    FloodAlert *flood = nullptr;

    double *fu_Surf = nullptr; /* Fraction of unfrozen landsurface */
    double *fu_Sub = nullptr; /* Fraction of unfrozen subsurface */
    _AccTemp *AccT_surf = nullptr;
    _AccTemp *AccT_sub = nullptr;
    double AccT_sub_max = 10;
    double AccT_sub_min = -10;
    double AccT_surf_max = 3;
    double AccT_surf_min = -3;
    
    double WatershedArea = 0.;
    // `ISFactor` and `windH` (top-level Model_Data) are unused fields:
    // they are never assigned, dereferenced, allocated or freed (the live
    // wind-height data flows through `Ele[i].windH` + `hot.windH`), so
    // they carry no nullptr default.
    double *ISFactor;        /* ISFactor is used to calculate ISMax from LAI */
    double *windH;        /* Height at which wind velocity is measured */
    _Lake *lake = nullptr;
    int NumLake = 0;
    double *QoutSurf = nullptr;

    /* QeleSurf/QeleSub are stored as one contiguous row-major
     * `double[NumEle*3]` block per array. Index convention matches
     * the MD_layout.hpp flat-3 idiom: at(i,j) ↔ `_flat[3*i + j]`.
     * Access via the inline `QeleSurfAt(i,j)` / `QeleSubAt(i,j)`
     * accessors below; hot-path call sites must use the accessor (no
     * bare `_flat[3*i + j]` indices). The PrintCtrl IO path uses the
     * flat overload InitIJ(...double *x_flat, int j, ...) — each
     * output file is 1 column out of 3 with stride 3 across rows; see
     * Model_Control.cpp. */
    double *QeleSurf_flat = nullptr; /* Overland Flux — flat NumEle*3 */
    double *QeleSub_flat = nullptr;  /* Subsurface Flux — flat NumEle*3 */
    //double ** FluxRiv;    /* River Segment Flux */
    double *QrivSurf = nullptr;        /* surface Flux between river and element */
    double *QrivSub = nullptr;        /* gw Flux between river and element */
    double *QrivDown = nullptr;
    double *QrivUp = nullptr;
    double *QsegSurf = nullptr;
    double *QsegSub = nullptr;

    double *QeleSurfTot = nullptr;
    double *QeleSubTot = nullptr;

    double *Qe2r_Surf = nullptr;
    double *Qe2r_Sub = nullptr;

    double *yEleWetFront = nullptr;        /* Weting Front */

    double *qElePrep = nullptr;        /* Precep. on each element */
    double *qEleETloss = nullptr;
    double *qEleNetPrep = nullptr;    /* Net precep. on each elment */
    double *qEleInfil = nullptr;    /* Variable infiltration rate */
    double *qEleExfil = nullptr;    /* Variable exfiltration rate */
    double *qEleRecharge = nullptr;    /* Recharge rate to GW */
    double *yEleSnowGrnd = nullptr;    /* Snow depth on ground element */
    double *yEleSnowCanopy = nullptr;    /* Snow depth on canopy element */
    double *yEleISmax = nullptr;    /* Maximum interception storage (liquid
                           * precep) */
    double *yEleISsnowmax = nullptr;    /* Maximum interception storage (snow) */
    double *qEleTF = nullptr;        /* Through Fall */

    double *yEleIS = nullptr;        /* Interception storage */
    double *yEleSnow = nullptr;        /* Snow depth on each element */
    double *yEleGW = nullptr;   // debug may not necessary
    double *yEleSurf = nullptr;   // debug may not necessary
    double *yEleUnsat = nullptr;   // debug may not necessary
//    double *yEleSM;   // Soil Moisture Ratio
    double *qEleETP = nullptr;    /* Potential ET  = qPotEvap * (1-VegFrac)+ qPotTran * VegFrac */
    double *qPotEvap = nullptr;   /* Potential Evaporation of Soil */
    double *qPotTran = nullptr;   /* Potential Transpiration of Vegetation */
    double *qEs = nullptr;    /* Evaporation from surface ponding */
    double *qEu = nullptr;    /* Evaporation from Unsat */
    double *qEg = nullptr;    /* Evaporation from GW */
    double *qTu = nullptr;    /* Transpiration from Unsat */
    double *qTg = nullptr;    /* Transpiration from GW */


    double *qEleE_IC = nullptr;    /* Evaporation from canopy interception */
    double *qEleEvapo = nullptr;    /* Evaporation from canopy interception */
    double *qEleTrans = nullptr;    /* Evaporation from canopy interception */

    double *iBeta = nullptr;
    double *qEleETA = nullptr;
    double *yRivStg = nullptr;   // debug may not necessary
    /* Lake variables */
    double *yLakeStg = nullptr;
    double *y2LakeArea = nullptr;
    double *QLakeSurf = nullptr;
    double *QLakeSub = nullptr;
    double *QLakeRivIn = nullptr;
    double *QLakeRivOut = nullptr;
    double *qLakeEvap = nullptr;
    double *qLakePrcp = nullptr;
    /* Per-edge scratch slots that avoid shared writes to per-lake
     * totals. Element->Lake surface/sub fluxes write to these slots
     * (size NumEle*3, indexed i*3+j) instead of the racy `QLakeSurf[ilake] += Q`
     * pattern; rhs_deterministic_gather() then sums them into
     * QLakeSurf/QLakeSub in a fixed order. */
    double *QeleSurf_lake = nullptr;
    double *QeleSub_lake = nullptr;
    /* Per-element scratch for lake-cell evap/prcp (NumEle sized).
     * Lake-cell elements write the pre-divided per-element
     * contribution; a gather in rhs_flux sums it to per-lake
     * qLakeEvap / qLakePrcp. That gather cannot be part of
     * rhs_deterministic_gather() because the lake clamp reads
     * qLakeEvap/qLakePrcp BEFORE rhs_deterministic_gather() is
     * called. */
    double *qEleEvapo_lake = nullptr;
    double *qElePrep_lake = nullptr;


    int NumSegmt;
    RiverSegement *RivSeg = nullptr;

    long ForcStartTime;

private:
    double *t_prcp = nullptr;
    double *t_temp = nullptr;
    double *t_rh = nullptr;
    double *t_sph = nullptr;
    double *t_wind = nullptr;
    double *t_rn = nullptr;
    double *t_vp = nullptr;
    double *t_lai = nullptr;
    double *t_mf = nullptr;
//    double *t_hc;  /* New defination: Height of Crop. void in temporary*/
public:
    /* Methods: */
    Model_Data();
    Model_Data(FileIn *f_in, FileOut *f_out);
    ~Model_Data();
    /* Model input/output */
    void loadinput();
    void initialize();
    /* Populate ElementHotData SoA static fields from
     * _Element AoS. Called from initialize() AFTER all element AoS
     * fields are loaded and BEFORE any RHS dispatch. Dynamic fields
     * (u_*) are seeded here too but get resynced by sync_hot_dynamic(i)
     * after each writer-method call during the RHS hot path. */
    void initialize_hot();
    /* Re-sync the four dynamic SoA fields (u_qi, u_qex,
     * u_effKH, u_satn) for element i from _Element AoS. Called after
     * Ele[i].updateElement(...) / updateLakeElement() / Flux_Infiltration() /
     * Flux_Recharge(). Inline so default-build emits no call overhead. */
    inline void sync_hot_dynamic(int i) {
        hot.u_qi[i]    = Ele[i].u_qi;
        hot.u_qex[i]   = Ele[i].u_qex;
        hot.u_effKH[i] = Ele[i].u_effKH;
        hot.u_satn[i]  = Ele[i].u_satn;
    }
    /* Inline accessors for the flattened
     * QeleSurf_flat / QeleSub_flat arrays. Index convention is
     * row-major: at(i,j) ↔ `_flat[3*i + j]`, matching the
     * MD_layout.hpp flat-3 idiom. Do not use bare
     * `*_flat[3*i + j]` indexing in the 4 hot-path TUs
     * (MD_ElementFlux.cpp / MD_f.cpp / MD_f_uncouple.cpp / MD_update.cpp)
     * — every read/write
     * MUST go through these accessors. Rationale:
     * (a) one source for the index expression — index-flip bugs
     * (3*j+i vs 3*i+j) are caught by checking the accessor once,
     * not 20 call sites; (b) future SIMD / NUMA tuning lands in one
     * place; (c) DEBUG bounds-check insertion point. The accessors
     * return references so they are usable on both LHS and RHS of
     * assignment without overhead in release builds. */
    inline double &QeleSurfAt(int i, int j) {
        return QeleSurf_flat[3*i + j];
    }
    inline double &QeleSubAt(int i, int j) {
        return QeleSub_flat[3*i + j];
    }
    inline double  QeleSurfAt(int i, int j) const {
        return QeleSurf_flat[3*i + j];
    }
    inline double  QeleSubAt(int i, int j) const {
        return QeleSub_flat[3*i + j];
    }
    void initializeLake();
    void initialize_output();
    void SetIC2Y(N_Vector udata1, N_Vector udata2, N_Vector udata3, N_Vector udata4, N_Vector udata5);
    void SetIC2Y(N_Vector udata);
    void LoadIC();
    /* screen print */
    void modelSummary(int end);
    int PrintInit(const char *fn, double t);
    
    void summary(N_Vector u1, N_Vector u2, N_Vector u3, N_Vector u4, N_Vector u5);
    void summary(N_Vector u);
    /* Output-time cache refresh — re-run
     * rhs_update + rhs_flux from Y(udata) so PCtrl-aliased output
     * buffers (QrivDown + siblings) reflect Y(tout) state, not the
     * RHS side-effect cache left by CV_NORMAL internal step at
     * t_internal != tout. Called from shud.cpp MainLoop between
     * summary(udata) and CS.ExportResults(t). */
    void recompute_for_output(N_Vector udata, double t);
    int ScreenPrint(double t, unsigned long it);
    int ScreenPrintu(double t, unsigned long it);
    /* methods in f function */
    /* The coupled path uses rhs_update / rhs_flux / rhs_apply declared
     * below. The methods here (f_loopET, f_loop1..5, f_applyDY_*,
     * f_applyDYi, f_updatei) serve the uncoupled path — called from
     * f.cpp's f_surf/f_unsat/f_gw/f_river/f_lake receivers. */
    void f_loopET(double t);
    void f_loop1(double t);
    void f_loop2(double t);
    void f_loop3(double t);
    void f_loop4(double t);
    void f_loop5(double t);

    void f_applyDY_surf(double * DY, double t);
    void f_applyDY_unsat(double * DY, double t);
    void f_applyDY_gw(double * DY, double t);
    void f_applyDY_river(double * DY, double t);
    void f_applyDYi(double * DY, double t, int flag);
    void f_updatei(double * Y, double * DY, double t, int flag);

    /* Coupled-path RHS: rhs_core() runs rhs_update (state update) ->
     * rhs_flux (flux computation) -> rhs_apply (DY accumulation),
     * switch-dispatched on the `ExecPolicy policy` argument.
     * See src/Model/MD_rhs_core.{cpp,hpp}. */
    void rhs_update(double * Y, double * DY, double t);
    void rhs_flux(double t);
    void rhs_apply(double * DY, double t);
    void rhs_core(double * Y, double * DY, double t, ExecPolicy policy);
    /* Unified deterministic gather called from rhs_flux.
     * Iterates the adjacency lists of MD_adjacency.hpp in their fixed
     * order to perform the segment->river/element, downstream and lake
     * sums, so results do not depend on the thread count. Body in
     * MD_rhs_core.cpp. */
    void rhs_deterministic_gather();
    
//    void updateWF(double dt);
    void CheckInputData();
    void InitFloodAlert(const char *fn);
    void updateRiverStage(N_Vector uY);
    void debugData();
    void debugData(const char *fn);
    void f_etFlux(int i, double t);
    void ET(double t, double tnext);
    void updateforcing(double t);
    double getArea();
private:
    void fillpits(int i);
    void tReadForcing(double t, int i);
    void ElementTable(const char *fn);
    void RiverTable(const char *fn);
    
    void LakeTable(const char *fn);
    int  LakeUniqueID();
    void LakeInitialize();
    void lake_readBathy(const char *fn);
    void lake_readIC(const char *fn);
    void lake_read_sp(const char *fn);
    
    /* Memory management: allocation and recycle */
    void malloc_Y();
    void malloc_EleRiv();
    void FreeData();
    
    /* put calibration file into the parameters */
    void copyCalib();
    void calibSoil();
    void calibGeol();
    void calibLandc();
    
    /* Check the input data */
    void CheckInput_forc();
    void CheckInput_mesh();
    void CheckInput_att();
    void CheckInput_soil();
    void CheckInput_geol();
    void CheckInput_landcover();
    
    /* Read input data */
    void read_calib(const char *fn);
    void read_para(const char *fn);
    void read_riv(const char *fn);
    void read_rivseg(const char *fn);
    void read_mesh(const char *fn);
    void read_cfgout(const char *fn);
    void setIO_ele(int x);
    void setIO_riv(int x);
    void setIO_lake(int x);
    
    void read_att(const char *fn);
    void read_soil(const char *fn);
    void read_geol(const char *fn);
    void read_lc(const char *fn);
    void read_forc_csv(const char *fn);
//    void read_rl(const char *fn);
    void read_lai(const char *fn);
    void read_mf(const char *fn);
    
    void read_ssEle(const char *fn);
    void read_bcEle1(const char *fn);
    void read_bcEle2(const char *fn);
    void read_bcRiv1(const char *fn);
    void read_bcRiv2(const char *fn);
    void read_bcLake1(const char *fn);
    void read_bcLake2(const char *fn);
//    void CorrectRiver(double eps);
    void rmSinks();
    
    /* Physical processes */
    void Flux_RiverDown(double t, int i);
    void applyBCSS(double *DY, int i);
    
    /* Methods for element calculation */
    void fun_Ele_sub(int i, double t);
    void fun_Ele_surface(int i, double t);
    void fun_Ele_Infiltraion(int i, double t);
    void fun_Ele_Recharge(int i, double t);
    void fun_Seg_surface(int iEle, int iRiv, int i);
    void fun_Seg_sub(int iEle, int iRiv, int i);
    void fun_Ele_lakeVertical(int i, double t);
    void fun_Ele_lakeHorizon(int i, double t);
    
    /* Functions */
    void TimeSpent();
    double WeirFlow_jtoi(double zi, double yi, double zj, double yj,
                    double zbank, double cwr, double width, double threshold);
    double updateArea();
};
#endif                /* Model_Data_hpp */

