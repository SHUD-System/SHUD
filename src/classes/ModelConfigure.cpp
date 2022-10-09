//  ModelConfigure.cpp
//
//  Created by Lele Shu on 9/18/18.
//  Copyright © 2018 Lele Shu. All rights reserved.
//

#include "ModelConfigure.hpp"
void Soil_Layer::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "soil_index");
    fprintf(fp, "%s\t", "infKsatV");
    fprintf(fp, "%s\t", "ThetaS");
    fprintf(fp, "%s\t", "ThetaR");
    fprintf(fp, "%s\t", "Alpha");
    fprintf(fp, "%s\t", "Beta");
    fprintf(fp, "%s\t", "hAreaF");
    fprintf(fp, "%s\t", "macKsatV");
    fprintf(fp, "%s\t", "infD");
}
void Soil_Layer::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
    fprintf(fp, "%g\t", infKsatV);
    fprintf(fp, "%g\t", ThetaS);
    fprintf(fp, "%g\t", ThetaR);
    fprintf(fp, "%g\t", Alpha);
    fprintf(fp, "%g\t", Beta);
    fprintf(fp, "%g\t", hAreaF);
    fprintf(fp, "%g\t", macKsatV);
    fprintf(fp, "%g\t", infD);
}
void Geol_Layer::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "geol_index");
    fprintf(fp, "%s\t", "KsatH");
    fprintf(fp, "%s\t", "KsatV");
    fprintf(fp, "%s\t", "SpYield");
    fprintf(fp, "%s\t", "geo_ThetaS");
    fprintf(fp, "%s\t", "geo_ThetaR");
    fprintf(fp, "%s\t", "geo_vAreaF");
    fprintf(fp, "%s\t", "macKsatH");
    fprintf(fp, "%s\t", "macD");
}
void Geol_Layer::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
    fprintf(fp, "%g\t", KsatH);
    fprintf(fp, "%g\t", KsatV);
    fprintf(fp, "%g\t", Sy);
    fprintf(fp, "%g\t", geo_ThetaS);
    fprintf(fp, "%g\t", geo_ThetaR);
    fprintf(fp, "%g\t", geo_vAreaF);
    fprintf(fp, "%g\t", macKsatH);
    fprintf(fp, "%g\t", macD);
}
void Landcover::printHeader(FILE *fp){
    fprintf(fp, "%s\t", "LC_index");
//    fprintf(fp, "%s\t", "LAImax");
    fprintf(fp, "%s\t", "VegFrac");
    fprintf(fp, "%s\t", "Albedo");
//    fprintf(fp, "%s\t", "Rs_ref");
//    fprintf(fp, "%s\t", "Rmin");
    fprintf(fp, "%s\t", "Rough");
    fprintf(fp, "%s\t", "RzD");
    fprintf(fp, "%s\t", "SoilDgrd");
    fprintf(fp, "%s\t", "ImpAF");
}
void Landcover::printInfo(FILE *fp){
    fprintf(fp, "%d\t", index);
//    fprintf(fp, "%g\t", LAImax);
    fprintf(fp, "%g\t", VegFrac);
    fprintf(fp, "%g\t", Albedo);
//    fprintf(fp, "%g\t", Rs_ref);
//    fprintf(fp, "%g\t", Rmin);
    fprintf(fp, "%g\t", Rough);
    fprintf(fp, "%g\t", RzD);
    fprintf(fp, "%g\t", SoilDgrd);
    fprintf(fp, "%g\t", ImpAF);
}
void Soil_Layer::applyCalib(calib_soil *g){
    infKsatV *= g->infKsatV;
    Alpha *= g->Alpha;
    Beta *= g->Beta;
    if(Beta < 1.1){
        printf("Warning: Beta value (%.g) less than 1.1 may lead in NAN value in the model, so default beta=1.1 is used.\n", Beta);
        Beta = 1.1;
//        myexit(ERRDATAIN);
    }
    hAreaF *= g->hAreaF;
    macKsatV *= g->macKsatV;
    infD *= g->infD;
}
void Soil_Layer::checkValue(){
    checkRange(Alpha, .05, 20., index - 1, "Alpha");
    checkRange(Beta, 1., 10., index - 1, "Beta");
    checkRange(infKsatV, 0., 1.0e3, index - 1, "infKsatV");
    checkRange(infD, 0., 10., index - 1, "infD");
    checkRange(hAreaF, 0., 100., index - 1, "hAreaF");
    checkRange(ThetaS, 0.001, 1., index - 1, "ThetaS");
}
void Geol_Layer::applyCalib(calib_geol *g){
    KsatH *= g->KsatH;
    KsatV *= g->KsatV;
    macD *= g->macD;
    macKsatH *= g->macKsatH;
    geo_vAreaF *= g->vAreaF;
    
    /* note: this one is different from above. */
    Sy = g->ThetaS * geo_ThetaS - g->ThetaR * geo_ThetaR;
}
void Geol_Layer::checkValue(){
    checkRange(Sy, 1.0e-4, 1., index - 1, "Specific Yield");
    checkRange(macD, 0., 10., index - 1, "macD");
    checkRange(KsatH, 0., 1.0e3, index - 1, "KsatH");
    checkRange(KsatV, 0., 1.0e3, index - 1, "KsatV");
    checkRange(geo_ThetaS, 0., 1., index - 1, "geo_ThetaS");
    checkRange(geo_vAreaF, 0., 1., index - 1, "geo_vAreaF");
}
void Landcover::applyCalib(calib_landcover *g){
//    LAImax  *= 1.0;
    VegFrac *= g->VegFrac;
    Albedo *= g->Albedo;
//    Rs_ref *= 1.0;
//    Rmin *= 1.0;
    Rough *= g->Rough;
    RzD *= g->RzD;
    SoilDgrd *= g->SoilDgd;
    ImpAF *= g->ImpAF;
}
void Landcover::checkValue(){
//    checkRange(LAImax, 0.0, 40., index - 1, "LAImax");
    checkRange(VegFrac, 0.0, 1.0, index - 1, "VegFrac");
    checkRange(Albedo, 0.0, 1.0, index - 1, "Albedo");
    checkRange(Rough, 0.0, 1.0, index - 1, "Rough");
    checkRange(SoilDgrd, 0.0, 1.0, index - 1, "SoilDgrd");
    checkRange(ImpAF, 0.0, 1.0, index - 1, "ImpAF");
    checkRange(RzD, 0.0, 10., index - 1, "RzD");
}
void globalCal::push(const char *optstr, double val){
    if (strcasecmp ("GEOL_KSATH", optstr) == 0)
        cgeol.KsatH = val;
    else if (strcasecmp ("GEOL_KSATV", optstr) == 0)
        cgeol.KsatV = val;
    else if (strcasecmp ("GEOL_KMACSATH", optstr) == 0)
        cgeol.macKsatH = val;
    else if (strcasecmp ("GEOL_DMAC", optstr) == 0)
        cgeol.macD = val;
    else if (strcasecmp ("GEOL_ThetaS", optstr) == 0)
        cgeol.ThetaS = val;
    else if (strcasecmp ("GEOL_ThetaR", optstr) == 0)
        cgeol.ThetaR = val;
    else if (strcasecmp ("GEOL_MACVF", optstr) == 0)
        cgeol.vAreaF = val;
    
    /*  Unsat soil layers */
    else if (strcasecmp ("SOIL_KINF", optstr) == 0)
        csoil.infKsatV = val;
    else if (strcasecmp ("SOIL_KMACSATV", optstr) == 0)
        csoil.macKsatV = val;
    else if (strcasecmp ("SOIL_DINF", optstr) == 0)
        csoil.infD = val;
    else if (strcasecmp ("SOIL_ALPHA", optstr) == 0)
        csoil.Alpha = val;
    else if (strcasecmp ("SOIL_BETA", optstr) == 0)
        csoil.Beta = val;
    else if (strcasecmp ("SOIL_MACHF", optstr) == 0)
        csoil.hAreaF = val;
    
    /*  Landcover layers */
    else if (strcasecmp ("LC_VEGFRAC", optstr) == 0)
        clandc.VegFrac = val;
    else if (strcasecmp ("LC_ALBEDO", optstr) == 0)
        clandc.Albedo = val;
    else if (strcasecmp ("LC_ROUGH", optstr) == 0)
        clandc.Rough = val;
    else if (strcasecmp ("LC_ISMAX", optstr) == 0)
        clandc.cISmax = val;
    else if (strcasecmp ("LC_DROOT", optstr) == 0)
        clandc.RzD = val;
    else if (strcasecmp ("LC_SoilDgd", optstr) == 0)
        clandc.SoilDgd = val;
    else if (strcasecmp ("LC_ImpAF", optstr) == 0)
        clandc.ImpAF = val;
    
    /*  Aquifer Depth */
    else if (strcasecmp ("AQ_DEPTH+", optstr) == 0)
        cAqD = val;
    
    /*  Forcing */
    else if (strcasecmp ("TS_PRCP", optstr) == 0)
        cPrep = val;
    else if (strcasecmp ("TS_SFCTMP+", optstr) == 0)
        cTemp = val;
    else if (strcasecmp ("TS_LAI", optstr) == 0)
        cLAItsd = val;
    else if (strcasecmp ("TS_MF", optstr) == 0)
        cMF = val;
    /* ET */
    else if (strcasecmp ("ET_IC", optstr) == 0)
        cE_ic = val;
    else if (strcasecmp ("ET_TR", optstr) == 0)
        cE_trans = val;
    else if (strcasecmp ("ET_SOIL", optstr) == 0)
        cE_Evapo = val;
    else if (strcasecmp ("ET_ETP", optstr) == 0)
        cETP = val;
    
    /*  Rivers */
    else if (strcasecmp ("RIV_ROUGH", optstr) == 0)
        criv.rivRough = val;
    else if (strcasecmp ("RIV_KH", optstr) == 0)
        criv.rivKsatH = val;
    else if (strcasecmp ("RIV_CWR", optstr) == 0)
        criv.rivCwr = val;
    else if (strcasecmp ("RIV_DPTH+", optstr) == 0)
        criv.rivDepth = val;
    else if (strcasecmp ("RIV_WDTH+", optstr) == 0)
        criv.rivWidth = val;
    else if (strcasecmp ("RIV_BSLOPE+", optstr) == 0)
        criv.rivBankSlope = val;
    else if (strcasecmp ("RIV_SINU", optstr) == 0)
        criv.rivSINU = val;
    else if (strcasecmp ("RIV_BEDTHICK", optstr) == 0)
        criv.rivBedThick = val;
    
    /* Frozen Soil (Permafrost) */
    else if (strcasecmp ("Fzn_submax", optstr) == 0)
        cfrozen.FT_sub_max = val;
    else if (strcasecmp ("Fzn_submin", optstr) == 0)
        cfrozen.FT_sub_min = val;
    else if (strcasecmp ("Fzn_subday", optstr) == 0)
        cfrozen.FT_sub_Day = val;
    
    else if (strcasecmp ("Fzn_surfmax", optstr) == 0)
        cfrozen.FT_surf_max = val;
    else if (strcasecmp ("Fzn_surfmin", optstr) == 0)
        cfrozen.FT_surf_min = val;
    else if (strcasecmp ("Fzn_surfday", optstr) == 0)
        cfrozen.FT_surf_Day = val;
    
    /* Initial Conditions */
    else if (strcasecmp ("IC_GW+", optstr) == 0)
        c_ic_gw = val;
    else if (strcasecmp ("IC_RIV+", optstr) == 0)
        c_ic_riv = val;
    /* Unrecognized Parameter Flag */
    else{
        printf
        ("\n  Parameter: %s cannot be recognized. Please see User's Manual for more details!\n",
         optstr);
        myexit(ERRFileIO);
    }//end ifelse
}
double globalCal::getValue(const char *optstr){
    double ret = 0.0;
    if (strcasecmp ("GEOL_KSATH", optstr) == 0)
        ret = cgeol.KsatH;
    else if (strcasecmp ("GEOL_KSATV", optstr) == 0)
        ret = cgeol.KsatV;
    else if (strcasecmp ("GEOL_KMACSATH", optstr) == 0)
        ret = cgeol.macKsatH;
    else if (strcasecmp ("GEOL_DMAC", optstr) == 0)
        ret = cgeol.macD;
    else if (strcasecmp ("GEOL_THETAS", optstr) == 0)
        ret = cgeol.ThetaS;
    else if (strcasecmp ("GEOL_THETAR", optstr) == 0)
        ret = cgeol.ThetaR;
    else if (strcasecmp ("GEOL_MACVF", optstr) == 0)
        ret = cgeol.vAreaF;
    
    /*  Unsat soil layers */
    else if (strcasecmp ("SOIL_KINF", optstr) == 0)
        ret = csoil.infKsatV;
    else if (strcasecmp ("SOIL_KMACSATV", optstr) == 0)
        ret = csoil.macKsatV;
    else if (strcasecmp ("SOIL_DINF", optstr) == 0)
        ret = csoil.infD;
    else if (strcasecmp ("SOIL_ALPHA", optstr) == 0)
        ret = csoil.Alpha;
    else if (strcasecmp ("SOIL_BETA", optstr) == 0)
        ret = csoil.Beta;
    else if (strcasecmp ("SOIL_MACHF", optstr) == 0)
        ret = csoil.hAreaF;
    
    /*  Landcover layers */
    else if (strcasecmp ("LC_VEGFRAC", optstr) == 0)
        ret = clandc.VegFrac;
    else if (strcasecmp ("LC_ALBEDO", optstr) == 0)
        ret = clandc.Albedo;
    else if (strcasecmp ("LC_ROUGH", optstr) == 0)
        ret = clandc.Rough;
    else if (strcasecmp ("LC_ISMAX", optstr) == 0)
        ret = clandc.cISmax;
    else if (strcasecmp ("LC_DROOT", optstr) == 0)
        ret = clandc.RzD;
    else if (strcasecmp ("LC_SoilDgd", optstr) == 0)
        ret = clandc.SoilDgd;
    else if (strcasecmp ("LC_ImpAF", optstr) == 0)
        ret = clandc.ImpAF;
    
    /*  Aquifer Depth */
    else if (strcasecmp ("AQ_DEPTH+", optstr) == 0)
        ret = cAqD;
    
    /*  Forcing */
    else if (strcasecmp ("TS_PRCP", optstr) == 0)
        ret = cPrep;
    else if (strcasecmp ("TS_SFCTMP+", optstr) == 0)
        ret = cTemp;
    else if (strcasecmp ("TS_LAI", optstr) == 0)
        ret = cLAItsd;
    else if (strcasecmp ("TS_MF", optstr) == 0)
        ret = cMF;
    /* ET */
    else if (strcasecmp ("ET_ETP", optstr) == 0)
        ret = cETP;
    else if (strcasecmp ("ET_IC", optstr) == 0)
        ret = cE_ic;
    else if (strcasecmp ("ET_TR", optstr) == 0)
        ret = cE_trans;
    else if (strcasecmp ("ET_SOIL", optstr) == 0)
        ret = cE_Evapo;
    
    /*  Rivers */
    else if (strcasecmp ("RIV_ROUGH", optstr) == 0)
        ret = criv.rivRough;
    else if (strcasecmp ("RIV_KH", optstr) == 0)
        ret = criv.rivKsatH;
    else if (strcasecmp ("RIV_CWR", optstr) == 0)
        ret = criv.rivCwr;
    else if (strcasecmp ("RIV_DPTH+", optstr) == 0)
        ret = criv.rivDepth;
    else if (strcasecmp ("RIV_BSLOPE", optstr) == 0)
        ret = criv.rivBankSlope;
    else if (strcasecmp ("RIV_WDTH+", optstr) == 0)
        ret = criv.rivWidth;
    else if (strcasecmp ("RIV_SINU", optstr) == 0)
        ret = criv.rivSINU;
    else if (strcasecmp ("RIV_BEDTHICK", optstr) == 0)
        ret = criv.rivBedThick;
    
    /* Frozen Soil (Permafrost) */
    else if (strcasecmp ("Fzn_submax", optstr) == 0)
        ret = cfrozen.FT_sub_max;
    else if (strcasecmp ("Fzn_submin", optstr) == 0)
        ret = cfrozen.FT_sub_min;
    else if (strcasecmp ("Fzn_subday", optstr) == 0)
        ret = cfrozen.FT_sub_Day;
    
    else if (strcasecmp ("Fzn_surfmax", optstr) == 0)
        ret = cfrozen.FT_surf_max;
    else if (strcasecmp ("Fzn_surfmin", optstr) == 0)
        ret = cfrozen.FT_surf_min;
    else if (strcasecmp ("Fzn_surfday", optstr) == 0)
        ret = cfrozen.FT_surf_Day;
    
    /* Initial condition */
    else if (strcasecmp ("IC_GW+", optstr) == 0)
        ret = c_ic_gw;
    else if (strcasecmp ("IC_RIV+", optstr) == 0)
        ret = c_ic_riv;

    /* Unrecognized Parameter Flag */
    else{
        printf
        ("\n  Parameter: %s cannot be recognized. Please see User's Manual for more details!\n",
         optstr);
        myexit(ERRFileIO);
    }//end ifelse
    return ret;
}
void globalCal::copy(const char **varname, int nv, double *x, int nx){
    if( nx < nv ){
        myexit(ERRCONSIS);
    }
    for(int i = 0; i < nv; i++){
        push(varname[i], x[i]);
    }
}
void globalCal::write(const char *fn){
    FILE *fp;
    fp = fopen (fn, "w");
    CheckFile(fp, fn);
    fprintf(fp, "#%s\n", "Geology Layers");
    fprintf(fp, "GEOL_KSATH\t%g\n", cgeol.KsatH);
    fprintf(fp, "GEOL_KSATV\t%g\n", cgeol.KsatV);
    fprintf(fp, "GEOL_KMACSATH\t%g\n", cgeol.macKsatH);
    fprintf(fp, "GEOL_DMAC\t%g\n", cgeol.macD);
    fprintf(fp, "GEOL_THETAS\t%g\n", cgeol.ThetaS);
    fprintf(fp, "GEOL_THETAR\t%g\n", cgeol.ThetaR);
    fprintf(fp, "GEOL_MACVF\t%g\n", cgeol.vAreaF);
    
    fprintf(fp, "#%s\n", "Soil Layers");
    fprintf(fp, "SOIL_KINF\t%g\n", csoil.infKsatV);
    fprintf(fp, "SOIL_KMACSATV\t%g\n", csoil.macKsatV);
    fprintf(fp, "SOIL_DINF\t%g\n", csoil.infD);
    fprintf(fp, "SOIL_ALPHA\t%g\n", csoil.Alpha);
    fprintf(fp, "SOIL_BETA\t%g\n", csoil.Beta);
    fprintf(fp, "SOIL_MACHF\t%g\n", csoil.hAreaF);
    fprintf(fp, "SOIL_DINF\t%g\n", csoil.infD);
    
    fprintf(fp, "#%s\n", "Land use");
    fprintf(fp, "LC_VEGFRAC\t%g\n", clandc.VegFrac);
    fprintf(fp, "LC_ALBEDO\t%g\n", clandc.Albedo);
    fprintf(fp, "LC_ROUGH\t%g\n", clandc.Rough);
    fprintf(fp, "LC_ISMAX\t%g\n", clandc.cISmax);
    fprintf(fp, "LC_DROOT\t%g\n", clandc.RzD);
    fprintf(fp, "LC_SoilDgd\t%g\n", clandc.SoilDgd);
    fprintf(fp, "LC_ImpAF\t%g\n", clandc.ImpAF);
    
    fprintf(fp, "#%s\n", "River stucture");
    fprintf(fp, "RIV_ROUGH\t%g\n", criv.rivRough);
    fprintf(fp, "RIV_KH\t%g\n", criv.rivKsatH);
    fprintf(fp, "RIV_DPTH+\t%g\n", criv.rivDepth);
    fprintf(fp, "RIV_WDTH+\t%g\n", criv.rivWidth);
    fprintf(fp, "RIV_SINU\t%g\n", criv.rivSINU);
    fprintf(fp, "RIV_CWR\t%g\n", criv.rivCwr);
    fprintf(fp, "RIV_BEDTHICK\t%g\n", criv.rivBedThick);
    
    fprintf(fp, "#%s\n", "Misc");
    fprintf(fp, "AQ_DEPTH+\t%g\n", cAqD);
    fprintf(fp, "TS_PRCP\t%g\n", cPrep);
    fprintf(fp, "TS_SFCTMP+\t%g\n", cTemp);
    fprintf(fp, "TS_LAI\t%g\n", cLAItsd);
    fprintf(fp, "TS_MF\t%g\n", cMF);
    fprintf(fp, "ET_ETP\t%g\n", cETP);
    fprintf(fp, "ET_IC\t%g\n", cE_ic);
    fprintf(fp, "ET_TR\t%g\n", cE_trans);
    fprintf(fp, "ET_SOIL\t%g\n", cE_Evapo);
    
    fprintf(fp, "#%s\n", "Frozen");
    fprintf(fp, "Fzn_submax\t%g\n", cfrozen.FT_sub_max);
    fprintf(fp, "Fzn_submin\t%g\n", cfrozen.FT_sub_min);
    fprintf(fp, "Fzn_subday\t%g\n", cfrozen.FT_sub_Day);
    
    fprintf(fp, "Fzn_surfmax\t%g\n", cfrozen.FT_surf_max);
    fprintf(fp, "Fzn_surfmin\t%g\n", cfrozen.FT_surf_min);
    fprintf(fp, "Fzn_surfday\t%g\n", cfrozen.FT_surf_Day);
    
    
    fprintf(fp, "#%s\n", "Intial Condition");
    fprintf(fp, "IC_GW+\t%g\n", c_ic_gw);
    fprintf(fp, "IC_RIV+\t%g\n", c_ic_riv);
    
    fclose(fp);
}
void globalCal::read(const char *fn){
    /*========= open *.calib file ==========*/
    FILE *fp;
    char str[MAXLEN], optstr[MAXLEN];
    fp = fopen (fn, "r");
    CheckFile(fp, fn);
    /* start reading calib_file */
    double val;
    while (fgets(str, MAXLEN, fp)) {
        if (str[0] == '#' || str[0] == '\n' || str[0] == '\0' || str[0] == ' '){
            continue;
        }
        sscanf (str, "%s %lf", optstr, &val);
        push(optstr, val);
    } // end of while gets()
    fclose (fp);
}
globalCal::globalCal(){
}

void globalCal::copy(globalCal *p){
    csoil.infKsatV = p->csoil.infKsatV;
    csoil.macKsatV = p->csoil.macKsatV;
    csoil.infD = p->csoil.infD;
    csoil.Alpha = p->csoil.Alpha;
    csoil.Beta = p->csoil.Beta;
    csoil.hAreaF = p->csoil.hAreaF;
    
    criv.rivRough = p->criv.rivRough;
    criv.rivBankSlope = p->criv.rivBankSlope;
    criv.rivCwr = p->criv.rivCwr;
    criv.rivKsatH = p->criv.rivKsatH;
    criv.rivDepth = p->criv.rivDepth;
    criv.rivWidth = p->criv.rivWidth;
    criv.rivSINU = p->criv.rivSINU;
    criv.rivBedThick = p->criv.rivBedThick;
    
    cgeol.KsatH = p->cgeol.KsatH;
    cgeol.KsatV = p->cgeol.KsatV;
    cgeol.macKsatH = p->cgeol.macKsatH;
    cgeol.macD = p->cgeol.macD;
    cgeol.ThetaS = p->cgeol.ThetaS;
    cgeol.vAreaF = p->cgeol.vAreaF;
    
    clandc.VegFrac = p->clandc.VegFrac;
    clandc.Albedo = p->clandc.Albedo;
    clandc.Rough = p->clandc.Rough;
    clandc.SoilDgd = p->clandc.SoilDgd;
    clandc.RzD = p->clandc.RzD;
    clandc.ImpAF = p->clandc.ImpAF;
    clandc.cISmax = p->clandc.cISmax;
    
    cAqD = p->cAqD ; // +
    cTemp = p->cTemp ; // +
    
    cPrep = p->cPrep ;
    cETP = p->cETP ;
    cE_ic = p->cE_ic ;
    cE_trans = p->cE_trans ;
    cE_Evapo = p->cE_Evapo ;
    cISmax = p->cISmax ;
    cLAItsd = p->cLAItsd ;
    cMF = p->cMF ;
}
