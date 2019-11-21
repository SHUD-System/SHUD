#include "Model_Data.hpp"
void Model_Data::CheckInputData(){
    CheckInput_att();
    CheckInput_mesh();
    CheckInput_geol();
    CheckInput_soil();
    CheckInput_landcover();
    CheckInput_forc();
}
void Model_Data::CheckInput_att(){
    int flag = 0;
    for(int i = 0; i< NumEle; i++){
        flag = flag | checkRange(Ele[i].iSoil, 1, NumSoil, i, "Soil");
        flag = flag | checkRange(Ele[i].iGeol, 1, NumGeol, i,  "Geol");
        flag = flag | checkRange(Ele[i].iLC, 1, NumLC,  i, "Landcover");
        flag = flag | checkRange(Ele[i].iForc, 1, NumForc,  i, "Forcing");
        flag = flag | checkRange(Ele[i].iMF, 1, NumMeltF,  i, "MeltFactor");
    }
    if(flag){
        myexit(ERRDATAIN);
    }
}
void Model_Data::CheckInput_mesh(){

}
void Model_Data::CheckInput_soil(){
    for(int i = 0; i < NumSoil; i++){
        Soil[i].checkValue();
    }
}
void Model_Data::CheckInput_geol(){
    for(int i = 0; i < NumGeol; i++){
        Geol[i].checkValue();
    }
}
void Model_Data::CheckInput_landcover(){
    for(int i = 0; i < NumLC; i++){
        LandC[i].checkValue();
    }
}
void Model_Data::CheckInput_forc(){
    /********************************************
    This function cannot find the abnormal values,
     but find possible error from WRONG UNIT,
     e.g. precipitation in mm/day, mm/hr;
     Temperature in Kelvin instead of Celcius
     *******************************************/
    for(int i = 0; i < NumForc; i++){
        tsd_weather[i].checkValue(i_prcp, 0., 0.4 * 24, "Prcp");
        /*******************************************
         Maximum hourly precipitation. 0.4 m/hr.
         ref: http://www.nws.noaa.gov/oh/hdsc/record_precip/record_precip_world.html
         ********************************************/
        tsd_weather[i].checkValue(i_temp, -70., 50., "Temp");
        tsd_weather[i].checkValue(i_rh, 0., 1.5, "RH"); /* Theoretical MAX rh is 1 (100%). */
        tsd_weather[i].checkValue(i_wind, 0., 50. * 86400., "Wind"); /*50m/s is 13-level Typhoon. */
        tsd_weather[i].checkValue(i_rn, 0., 1360. * 86400, "Radiation");
        /* max = solar constant. 1360 */
    }
}

