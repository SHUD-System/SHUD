#include "TabularData.hpp"

TabularData::TabularData(){
    
}

TabularData::~TabularData(){
    reset();
}
void TabularData::reset(){
    if(nrow > 0){
        for(int i = 0; i < nrow; i++){
            delete[] x[i];
        }
        delete[] x;
    }
    nrow = 0;
    ncol = 0;
    header[0] = '\0';
}
int TabularData::read(FILE * fp){
    int tmp;
    read(fp, &tmp);
    return nrow;
}

int TabularData::read(FILE * fp, int *NumCol){
    reset();
    char str[MAXLEN];
    char * ps;
    fgets(str, MAXLEN, fp); // Dimension of the table
    sscanf(str, "%d %d", &nrow, &ncol);
    *NumCol = ncol;
    if(ncol < 1){
        fprintf(stderr, "WARNING: Number of columns is %d\n", ncol);
    }
    x = new double*[nrow];
    for(int i = 0; i < nrow; i++){
        x[i] = new double[ncol];
    }
    fgets(header, MAXLEN, fp); // the header
    for (int i=0; i < nrow && fgets(str, MAXLEN, fp) ; i++){
        ps = str;
        for(int j = 0; j < ncol; j++){
            x[i][j] = strtold(ps, &ps);
        }
    }
    //    for (int i=0; i < nrow; i++){
    //        for(int j = 0; j < ncol; j++){
    //            printf("%e \t ", x[i][j]);
    //        }
    //        printf("\n");
    //    }
    return nrow;
}
int TabularData::read(const char *fn){
    FILE * fp =  fopen(fn, "r");
    CheckFile(fp, fn);
    read(fp);
    fclose(fp);
    return nrow;
}
