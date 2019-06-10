#include <stdlib.h>
#include <stdio.h>
#include <string.h>
//#include <libxml/tree.h>
//#include <libxml/parser.h>
#include <assert.h>
#include <ctype.h>
#include "parameters.h"


static FILE * parmsin;
static FILE * parmsout;
static FILE * matsout;

static System_Parameters currParams = {0};
static char linebuffer[5000];
float accumulated[TRANSMATDATASIZE];

#define MIRSTRENGTHTODEFECT(x) 1.0 - x
#define DEFECTTOMIRSTRENGTH(x) 1.0 - x
#define PMRTOPPROB(x) x/(1.0+x)
#define PPROBTOPMR(x) x/(1.0-x)


int isendoffile(){
    return feof(parmsin);
}

void writeSystemParameters(const System_Parameters p, FILE * out){
    fprintf(out,"%d\t%d\t%d\t%d\t%d\t%d\t",p.nMess,p.nMicro,p.c_Parms.d_lo,p.c_Parms.d_hi,p.l_Parms.d_lo,p.l_Parms.d_hi);
    fprintf(out,"%f\t%f\t%f\t%f\t%f\t%f\t%f\n",p.c_Parms.pnz,MIRSTRENGTHTODEFECT(p.l_Parms.miRStrength),PPROBTOPMR(p.c_Parms.pprob),
        p.l_Parms.ptarg,p.s_Parms.rho,p.s_Parms.Cnoise0,p.s_Parms.Cnoise1);
}

//A rather unsafe way to check for decimal numbers - be careful, it doesn't even ensure that
//there's only one decimal pt!!!!
int isADecimalNumberNonScientific(const char * str){
    const char * p =str;
    int until = strlen(str);
    int isaDecimal = 1;
    for(;until > 0; until--,p++){
        if ( (*p != '.') && (!isdigit(*p)) ){
            isaDecimal = 0;
            break;
        }
    }
    return isaDecimal;
}

#define CHECKNUMBER(str,linebuf) if(!isADecimalNumberNonScientific(str)){       \
        fprintf(stderr,"Line is not in allowed format:\n%s",linebuf);   \
        exit(-200);                                                     \
    }                                                                   \


//unsafe but practical
System_Parameters textReadSystemParameters(){
    System_Parameters toRet = {0};
    char linebuf[5000];
    char linebufcopy[5000];
    char * str;
    char *p = NULL;
    p = fgets(linebuf,5000,parmsin);
    if (NULL == p) return toRet;
    strcpy(linebufcopy,linebuf);

    str = strtok(linebufcopy,"\t");
    CHECKNUMBER(str,linebuf); toRet.nMess = atoi(str); str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.nMicro = atoi(str); str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.c_Parms.d_lo = atoi(str); str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.c_Parms.d_hi = atoi(str); str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.l_Parms.d_lo = atoi(str);  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.l_Parms.d_hi = atoi(str);  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.c_Parms.pnz = strtof(str,NULL);  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.l_Parms.miRStrength = DEFECTTOMIRSTRENGTH(strtof(str,NULL));  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.c_Parms.pprob = PMRTOPPROB(strtof(str,NULL));  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.l_Parms.ptarg = strtof(str,NULL);  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.s_Parms.rho = strtof(str,NULL);  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.s_Parms.Cnoise0 = strtof(str,NULL);  str = strtok(NULL,"\t");
    CHECKNUMBER(str,linebuf); toRet.s_Parms.Cnoise1 = strtof(str,NULL);  str = strtok(NULL,"\t");
    return toRet;
}



void openInputOutputFiles(const char * fnameIn,const char * fnameOut){
    parmsin = fopen(fnameIn,"r");
    //eliminates header
    if (NULL == fgets(linebuffer,5000,parmsin)){
        fprintf(stderr, "invalid file\n");
        parmsin = NULL;
        return;
    }
    char tmp[1000];
    //auto null-terminates
    snprintf(tmp,1000,"%s.params",fnameOut);
    parmsout = fopen(tmp,"w");
    if(!parmsout){
        fprintf(stderr,"Failed to open params output file: %s \n",tmp);
        exit(-131);
    }

    snprintf(tmp,1000,"%s.matrices",fnameOut);
    matsout = fopen(tmp,"w");
    if(!matsout){
        fprintf(stderr,"Failed to open matrices output file: %s \n",tmp);
        exit(-131);
    }
    //printf("%s\n",linebuffer);
}

System_Parameters * nextParameters(){
    currParams = textReadSystemParameters();
    if(currParams.nMess == 0){
        return NULL;
    }else if(currParams.nMess != NMESS){
        fprintf(stderr,"system unable to accept any value for nMess other than %d, current value is: %d\n",NMESS,currParams.nMess);
        exit(-800);
    }
    //writeSystemParameters(currParams,stdout);
    return &currParams;
}

void closeInputOutputFiles(){
    fclose(parmsin);
    fclose(parmsout);
    fclose(matsout);
}




void accumulateTransmats(const int * h_tmats, int nMats,int numVecsPerOne){

    int i,j;
    for(i=0;i<TRANSMATDATASIZE;i++)
        accumulated[i]=0.0;
    const int *currItem = h_tmats;
    for(i=0;i<nMats;i++){
        for(j=0;j<TRANSMATDATASIZE;j++){
            accumulated[j] = accumulated[j] + (float) *currItem++;
        }
    }
    float rowSum = (float) (nMats * numVecsPerOne);
    //divide by sum
    for(i=0;i<TRANSMATDATASIZE;i++){
        accumulated[i] /= rowSum;
    }

}

void writeData(const int * h_tmats,int nMats){
    writeSystemParameters(currParams,parmsout);
    accumulateTransmats(h_tmats,nMats,NVPO);
    int i;
    float * currElem = accumulated;
    for (i=1;i<TRANSMATDATASIZE;i++,currElem++){
        fprintf(matsout,"%f\t",*currElem);
    }
    fprintf(matsout,"%f\n",*currElem);
}
