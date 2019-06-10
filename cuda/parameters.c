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

////////////////////////////////////////////////////////////////////////////////
// The subroutines below can parse an XML file with parameters (see "config.xml")
//to test a single set of parameters if desired
////////////////////////////////////////////////////////////////////////////////

// static int getXMLint(xmlDocPtr doc, xmlNodePtr cur, const char * name) {

// 	xmlChar *key;
// 	int toRet;
// 	int found = 0;
// 	while (cur != NULL) {
// 		if ((!xmlStrcmp(cur->name, (const xmlChar* ) name))) {
// 			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
			
// 			toRet = strtol((const char *)key, NULL, 10);
// 			//printf("toRet=%d\n",toRet);
// 			xmlFree(key);
// 			found = 1;
// 		}
// 		//printf("%s\n",cur->name);
// 		cur = cur->next;
// 	}
// 	if(!found){
// 		printf("not found  %s\n",name);
// 		exit(-50);
// 	}
// 	assert(found);
// 	return toRet;
// }


// static double getXMLnumeric(xmlDocPtr doc, xmlNodePtr cur, const char * name) {

// 	xmlChar *key;
// 	double toRet;
// 	int found = 0;
// 	while (cur != NULL) {
// 		if ((!xmlStrcmp(cur->name, (const xmlChar* ) name))) {
// 			key = xmlNodeListGetString(doc, cur->xmlChildrenNode, 1);
// 			//printf("%s\n",key);

// 			toRet = strtod((const char *)key, NULL);
// 			xmlFree(key);
// 			found = 1;
// 		}
// 		cur = cur->next;
// 	}
// 	assert(found);
// 	return toRet;
// }


// static L_Array_Init_Parms parseLParams(xmlDocPtr doc, xmlNodePtr cur) {
// 	// globalDims.nMess = globalDims.nMicro = 0;
// 	// int foundnMess, foundnMicro;
// 	// foundnMicro = foundnMess = 0;
// 	cur = cur->xmlChildrenNode;
// 	L_Array_Init_Parms toRet = {0};
// 	//make sure it matches the right name
// 	int foundLow,foundHi,foundptarg,foundstrength;
// 	foundLow=foundHi=foundptarg=foundstrength=0;
// 	while (cur != NULL) {
// 	if(!foundLow)
// 		{toRet.d_lo = getXMLint(doc, cur, "d_lo");foundLow=1;}
// 		if(!foundHi)
// 		{toRet.d_hi = getXMLint(doc, cur, "d_hi");foundHi=1;}
// 		if(!foundptarg)
// 		{toRet.ptarg = getXMLnumeric(doc,cur,"ptarg");foundptarg=1;}
// 		if(!foundstrength)
// 		{toRet.miRStrength = getXMLnumeric(doc,cur,"miRStrength");foundstrength=1;}
// 		cur = cur->next;
// 	}
// 		return toRet;
// }


// static C_Array_Init_Parms parseCParams(xmlDocPtr doc, xmlNodePtr cur) {
// 	// globalDims.nMess = globalDims.nMicro = 0;
// 	// int foundnMess, foundnMicro;
// 	// foundnMicro = foundnMess = 0;
// 	cur = cur->xmlChildrenNode;
// 	C_Array_Init_Parms toRet;
// 	//make sure it matches the right name
// 	int foundIn,foundLow,foundHi,foundPnz,foundPProb;
// 	foundIn=foundLow=foundHi=foundPnz=foundPProb=0;
// 	while (cur != NULL) {
// 		if(!foundIn)
// 		{toRet.byIndegree = getXMLint(doc, cur, "byIndegree");foundIn=1;}
// 		if(!foundLow)
// 		{toRet.d_lo = getXMLint(doc, cur, "d_lo");foundLow=1;}
// 		if(!foundHi)
// 		{toRet.d_hi = getXMLint(doc, cur, "d_hi");foundHi=1;}
// 		if(!foundPnz)
// 		{toRet.pnz = getXMLnumeric(doc,cur,"pnz");foundPnz=1;}
// 		if(!foundPProb)
// 			{toRet.pprob = getXMLnumeric(doc,cur,"pprob");foundPProb=1;}
// 		cur = cur->next;
// 	}
// 		return toRet;
// }


// static Simulation_Params parseSimParams(xmlDocPtr doc, xmlNodePtr cur) {
// 	cur = cur->xmlChildrenNode;
// 	Simulation_Params toRet;
// 	//make sure it matches the right name
// 	int foundrho,foundc0,foundc1,foundl0,foundl1;
// 	foundrho=foundc0=foundc1=foundl0=foundl1=0;
// 	while (cur != NULL) {
// 		if(!foundrho)
// 		{toRet.rho = getXMLnumeric(doc, cur, "rho");foundrho=1;}
// 		if(!foundc0)
// 		{toRet.Cnoise0 = getXMLnumeric(doc, cur, "cnoise0");foundc0=1;}
// 		if(!foundc1)
// 		{toRet.Cnoise1 = getXMLnumeric(doc, cur, "cnoise1");foundc1=1;}
// 		if(!foundl0)
// 		{toRet.Lnoise0 = getXMLnumeric(doc,cur,"lnoise0");foundl0=1;}
// 		if(!foundl1)
// 			{toRet.Lnoise1 = getXMLnumeric(doc,cur,"lnoise1");foundl1=1;}
// 		cur = cur->next;
// 	}
// 		return toRet;
// }



// void initializeParametersXML( const char * cwd,System_Parameters * sys_params){
// 	xmlDocPtr doc;
// 	xmlNodePtr cur;
// 	char str[400];
// 	sprintf(str, "%s/config.xml", cwd);
// 	doc = xmlParseFile(str);
//     sys_params->c_Parms.d_lo = sys_params->c_Parms.d_hi = sys_params->c_Parms.pnz = sys_params->c_Parms.pprob = sys_params->c_Parms.byIndegree = -1;

// 	if (doc == NULL) {
// 		fprintf(stderr, "config.xml parsing failed");
// 		exit(-85);
// 	}

// 	cur = xmlDocGetRootElement(doc);

// 	if (cur == NULL) {
// 		fprintf(stderr, "empty document\n");
// 		xmlFreeDoc(doc);
// 		exit(-84);
// 	}

// 	if (xmlStrcmp(cur->name, (const xmlChar *) "sysparams")) {
// 		fprintf(stderr, "Root node failure - should be \"sysparams\", read as \"%s\"\n", cur->name);
// 		xmlFreeDoc(doc);
// 		exit(-86);
// 	}
// 	printf("%s:%d - initializeParameters called\n",__FILE__,__LINE__);
// 		cur = cur->xmlChildrenNode;
// 		while (cur != NULL) {
// 			if ((!xmlStrcmp(cur->name, (const xmlChar *) "c_params")))
// 				sys_params->c_Parms = parseCParams(doc, cur);
// 			if ((!xmlStrcmp(cur->name, (const xmlChar *) "l_params")))
// 				sys_params->l_Parms = parseLParams(doc, cur);
// 			if ((!xmlStrcmp(cur->name, (const xmlChar *) "sim_params")))
// 				sys_params->s_Parms = parseSimParams(doc, cur);
			
// 			cur = cur->next;
			
// 		}
// 		//Check validity
// 		// check_parameters_validity();
// 		xmlFreeDoc(doc);

// }

