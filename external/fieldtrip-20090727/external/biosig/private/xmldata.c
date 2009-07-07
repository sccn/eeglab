/*
//   LICENSE:
//  ========= 
//
//   Copyright (C) Peter Rydesäter 2002, Mitthögskolan, SWEDEN
//
//   This program is free software; you can redistribute it and/or
//   modify it under the terms of the GNU General Public License
//   as published by the Free Software Foundation; either version 2
//   of the License, or (at your option) any later version.
//
//   This program is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with this program; if not, write to the Free Software
//   Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
//   
//
//   Please, send me an e-mail if you use this and/or improve it. 
//   If you think this is usful fore your work, please remember me and my
//   hard work with it!
//
//   Peter.Rydesater@mh.se
//
//   Peter Rydesäter
//   Mitthögskolan (Mid Sweden University)
//   SE-831 25 ÖSTERSUND
//   SWEDEN
//
//
*/


#include "mex.h"
#include <string.h>

mxArray **globalplhs=NULL;

unsigned short int *matstr={0};
int matlen=0;
int matpos=0;
#define STREND (-1)

mxArray *myCreateScalarDouble(double val){
    mxArray *ptr=mxCreateDoubleMatrix(1,1,0);
    mxGetPr(ptr)[0]=val;
    return ptr;
}


int iseof()
{
    return (matpos>=matlen || matpos<0);
}

int nextch()
{
    if(iseof())
	return STREND;
    else
	return matstr[matpos];
}

int readch()
{
    if(iseof())
	return STREND;
    else
	return matstr[matpos++];
}

void unreadch(int ch)
{
    if(matpos>0)
	matpos--;
}

void skipblank()
{
    while(!iseof()){
	int ch=readch();
	if(!isspace(ch)){
	    unreadch(ch);
	    break;
	}
    }
}

void readlabel(char *str){
    int n=0;
    while(!iseof()){
	int ch=readch();
	if(ch=='=' || isspace(ch) || ch=='>' || ch=='/' || iseof()){
	    unreadch(ch);
	    break;
	}
	if(!isalnum(ch)){
	    ch='_';    
	    if(n==0) str[n++]='x';
	}	    
	str[n++]=ch;
    }
    str[n]=0;
}

void readcom(char *str,int n){
    while(!iseof()){
	int ch=readch();
	str[n++]=ch;
	if(ch=='>')
	    break;	
    }
    str[n]=0;
}

void readtext(char *str){
    int n=0;
    skipblank();                 /* remove blank; */
    while(!iseof()){
	if(nextch()=='<')
	    break;
	str[n++]=readch();
    }
    while(n>0 && isspace(str[n-1])) /* remove blank; */
	n--;
    str[n]=0;
}

void readparvalue(char *str){
    int n=0;
    skipblank();
    if(nextch()=='"'){
	readch();
	while(!iseof()){
	    int ch=readch();
	    if(ch=='"')
		break;
	    str[n++]=ch;
	}
    }else{
	while(!iseof()){
	    int ch=nextch();
	    if(isspace(ch) || ch=='/' || ch=='>'){
		break;
	    }
	    str[n++]=readch();
	}
    }
    str[n]=0;
}
/*////////////////////////

////////////////////////*/

void mktype(int type){
    globalplhs[2]=myCreateScalarDouble((double)type);
}

void mktagend(char *str){
    int dims[]={0,0};
    globalplhs[0]=mxCreateString(str);
    globalplhs[1]=mxCreateString("");
}

void mktagstart(char *str){
    globalplhs[0]=mxCreateString(str);
    globalplhs[1]=mxCreateStructMatrix(1,1,0,NULL);
}

void mktagpar(char *str){
    mxAddField(globalplhs[1],str);
}

void mktagparvalue(char *str){
    mxSetFieldByNumber(globalplhs[1],0,mxGetNumberOfFields(globalplhs[1])-1,mxCreateString(str));
}

void mktext(char *str){
    if(str[0]==0)
	return;
    globalplhs[0]=mxCreateString("");
    globalplhs[1]=mxCreateString(str);
    globalplhs[2]=myCreateScalarDouble(0);
}

/*////////////////////////*/

void readmk_text(){
    char buff[100000];
    readtext(buff);
    mktext(buff);
    return;
    }

void readmk_tag(){
    char buff[100000];
    int tagtype=1;
    while(readch()!='<')
	if(iseof()) return;
    skipblank();
    if(nextch()=='/'){
	tagtype=3;
	readch();
	skipblank();
    }
    if(nextch()=='!'){
	buff[0]='<';
	readcom(buff,1);
	mktext(buff);
	return;
    }
    readlabel(buff);
    if(tagtype==3)
	mktagend(buff);	
    else
	mktagstart(buff);
    while(!iseof()){
	skipblank();
	if(nextch()=='/'){
	    readch();
	    tagtype=2;
	}
	skipblank();
	if(nextch()=='>'){
	    readch();
	    mktype(tagtype);
	    break;
	}
	readlabel(buff);
	mktagpar(buff);
	skipblank();
	if(nextch()=='='){
	    readch();
	    readparvalue(buff);
	    mktagparvalue(buff);
	}
    }
}

void mexFunction(int nlhs,mxArray *plhs[],int nrhs, const mxArray *prhs[]){
    globalplhs=plhs;
    matpos=0;
    if(nrhs<2)
	return;
    if(mxGetClassID(prhs[0])!=mxCHAR_CLASS)
	return;
    matstr=(unsigned short int *)mxGetPr(prhs[0]);
    matlen=mxGetNumberOfElements(prhs[0]);
    matpos=mxGetScalar(prhs[1])-1;
    if(nlhs<3)
	return;
    if(!iseof()){
	readmk_text();
	if(globalplhs[0]==NULL)
	    readmk_tag();
    }
    if(nlhs>0 && globalplhs[0]==NULL)
	globalplhs[0]=myCreateScalarDouble(0);
    if(nlhs>1 && globalplhs[1]==NULL)
	globalplhs[1]=myCreateScalarDouble(0);
    if(nlhs>2 && globalplhs[2]==NULL)
	globalplhs[2]=myCreateScalarDouble(0);
    if(nlhs>3 && globalplhs[3]==NULL)
	globalplhs[3]=myCreateScalarDouble((matpos+1)*(!iseof()));
}
