//=================================================================================
//  mexSLIC.c
//  Superpixel Segmentation with the DBSCAN algorithm
//  Author: 2016-10-24, jlfeng

//=================================================================================
/*Redistribution and use in source and binary forms, with or without
 modification, are permitted provided that the following conditions are met
 
 * Redistributions of source code must retain the above copyright
 notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright
 notice, this list of conditions and the following disclaimer in the
 documentation and/or other materials provided with the distribution.
 * Neither the name of EPFL nor the
 names of its contributors may be used to endorse or promote products
 derived from this software without specific prior written permission.
 
 THIS SOFTWARE IS PROVIDED BY THE REGENTS AND CONTRIBUTORS ``AS IS'' AND ANY
 EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
 WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
 DISCLAIMED. IN NO EVENT SHALL THE REGENTS AND CONTRIBUTORS BE LIABLE FOR ANY
 DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
 (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
 ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */
 
#include<mex.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <string.h>  

#define MX_MAX(a,b) ((a) > (b) ? (a) : (b)) 
#define MX_MIN(a,b) ((a) < (b) ? (a) : (b)) 

_inline double GetEuclidDist(double *pf1, double*pf2,  int nd)
{
    double dist=0;
    int ii;
    for (ii=0; ii<nd;ii++)
    {
        dist+=pow(pf1[ii]-pf2[ii],2);
    }
    return dist;
}

void NeighborSearch(double *ptrImgData, mwSize *dims, int searchRange,
    int pp, int ppn, double thrDist, int *listNbrPix, int *numNbrPix)
{
     int ppnn0, ppnn;
     int imgHeight=dims[0], imgWidth=dims[1], imgDepth=dims[2];
     int imgSize=imgHeight*imgWidth;
     double *ptrDataP0=ptrImgData+pp*imgDepth;   
     double *ptrDataPn=ptrImgData+ppn*imgDepth;  
     double *ptrDataPnn;     
     double distTemp;
     int dpx, dpy;
     
     for (dpy=-searchRange; dpy<=searchRange;dpy++)
     {
         ppnn0=ppn+dpy*imgWidth;        
         for (dpx=-searchRange; dpx<=searchRange;dpx++)
         {
             ppnn=ppnn0+dpx;
             if (ppnn>=0 && ppnn<imgSize)
             {
                 ptrDataPnn=ptrImgData+ppnn*imgDepth;  
                 distTemp=0.6*GetEuclidDist(ptrDataPnn,ptrDataP0,imgDepth)
                    +0.4*GetEuclidDist(ptrDataPnn,ptrDataPn,imgDepth);
                 if (distTemp<thrDist)
                 {
                     listNbrPix[*numNbrPix]=ppnn;
                     *numNbrPix=*numNbrPix+1;
                 }
             }
         }
     }
}

void ExecuteDBSCAN(double *ptrImgData, mwSize *dims, double thrDist,  
    int maxSpSize, int *spLabel, int *numSpRet)
{
    int imgHeight=dims[0], imgWidth=dims[1], imgDepth=dims[2];
    int imgSize=imgHeight*imgWidth;
    int pr, pc, pp, prn, pcn, ppn;
    int numSp, numSpPix;
    double dist;
    int *listNbrPix, idxNbrPix, numNbrPix;
    int *listCadPix, idxCadPix, numCadPix;
    int searchRange=1;
    int nbrSize=(searchRange+1)*(searchRange+1);
    int minNbrPix=3;
    int labelPn;
    int deltR[4]={-1,1,0,0};
    int deltC[4]={0,0,-1,1};
    int ii;
    double minDist;
     
    listCadPix=mxMalloc(sizeof(int)*imgSize);
    memset(listCadPix, 0, sizeof(int)*imgSize);
    listNbrPix=mxMalloc(sizeof(int)*nbrSize*2);
    memset(listNbrPix, 0, sizeof(int)*nbrSize*2);

    //Set the initial value of labels
    for (pp=0;pp<imgSize;pp++)
    {
        spLabel[pp]=-1;
    }    
    numSp=0;
    
    for (pr=0;pr<imgHeight;pr++)
    {
        for (pc=0;pc<imgWidth; pc++)
        {
            pp=pr*imgWidth+pc;
            if (-1==spLabel[pp])// find a new start point
            {                
                listCadPix[0]=pp;
                numCadPix=1;
                NeighborSearch(ptrImgData, dims, searchRange, pp, pp, thrDist, listCadPix, &numCadPix);
                if (numCadPix<minNbrPix)
                {
                    spLabel[pp]=0; //Mark as visited
                    continue;
                }
                else
                {
                    numSp+=1;//Start a new superpixel
                    spLabel[pp]=numSp;
                    numSpPix=1;
                    idxCadPix=1;
                    while (idxCadPix<numCadPix && numSpPix<maxSpSize)
                    {
                        ppn=listCadPix[idxCadPix];
                        labelPn=spLabel[ppn];                        
                        if (labelPn<1) //Have not been assigned to a superpixel
                        {
                            spLabel[ppn]=numSp;
                            numSpPix++;
                            if (labelPn<0)//Have Not Searched the Neighborhood
                            {
                                numNbrPix=0;
                                NeighborSearch(ptrImgData, dims,searchRange, pp, ppn , thrDist, listNbrPix, &numNbrPix);
                                if (numNbrPix>=minNbrPix)
                                {
                                    memcpy(listCadPix+numCadPix,listNbrPix,numNbrPix*sizeof(int));
                                    numCadPix+=numNbrPix;
                                }
                            }
                        }                              
                        idxCadPix++;
                    }// while (idxCadPix<numCadPix)                     
                }// if (numCadPix<minNbrPix)
            }// if (spLabel[pp]<0)
        }// for (pc=0;pc<imgWidth; pc++)
    } //for (pr=0;pr<imgHeight;pr++) 
     *numSpRet=numSp;
 
    //Clear noise pixels    
    for (pr=0;pr<imgHeight;pr++)
    {
        for (pc=0;pc<imgWidth; pc++)
        {
            pp=pr*imgWidth+pc;
            if (0==spLabel[pp])
            {
                minDist=DBL_MAX;
                for (ii=0;ii<4;ii++)
                {
                    prn=pr+deltR[ii];
                    pcn=pc+deltC[ii];
                    if (prn>-1 && prn<imgHeight && pcn>-1 && pcn<imgWidth)
                    {
                        ppn=prn*imgWidth+pcn;
                        if (spLabel[ppn]>0)
                        {
                            dist=GetEuclidDist(ptrImgData+pp*imgDepth,
                                ptrImgData+ppn*imgDepth,imgDepth);
                            if (dist<minDist)
                            {
                                 minDist=dist;
                                 spLabel[pp]=spLabel[ppn];
                            }
                        } 
                    }
                }
            }
        }
    }    
    mxFree(listCadPix);
    mxFree(listNbrPix);
}

void SuperpixelRelabeling(int *spLabel, mwSize *dims, int numSpIn, int *spLabelC, int *numSpOut)
{
    int height=dims[0], width=dims[1], imgSize=width*height;
    int thrSpSize=imgSize/(numSpIn*2);
    int idxPixel, idxPixelAnchor;
    int px,py, pxn, pyn, idxn;    
    int dx[4]={-1,0,1,0}; // 4-connection neighborhood
    int dy[4]={0,-1,0,1};
    int numSp=0;
    int adjLabel=0;
    int *vecIdxPixelInSp=mxMalloc(2*imgSize*sizeof(int));
    int numPixelInSp, idxPixelInSp;
    for (idxPixel=0;idxPixel<imgSize;idxPixel++)
    {
        spLabelC[idxPixel]=-1;
    }
    
    idxPixelAnchor=0;
    for (py=0;py<height;py++)
    {
        for (px=0;px<width;px++)
        {
            if (spLabelC[idxPixelAnchor]<0)// find a new superpixel
            {
                spLabelC[idxPixelAnchor]=numSp;
                vecIdxPixelInSp[0]=px;
                vecIdxPixelInSp[1]=py;
                
                for (idxn=0;idxn<4;idxn++)// search the neighboring superpixel
                {
                    pxn=px+dx[idxn];
                    pyn=py+dy[idxn];
                    if ((pxn>-1 && pxn<width) &&(pyn>-1 && pyn<height))
                    {
                        idxPixel=pyn*width+pxn;
                        if (spLabelC[idxPixel]>-1)
                        {
                            adjLabel=spLabelC[idxPixel];
                        }
                    }
                }
                // Search pixels of the same superpixel
                numPixelInSp=1;
                idxPixelInSp=0;
                while (idxPixelInSp<numPixelInSp)
                {
                    for (idxn=0;idxn<4;idxn++)
                    {
                        pxn=vecIdxPixelInSp[idxPixelInSp*2]+dx[idxn];
                        pyn=vecIdxPixelInSp[idxPixelInSp*2+1]+dy[idxn];
                        if ((pxn>-1 && pxn<width) &&(pyn>-1 && pyn<height))
                        {
                            idxPixel=pyn*width+pxn;
                            if (spLabelC[idxPixel]<0 && spLabel[idxPixel]==spLabel[idxPixelAnchor])
                            {
                                vecIdxPixelInSp[numPixelInSp*2]=pxn;
                                vecIdxPixelInSp[numPixelInSp*2+1]=pyn;
                                spLabelC[idxPixel]=numSp;
                                numPixelInSp++;
                            }
                        }
                    }
                    idxPixelInSp++;
                }
                
                if (numPixelInSp<thrSpSize)
                {
                    for (idxPixelInSp=0;idxPixelInSp<numPixelInSp;idxPixelInSp++)
                    {
                        idxPixel=vecIdxPixelInSp[idxPixelInSp*2+1]*width+vecIdxPixelInSp[idxPixelInSp*2];
                        spLabelC[idxPixel]=adjLabel;
                    }
                }
                else
                {
                    numSp++;
                }                
            }
             idxPixelAnchor++;
        }
    }
    
    *numSpOut=numSp;
    mxFree(vecIdxPixelInSp);
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    int imgWidth, imgHeight, imgDepth,imgSize;
    double *ptrInputData, *ptrImgData,*ptrPixelData;
    int *spLabel, *spLabelM;
    int numSpIn, numSpOut, numSpMid;
    int meanSpSize, minSpSize, maxSpSize;
    double thrDist;
    int px,py,pd, idxPixel;
    mwSize numdims, *dims, imgDims[3];
    int *spLabelOut, *ptrNumSpOut;
    int ii;
    
    //Check input parameter
    if (nrhs<1)
    {
        mexErrMsgTxt("No input img.");
    }
    else if (nrhs>3)
    {
        mexErrMsgTxt("Too many input arguments.");
    }
    
    // Get data from input augments    
    numdims=mxGetNumberOfDimensions(prhs[0]) ;
    dims=mxGetDimensions(prhs[0]);
    ptrInputData=(double *)mxGetData(prhs[0]);
    imgHeight=dims[0];
    imgWidth=dims[1];
    imgDepth=(numdims==2)?1:dims[2];
    imgSize=imgHeight*imgWidth;
    numSpIn=(int)mxGetScalar(prhs[1]);
    thrDist=mxGetScalar(prhs[2]);
    
     // Allocate memory for temporary data
    ptrImgData=mxMalloc(sizeof(double)*imgSize*imgDepth);
    spLabel=mxMalloc(sizeof(int)*imgSize);
    spLabelM=mxMalloc(sizeof(int)*imgSize); 
    memset(ptrImgData,0,sizeof(double)*imgSize*imgDepth);
    memset(spLabel,0,sizeof(int)*imgSize);
    memset(spLabelM,0,sizeof(int)*imgSize);    
    
    meanSpSize=imgSize/numSpIn;
    minSpSize=meanSpSize/2;
    maxSpSize=meanSpSize+minSpSize;
    if (meanSpSize<10)
    {
        mexErrMsgTxt("Mean size of superpixel is too small.");
    }
    
    //img data copy
    ii=0;
    for (px=0;px<imgWidth;px++)        
     {
         for (py=0;py<imgHeight;py++)    
         {
             idxPixel=py*imgWidth+px;
             ptrPixelData=ptrImgData+idxPixel*imgDepth;
             for (pd=0;pd<imgDepth;pd++)
             {
                 ptrPixelData[pd]=ptrInputData[ii+pd*imgSize];
             }             
             ii++;
         }
     }     
     
    imgDims[0]=imgHeight;
    imgDims[1]=imgWidth;
    imgDims[2]=imgDepth;
    //Perform DBSCAN algorithm
    ExecuteDBSCAN(ptrImgData, imgDims, thrDist, maxSpSize,spLabel, &numSpMid);
    // Merging small superpixels
    SuperpixelRelabeling(spLabel, imgDims, numSpIn, spLabelM, &numSpOut);
    
    // Assign output augments
    plhs[0] = mxCreateNumericMatrix(imgHeight,imgWidth,mxINT32_CLASS,mxREAL);
    spLabelOut=(int *)mxGetData(plhs[0]);
    ii=0;
    for (px=0;px<imgWidth;px++)
    {
        for (py=0;py<imgHeight;py++)
        {
             idxPixel=py*imgWidth+px;
             spLabelOut[ii]=spLabelM[idxPixel];
             ii++;
        }
    }
    plhs[1] = mxCreateNumericMatrix(1,1,mxINT32_CLASS,mxREAL);
    ptrNumSpOut = (int*)mxGetData(plhs[1]);
    *ptrNumSpOut=numSpOut;
    
    // Deallocate memory
    mxFree(ptrImgData);
    mxFree(spLabel);
    mxFree(spLabelM);
}

