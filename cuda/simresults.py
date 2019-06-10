#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  9 15:38:41 2019

@author: russell
"""

import numpy as np

def getFileData(prefix):
    paramsFilename = prefix + ".params"
    matsFilename = prefix + ".matrices"
    params = np.loadtxt(paramsFilename)
    mats = np.loadtxt(matsFilename)
    nmatRows = mats.shape[0]
    nparamRows = params.shape[0]
    if nmatRows !=  nparamRows:
        raise AttributeError("Number of params rows {0} ne number of matrices rows {1}!!".format(nmatRows,nparamRows))
    #try to get appropriate matrix size
    matSize = np.sqrt(mats.shape[1])
    if not np.isclose(matSize,np.rint(matSize)):
        raise AttributeError("Not square! sqrt = {0}".format(matSize))
    matSize = int(np.rint(matSize))
    mats = mats.reshape((nmatRows,matSize,matSize))
    #need to transpose the latter axes
    mats = np.transpose(mats,(0,2,1))
    return params,mats

def getSD(mat):
    w,v = np.linalg.eig(mat)
    if np.all(np.isclose(w[0],1.0)):
        statd = v[:,0]
        if not np.all(np.isclose(np.imag(statd),0.0)):
            raise Exception("statdist vector has nonzero complex part!")
        statd = np.real(statd)
        statd = statd / np.sum(statd,keepdims=True)
        return statd
    raise Exception("First eigenvalue is not one!")

def meanandsd(s):
    mn = np.arange(len(s))
    mn2 = mn * mn
    meanval = np.dot(s,mn)
    second = np.dot(s,mn2)
    mysd=np.sqrt(second - (meanval*meanval))
    return meanval,mysd

def condent(mat):
    mlnm=mat*np.ma.log(mat).filled(0)
    return -1*np.sum(mlnm,0)

def entropyrate(stat,ce):
    return np.dot(stat,ce)

def processmatrices(mats):
    statds=np.zeros((mats.shape[0],mats.shape[1]))
    condents=np.zeros((mats.shape[0],mats.shape[1]))
    erates=np.zeros(mats.shape[0])
    means=np.zeros(mats.shape[0])
    stds=np.zeros(mats.shape[0])
    print(statds.shape)
    print(mats.shape)
    for i in range(mats.shape[0]):
        mat=mats[i]
        try:
            stationary = getSD(mat)
            conde = condent(mat)
            erate = entropyrate(stationary,conde)
            m,s = meanandsd(stationary)
            statds[i,:]=stationary
            condents[i,:] = conde
            means[i]=m
            stds[i] = s
            erates[i]= erate
        except Exception as e:
            print(str(e))
    return np.hstack((erates[:,None],means[:,None],stds[:,None],condents))
    
header="Mess\tnMicro\ttf_deg_low\ttf_deg_high\tm_deg_low\tm_deg_high\ttf_p_nz\tdefect\tpmr\tm_p_nz\trho\tp_zero\tp_one\tentropy_rate\texpression_mean\texpression_std"
def main(prefix):
    params,matrices=getFileData(prefix)
    processed=processmatrices(matrices)
    outdata = np.hstack((params,processed))
    condentcols="\t".join(["condent_{0}".format(c) for c in np.arange(matrices.shape[1])])
    hdr = header + "\t" + condentcols
    np.savetxt(prefix+".processed.txt",outdata,delimiter='\t',header=hdr,fmt='%1.6f')