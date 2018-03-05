#!/usr/bin/python
# -*- coding: cp1252 -*-
import numpy as np
class Trace:
    def __init__(self, nbvar=1):
        self.nbvar=nbvar;
        self.Data=np.zeros(nbvar);
        self.Cost=[1e-15];
        self.DeltaTb=[];#Difference between mod and obs Tb
        self.Bias=[];#To correct the temperature bias of RT model
class StoreArray:
    def __init__(self, nbx=1, nby=1, nbvar=1):
        self.Data=np.zeros((nbx, nby, nbvar));