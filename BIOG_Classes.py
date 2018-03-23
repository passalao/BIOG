#!/usr/bin/python
# -*- coding: cp1252 -*-
import numpy as np
import sys

class Trace:
    def __init__(self, nbvar=1):
        self.nbvar=nbvar;
        self.Data=np.zeros(nbvar);
        self.Cost=[sys.float_info.min];
        self.DeltaTb=[];#Difference between mod and obs Tb
        self.Bias=[sys.float_info.min];#To correct the temperature bias of RT model
class StoreArray:
    def __init__(self, nbx=1, nby=1, nbvar=1):
        self.Data=np.zeros((nbx, nby, nbvar));

class Site:
    def __init__(self, Names=[], Coords=[[],[]]):
        self.Names=Names#[]*nbsites;
        self.Coords=Coords#np.zeros((nbsites,2))