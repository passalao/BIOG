#!/usr/bin/python
# -*- coding: cp1252 -*-
import numpy as np
class Trace:
    def __init__(self, nbvar=1):
        self.nbvar=nbvar;
        self.Data=np.zeros(nbvar);
        self.Cost=[1e-6];
