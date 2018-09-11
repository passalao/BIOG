import scipy.optimize as opt
import numpy as np
def Func(t):
    x=t[0]
    y=t[1]
    return x**2+y**2

def Cons(t):
    x=t[0]
    y=t[1]
    return x+y-1

def JacFunc(t):
    x=t[0]
    y=t[1]
    Jac=np.array([2*x,2*y])
    return Jac

def JacCons(t):
    x=t[0]
    y=t[1]
    Jac=np.array([1,1])
    return Jac

def FunWithCons(t):
    x=t[0]
    y=t[1]
    return Func(t)+Lambda*Cons(t)

def JacFunWithCons(t):
    x=t[0]
    y=t[1]
    return JacFunc(t)+Lambda*JacCons(t)

X0=np.array([10,3])
#BestValue = opt.fmin_l_bfgs_b(FunWithCons, X0, fprime=JacFunWithCons)
#BestValue = opt.minimize(Func, X0, jac=JacFun, constraints={"eq", Cons, JacCons})
BestValue = opt.fmin_slsqp(Func, X0, fprime=JacFunc, f_eqcons=Cons, fprime_eqcons=JacCons)

print(BestValue)

