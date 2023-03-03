from sympy import *
import numpy as np

def round_expr(expr, num_digits):
    return expr.xreplace({n : round(n, num_digits) for n in expr.atoms(Number)})

def normalize(f_samp,L_freq):
    return (np.array(L_freq)/f_samp)*2*np.pi

def to_analog(L_w):
    return np.tan(np.array(L_w)/2)

def Bp_to_lp(Op1,Op2,Os1,Os2):
    B = Op2-Op1
    Oo = np.sqrt(Op2*Op1)

    def transfer_fn(s):
        return (s**2+Oo**2)/(B*s)
    def freq_transform(Omega):
        return (transfer_fn(Omega*1j)/(1j)).real
    
    return np.vectorize(freq_transform)([Op1,Op2,Os1,Os2]),B,Oo,transfer_fn

def sb_to_lp(Op1,Op2,Os1,Os2):
    B = Op2-Op1
    Oo = np.sqrt(Op2*Op1)

    def transfer_fn(s):
        return (B*s)/(s**2+Oo**2)
    def freq_transform(Omega):
        return (transfer_fn(Omega*1j)/(1j)).real
    
    return np.vectorize(freq_transform)([Op1,Op2,Os1,Os2]),B,Oo,transfer_fn
