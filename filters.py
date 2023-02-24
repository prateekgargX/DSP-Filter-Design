import numpy as np

def Butterworth(delta1,delta2,omega_p,omega_s):
  def params(delta1,delta2,omega_p,omega_s):
    d1 = 1/(1-delta1)**2 - 1
    d2 = 1/delta2**2 - 1

    # here we take "optimal" N
    N = np.ceil(0.5*(np.log(d2)-np.log(d1))/(np.log(omega_s)-np.log(omega_p)))

    # Upper and lower bounds on omega_c
    U_omega_c = omega_s/(d2**(1/(2*N)))
    L_omega_c = omega_p/(d1**(1/(2*N)))

    # in principle L_omega_c < omega_c < U_omega_c, omega_c = np.random.uniform(L_omega_c,U_omega_c)
    # but for reproducibilty,
    omega_c = (L_omega_c+U_omega_c)/2 
    return int(N),omega_c

  N,omega_c = params(delta1,delta2,omega_p,omega_s)
  

  theta =  ((2*np.linspace(0,N-1,N)+1)/(2*N)+0.5)*np.pi*(1j)
  poles = omega_c*np.exp(theta)

  def systemfunc(x):
    deno = np.prod(x-poles) 
    return (omega_c**N)/deno

  return systemfunc,N,omega_c,poles 

