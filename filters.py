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

def Chebyshev(delta1,delta2,omega_p,omega_s):
  def params(delta1,delta2,omega_p,omega_s):
    d1 = 1/(1-delta1)**2 - 1
    d2 = 1/delta2**2 - 1 
    # Here we take "optimal" epsilon
    eps = np.sqrt(d1)
    # Here we take "optimal" N
    N = np.ceil(np.arccosh(np.sqrt(d2/d1))/np.arccosh(omega_s/omega_p))

    return int(N),eps

  def C(n,x):
    if n==0 : return 1
    elif n==1 : return x
    else : return 2*x*C(n-1,x)-C(n-2,x)

  N,eps = params(delta1,delta2,omega_p,omega_s)

  Ak = ((2*(np.linspace(0,2*N-1,2*N))+1)*np.pi)/(2*N)
  Bk = (np.arcsinh(1/eps))/N
  all_poles = omega_p*(np.sin(Ak)*np.sinh(Bk)+(1j)*np.cos(Ak)*np.cosh(Bk))
  poles = all_poles[all_poles.real<0]

  def systemfunc(x):
    nume = np.prod(-poles)*(1/np.sqrt(1+eps**2) if N%2==0 else 1)
    deno = np.prod(x-poles)
    return nume/deno

  return systemfunc,N,eps,poles
  
  
