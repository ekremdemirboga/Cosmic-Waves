from re import L
import numpy as np


def diff1_o4(u, dx):
    """
    First finite difference derivative of 4th order
    with one-sided differences at the boundaries
    u   : a vector
    dx  : the cell size
    du  : derivative of vector u    
    """
    du = np.zeros_like(u)
    # sided differences at the left boundary
    du[0] = -25.*u[0]+48.*u[1]-36.*u[2]+16.*u[3]-3.*u[4]
    du[1] = -3.*u[0]-10.*u[1]+18.*u[2] - 6.*u[3] + u[4]


    # sided differences at the right boundary
    du[-1] = 25.*u[-1]-48.*u[-2]+36.*u[-3]-16.*u[-4]+3.*u[-5]
    du[-2] = 3.*u[-1]+10.*u[-2]-18.*u[-3] + 6.*u[-4] - u[-5]


    # centered differences in the bulk
    du[2:-2] = -u[4:]+8.*u[3:-1]-8.*u[1:-3]+u[:-4]
    # divide by cell size
    du /= 12.*dx
    return du


def diff2_o4(u, dx):
    """
    Second finite difference derivative of 4th order
    with one-sided differences at the boundaries
    u   : a vector
    hg  : the cell size
    d2u : second derivative of vector u
    """
    d2u = np.zeros_like(u)
    # one-sided differences at the boundaries
    d2u[0] = 45.*u[0]-154.*u[1]+214.*u[2]-156.*u[3]+61.*u[4]-10.*u[5]
    d2u[1] = 10.*u[0] - 15.*u[1] - 4.*u[2] + 14.*u[3] - 6.*u[4] + u[5]

    d2u[-1] = 45.*u[-1]-154.*u[-2]+214.*u[-3]-156.*u[-4]+61.*u[-5]-10.*u[-6]
    d2u[-2] = 10.*u[-1] - 15.*u[-2] - 4.*u[-3] + 14.*u[-4] - 6.*u[-5] + u[-6]

    # centered differences in the bulk
    d2u[2:-2] = -u[4:]+16.*u[3:-1]-30.*u[2:-2]+16.*u[1:-3]-u[:-4]
    # divide by cell size
    d2u /= 12.*dx**2
    return d2u

def rk4(cur, coefs, dt, dx, tds):
    """
    Runge--Kutta time integration
    cur   : a 2 x N matrix containing the function and its time derivative
    coefs : a 5 x N matrix containing the coefficients for the wave equation
    dt    : the time grid size
    dx    : the spatial grid size
    tds   : downsampling factor in time
    returns nex, the evolved function and its time derivative, same format as cur
    """
    nex = np.zeros_like(cur)
    nex[:, :] = cur
    # Runge-Kutta RK4
    for i in range(tds):
        k1 = dt * wave1drhs(nex,        coefs, dx)
        k2 = dt * wave1drhs(nex+0.5*k1, coefs, dx)
        k3 = dt * wave1drhs(nex+0.5*k2, coefs, dx)
        k4 = dt * wave1drhs(nex+k3,     coefs, dx)
        nex += (k1+2.*(k2+k3)+k4)/6.
    return nex

def gaussian(x, x0, sigma):
    """
    Normalized Gaussian function
    x     : grid, or time
    x0    : center of gaussian
    sigma : width of gaussian
    """
    return np.exp(-(x-x0)**2/sigma**2)/(np.sqrt(2.*np.pi)*sigma)


def gaussian_spti(x, t, x0, t0, sigma):
    """
    Normalized Gaussian function in space and time (for Green function source)
    x     : spatial grid
    t     : current time
    x0    : spatial center of gaussian
    t0    : time center of gaussian
    sigma : width of gaussian
    """
    return gaussian(t, t0, sigma)*gaussian(x, x0, sigma)

def set_coefs(r, ell, re=1, rc=10, L=1,coord = 'Parikh',metric = 'dS'):
    """
    c0 : A_tr
    c1 : A_rr
    c2 : B_t
    c3 : B_r
    c4 : C
    """
    coefs = np.zeros((5, len(r)))

    # Coefficients for Schwarzschild-de Sitter wave equation in hyperboloidal coordinates
    if metric == 'dS':
        if (coord == 'P'):
            coefs[0, :] = 2*r/L
            coefs[1, :] = 1-r**2/L**2
            coefs[2, :] = -1/L
            coefs[3, :] = -2*r/L**2
            coefs[4, :] = -2/L**2 #+ ell**2/r**2
        elif (coord == 'E'):
            coefs[0, :] = 2*r**2/(L**2+r**2)
            coefs[1, :] = -1 + 2*L**2/(L**2+r**2)
            coefs[2, :] = -2*r/(L**2+r**2)
            coefs[3, :] = -2*r/(L**2+r**2)
            coefs[4, :] = (ell**2*L**2 -2*r**2)/(L**2*r**2 + r**4)
    elif metric == 'SdS':
        ldS2 = re**2 + re*rc + rc**2
        if (coord == 'choice1'):
            ## old ones
            coefs[0, :] = -0.5*((2*r - rc - re)*(rc - re)*(r + rc + re))/(ldS2*r)
            coefs[1, :] = -0.25*((r - rc)*(r - re)*(rc - re)**2*(r + rc + re)**2)/(ldS2**2*r**2)
            coefs[2, :] = -0.5*((rc - re)*(r + rc + re))/(ldS2*r)
            coefs[3, :] = ((rc - re)**2*(r + rc + re)*(-2*r**3 + rc*re*(rc + re)))/(4.*ldS2**2*r**3)
            coefs[4, :] = -0.25*((rc - re)**2*(r + rc + re)*(ell**2*ldS2*r - 2*r**3 + rc*re*(rc + re)))/(ldS2**2*r**4)



    # # Coefficients for the time independent RWZ equation
    # coefs[0, :] = ((-2 + x)*(2 - 9*x + 5*x**2))/(8.*(-5 + 3*x))
    # coefs[1, :] = -((-2 + x)**2*(-1 + x)**2*x)/(32.*(-5 + 3*x))
    # coefs[2, :] = ((-5 + x)*(-1 + x))/(4.*(1 + x)*(-5 + 3*x))
    # coefs[3, :] = -((-2 + x)**2*(-1 + x)*(-1 + 3*x +
    #                                       2*x**2))/(32.*(1 + x)*(-5 + 3*x))
    # # For scalar perturbations
    # # coefs[4,:] = ((-2 + x)**2*(1 + ell*(1 + ell) - \
    # #                            x + ell*(1 + ell)*x))/(16.*(1+x)**2*(-5+3*x))
    # # For gravitational perturbations
    # coefs[4, :] = ((-2 + x)**2*(-3 + ell*(1 + ell) +
    #                             3*x + ell*(1 + ell)*x))/(16.*(1+x)**2*(-5+3*x))
    return coefs


def wave1drhs(u, coefs, dx):
    """
    Evaluate the RHS for the wave equation with time independent but otherwise
    general coefficients and a general source term
    The equation is written of the form
    rhs1 = u_t
    rhs2 = c0 u_tr + c1 u_rr + c2 u_t + c3 u_r + c4 u + s
    u     : 2 x N matrix representing the unknown and its time derivative (u,u_t)
    coefs : 5 x N matrix for the coefficients of the general wave equation
    dx    : the grid cell size
    returns the right-hand side of the wave equation
    """
    rhs = np.zeros_like(u)
    # Time derivative.
    rhs[0, :] = u[1, :]
    # Wave equation
    rhs[1, :] = coefs[0, :]*diff1_o4(u[1, :], dx) + coefs[1, :]*diff2_o4(u[0, :], dx)\
        + coefs[2, :]*u[1, :] + coefs[3, :]*diff1_o4(u[0, :], dx) \
        + coefs[4, :]*u[0, :]
    return rhs


def wave1dsolve(initial, coefs, dx, dt, tds, nsteps):
    """
    Solve the 1D wave equation with general coefficients:
    initial : the initial data, a 2 x N array (function and time derivative)
    coefs   : the coefficients of the wave equation
    dx      : the spatial grid size
    dt      : the time grid size
    nsteps  : the number of time steps to make
    returns data, the solution in the form of an nsteps x 2 x N array
    """
    data = np.zeros(((nsteps+1,)+initial.shape))
    data[0, :, :] = initial
    for i in range(1, nsteps+1):
        data[i, :, :] = rk4(data[i-1, :, :], coefs, dt, dx, tds)
        print("wave1dsolve(): step "+str(i))
    print("Finished evolution at t = "+str(dt*nsteps*tds))
    return data
