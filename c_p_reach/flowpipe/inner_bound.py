import numpy as np
import picos
import control
import itertools
import scipy

def omega_solve_control_gain(omega1, omega2, omega3):
    # A = np.zeros((3,3))
    A =  -np.array([[0, -omega3, omega2],
                    [omega3, 0, -omega1],
                    [-omega2, omega1, 0]])
    B = np.array([[1, 0, 0],
                  [0, 1, 0],
                  [0, 0, 1]]) # control alpha1,2,3
    Q = 10*np.eye(3)  # penalize state
    R = 1*np.eye(3)  # penalize input
    K, _, _ = control.lqr(A, B, Q, R) 
    K = -K # rescale K, set negative feedback sign
    BK = B@K
    return B, K, BK , A+B@K

def omegaLMIs(alpha, A, B, verbosity=0):
    prob = picos.Problem()
    P = picos.SymmetricVariable('P', (3, 3))
    mu1 = picos.RealVariable('mu_1')
    for Ai in A:
        block_eq1 = picos.block([
                [Ai.T*P + P*Ai + alpha*P, P.T*B],
                [B.T*P, -alpha*mu1*np.eye(3)]])
        prob.add_constraint(block_eq1 << 0)
    prob.add_constraint(P >> 1)
    prob.add_constraint(mu1 >> 1e-5)
    prob.set_objective('min', mu1)
    try:
        prob.solve(solver="cvxopt", options={'verbosity': verbosity})
        cost = mu1.value
    except Exception as e:
        print(e)
        cost = -1
    return {
        'cost': cost,
        'prob': prob,
        'mu1': mu1.value,
        'P': np.round(np.array(P.value), 3),
        'alpha':alpha
        }

def find_omega_invariant_set(omega1, omega2, omega3, verbosity=0):
    iterables =[omega1, omega2, omega3]
    omega = []
    for m in itertools.product(*iterables):
        m = np.array(m)
        omega.append(m)
    
    A = []
    eig = []
    for omegai in omega:
        omega1 = omegai[0]
        omega2 = omegai[1]
        omega3 = omegai[2]
        B, K, BK, Ai = omega_solve_control_gain(omega1, omega2, omega3)
        A.append(Ai)
        eig.append(np.linalg.eig(Ai)[0])
    
    # we use fmin to solve a line search problem in alpha for minimum gamma
    if verbosity > 0:
        print('line search')
    
    # we perform a line search over alpha to find the largest convergence rate possible
    alpha_1 = -np.real(np.max(eig)) # smallest magnitude value from eig-value, and range has to be positive
    alpha_opt = scipy.optimize.fminbound(lambda alpha: omegaLMIs(alpha, A, B, verbosity=verbosity)['cost'], x1=1e-5, x2=alpha_1, disp=True if verbosity > 0 else False)
    
    # if the alpha optimization fail, pick a fixed value for alpha.
    sol = omegaLMIs(alpha_opt, A, B)
    prob = sol['prob']
    if prob.status == 'optimal':
        P = prob.variables['P'].value
        mu1 =  prob.variables['mu_1'].value
        if verbosity > 0:
            print(sol)
    else:
        raise RuntimeError('Optimization failed')
        
    return sol

def omega_invariant_set_points(sol, t, w1_norm, beta): # w1_norm: distrubance in alpha  
    P = sol['P']
    val = np.real(beta*np.exp(-sol['alpha']*t) + (sol['mu1']*w1_norm**2)*(1-np.exp(-sol['alpha']*t))) # V(t)
    #print('val', val)
    
    # 1 = xT(P/V(t))x, equation for the ellipse
    evals, evects = np.linalg.eig(P/val)
    radii = 1/np.sqrt(evals)
    R = evects@np.diag(radii)
    R = np.real(R)
    
    # draw sphere
    points = []
    n = 25
    for u in np.linspace(0, 2*np.pi, n):
        for v in np.linspace(0, 2*np.pi, 2*n):
            points.append([np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)])
    for v in np.linspace(0, 2*np.pi, 2*n):
        for u in np.linspace(0, 2*np.pi, n):
            points.append([np.cos(u)*np.sin(v), np.sin(u)*np.sin(v), np.cos(v)])
    points = np.array(points).T
    return R@points, val

def omega_bound(omega1, omega2, omega3, a_dist, beta):
    sol = find_omega_invariant_set(omega1, omega2, omega3)
    points, _ = omega_invariant_set_points(sol, 20, a_dist, beta) 
    max_omega1 = points[0,:].max()
    min_omega1 = abs(points[0,:].min())
    max_omega2 = points[1,:].max()
    min_omega2 = abs(points[1,:].min())
    max_omega3 = points[2,:].max()
    min_omega3 = abs(points[2,:].min())
    bound = np.linalg.norm(np.array([max_omega1,min_omega1,max_omega2,min_omega2,max_omega3,min_omega3]), np.inf)
    return bound