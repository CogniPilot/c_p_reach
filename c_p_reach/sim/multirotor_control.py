import numpy as np
import casadi as ca
import control
from c_p_reach.lie.SE23 import *
from .multirotor_ref_traj import f_ref
from scipy import signal
import matplotlib.pyplot as plt
from c_p_reach.flowpipe.outer_bound import inv_bound

def se23_solve_control(ax,ay,az,omega1,omega2,omega3):
    A = -ca.DM(SE23Dcm.ad_matrix(np.array([0,0,0,ax,ay,az,omega1,omega2,omega3]))+SE23Dcm.adC_matrix())
    B = np.array([[0,0,0,0], # vx
                  [0,0,0,0], # vy
                  [0,0,0,0], # vz
                  [0,0,0,0], # ax
                  [0,0,0,0], # ay
                  [1,0,0,0], # az
                  [0,1,0,0], # omega1
                  [0,0,1,0], # omega2
                  [0,0,0,1]]) # omega3 # control omega1,2,3, and az
    Q = 10*np.eye(9)  # penalize state
    R = 1*np.eye(4)  # penalize input
    K, _, _ = control.lqr(A, B, Q, R) 
    K = -K # rescale K, set negative feedback sign
    BK = B@K
    return B, K, BK , A+B@K

def control_law(B, K, e):
    e = SE23Dcm.vee(e)
    u = SE23Dcm.diff_correction_inv(e)@B@K@e # controller input
    return u

def compute_control(t, y_vect, ref, freq_d, w1_mag, w2_mag, dist): # w1_mag: acceleration, w2_mag: omega
    # reference data (from planner, function of time)
    # x_opt is t coeff of poly
    T_opt = ref['T']
    T = np.cumsum(T_opt)
    for i in range(T.shape[0]):
        if i==0 and t <= T[i]:
            x_opt = np.flip(ref['poly_x'][0:8])
            y_opt = np.flip(ref['poly_y'][0:8])
            z_opt = np.flip(ref['poly_z'][0:8])
            beta = t/T_opt[0]
            break
        elif T[i-1] < t <= T[i]:
            x_opt = np.flip(ref['poly_x'][8*(i):8*(i+1)])
            y_opt = np.flip(ref['poly_y'][8*(i):8*(i+1)])
            z_opt = np.flip(ref['poly_z'][8*(i):8*(i+1)])
            beta = (t-np.sum(T_opt[:i]))/T_opt[i]

    vx_opt = np.polyder(x_opt)
    vy_opt = np.polyder(y_opt)
    vz_opt = np.polyder(z_opt)
    ax_opt = np.polyder(vx_opt)
    ay_opt = np.polyder(vy_opt)
    az_opt = np.polyder(vz_opt)
    jx_opt = np.polyder(ax_opt)
    jy_opt = np.polyder(ay_opt)
    jz_opt = np.polyder(az_opt)
    sx_opt = np.polyder(jx_opt)
    sy_opt = np.polyder(jy_opt)
    sz_opt = np.polyder(jz_opt)

    # reference input at time t
    # world frame
    r_vx = np.polyval(vx_opt, beta)
    r_vy = np.polyval(vy_opt, beta)
    r_vz = np.polyval(vz_opt, beta)
    r_ax = np.polyval(ax_opt, beta)
    r_ay = np.polyval(ay_opt, beta)
    r_az = np.polyval(az_opt, beta)
    r_jx = np.polyval(jx_opt, beta)
    r_jy = np.polyval(jy_opt, beta)
    r_jz = np.polyval(jz_opt, beta)
    r_sx = np.polyval(sx_opt, beta)
    r_sy = np.polyval(sy_opt, beta)
    r_sz = np.polyval(sz_opt, beta)

    # print(r_vy)

    # body frame
    ref = f_ref(0, 0, 0, [r_vx, r_vy, r_vz], [r_ax, r_ay, r_az], [r_jx, r_jy, r_jz], [r_sx, r_sy, r_sz], 1, 9.8, 1, 1, 1, 0)
    R = np.array(ref[1])
    R_eb = np.linalg.inv(R) # earth frame to body frame
    ae = np.array([r_ax,r_ay,r_az])
    ab = R@ae
    omega = ref[2]
    omega = np.array(omega).reshape(3,)
    r_omega1 = omega[0]
    r_omega2 = omega[1]
    r_omega3 = omega[2]
    
    # initial states of vehicle and reference and error
    e = SE23Dcm.wedge(np.array([y_vect[0], y_vect[1], y_vect[2], y_vect[3], y_vect[4], y_vect[5], y_vect[6], y_vect[7], y_vect[8]])) # log error
    
    B, K, _, _ = se23_solve_control(0, 0, 9.8, 0, 0, 0) # time-invariant at hover

    vr = np.array([0,0,0,ab[0],ab[1],ab[2],r_omega1,r_omega2,r_omega3])

    # disturbance
    if dist == 'sine':
        phi = 0.8
        phi2 = 0.5
        wax = np.cos(2*np.pi*freq_d*t+phi)*w1_mag
        way = np.sin(2*np.pi*freq_d*t+phi)*w1_mag/np.sqrt(2)
        waz = np.sin(2*np.pi*freq_d*t+phi)*w1_mag
        womega1 = np.cos(2*np.pi*freq_d*t+phi2)*w2_mag
        womega2 = np.sin(2*np.pi*freq_d*t+phi2)*w2_mag
        womega3 = np.sin(2*np.pi*freq_d*t+phi2)*w2_mag
    elif dist  == 'square':
        # phi = 1
        # phi2 = 0.8
        # wax = np.cos(2*np.pi*freq_d*t+phi)*w1_mag
        # way = np.sin(2*np.pi*freq_d*t+phi)*w1_mag/np.sqrt(2)
        # waz = way
        # womega1 = np.cos(2*np.pi*freq_d*t+phi2)*w2_mag
        # womega2 = np.sin(2*np.pi*freq_d*t+phi2)*w2_mag
        # womega3 = 0
        wax = signal.square(2*np.pi*freq_d*t+np.pi)*w1_mag/np.sqrt(2)
        way = signal.square(2*np.pi*freq_d*t)*w1_mag/2
        waz = signal.square(2*np.pi*freq_d*t)*w1_mag/2
        womega1 = signal.square(2*np.pi*freq_d*t+np.pi)*w2_mag/2
        womega2 = signal.square(2*np.pi*freq_d*t)*w2_mag/2
        womega3 = 0
    w = np.array([0,0,0,wax,way,waz,womega1,womega2,womega3])
    # print(w)
    
    # control law applied to log-linear error
    u = np.array(ca.DM(control_law(B, K, e))).reshape(9,)
        
    # log error dynamics
    U = ca.DM(SE23Dcm.diff_correction(SE23Dcm.vee(e)))
    # these dynamics don't hold exactly unless you can move sideways
    A = -ca.DM(SE23Dcm.ad_matrix(vr)+SE23Dcm.adC_matrix())
    # e_dot = (A+ B@K)@ca.DM(SE23Dcm.vee(e)) + U@w # vector form
    e_dot = (A)@ca.DM(SE23Dcm.vee(e)) + U@(u+w)
    e_dot = np.array(e_dot).reshape(9,)

    return [
            # log error
            e_dot[0],
            e_dot[1],
            e_dot[2],
            e_dot[3],
            e_dot[4],
            e_dot[5],
            e_dot[6],
            e_dot[7],
            e_dot[8],
        ]

def simulate_rover(ref, freq_d, w1, w2, x0, y0, z0, vx0, vy0, vz0, theta1_0, theta2_0, theta3_0, dist, plot=False):
    t = np.arange(0,np.sum(ref['T']),0.01)
    X0 = SE23Dcm.matrix(np.array([x0, y0, z0, vx0, vy0, vz0, theta1_0, theta2_0, theta3_0]))  # initial vehicle position in SE2(3)
    X0_r = SE23Dcm.matrix(np.array([0,0,0,0,0,0,0,0,0]))  # initial reference position in SE2(3)
    e0 = SE23Dcm.log(SE23Dcm.inv(X0)@X0_r)  # initial log of error log(X^-1Xr)
    x0 = [ca.DM(SE23Dcm.vee(e0))][0] # intial state for system in vector form
    y0= np.array(x0).reshape(9,)

    import scipy.integrate
    res = scipy.integrate.solve_ivp(
        fun=compute_control,
        t_span=[t[0], t[-1]], t_eval=t,
        y0=y0, args=[ref, freq_d, w1, w2, dist], rtol=1e-6, atol=1e-9)
    return res

def compute_exp_log_err(rx, ry, rz, rvx, rvy, rvz, rtheta1, rtheta2, rtheta3, ex, ey, ez, evx, evy, evz, etheta1, etheta2, etheta3):
    zeta = SE23Dcm.wedge(np.array([ex, ey, ez, evx, evy, evz, etheta1, etheta2, etheta3]))
    eta = SE23Dcm.exp(zeta)
    eta_inv = ca.DM(SE23Dcm.inv(eta))
    Xr = ca.DM(SE23Dcm.matrix(np.array([rx, ry, rz, rvx, rvy, rvz, rtheta1, rtheta2, rtheta3])))
    X = Xr@eta_inv
    x_vect = ca.DM(SE23Dcm.vector(X))
    x_vect = np.array(x_vect).reshape(9,)
    return x_vect

def plot_sim(ref, abound, omegabound, flowpipes, num_pipes, axis):
    freq = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]

    fig = plt.figure(figsize=(15,15))
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.2, top=0.95)

    T0_list = ref['T']
    T = np.cumsum(T0_list)
    
    label_added =False
    for f in freq:
        res = simulate_rover(ref, f, abound, omegabound, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 'sine')
        t = res['t']
    
        y_vect = res['y']
        ex, ey, ez, evx, evy, evz, etheta1, etheta2, etheta3 = [y_vect[i, :] for i in range(len(y_vect))]
        exp_log_err = np.zeros((9,len(t)))

        for j in range(len(t)):
            for i in range(T.shape[0]):
                if i==0 and t[j] <= T[i]:
                    xr = np.flip(ref['poly_x'][0:8])
                    yr = np.flip(ref['poly_y'][0:8])
                    zr = np.flip(ref['poly_z'][0:8])
                    beta = t[j]/T0_list[0]
                    break
                elif T[i-1] < t[j] <= T[i]:
                    xr = np.flip(ref['poly_x'][8*(i):8*(i+1)])
                    yr = np.flip(ref['poly_y'][8*(i):8*(i+1)])
                    zr = np.flip(ref['poly_z'][8*(i):8*(i+1)])
                    beta = (t[j]-np.sum(T0_list[:i]))/T0_list[i]
            vx_opt = np.polyder(xr)
            vy_opt = np.polyder(yr)
            vz_opt = np.polyder(zr)
            ax_opt = np.polyder(vx_opt)
            ay_opt = np.polyder(vy_opt)
            az_opt = np.polyder(vz_opt)
            jx_opt = np.polyder(ax_opt)
            jy_opt = np.polyder(ay_opt)
            jz_opt = np.polyder(az_opt)
            sx_opt = np.polyder(jx_opt)
            sy_opt = np.polyder(jy_opt)
            sz_opt = np.polyder(jz_opt)

            # reference input at time t
            # world frame
            r_x = np.polyval(xr, beta)
            r_y = np.polyval(yr, beta)
            r_z = np.polyval(zr, beta)
            r_vx = np.polyval(vx_opt, beta)
            r_vy = np.polyval(vy_opt, beta)
            r_vz = np.polyval(vz_opt, beta)
            r_ax = np.polyval(ax_opt, beta)
            r_ay = np.polyval(ay_opt, beta)
            r_az = np.polyval(az_opt, beta)
            r_jx = np.polyval(jx_opt, beta)
            r_jy = np.polyval(jy_opt, beta)
            r_jz = np.polyval(jz_opt, beta)
            r_sx = np.polyval(sx_opt, beta)
            r_sy = np.polyval(sy_opt, beta)
            r_sz = np.polyval(sz_opt, beta)
            ref_v = f_ref(0, 0, 0, [r_vx, r_vy, r_vz], [r_ax, r_ay, r_az], [r_jx, r_jy, r_jz], [r_sx, r_sy, r_sz], 1, 9.8, 1, 1, 1, 0)
            R = ref_v[1]
            theta = ca.DM(Euler.from_dcm(R))
            theta = np.array(theta).reshape(3,)
            r_theta1 = theta[0]
            r_theta2 = theta[1]
            r_theta3 = theta[2]
            omega = ref_v[2]
            omega = np.array(omega).reshape(3,)
            exp_log_err[:,j] = np.array([compute_exp_log_err(r_x, r_y, r_z, r_vx, r_vy, r_vz, r_theta1, r_theta2, r_theta3,
                                                             ex[j], ey[j], ez[j], evx[j], evy[j], evz[j], etheta1[j], etheta2[j], etheta3[j])])
            # exp_log_errsq[:,j] = np.array([compute_exp_log_err(r_x, r_y, r_z, r_vx, r_vy, r_vz, r_theta1, r_theta2, r_theta3,
            #                                                 exsq[j], eysq[j], ezsq[j], evxsq[j], evysq[j], evzsq[j], etheta1sq[j], etheta2sq[j], etheta3sq[j])])
        if axis == 'xy':
            if not label_added:
                plt.plot(exp_log_err[0,:], exp_log_err[1,:], 'g', label='sim',linewidth=0.7)
                label_added = True
            else:
                plt.plot(exp_log_err[0,:], exp_log_err[1,:], 'g',linewidth=0.7)
        elif axis == 'xz':
            if not label_added:
                plt.plot(exp_log_err[0,:], exp_log_err[2,:], 'g', label='sim',linewidth=0.7)
                label_added = True
            else:
                plt.plot(exp_log_err[0,:], exp_log_err[2,:], 'g',linewidth=0.7)
            
    if axis == 'xy':
        plt.plot(ref['x'], ref['y'], 'r--', label='ref')
        plt.ylabel('y')
    elif axis == 'xz':
        plt.plot(ref['x'], ref['z'], 'r--', label='ref')
        plt.ylabel('z')
    plt.xlabel('x')
    plt.grid(True)


    # h_nom = plt.plot(nom[:,0], nom[:,1], color='k', linestyle='-')
    for facet in range(num_pipes):
        hs_ch_LMI = plt.plot(flowpipes[facet][:,0], flowpipes[facet][:,1], color='c', linestyle='--')

    # plt.axis('equal')
    plt.title('Flow Pipes')
    # plt.xlabel('x')
    # plt.ylabel('z')
    lgd = plt.legend(loc=2, prop={'size': 18})
    ax = lgd.axes
    handles, labels = ax.get_legend_handles_labels()
    handles.append(hs_ch_LMI[0])
    labels.append('Flow Pipes')
    lgd._legend_box = None
    lgd._init_legend_box(handles, labels)
    lgd._set_loc(lgd._loc)
    lgd.set_title(lgd.get_title().get_text())

    return 

def plot_timehis(sol_LMI, ref, abound, omegabound, n_time, ebeta):
    fig = plt.figure(figsize=(15,15))
    plt.rcParams.update({'font.size': 18})
    fig.subplots_adjust(hspace=0.2, top=0.95)

    # calculte bound along time (small disturbance case)
    T0_list = ref['T']
    T = np.cumsum(T0_list)
    t_vect = np.linspace(1e-5,T[-1],n_time)
    invbound = np.zeros((6,n_time))
    for j in range(0,n_time):
        for i in range(T.shape[0]):
            if i==0 and t_vect[j] <= T[i]:
                xr = np.flip(ref['poly_x'][0:8])
                yr = np.flip(ref['poly_y'][0:8])
                zr = np.flip(ref['poly_z'][0:8])
                beta = t_vect[j]/T0_list[0]
                break
            elif T[i-1] < t_vect[j] <= T[i]:
                xr = np.flip(ref['poly_x'][8*(i):8*(i+1)])
                yr = np.flip(ref['poly_y'][8*(i):8*(i+1)])
                zr = np.flip(ref['poly_z'][8*(i):8*(i+1)])
                beta = (t_vect[j]-np.sum(T0_list[:i]))/T0_list[i]
        rx = np.polyval(xr, beta)
        ry = np.polyval(yr, beta)
        rz = np.polyval(zr, beta)
        ib = inv_bound(sol_LMI, t_vect[i], abound, omegabound, ebeta)
        ib[0] = rx + ib[0]
        ib[1] = ry + ib[1]
        ib[2] = rz + ib[2]
        ib[3] = rx + ib[3]
        ib[4] = ry + ib[4]
        ib[5] = rz + ib[5]
        invbound[:,j] = ib

    freq = [0.01] #[0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09]
    label_added =False
    for f in freq:
        res = simulate_rover(ref, f, abound, omegabound, 0, 0, 0, 0, 0, 0, 0, 0, 0, 'sine')
        t = res['t']
    
        y_vect = res['y']
        ex, ey, ez, evx, evy, evz, etheta1, etheta2, etheta3 = [y_vect[i, :] for i in range(len(y_vect))]
        exp_log_err = np.zeros((9,len(t)))

        for j in range(len(t)):
            for i in range(T.shape[0]):
                if i==0 and t[j] <= T[i]:
                    xr = np.flip(ref['poly_x'][0:8])
                    yr = np.flip(ref['poly_y'][0:8])
                    zr = np.flip(ref['poly_z'][0:8])
                    beta = t[j]/T0_list[0]
                    break
                elif T[i-1] < t[j] <= T[i]:
                    xr = np.flip(ref['poly_x'][8*(i):8*(i+1)])
                    yr = np.flip(ref['poly_y'][8*(i):8*(i+1)])
                    zr = np.flip(ref['poly_z'][8*(i):8*(i+1)])
                    beta = (t[j]-np.sum(T0_list[:i]))/T0_list[i]

            vx_opt = np.polyder(xr)
            vy_opt = np.polyder(yr)
            vz_opt = np.polyder(zr)
            ax_opt = np.polyder(vx_opt)
            ay_opt = np.polyder(vy_opt)
            az_opt = np.polyder(vz_opt)
            jx_opt = np.polyder(ax_opt)
            jy_opt = np.polyder(ay_opt)
            jz_opt = np.polyder(az_opt)
            sx_opt = np.polyder(jx_opt)
            sy_opt = np.polyder(jy_opt)
            sz_opt = np.polyder(jz_opt)

            # reference input at time t
            # world frame
            r_x = np.polyval(xr, beta)
            r_y = np.polyval(yr, beta)
            r_z = np.polyval(zr, beta)
            r_vx = np.polyval(vx_opt, beta)
            r_vy = np.polyval(vy_opt, beta)
            r_vz = np.polyval(vz_opt, beta)
            r_ax = np.polyval(ax_opt, beta)
            r_ay = np.polyval(ay_opt, beta)
            r_az = np.polyval(az_opt, beta)
            r_jx = np.polyval(jx_opt, beta)
            r_jy = np.polyval(jy_opt, beta)
            r_jz = np.polyval(jz_opt, beta)
            r_sx = np.polyval(sx_opt, beta)
            r_sy = np.polyval(sy_opt, beta)
            r_sz = np.polyval(sz_opt, beta)
            ref_v = f_ref(0, 0, 0, [r_vx, r_vy, r_vz], [r_ax, r_ay, r_az], [r_jx, r_jy, r_jz], [r_sx, r_sy, r_sz], 1, 9.8, 1, 1, 1, 0)
            R = ref_v[1]
            theta = ca.DM(Euler.from_dcm(R))
            theta = np.array(theta).reshape(3,)
            r_theta1 = theta[0]
            r_theta2 = theta[1]
            r_theta3 = theta[2]
            omega = ref_v[2]
            omega = np.array(omega).reshape(3,)
            exp_log_err[:,j] = np.array([compute_exp_log_err(r_x, r_y, r_z, r_vx, r_vy, r_vz, r_theta1, r_theta2, r_theta3,
                                                            ex[j], ey[j], ez[j], evx[j], evy[j], evz[j], etheta1[j], etheta2[j], etheta3[j])])
            # exp_log_errsq[:,j] = np.array([compute_exp_log_err(r_x, r_y, r_z, r_vx, r_vy, r_vz, r_theta1, r_theta2, r_theta3,
            #                                                 exsq[j], eysq[j], ezsq[j], evxsq[j], evysq[j], evzsq[j], etheta1sq[j], etheta2sq[j], etheta3sq[j])])
    
        if not label_added:
            plt.plot(t, exp_log_err[2,:], 'g', label='sim z',linewidth=0.5)
            label_added = True
        else:
            plt.plot(t, exp_log_err[2,:], 'g',linewidth=0.5)

    t_vect = np.linspace(1e-5,np.cumsum(T0_list)[-1],n_time)
    plt.plot(ref['t'], ref['z'], 'r', label='ref z')
    plt.plot(t_vect, invbound[2,:], 'c', label='LMI')
    plt.plot(t_vect, invbound[5,:], 'c')
    plt.xlabel('t, sec')
    plt.ylabel('z, m')
    plt.grid(True)
    plt.legend(loc=2)
