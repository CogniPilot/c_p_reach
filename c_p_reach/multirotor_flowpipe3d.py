import numpy as np
import matplotlib.pyplot as plt
import casadi as ca
import os

from c_p_reach.lie.SE23 import *
from c_p_reach.flowpipe.inner_bound import *
from c_p_reach.flowpipe.outer_bound import *
from c_p_reach.flowpipe.flowpipe import *
from c_p_reach.sim.multirotor_control import *
from c_p_reach.sim.multirotor_plan import *

import argparse
import pkg_resources

def run():
    version = pkg_resources.get_distribution('c_p_reach').version
    parser = argparse.ArgumentParser(f'multirotor_flowpipe3d {version:s}')
    parser.add_argument("--trajectory", default="traj.json")
    parser.add_argument("--position_target", default=2.0)
    parser.add_argument("--attitude_target", default=0.5)
    parser.add_argument("--gyro_noise", default=0.01)
    parser.add_argument("--plots", action='store_true')
    args = parser.parse_args()

    do_plots = args.plots
    os.makedirs("fig", exist_ok=True)

    print('\n==============================================================')
    print('C-P Reach  version:', version)
    print('==============================================================')

    print('beginning reachability analysis for multirotor')
    n_legs = 10
    poly_deg = 7
    min_deriv = 4  # min snap
    bc_deriv = 4

    # Set disturbance here
    w1 = 0.0 # disturbance for translational (impact a)  thrust disturbance for outer loop
    w2 = args.gyro_noise # disturbance for angular (impact alpha)  inner loop angular disturbance  BKd

    print('reading trajectory')
    bc = np.array(
            [  # boundary conditions
                [
                    [0, 0, 0],
                    [1, 0, 0],
                    [1, 1, 1],
                    [2, 1, 1],
                    [2, 2, 1],
                    [1, 2, 0],
                    [0, 2, 0],
                    [-1, 2, 0],
                    [-2,2,0],
                    [-2,1,0],
                    [-2,0,0]
                ],  # pos
                [
                    [0, 0, 0],
                    [0.3, 0, 0],
                    [0, 0.3, 0.3],
                    [0.3, 0, 0],
                    [0, 0.3, 0],
                    [-0.3, 0, 0],
                    [-0.3, 0, -0.3],
                    [-0.3, 0, 0],
                    [-0.3, 0, 0],
                    [0, -0.3, 0],
                    [0, 0, 0]
                ],  # vel
                [
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]
                ],  # acc
                [
                    [0, 0, 0], 
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0],
                    [0, 0, 0]
                ],
            ]  # jerk
        )
    k_time = 1e5

    print('finding cost function')
    cost = find_cost_function(
        poly_deg=poly_deg,
        min_deriv=min_deriv,
        rows_free=[],
        n_legs=n_legs,
        bc_deriv=bc_deriv,
    )

    print('planning trajectory')
    ref = planner(bc, cost, n_legs, poly_deg, k_time)
    ax = [np.max(ref['ax'])]
    ay = [np.max(ref['ay'])]
    az = [-np.min(ref['az'])+9.8]
    omega1 = [np.max(ref['omega1'])]
    omega2 = [np.max(ref['omega2'])]
    omega3 = [np.max(ref['omega3'])]

    if do_plots:
        fig = plt.figure(figsize=(8,8))
        traj_x = compute_trajectory(ref['x'], ref['T'], poly_deg=poly_deg, deriv=0) # points
        traj_y = compute_trajectory(ref['y'], ref['T'], poly_deg=poly_deg, deriv=0)
        traj_z = compute_trajectory(ref['z'], ref['T'], poly_deg=poly_deg, deriv=0)
        axis = fig.add_subplot(111, projection="3d")
        axis.plot(ref["x"], ref["y"], ref["z"]);
        axis.set_xlabel('x, m', labelpad=10)
        axis.set_ylabel('y, m', labelpad=12)
        axis.set_zlabel('z, m', rotation=90, labelpad=8)
        axis.set_title('Reference Trajectory')
        plt.axis('auto')
        plt.tight_layout()

        plt.savefig('fig/reference_trajectory.png')
        plt.close()

    print('solving LMI for dynamics')
    # solve LMI
    # J omega_dot = k_rollrate * (roll_rate_ref - (roll_rate + gyro_disturbance))
    # J omega_dot = k_rollrate * (roll_rate_ref - roll_rate) + k_rollrate * gyro_disturbance

    sol = find_omega_invariant_set(omega1, omega2, omega3) 
    mu_inner = sol['mu1']
    #print(mu_inner)

    # Initial condition
    P = sol['P']
    e0 = np.array([0,0,0]) # initial error
    beta = (e0.T@P@e0) # initial Lyapnov value

    # find bound
    omegabound = omega_bound(omega1, omega2, omega3, w2, beta) # result for inner bound
    #print(omegabound)

    print('solving LMI for kinematics')
    # solve LMI
    sol_LMI = find_se23_invariant_set(ax, ay, az, omega1, omega2, omega3)
    mu_outer = sol_LMI['mu3']
    #print(mu_outer)

    # Initial condition
    e = np.array([0,0,0,0,0,0,0,0,0]) # initial error in Lie group (nonlinear)

    # transfer initial error to Lie algebra (linear)
    e0 = ca.DM(SE23Dcm.vee(SE23Dcm.log(SE23Dcm.matrix(e))))
    e0 = np.array([e0]).reshape(9,)
    ebeta = e0.T@sol_LMI['P']@e0

    print('finding invariant set')

    # find invairant set points in Lie algebra (linear)
    points, val = se23_invariant_set_points(sol_LMI, 20, w1, omegabound, ebeta)
    points_theta, val = se23_invariant_set_points_theta(sol_LMI, 20, w1, omegabound, ebeta)

    # map invariant set points to Lie group (nonlinear)
    inv_points = exp_map(points, points_theta)

    # currently assumes
    # x_dot = Ax + BKx + d
    # want to change to 
    # x_dot = Ax + BK(x + d) = (A + BK)x + BKd

    # BK = 3 for inner loop
    # x_inf  < mu d_inf
    BK = 3

    mu_total = (mu_outer*mu_inner)*BK
    gyro_noise_req_for_attitude_target = args.attitude_target/mu_total
    gyro_noise_req_for_position_target = args.position_target/mu_total
    print('\n\n==============================================================')
    print('RESULTS')
    print('==============================================================')
    print('mu_total', mu_total)
    print('gyro noise req for attitude target:', gyro_noise_req_for_attitude_target)
    print('gyro noise req for position target:', gyro_noise_req_for_position_target)

    if do_plots:
        print('plotting 2D invariant sets')
        plt.figure(figsize=(14,7))
        plt.rcParams.update({'font.size': 12})
        ax1 = plt.subplot(121)
        ax1.plot(points[0, :], points[1, :], 'g', label='with Dynamic Inversion')
        # ax.plot(pointscl[0, :], pointscl[1, :], 'b', linewidth=0.5, label='without Dynamic Inversion')
        ax1.set_xlabel('$\\zeta_x$, m')
        ax1.set_ylabel('$\\zeta_y$, m')
        # ax.plot(xopt.x,xopt.y,'ro')
        plt.axis('equal')
        plt.grid(True)
        # plt.legend(loc=1)
        ax2 = plt.subplot(122)
        ax2.plot(inv_points[0, :-1], inv_points[1, :-1], 'g', label='with Dynamic Inversion')
        # ax2.plot(inv_pointscl[0, :-1], inv_pointscl[1, :-1], 'b', linewidth=0.5, label='without Dynamic Inversion')
        # ax2.plot(e[0],e[1],'ro')
        ax2.set_xlabel('$\\eta_x$, m')
        ax2.set_ylabel('$\\eta_y$, m')
        plt.grid(True)
        # plt.legend(loc=1)
        plt.axis('equal')
        plt.tight_layout()
        ax1.set_title('Invariant Set in Lie Algebra', fontsize=20)
        ax2.set_title('Invariant Set in Lie Group', fontsize=20)
        plt.savefig('fig/Invariant_l.png')
        plt.close()

    if do_plots:
        print('plotting 3d invariant sets')
        plt.figure(figsize=(14,7))
        ax1 = plt.subplot(121, projection='3d', proj_type='ortho', elev=40, azim=20)
        # ax.plot3D(e0[0], e0[1], e0[2], 'ro');
        ax1.plot3D(points[0, :], points[1, :], points[2, :],'g', label='with Dynamic Inversion')
        # ax.plot3D(pointscl[0, :], pointscl[1, :], pointscl[2, :],'b', linewidth=0.5, label='without Dynamic Inversion')
        ax1.set_xlabel('$\\zeta_x$, m')
        ax1.set_ylabel('$\\zeta_y$, m')
        ax1.set_zlabel('$\\zeta_z$, rad', labelpad=1)
        ax1.set_title('Invariant Set in Lie Algebra', fontsize=20)
        # plt.subplots_adjust(left=9, right=10, top=0.5, bottom=0.08)
        # plt.tight_layout()
        plt.axis('auto')
        plt.tight_layout(pad=0.4, w_pad=0.5, h_pad=1.0)
        # plt.legend(loc=1)
        ax2 = plt.subplot(122, projection='3d', proj_type='ortho', elev=40, azim=20)
        # ax2.plot3D(e[0], e[1], e[2], 'ro');
        ax2.plot3D(inv_points[0, :], inv_points[1, :], inv_points[2, :], 'g', label='with Dynamic Inversion')
        # ax2.plot3D(inv_pointscl[0, :], inv_pointscl[1, :], inv_pointscl[2, :], 'b', linewidth=0.5, label='without Dynamic Inversion')
        # ax2.plot3D(e[0], e[1], e[2], 'ro');
        ax2.set_xlabel('$\\eta_x$, m')
        ax2.set_ylabel('$\\eta_y$, m')
        ax2.set_zlabel('$\\eta_z$, rad')
        ax2.set_title('Invariant Set in Lie Group', fontsize=20)
        plt.axis('auto')
        plt.subplots_adjust(left=0.45, right=1, top=0.5, bottom=0.08)
        # plt.legend(loc=1)
        plt.tight_layout()
        plt.savefig('fig/Invariant3d_l.png')
        plt.close()


    if do_plots:
        print('calculating interval hull')
        # Calculate convex hull for flow pipes
        n = 30 # number of flow pipes
        flowpipes_traj, intervalhull_traj, nom_traj, t_vect = flowpipes(ref, n, ebeta, w1, omegabound, sol_LMI, 'xy')
        plt.savefig('fig/interval_hull.png')
        plt.close()

        print('plotting flow pipes')
        plot_flowpipes(nom_traj, flowpipes_traj, n, 'xy')
        plt.savefig('fig/flow_pipes.png')
        plt.close()

        # print('plotting sim')
        # plot_sim(ref, w1, omegabound, flowpipes_traj, n, 'xy')

        # print('plotting time history')
        # plot_timehis(sol_LMI, ref, w1, w2, 40, ebeta)
