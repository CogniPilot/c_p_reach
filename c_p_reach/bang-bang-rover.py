import numpy as np
import matplotlib.pyplot as plt
import casadi as ca
import json
import os

from c_p_reach.lie.SE23 import *
from c_p_reach.flowpipe.inner_bound import *
from c_p_reach.flowpipe.outer_bound import *
from c_p_reach.flowpipe.flowpipe import *
from c_p_reach.sim.multirotor_control import *
from c_p_reach.sim.multirotor_plan import *

import sympy
import sympy.physics.mechanics as me

import argparse
import pkg_resources

def run():
    version = pkg_resources.get_distribution('c_p_reach').version
    parser = argparse.ArgumentParser(f'multirotor_flowpipe3d {version:s}')
    parser.add_argument("--rover_exploit", action = 'store_true')
    parser.add_argument("--rover_specs", default="rover_specs.json")
    parser.add_argument("--trajectory", default="traj.json")
    parser.add_argument("--position_target", default=2.0)
    parser.add_argument("--attitude_target", default=0.5)
    parser.add_argument("--gyro_noise", default=0.01)
    parser.add_argument("--plots", action='store_true')
    args = parser.parse_args()

    do_plots = args.plots
    do_rover = args.rover_exploit
    os.makedirs("fig", exist_ok=True)

    print('\n==============================================================')
    print('C-P Reach  version:', version)
    print('==============================================================')

    # Rover Analysis
    if do_rover:
        print('Beginning Reachability Analysis for Rover')
        # Retrieve rover specs from json file
        f = open(args.rover_specs)
        rover = json.load(f)

        # Steering Exploit (Reachable set with bang-bang controller and EMI interference with Magnometer)
        pi= math.pi
        epsilon = rover['eps_controller'] * pi/180 # Bang-bang controller value.
        disturbance = rover['emi_disturbance']*pi/180 # Magnetometer Disturbance.
        y_0 = rover['y_0'] # Initial y position
        x_0 = rover['x_0'] # Initial x position
        y_f = rover['y_f'] # Final y position
        x_f = rover['x_f'] # Final x position
        step_size = 0.01


        X = np.arange(x_0,0,step_size)
        Y = np.zeros_like(X)
        Y2 = np.zeros_like(X)
        y = y_0
        y2 = y_0
        
        for i,x in enumerate(X):
            theta_h = math.atan2(y-y_f,x-x_f) + disturbance
            dy_dx = math.tan(math.atan2(y-y_f,x-x_f)+ epsilon - disturbance)
            dy_dx_2 = math.tan(math.atan2(y2-y_f,x-x_f) -epsilon - disturbance)
            Y[i] = y
            Y2[i] = y2
            y = y + dy_dx*step_size
            y2 = y2 + dy_dx_2*step_size
            
            #print(f"x: {x:.2f}, y: {y:.2f}, theta_h: {theta_h*180/pi:.2f}, theta: {math.atan(dy_dx)*180/pi:.2f}, dy_dx: {dy_dx:.3f}")

        plt.plot(X,Y,zorder = 4,color="blue")
        plt.plot(X,Y2,zorder = 4, color="blue")
        plt.fill_between(X, Y, Y2, alpha=0.3, color="blue")
        #plt.fill_between(X, Y2, 0, alpha=0.3, color="blue")
        plt.scatter(x_0, y_0, color="black",zorder=5)
        plt.annotate('Rover', (x_0+0.8, y_0-0.6), textcoords="offset points", xytext=(10,10), ha='center', zorder=10)
        plt.scatter(x_f, y_f, color="black",zorder=5)
        plt.annotate('Goal', (x_f-1, y_f-0.6), textcoords="offset points", xytext=(10,10), ha='center', zorder=10)
        ax = plt.subplot(1,1,1)
        
        # ax.plot(x_0,y_0,'bo')
        ax.set_xlabel('X-position of Rover (m)')
        ax.set_ylabel('Y-position of Rover (m)')
        ax.set_axisbelow(True)
        plt.grid(True)
        # plt.legend(loc=1)
        plt.axis('equal')
        #plt.tight_layout()
        ax.set_title(f'Possible Locations for Rover With a Disturbance of {disturbance*180/pi:.2f}$\\degree$', fontsize=14)
        plt.savefig("fig/rover_reachable_positions.png")
        plt.close()
        

        # Roll-over explot
        m, g, v, r = sympy.symbols('m, g, v, r', real=True)
        theta_t, theta_f, theta, l = sympy.symbols('theta_t, theta_f, theta, l',real=True)
        frame_e = me.ReferenceFrame('e') # earth local level frame
        frame_t = frame_e.orientnew('t', 'Axis', (frame_e.z, theta_t)) # terrain frame
        frame_b = frame_t.orientnew('b', 'Axis', (frame_e.z, theta)) # terrain frame
        frame_f = frame_b.orientnew('f', 'Axis', (frame_e.z, theta_f)) # body frame
        # position vector from tire rotation point to center of mass
        r_ap = l*frame_f.x
        W = -m*g*frame_e.y
        Fc = -(m*v**2/r)*frame_t.x
        F = W + Fc
        M = r_ap.cross(F)

        eq1 = (M.to_matrix(frame_e)[2]/(g*l*m)).simplify()
        eq2 = eq1.subs(theta, 0)
        print(eq1, eq2)
        v_expr = sympy.solve(eq1.subs(theta, 0), v)[1]

        def plot_roll_over_analysis(theta_f_val,rad):
            f_v_expr = sympy.lambdify([g, r, theta_f, theta_t], [v_expr])
            
            plt.figure()
            for r_val in [rad]:
                theta_t_vals = np.linspace(0, np.pi/2-theta_f_val, 1000)
                v_vals = f_v_expr(g=9.8, r=r_val, theta_f=theta_f_val,theta_t=theta_t_vals)[0]
                plt.plot(np.rad2deg(theta_t_vals), v_vals, label='r={:4.0f} m'.format(r_val))
            #plt.plot(9.17,10,'ro',label='unsuccessful')
            #plt.plot(9.17,15,'go',label='successful')
            plt.fill_between(np.rad2deg(theta_t_vals), v_vals, 0, alpha=0.3, color="green",label='safe')
            plt.fill_between(np.rad2deg(theta_t_vals), v_vals, 50, alpha=0.3, color="red",label='unsafe')
            #plt.fill_between(X, Y2, 0, alpha=0.3, color="blue")
            plt.grid()
            plt.xlabel('Terrain Angle [deg]')
            plt.ylabel('Velocity for Roll Over [m/s]')
            plt.title('Roll Over Velocity $\\theta_f$ = {:0.1f} deg'.format(np.rad2deg(theta_f_val)))
            plt.legend()
            plt.ylim(bottom=0,top=17)

            # Individual points
            terrain_angles = [ 0. 10. 20. 30. 40. 50. 60.]
            velocities = [ 4  6  8 10 12 14]
            results = [['S', 'S', 'S', 'S', 'S', 'U'], ['S', 'S', 'S', 'S', 'U', 'U'], ['S', 'S', 'S', 'S', 'U', 'U'], ['S', 'S', 'S', 'S', 'U', 'U'], ['S', 'S', 'S', 'U', 'U', 'U'], ['S', 'S', 'U', 'U', 'U', 'U'], ['S', 'U', 'U', 'U', 'U', 'U']]
            
            for t,t_angle in enumerate(terrain_angles):
                for v,vel in enumerate(velocities):
                    if results[t,v] == "S":
                        sim_color = "green"
                    else:
                        sim_color = "red"
                plt.scatter(t_angle, vel, color=sim_color, zorder=5)

            # Save plot
            plt.savefig('fig/roll_over_velocity.png',format='png')
        
        #lx = 1/2 rover width
        lx = rover['width']/2#0.105
        #ly = COM to ground
        ly = rover['COM_height']#0.06
        # Turn radius
        rad=rover['turn_radius']
        plot_roll_over_analysis(theta_f_val=np.arctan(ly/lx), rad=rad)

                

    # Mulirotor Analysis
    else:
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
