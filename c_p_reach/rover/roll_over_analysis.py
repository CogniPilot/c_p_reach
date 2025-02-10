import pandas as pd
import sympy
import sympy.physics.mechanics as me
import numpy as np 
from matplotlib import pyplot as plt

def analyze(rover):
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

        # # Individual points
        # sim_points
        # for sim_p in sim_points:
        #     if sim_p[2] == 1:
        #         sim_color = "red"
        #     else:
        #         sim_color = "green"
        #plt.scatter(sim_p[0], sim_p[1], color=sim_color, zorder=5)

        # Save plot
        plt.savefig('fig/roll_over_velocity.png',format='png')

    #lx = 1/2 rover width
    lx = rover['width']/2#0.105
    #ly = COM to ground
    ly = rover['COM_height']#0.06
    # Turn radius
    rad=rover['turn_radius']
    plot_roll_over_analysis(theta_f_val=np.arctan(ly/lx), rad=rad)