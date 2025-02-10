import math
import numpy as np
from matplotlib import pyplot as plt
def analyze(rover):
    pi= math.pi
    disturbance = rover['emi_disturbance'] # Magnetometer Disturbance.
    y_0 = rover['y_0'] # Initial y position
    x_0 = rover['x_0'] # Initial x position
    y_f = rover['y_f'] # Final y position
    x_f = rover['x_f'] # Final x position
    turn_radius = rover['turn_radius']

    # 1 - Positive disturbance
    # 2 - Negative disturbance     

    # Initial direction
    psi_0 = math.atan2(y_f-y_0, x_f-x_0)

    # Phase 1 (go straight 7m)
    x_1 = 7*math.cos(psi_0) + x_0
    y_1 = 7*math.sin(psi_0) + y_0

    X1 = [x_0, x_1]
    Y1 = [y_0, y_1]
    X2 = [x_0, x_1]
    Y2 = [y_0, y_1]

    print(f"psi_0: {psi_0*180/pi}, x_1: {x_1}, y_1: {y_1}")

    plt.plot(X1,Y1, color="blue")
    plt.plot(X2,Y2, color="blue")
    

    # Phase 2 (turn)
    # Locate center point of arc, based on turn radius.
    gamma = (90-psi_0*180/pi)*pi/180
    center_point = [x_1 + turn_radius*math.cos(gamma), y_1-turn_radius*math.sin(gamma)]
    #print(f"center point: {center_point}")

    psi1_2 = psi_0 + (disturbance)*pi/180
    psi2_2 = psi_0 - (disturbance)*pi/180
    
    #print(f"psi1_2 = {psi1_2*180/pi} psi2_2 = {psi2_2*180/pi}")

    # First trajectory leaves early.
    psi1 = np.linspace(psi_0 + 90*pi/180, psi1_2, 100)
    # Second trajectory leaves late.
    psi2 = np.linspace(psi_0 + 90*pi/180, psi2_2, 100)
    
    # First and second trajectory
    X1 = turn_radius*np.cos(psi1) + center_point[0]
    Y1 = turn_radius*np.sin(psi1) + center_point[1]
    X2 = turn_radius*np.cos(psi2) + center_point[0]
    Y2 = turn_radius*np.sin(psi2) + center_point[1]

    #plt.scatter(center_point[0],center_point[1])

    # plt.xlim(-10, 20)
    # plt.ylim(-20, 10)
    plt.plot(X2,Y2, color="blue")
    plt.plot(X1,Y1, color="blue")


    # Obtain end points.
    X_end_bound = []
    Y_end_bound = []
    for idx in range(len(psi2)):
        if psi2[idx] < psi1_2:
            X_end_bound.append(7*np.sin(psi2[idx]) + X2[idx])
            Y_end_bound.append(-7*np.cos(psi2[idx]) + Y2[idx])

    plt.plot(X_end_bound, Y_end_bound, color="blue")

    
    # End of phase 2
    x1_2 = X1[-1]
    y1_2 = Y1[-1]
    x2_2 = X2[-1]
    y2_2 = Y2[-1]

    # min and max end points
    x1_end = X_end_bound[0]
    x2_end = X_end_bound[-1]
    y1_end = Y_end_bound[0]
    y2_end = Y_end_bound[-1]

    X1_leg = [x1_2, x1_end]
    Y1_leg = [y1_2, y1_end]
    X2_leg = [x2_2, x2_end]
    Y2_leg = [y2_2, y2_end]

    plt.plot(X1_leg,Y1_leg, color="blue")
    plt.plot(X2_leg,Y2_leg, color="blue", label= "Abstract Sim")
    
    # num_points = len(Y_end_bound)
    # fill_region_x = np.linspace(x1_end, x2_end, num_points)
    # fill_region_y1 = np.linspace(y1_2, y1_end, num_points)

    #print(fill_region_x[0],Y_end_bound[0], fill_region_y1[0])

    # # Gazebo Points
    # data = pd.read_csv('fig/points.csv')
    # plt.plot(data['X'], data['Y'],zorder=10,color="orange")

    # Casadi Points
    # data2 = pd.read_csv('fig/casadi_points.csv')
    # plt.plot(-data2['Y'], data2['X'],zorder=10,color="orange", label = "Low Fidelity Sim")
    # plt.legend(fontsize=14,loc='best')
    

    

    #plt.fill_between(fill_region_x, Y_end_bound,fill_region_y1, alpha=0.3, color="blue")

    



    # # Phase 3, get end points after turn. (min and max)
    # x1_2 = X1[-1]
    # y1_2 = Y1[-1]
    # x2_2 = X2[-1]
    # y2_2 = Y2[-1]
    
    # x1_f = 7*math.sin(psi1_2) + x1_2
    # y1_f = -7*math.cos(psi1_2) + y1_2
    # x2_f = 7*math.sin(psi2_2) + x2_2
    # y2_f = -7*math.cos(psi2_2) + y2_2

    # # print(x1_2, y1_2)
    # # print(psi1_2*180/pi)
    # # print(x1_f, y1_f)

    # plt.scatter(x1_f, y1_f)
    # plt.scatter(x2_f, y2_f)

    

    # for idx in range(len(phis1))
    #     x_f = 7*math.sin(phis1[idx]) + X1[idx]
    #     y_f = -7*math.cos(phis1[idx]) + Y1[idx]
    
    # X_end_bound = turn_radius*np.cos(psi_range) + center_point[0]
    # Y_end_bound = turn_radius*np.sin(theta0) + center_point[1]
    #plt.scatter(x2_f, y2_f)
    
        

    # - disturbance

        
        #print(f"x: {x:.2f}, y: {y:.2f}, theta_h: {theta_h*180/pi:.2f}, theta: {math.atan(dy_dx)*180/pi:.2f}, dy_dx: {dy_dx:.3f}")

    # plt.plot(X,Y,zorder = 4,color="blue")
    # plt.plot(X,Y2,zorder = 4, color="blue")
    # plt.fill_between(X, Y, Y2, alpha=0.3, color="blue")
    # #plt.fill_between(X, Y2, 0, alpha=0.3, color="blue")
    plt.scatter(x_0, y_0, color="black",zorder=5)
    plt.annotate('Rover Starting Position', (x_0+3.5, y_0-0.6), textcoords="offset points", xytext=(10,10), ha='center', zorder=10)
    # plt.scatter(x_f, y_f, color="black",zorder=5)
    # plt.annotate('Goal', (x_f-1, y_f-0.6), textcoords="offset points", xytext=(10,10), ha='center', zorder=10)
    ax = plt.subplot(1,1,1)
    
    # # ax.plot(x_0,y_0,'bo')
    ax.set_xlabel('X-position of Rover (m)')
    ax.set_ylabel('Y-position of Rover (m)')
    ax.set_axisbelow(True)
    plt.grid(True)
    # # plt.legend(loc=1)
    plt.axis('equal')
    #plt.tight_layout()
    ax.set_title(f'Possible Trajectories for Rover With a Disturbance of {disturbance:.0f}$\\degree$', fontsize=14)
    plt.savefig("fig/rover_reachable_positions.png")
    plt.close()