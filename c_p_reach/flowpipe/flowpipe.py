import numpy as np
from pytope import Polytope
import casadi as ca
import matplotlib.pyplot as plt
from c_p_reach.lie.se3 import *
from .IntervalHull import qhull2D, minBoundingRect
from .outer_bound import se23_invariant_set_points, se23_invariant_set_points_theta, exp_map

def rotate_point(point, angle):
    new_point = np.array([point[0] * np.cos(angle) - point[1] * np.sin(angle),
                       point[0] * np.sin(angle) + point[1] * np.cos(angle)])
    return new_point 

def flowpipes(ref, n, beta, w1, omegabound, sol, axis):

    x_r = ref['x']
    y_r = ref['y']
    z_r = ref['z']
    
    #####NEED to change this if wants to show different axis#####
    if axis == 'xy':
        nom = np.array([x_r,y_r]).T
    elif axis == 'xz':
        nom = np.array([x_r,z_r]).T

    flowpipes = []
    intervalhull = []
    t_vect = []
    
    step0 = int(len(x_r)/n)
    
    a = 0    
    for i in range(n):
        if i < len(x_r)%n:
            steps = step0 + 1
        else:
            steps = step0
        b = a + steps
        if i == n-1:
            nom_i = nom[a:len(x_r)+1,:]
            if nom_i.shape[0] < 3:
                nom_i = np.vstack((nom_i, np.array(nom[-1,:])+0.01))
        else:
            nom_i = nom[a:b+1,:]
        # Get interval hull
        hull_points = qhull2D(nom_i)
        (rot_angle, area, width, height, center_point, corner_points) = minBoundingRect(hull_points)
            
        t = 0.05*a
        t_vect.append(t)
        
        if t == 0:
            t = 1e-3
        points, val1 = se23_invariant_set_points(sol, t, w1, omegabound, beta) # invariant set at t0 in that time interval
        points2, val2 = se23_invariant_set_points(sol, 0.05*b, w1, omegabound, beta) # invariant set at t final in that time interval
        points_theta, _ = se23_invariant_set_points_theta(sol, t, w1, omegabound, beta)
        
        if val2 > val1: 
            points = points2
            points_theta, _ = se23_invariant_set_points_theta(sol, 0.05*b, w1, omegabound, beta)

        inv_points = exp_map(points, points_theta)

        ######NEED TO change this if want to show other (0 for delete x, 1 for y, 2 for z)#######
        if axis == 'xy':
            inv_points = np.delete(inv_points,2,0) # we want to show x-z, delete y
        elif axis == 'xz': 
            inv_points = np.delete(inv_points,1,0) # we want to show x-z, delete y

        inv_set = [[],[]]
        ang = np.linspace(0, np.pi, 10)
        for theta in ang:
            inv_set1 = rotate_point(inv_points, theta) # it only gives you x and y
            inv_set = np.append(inv_set, inv_set1, axis = 1) 
            
        P2 = Polytope(inv_set.T) 
        
        # minkowski sum
        P1 = Polytope(corner_points) # interval hull
        
        P = P1 + P2 # sum

        p1_vertices = P1.V
        p_vertices = P.V

        p_vertices = np.append(p_vertices, p_vertices[0].reshape(1,2), axis = 0) # add the first point to last, or the flow pipes will miss one line
        
        # create list for flow pipes and interval hull
        flowpipes.append(p_vertices)
        intervalhull.append(P1.V)
        
        a = b
    return flowpipes, intervalhull, nom, t_vect

def plot_flowpipes(nom, flowpipes, n, axis):
    # flow pipes
    plt.figure(figsize=(15,15))
    h_nom = plt.plot(nom[:,0], nom[:,1], color='k', linestyle='-')
    for facet in range(n):
        hs_ch_LMI = plt.plot(flowpipes[facet][:,0], flowpipes[facet][:,1], color='c', linestyle='--')

    # plt.axis('equal')
    plt.title('Flow Pipes')
    if axis == 'xy':
        plt.ylabel('y')
    elif axis == 'xz':
        plt.ylabel('z')
    plt.xlabel('x')
    lgd = plt.legend(loc=2, prop={'size': 18})
    ax = lgd.axes
    handles, labels = ax.get_legend_handles_labels()
    handles.append(h_nom[0])
    labels.append('Reference Trajectory')
    handles.append(hs_ch_LMI[0])
    labels.append('Flow Pipes')
    lgd._legend_box = None
    lgd._init_legend_box(handles, labels)
    lgd._set_loc(lgd._loc)
    lgd.set_title(lgd.get_title().get_text())
