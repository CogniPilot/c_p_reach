import casadi as ca
from .so3 import *
import math

from c_p_reach.lie.util import series_dict
from c_p_reach.lie.matrix_lie_group import MatrixLieGroup
from c_p_reach.lie.so3 import Quat, Euler, Mrp, Dcm


class _SE23(MatrixLieGroup):
    def __init__(self, SO3=None):
        if SO3 == None:
            self.SO3 = Dcm
        else:
            self.SO3 = SO3
        super().__init__(
            group_params= 6 + self.SO3.group_params,
            algebra_params= 6 + self.SO3.algebra_params,
            group_shape=(5, 5),
        )

    def matrix(self, v): # SE3 Lie Group matrix from 6x1 vector
        theta = v[6:9]
        R = Dcm.from_euler(theta)
        p = ca.SX(3, 2)
        p[0, 0] = v[3]
        p[1, 0] = v[4]
        p[2, 0] = v[5]
        p[0, 1] = v[0]
        p[1, 1] = v[1]
        p[2, 1] = v[2]
        horz = ca.horzcat(R, p)
        lastRow1 = ca.SX([0, 0, 0, 1, 0]).T
        lastRow2 = ca.SX([0, 0, 0, 0, 1]).T
        return ca.vertcat(horz, lastRow1, lastRow2)
    
    def vector(self, X): # 9x1 vector from SE3 Lie Group matrix
        v = ca.SX(9,1)
        R = X[0:3, 0:3]
        theta = Euler.from_dcm(R)
        p = X[0:3,4]
        v = X[0:3,3]
        x = ca.vertcat(p,v,theta)
        return x

    def adC_matrix(self):
        adC = ca.SX(9,9)
        adC[0, 3] = 1
        adC[1, 4] = 1
        adC[2, 5] = 1

        return adC


    def ad_matrix(self, v):
        """
        takes 9x1 lie algebra
        input vee operator [x,y,z,vx,vy,vz,theta1,theta2,theta3]
        """
        x = v[0]
        y = v[1]
        z = v[2]
        vx = v[3]
        vy = v[4]
        vz = v[5]
        theta1 = v[6]
        theta2 = v[7]
        theta3 = v[8]

        ad_se3 = ca.SX(9, 9)
        ad_se3[0, 1] = -theta3
        ad_se3[0, 2] = theta2
        ad_se3[1, 0] = theta3
        ad_se3[1, 2] = -theta1
        ad_se3[2, 0] = -theta2
        ad_se3[2, 1] = theta1
        ad_se3[0, 7] = -z
        ad_se3[0, 8] = y
        ad_se3[1, 6] = z
        ad_se3[1, 8] = -x
        ad_se3[2, 6] = -y
        ad_se3[2, 7] = x
        ad_se3[3, 4] = -theta3
        ad_se3[3, 5] = theta2
        ad_se3[4, 3] = theta3
        ad_se3[4, 5] = -theta1
        ad_se3[5, 3] = -theta2
        ad_se3[5, 4] = theta1
        ad_se3[3, 7] = -vz
        ad_se3[3, 8] = vy
        ad_se3[4, 6] = vz
        ad_se3[4, 8] = -vx
        ad_se3[5, 6] = -vy
        ad_se3[5, 7] = vx
        ad_se3[6, 7] = -theta3
        ad_se3[6, 8] = theta2
        ad_se3[7, 6] = theta3
        ad_se3[7, 8] = -theta1
        ad_se3[8, 6] = -theta2
        ad_se3[8, 7] = theta1
        return ad_se3

    def vee(self, X):
        """
        This takes in an element of the SE3 Lie Group (Wedge Form) and returns the se3 Lie Algebra elements
        """
        v = ca.SX(9, 1)
        v[0, 0] = X[0, 4]  # x
        v[1, 0] = X[1, 4]  # y
        v[2, 0] = X[2, 4]  # z
        v[3, 0] = X[0, 3]  # vx
        v[4, 0] = X[1, 3]  # vy
        v[5, 0] = X[2, 3]  # vz
        v[6, 0] = X[2, 1]  # theta0
        v[7, 0] = X[0, 2]  # theta1
        v[8, 0] = X[1, 0]  # theta2
        return v

    def wedge(self, v):
        """
        This takes in an element of the se3 Lie Algebra and returns the se3 Lie Algebra matrix

        v: [x,y,z,theta0,theta1,theta2]
        """
        X = ca.SX.zeros(5, 5)
        X[0, 4] = v[0]
        X[1, 4] = v[1]
        X[2, 4] = v[2]
        X[0, 3] = v[3]
        X[1, 3] = v[4]
        X[2, 3] = v[5]
        X[:3, :3] = Dcm.wedge(v[6:9])
        return X

    def exp(self, v):  # accept input in wedge operator form
        v = self.vee(v)
        # v = [x,y,z,theta1,theta2,theta3]
        v_so3 = v[
            6:9
        ]  # grab only rotation terms for so3 uses ##corrected to v_so3 = v[3:6]
        X_so3 = Dcm.wedge(v_so3)  # wedge operator for so3
        theta = ca.norm_2(
            Dcm.vee(X_so3)
        )  # theta term using norm for sqrt(theta1**2+theta2**2+theta3**2)

        # translational components u
        u = ca.SX(3, 2)
        u[0, 0] = v[3]
        u[1, 0] = v[4]
        u[2, 0] = v[5]
        u[0, 1] = v[0]
        u[1, 1] = v[1]
        u[2, 1] = v[2]

        R = Dcm.exp(
            v_so3
        )  #'Dcm' for direction cosine matrix representation of so3 LieGroup Rotational

        A = series_dict["sin(x)/x"](theta)
        B = series_dict["(1 - cos(x))/x^2"](theta)
        C = (1 - A) / theta**2

        V = ca.SX.eye(3) + B * X_so3 + C * X_so3 @ X_so3

        horz = ca.horzcat(R, ca.mtimes(V, u))

        lastRow1 = ca.SX([0, 0, 0, 1, 0]).T
        lastRow2 = ca.SX([0, 0, 0, 0, 1]).T

        return ca.vertcat(horz, lastRow1, lastRow2)

    def identity(self):
        return ca.SX.eye(5)

    def product(self, a, b):
        self.check_group_shape(a)
        self.check_group_shape(b)
        return a @ b

    def inv(self, a):  # input a matrix of SX form from casadi
        Rt = a[0:3, 0:3].T
        vp = a[0:3, 3:5]
        horz = ca.horzcat(Rt, -ca.mtimes(Rt, vp))

        lastRow1 = ca.SX([0, 0, 0, 1, 0]).T
        lastRow2 = ca.SX([0, 0, 0, 0, 1]).T
        
        return ca.vertcat(horz, lastRow1, lastRow2)

    def log(self, G):
        R = G[:3, :3]
        theta = ca.arccos((ca.trace(R) - 1) / 2)
        wSkew = Dcm.wedge(Dcm.log(R))
        A = series_dict["sin(x)/x"](theta)
        B = series_dict["(1 - cos(x))/x^2"](theta)
        V_inv = (
            ca.SX.eye(3)
            - wSkew / 2
            + (1 / theta**2) * (1 - A / (2 * B)) * wSkew @ wSkew
        )

        t = ca.SX(3, 2)
        t[0, :] = G[0, 3:5]
        t[1, :] = G[1, 3:5]
        t[2, :] = G[2, 3:5]

        uInv = V_inv @ t
        horz2 = ca.horzcat(wSkew, uInv)
        lastRow2 = ca.SX([0, 0, 0, 0, 0]).T
        return ca.vertcat(horz2, lastRow2, lastRow2)

    def diff_correction(self, v):  # U Matrix for se3 with input vee operator
        return ca.inv(self.diff_correction_inv(v))

    def diff_correction_inv(self, v):  # U_inv of se3 input vee operator

        ad = self.ad_matrix(v)
        ad_zeta_k = ca.SX_eye(9)
        u_inv = ca.SX.eye(9)

        for k in range(1, 10):
            ad_zeta_k = ad_zeta_k@ad
            u_inv += ad_zeta_k/math.factorial(k+1)

        return u_inv

        # u_inv = ca.SX(6, 6)
        # u1 = c2*(-v[4]**2 - v[5]**2) + 1
        # u2 = -c1*v[5] + c2*v[3]*v[4]
        # u3 = c1 * v[4] + c2 * v[3]*v[5]
        # u4 = c2 * (-2*v[4]*v[1]-2*v[5]*v[2])
        # u5 = -c1 * v[2] + c2*(v[4]*v[0]+v[3]*v[1])
        # u6 = c1 * v[1] + c2*(v[3]*v[2]+v[5]*v[0])
        # uInvR1 = ca.vertcat(u1,u2,u3,u4,u5,u6)

        # u1 = c1 * v[5] + c2 * v[3] * v[4]
        # u2 = c2 *(-v[3]**2 - v[5]**2)+1
        # u3 = -c1*v[3] + c2 * v[4]*v[5]
        # u4 = c1 * v[2] + c2 * (v[3]*v[1]+v[4]*v[0])
        # u5 = c2* (-2*v[3] * v[0] -2*v[5]*v[2])
        # u6 = -c1 * v[0] + c2 * (v[4]*v[2] + v[5] *v[1])
        # uInvR2 = ca.vertcat(u1,u2,u3,u4,u5,u6)

        # u1 = -c1 * v[4] + c2 * v[3] * v[5]
        # u2 = c1 * v[3] + c2 * v[4] * v[5]
        # u3 = c1 * (-v[3] **2  - v[4]**2) +1
        # u4 = -c1 * v[1] + c2 * (v[3]*v[2] + v[5]*v[0])
        # u5 = c1 * v[0] + c2 * (v[4]*v[2] + v[5] *v[1])
        # u6 = c2 * (-2*v[3]*v[0] - 2*v[4] *v[1])
        # uInvR3 = ca.vertcat(u1,u2,u3,u4,u5,u6)

        # u1 = 0
        # u2 = 0
        # u3 = 0
        # u4 = c2 * (- v[4]**2 - v[5]**2) +1
        # u5 = -c1*v[5] + c2*v[3]*v[4]
        # u6 = c1 * v[4] + c2 * v[3] * v[5]
        # uInvR4 = ca.vertcat(u1,u2,u3,u4,u5,u6)

        # u1 = 0
        # u2 = 0
        # u3 = 0
        # u4 = c1 * v[5] + c2 * v[3] * v[4]
        # u5 = c2 * (-v[3]**2 - v[5]**2) +1
        # u6 = -c1 * v[3] + c2 * v[4] *v[5]
        # uInvR5 = ca.vertcat(u1,u2,u3,u4,u5,u6)

        # u1 = 0
        # u2 = 0
        # u3 = 0
        # u4 = -c1 * v[4] + c2 * v[3] * v[5]
        # u5 = c1 * v[3] + c2 * v[4] * v[5]
        # u6 = c2 * (-v[3] **2 - v[4]**2)+1
        # uInvR6 = ca.vertcat(u1,u2,u3,u4,u5,u6)

        # u_inv = ca.transpose(ca.horzcat(uInvR1,uInvR2,uInvR3,uInvR4, uInvR5, uInvR6))
        # return u_inv

    # verify this with series solution

    # https://github.com/jgoppert/pyecca/blob/master/pyecce/estimation/attitude/algorithms/mrp.py
    # Use this to try to get casadi to draw a plot for this
    # line 112-116 help for drawing plots

    # https://github.com/jgoppert/pyecca/blob/master/pyecca/estimation/attitude/algorithms/common.py
    # This import to import needed casadi command

    # New notes (Oct 18 22)
    ## sympy.cse(f) to find common self expression using sympy to clean up the casadi plot
    # cse_def, cse_expr = sympy.cse(f)

    # Oct 25 Update
    # work on updating SE2 umatrix to casadi
    # get SE2 U matrix
    # u matrix can be found through casadi using inverse function for casadi


SE23Dcm = _SE23(Dcm)
SE23Euler = _SE23(Euler)
SE23Quat = _SE23(Quat)
SE23Mrp = _SE23(Mrp)
