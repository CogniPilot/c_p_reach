
import casadi as ca

class Integrator:

    def __init__(self):

        # declare states
        x = ca.SX.sym('x');
        y = ca.SX.sym('y');

        # declare state vector
        self.x = ca.vertcat(
            x,
            y);
        
        # declare state derivative equations
        der_x = 1;
        der_y = x;

        # declare state derivative vector
        self.x_dot = ca.vertcat(
            der_x,
            der_y);
        self.ode = ca.Function('ode', [self.x], [self.x_dot])



