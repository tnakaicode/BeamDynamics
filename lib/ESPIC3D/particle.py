import numpy as np
import math
from solvers.particleUtils import velocityToMomentum, momentumToVelocity
#import scipy.constants


class particle(object):
    def __init__(self, *args):
        self.mass = args[0]
        self.charge = args[1]
        self.position = np.array(args[2])
        self.velocity = np.array(args[3])
        self.chargeOverMass = self.charge / self.mass

    # for a test, do just magnetic field and check circle

    # Boris push
    # need to use momentum for relativistic
    # google relativistic boris push

    # E and B are supposed to be at particle position
    def push(self, dt, E, B=None):
        if B:
            t = 0.5 * dt * self.chargeOverMass * B
            s = 2.0 * t / (1.0 + np.dot(t, t))

            V_minus = self.velocity + 0.5 * dt * self.chargeOverMass * E
            V_prime = V_minus + np.cross(V_minus, t)
            V_plus = V_minus + np.cross(V_prime, s)
            self.velocity = V_plus + 0.5 * dt * self.chargeOverMass * E
        else:
            self.velocity = self.velocity + dt * self.chargeOverMass * E

        self.position = self.position + dt * self.velocity

    def push_relativistic(self, dt, E, B=None):
        if B:
            blah = 0
            #t             = 0.5*dt*self.chargeOverMass*B
            #s             = 2.0*t/(1.0 + np.dot(t,t))

            #V_minus       = self.velocity + 0.5*dt*self.chargeOverMass*E
            #V_prime       = V_minus       + np.cross(V_minus,t)
            #V_plus        = V_minus       + np.cross(V_prime,s)
            #self.velocity = V_plus        + 0.5*dt*self.chargeOverMass*E
        else:
            momentum = velocityToMomentum(self.mass, self.velocity)
            momentum = momentum + dt * self.charge * E
            self.velocity = momentumToVelocity(self.mass, momentum)

        self.position = self.position + dt * self.velocity
