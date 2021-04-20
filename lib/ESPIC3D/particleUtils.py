import numpy as np
import math
import scipy.constants


def gamma(velocity):
    VoverC = velocity / scipy.constants.speed_of_light

    return pow(1.0 - np.dot(VoverC, VoverC), -0.5)


def velocityToMomentum(mass, velocity):
    return gamma(velocity) * mass * velocity


def momentumToVelocity(mass, momentum):
    momentumOverMass = momentum / mass
    momentumOverMassOverC = momentumOverMass / scipy.constants.speed_of_light

    return momentumOverMass * pow(1.0 + np.dot(momentumOverMassOverC, momentumOverMassOverC), -0.5)
