# output
#   columns (in order): time, CG height, pitch, roll,
#                       wheel loads (FL,FR,RL,RR), damper travel (same
#                       order)
#   time increases down rows

# data_in
# columns (in order): Time [s], Ax [G], Ay [G], Az [G], r1 [in], r2 [in]
#                     r3 [in], r4 [in]

# THE TIME STEP MUST BE CONSTANT!

# Linearized Transient Suspension Analysis
# By Thomas Bennett
# 10/10/2020

# Key Assumptions:
# - Linear suspension behavior (Springs, Dampers, IR's, etc)
# - Neglegible movement of CG's (sprung and unsprung)
# - No "semi-sprung" modeling
# - Neglegible roll center movement
# - Sprung mass MOI estimated using principal MOI's for a hollow cylinder
# - Damping coefficients estimated (coefficients can be changed or derived
#       from real setups later)
# - Equal weight distribution (also may be corrected later)

import math
import numpy as np
from scipy import integrate
import matplotlib.pyplot as plt
import time

start_time = time.time()

data_in = np.genfromtxt('Miata_run.csv', delimiter=',')
# hi

# Vehicle Properties

class Weight_Struct(object):
    def __init__(self, distribution=.0, vehicle=.0, driver=.0, unsprung_LF=.0, unsprung_RF=.0, unsprung_LR=.0,
                 unsprung_RR=.0):
        self.distribution = distribution
        self.vehicle = vehicle
        self.driver = driver
        self.unsprung_LF = unsprung_LF
        self.unsprung_RF = unsprung_RF
        self.unsprung_LR = unsprung_LR
        self.unsprung_RR = unsprung_RR
        self.vehicle_total = self.vehicle + self.driver


W = Weight_Struct

# Sprung Mass:  Pg 671 / 354 pdf
W.distribution = 0.54  # Weight distribution longitudinal
W.vehicle = 2167.1  # Weigh of Vehicle (lbs) (Unless specified otherwise)
W.driver = 0  # Weight of driver (lbs) (Unless specified otherwise)
W.vehicle_total = W.vehicle + W.driver  # Combined weight of vehicle and driver (lbs)
W.unsprung_LF = 59.525  # Left Front unsprung weight Pg 671 book / Pg 354 pdf Verify from 2019 (Wheel, upright,
# caliper, rotor, hub, etc and half of control arm,bell crank, pull rod)
W.unsprung_RF = 59.525  # Right Front unsprung weight Pg 671 book / Pg 354 pdf
W.unsprung_LR = 54.013  # Left Rear unsprung weight Pg 671 book / Pg 354 pdf Verify from 2019 (Wheel, upright,
# caliper, rotor, hub, tulip, etc and half of control arm,bell crank, pull rod)
W.unsprung_RR = 54.013  # Right Rear unsprung weight Pg 671 book / Pg 354 pdf
W.unsprung_F = W.unsprung_LF + W.unsprung_RF  # W_UF - Unsprung mass front of car
W.unsprung_R = W.unsprung_LR + W.unsprung_RR  # W_UR - Unsprung mass rear of car
W.sprung = W.vehicle_total - (W.unsprung_F + W.unsprung_R)  # Sprung mass of vehicle Pg 671 book / Pg 354 pdf
W.F_axle = W.vehicle_total * W.distribution  # Front axle weight (lbs) - with driver
W.R_axle = W.vehicle_total * (1 - W.distribution)  # Rear axle weight (lbs) - with driver

# NEW TO V7 #
W.spr_dis = (W.vehicle_total * W.distribution - W.unsprung_F) / W.sprung  # Long. weight distribution of sprung mass


class Dimensions:
    def __init__(self, track_F=0, track_R=0, wheelbase=0, IR_front=0, IR_rear=0, IR_change_rate_front=0,
                 IR_change_rate_rear=0):
        self.track_F = track_F
        self.track_R = track_R
        self.wheelbase = wheelbase
        self.Installation_Ratio_Front = IR_front
        self.Installation_Ratio_Rear = IR_rear
        self.Installation_Ratio_Change_Rate_Front = IR_change_rate_front
        self.Installation_Ratio_Change_Rate_Rear = IR_change_rate_rear


D = Dimensions

D.track_F = 1.525 / .0254 / 12  # Front Track - Convert inches to feet (ft)
D.track_R = 1.525 / .0254 / 12  # Rear Track - Convert inches to feet (ft)
D.wheelbase = 2.310 / .0254 / 12  # CAD Wheelbase - Convert inches to feet (ft)
D.a = D.wheelbase * W.distribution  # Distance from Front Axle to CG            - NOTE: May change how i calculate
# this....
D.b = D.wheelbase * (1 - W.distribution)  # Distance from Rear Axle to CG
D.b_s = ((W.vehicle_total * D.b) - (W.unsprung_F * D.wheelbase)) / W.sprung  # Longitudinal location of sprung mass
# CG Pg 673 / Pg 355
D.a_s = D.wheelbase - D.b_s  # Distance from front axle to Sprung mass CG location Pg 673 / Pg 355
D.Installation_Ratio_Front = 0.904  # 0.8913
D.Installation_Ratio_Rear = 0.8671  # 0.8543        0.87
D.Installation_Ratio_Change_Rate_Front = 0  # Look into  - From CJ 0
D.Installation_Ratio_Change_Rate_Rear = 0  # Look into  - From CJ 0


class Height_Struct:
    def __init__(self, CG=0, unsprung_F=0, unsprung_R=0, RC_front=0, RC_Rear=0, ride_height_F=0, ride_height_R=0):
        self.CG = CG
        self.unsprung_F = unsprung_F
        self.unsprung_R = unsprung_R
        self.RC_Front = RC_front
        self.RC_Rear = RC_Rear
        self.Ride_Height_F = ride_height_F
        self.Ride_Height_R = ride_height_R


H = Height_Struct  # look over unsprung cg height and sprung....

H.CG = 16.5 / 12  # Center gravity - Converted to ft - (NOTE: Need to re-measure our CG)
H.unsprung_F = 6 / 12  # z_WF - Front Height of unsprung mass (ft) - (Assumed at front wheel center) Pg 672 book /
# Pg 355 pdf
H.unsprung_R = 6 / 12  # z_WR - Rear Height of unsprung mass (ft) - (Assumed at  wheel center) Pg 672 book / Pg
# 355 pdf
H.RC_Front = 0 / 12  # Converted to (ft) Front Roll Center height z_RF Note: Make this into an array to analysis a
# variety of different roll centers???
H.RC_Rear = 0 / 12  # Converted to (ft) Rear Roll Center height z_RR Note: Make this into an array to analysis a
# variety of different roll centers???

# Note: Below I moved this into the for loop because it need to iterate and calculate the new heights
H.Roll_Axis_Inclination = math.asin((H.RC_Rear - H.RC_Front) / D.wheelbase)  # (Radians) - Roll Axis angle
H.Roll_Axis_Inclination_Degree = H.Roll_Axis_Inclination * (180 / math.pi)  # In Degrees - Roll Axis angle
H.sprung = ((W.vehicle_total / W.sprung) * H.CG) - ((W.unsprung_F / W.sprung) * H.unsprung_F) - (
        (W.unsprung_R / W.sprung) * H.unsprung_R)  # (ft) Height of sprung mass from ground (h_s)

H.SprungToRollAxis = (H.sprung - (H.RC_Front + (H.RC_Rear - H.RC_Front) * D.a_s / D.wheelbase)) * math.cos(
    math.atan((H.RC_Rear - H.RC_Front) / D.wheelbase))
H.RollCouple = H.SprungToRollAxis  # h2 - Sprung weight (lb) to perpendicular distance to NRA (ft)

H.Ride_Height_F = 1  # Front ride height (in)
H.Ride_Height_R = 1.045  # Rear ride height (in)

# Tire Properties

K_SidewallHeight = 3.0 / 12  # ft  Verified from actual @ 10.4 psi
K_UnloadedRadius = 10.2 / 12  # ft  Verified from actual @ 10.4 psi
K_Pressure = 12  # 11.9psi  Chosen because Max AligningTorque per CorneringStiffness from tire data)
K_Tire = 13869.5  # lb/ft  Post break-in dynamic spring rate @ 12 psi from tire data)
K_Deflection = (W.vehicle_total / 4) / K_Tire  # in  Assuming equal weight distribution
StaticLoadWheelCenterheight = K_UnloadedRadius * 12 - K_Deflection  # in  Assuming equal weight distribution
Tire_Damping = 34.26088  # [lbf-s/ft] Assuming same damping coefficient for each tire

#######################################################################################################################
# Code unique to the transient analysis begins here

G = 32.185  # ft/s/s, Gravitational Acceleration


class Inertia:
    def __init__(self, xx=0, xxRC=0, yy=0):
        self.xx = xx
        self.xx_RC = xxRC
        self.yy = yy


I = Inertia

cyl_do = 1128 / 25.4 / 12  # ft
cyl_di = cyl_do / 2 / 12  # ft
cyl_h = 2917 / 25.4 * (5 / 6) / 12  # ft
I.xx = (cyl_do ** 2 + cyl_di ** 2) / 8 * W.sprung / 32.185  # slug*ft^2
I.xx_RC = I.xx + W.sprung / 32.185 * H.SprungToRollAxis ** 2  # slug*ft^2
I.yy = (3 * cyl_do ** 2 + 3 * cyl_di ** 2 + 4 * cyl_h ** 2) / 48 * W.sprung / 32.185  # slug*ft^2

# Inertia from CAD

I.xx = W.sprung / G * (5.2493 ** 2 + 3.08399 ** 2) / 12  # slug*ft^2
I.xx_RC = I.xx + W.sprung / G * H.SprungToRollAxis ** 2  # slug*ft^2
I.yy = W.sprung / G * (12.1063 ** 2 + 3.08399 ** 2) / 12  # slug*ft^2

# Suspension Parameters
K_ARB_F = 28764.92  # lbf*ft/rad
K_ARB_R = 3540.2983  # lbf*ft/rad
SpringRateF = 291.22  # lbf/in
SpringRateR = 439.68  # lbf/in
DampRateLF = 37.115  # lbf-s/in 88.8 maximum from Ohlins
DampRateRF = 37.115  # lbf-s/in 88.8 maximum from Ohlins
DampRateLR = 39.971  # lbf-s/in 88.8 maximum from Ohlins
DampRateRR = 39.971  # lbf-s/in 88.8 maximum from Ohlins


def dampFuncF(vel):
    if math.fabs(vel) <= 1 / 12:
        return 37.115
    else:
        return 15


def dampFuncR(vel):
    if math.fabs(vel) <= 1 / 12:
        return 39.971
    else:
        return 16


IR_F = 1
IR_R = 1
wheelRateF = SpringRateF * IR_F ** 2 * 12  # lbf/ft
wheelRateR = SpringRateR * IR_R ** 2 * 12  # lbf/ft
wheelDampLF = DampRateLF * IR_F ** 2 * 12  # lbf*s/ft
wheelDampRF = DampRateRF * IR_F ** 2 * 12  # lbf*s/ft
wheelDampLR = DampRateLR * IR_R ** 2 * 12  # lbf*s/ft
wheelDampRR = DampRateRR * IR_R ** 2 * 12  # lbf*s/ft
# Is this legit?
wheelDampF = (wheelDampLF + wheelDampRF) / 2
wheelDampR = (wheelDampLR + wheelDampRR) / 2

# The system will be simulated using the linear state variable method. The
# system will be put into the form "dx/dt = A*x+B*u" where "x" is the state
# vector and "u" is the input vector. "A" and "B" are constant matrices.

# The entries of the "A" and "B" matrices are calculated below.

# Pitch Dynamics Coefficients
a_p = 2 * D.wheelbase * (W.spr_dis * wheelRateR - (1 - W.spr_dis) * wheelRateF) / I.yy
b_p = 2 * D.wheelbase * (W.spr_dis * wheelDampR - (1 - W.spr_dis) * wheelDampF) / I.yy
c_p = -2 * (D.wheelbase ** 2) * (((1 - W.spr_dis) ** 2) * wheelRateF + (W.spr_dis ** 2) * wheelRateR) / I.yy
d_p = -2 * (D.wheelbase ** 2) * (((1 - W.spr_dis) ** 2) * wheelDampF + (W.spr_dis ** 2) * wheelDampR) / I.yy
e_p = D.wheelbase * (1 - W.spr_dis) * wheelRateF / I.yy
f_p = D.wheelbase * (1 - W.spr_dis) * wheelDampF / I.yy
g_p = -D.wheelbase * W.spr_dis * wheelRateR / I.yy
h_p = -D.wheelbase * W.spr_dis * wheelDampR / I.yy
Ax_p = W.sprung * H.CG / I.yy

# Roll Dynamics Coefficients
a_r = -((D.track_F ** 2 * wheelRateF + D.track_R ** 2 * wheelRateR) / 2 + K_ARB_F + K_ARB_R) / I.xx_RC
b_r = -(D.track_F ** 2 * wheelDampF + D.track_R ** 2 * wheelDampR) / 2 / I.xx_RC
c_r = (D.track_F * wheelRateF / 2 + K_ARB_F / D.track_F) / I.xx_RC
d_r = D.track_F * wheelDampF / 2 / I.xx_RC
e_r = (D.track_R * wheelRateR / 2 + K_ARB_R / D.track_R) / I.xx_RC
f_r = D.track_R * wheelDampR / 2 / I.xx_RC
Ay_r = -W.sprung * H.SprungToRollAxis / I.xx_RC

# Height Dynamics Coefficients
a_z = -2 * (wheelRateF + wheelRateR) / (W.sprung / G)
b_z = -2 * (wheelDampF + wheelDampR) / (W.sprung / G)
c_z = 2 * D.wheelbase * (W.spr_dis * wheelRateR - (1 - W.spr_dis) * wheelRateF) / (W.sprung / G)
d_z = 2 * D.wheelbase * (W.spr_dis * wheelDampR - (1 - W.spr_dis) * wheelDampF) / (W.sprung / G)
e_z = wheelRateF / (W.sprung / G)
f_z = wheelDampF / (W.sprung / G)
g_z = wheelRateR / (W.sprung / G)
h_z = wheelDampR / (W.sprung / G)
Az_z = -W.sprung / (W.sprung / G)

# Wheel 1 Dynamics Coefficients
a_u1 = wheelRateF / (W.unsprung_LF / G)
b_u1 = wheelDampLF / (W.unsprung_LF / G)
c_u1 = wheelRateF / (W.unsprung_LF / G) * D.wheelbase * (1 - W.spr_dis)
d_u1 = wheelDampLF / (W.unsprung_LF / G) * D.wheelbase * (1 - W.spr_dis)
e_u1 = (wheelRateF * D.track_F / 2 + K_ARB_F / D.track_F) / (W.unsprung_LF / G)
f_u1 = wheelDampLF / (W.unsprung_LF / G) * D.track_F / 2
g_u1 = -(wheelRateF + K_ARB_F / D.track_F ** 2 + K_Tire) / (W.unsprung_LF / G)
h_u1 = -(wheelDampLF + Tire_Damping) / (W.unsprung_LF / G)
k_u1 = K_ARB_F / D.track_F ** 2 / (W.unsprung_LF / G)
r_u1 = K_Tire / (W.unsprung_LF / G)
drdt_u1 = Tire_Damping / (W.unsprung_LF / G)
Ay_u1 = -(W.sprung * H.RC_Front / 2 + W.unsprung_F * H.unsprung_F) / D.track_F / (W.unsprung_LF / G)
Az_u1 = -W.unsprung_LF / (W.unsprung_LF / G)

# Wheel 2 Dynamics Coefficients
a_u2 = wheelRateF / (W.unsprung_RF / G)
b_u2 = wheelDampRF / (W.unsprung_RF / G)
c_u2 = wheelRateF / (W.unsprung_RF / G) * D.wheelbase * (1 - W.spr_dis)
d_u2 = wheelDampRF / (W.unsprung_RF / G) * D.wheelbase * (1 - W.spr_dis)
e_u2 = -(wheelRateF * D.track_F / 2 + K_ARB_F / D.track_F) / (W.unsprung_RF / G)
f_u2 = -wheelDampRF / (W.unsprung_RF / G) * D.track_F / 2
g_u2 = K_ARB_F / D.track_F ** 2 / (W.unsprung_RF / G)
h_u2 = -(wheelRateF + K_ARB_F / D.track_F ** 2 + K_Tire) / (W.unsprung_RF / G)
k_u2 = -(wheelDampRF + Tire_Damping) / (W.unsprung_RF / G)
r_u2 = K_Tire / (W.unsprung_RF / G)
drdt_u2 = Tire_Damping / (W.unsprung_RF / G)
Ay_u2 = (W.sprung * H.RC_Front / 2 + W.unsprung_F * H.unsprung_F) / D.track_F / (W.unsprung_RF / G)
Az_u2 = -W.unsprung_RF / (W.unsprung_RF / G)

# Wheel 3 Dynamics Coefficients
a_u3 = wheelRateR / (W.unsprung_LR / G)
b_u3 = wheelDampLR / (W.unsprung_LR / G)
c_u3 = -wheelRateR / (W.unsprung_LR / G) * D.wheelbase * W.spr_dis
d_u3 = -wheelDampLR / (W.unsprung_LR / G) * D.wheelbase * W.spr_dis
e_u3 = (wheelRateR * D.track_R / 2 + K_ARB_R / D.track_R) / (W.unsprung_LR / G)
f_u3 = wheelDampLR / (W.unsprung_LR / G) * D.track_R / 2
g_u3 = -(wheelRateR + K_ARB_R / D.track_R ** 2 + K_Tire) / (W.unsprung_LR / G)
h_u3 = -(wheelDampLR + Tire_Damping) / (W.unsprung_LR / G)
k_u3 = K_ARB_R / D.track_R ** 2 / (W.unsprung_LR / G)
r_u3 = K_Tire / (W.unsprung_LR / G)
drdt_u3 = Tire_Damping / (W.unsprung_LR / G)
Ay_u3 = -(W.sprung * H.RC_Rear / 2 + W.unsprung_R * H.unsprung_R) / D.track_R / (W.unsprung_LR / G)
Az_u3 = -W.unsprung_LR / (W.unsprung_LR / G)

# Wheel 4 Dynamics Coefficients
a_u4 = wheelRateR / (W.unsprung_RR / G)
b_u4 = wheelDampRR / (W.unsprung_RR / G)
c_u4 = -wheelRateR / (W.unsprung_RR / G) * D.wheelbase * W.spr_dis
d_u4 = -wheelDampRR / (W.unsprung_RR / G) * D.wheelbase * W.spr_dis
e_u4 = -(wheelRateR * D.track_R / 2 + K_ARB_R / D.track_R) / (W.unsprung_RR / G)
f_u4 = -wheelDampRR / (W.unsprung_RR / G) * D.track_R / 2
g_u4 = K_ARB_R / D.track_R ** 2 / (W.unsprung_RR / G)
h_u4 = -(wheelRateR + K_ARB_R / D.track_R ** 2 + K_Tire) / (W.unsprung_RR / G)
k_u4 = -(wheelDampRR + Tire_Damping) / (W.unsprung_RR / G)
r_u4 = K_Tire / (W.unsprung_RR / G)
drdt_u4 = Tire_Damping / (W.unsprung_RR / G)
Ay_u4 = (W.sprung * H.RC_Rear / 2 + W.unsprung_R * H.unsprung_R) / D.track_R / (W.unsprung_RR / G)
Az_u4 = -W.unsprung_RR / (W.unsprung_RR / G)

# Definition of the "A" matrix
A = np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [a_z, b_z, c_z, d_z, 0, 0, e_z, f_z, e_z, f_z, g_z, h_z, g_z, h_z],
              [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
              [a_p, b_p, c_p, d_p, 0, 0, e_p, f_p, e_p, f_p, g_p, h_p, g_p, h_p],
              [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, a_r, b_r, c_r, d_r, -c_r, -d_r, e_r, f_r, -e_r, -f_r],
              [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
              [a_u1, b_u1, c_u1, d_u1, e_u1, f_u1, g_u1, h_u1, k_u1, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
              [a_u2, b_u2, c_u2, d_u2, e_u2, f_u2, g_u2, 0, h_u2, k_u2, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
              [a_u3, b_u3, c_u3, d_u3, e_u3, f_u3, 0, 0, 0, 0, g_u3, h_u3, k_u3, 0],
              [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
              [a_u4, b_u4, c_u4, d_u4, e_u4, f_u4, 0, 0, 0, 0, g_u4, 0, h_u4, k_u4]])


def AUpdate(velLF, velRF, velLR, velRR):
    newLF = dampFuncF(velLF)
    newRF = dampFuncF(velRF)
    newLR = dampFuncR(velLR)
    newRR = dampFuncR(velRR)

    new_wheelDampLF = newLF * IR_F ** 2 * 12
    new_b_u1 = new_wheelDampLF / (W.unsprung_LF / G)
    new_d_u1 = new_wheelDampLF / (W.unsprung_LF / G) * D.wheelbase * (1 - W.spr_dis)
    new_f_u1 = new_wheelDampLF / (W.unsprung_LF / G) * D.track_F / 2
    new_h_u1 = -(new_wheelDampLF + Tire_Damping) / (W.unsprung_LF / G)

    new_wheelDampRF = newRF * IR_F ** 2 * 12
    new_b_u2 = new_wheelDampRF / (W.unsprung_RF / G)
    new_d_u2 = new_wheelDampRF / (W.unsprung_RF / G) * D.wheelbase * (1 - W.spr_dis)
    new_f_u2 = - new_wheelDampRF / (W.unsprung_RF / G) * D.track_F / 2
    new_k_u2 = -(new_wheelDampRF + Tire_Damping) / (W.unsprung_RF / G)

    new_wheelDampLR = newLR * IR_R ** 2 * 12
    new_b_u3 = new_wheelDampLR / (W.unsprung_LR / G)
    new_d_u3 = - new_wheelDampLR / (W.unsprung_LR / G) * D.wheelbase * W.spr_dis
    new_f_u3 = new_wheelDampLR / (W.unsprung_LR / G) * D.track_R / 2
    new_h_u3 = -(new_wheelDampLR + Tire_Damping) / (W.unsprung_LR / G)

    new_wheelDampRR = newRR * IR_R ** 2 * 12
    new_b_u4 = new_wheelDampRR / (W.unsprung_RR / G)
    new_d_u4 = - new_wheelDampRR / (W.unsprung_RR / G) * D.wheelbase * W.spr_dis
    new_f_u4 = - new_wheelDampRR / (W.unsprung_RR / G) * D.track_R / 2
    new_k_u4 = -(new_wheelDampRR + Tire_Damping) / (W.unsprung_RR / G)

    return np.array([[0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [a_z, b_z, c_z, d_z, 0, 0, e_z, f_z, e_z, f_z, g_z, h_z, g_z, h_z],
                     [0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
                     [a_p, b_p, c_p, d_p, 0, 0, e_p, f_p, e_p, f_p, g_p, h_p, g_p, h_p],
                     [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, a_r, b_r, c_r, d_r, -c_r, -d_r, e_r, f_r, -e_r, -f_r],
                     [0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0],
                     [a_u1, new_b_u1, c_u1, new_d_u1, e_u1, new_f_u1, g_u1, new_h_u1, k_u1, 0, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0],
                     [a_u2, new_b_u2, c_u2, new_d_u2, e_u2, new_f_u2, g_u2, 0, h_u2, new_k_u2, 0, 0, 0, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0],
                     [a_u3, new_b_u3, c_u3, new_d_u3, e_u3, new_f_u3, 0, 0, 0, 0, g_u3, new_h_u3, k_u3, 0],
                     [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1],
                     [a_u4, new_b_u4, c_u4, new_d_u4, e_u4, new_f_u4, 0, 0, 0, 0, g_u4, 0, h_u4, new_k_u4]])


# Definition of the "B" matrix

B = np.array([[0, 0, 0, 0, 0, 0, 0],
              [0, 0, Az_z, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [Ax_p, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, Ay_r, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, Ay_u1, Az_u1, r_u1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, Ay_u2, Az_u2, 0, r_u2, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, Ay_u3, Az_u3, 0, 0, r_u3, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, Ay_u4, Az_u4, 0, 0, 0, r_u4]])

E = np.array([[0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, drdt_u1, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, drdt_u2, 0, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, drdt_u3, 0],
              [0, 0, 0, 0, 0, 0, 0],
              [0, 0, 0, 0, 0, 0, drdt_u4]])

B_bar = np.dot(A, E) + B  # Effective "B" Matrix. See Linearized Dynamic Suspension Model page 6.

# Simulation setup
# tspan = np.linspace(0, data_in[0, 0], num_steps)  # Time definition; ensures constant time step
tspan = data_in[1:, 6]  # Alternate time vector definition.
tstep = data_in[1:, 6][1] - data_in[1:, 6][0]
x0 = np.zeros(14)  # Initial state vector values (set to zero)

# System inputs [Ax Ay Az r1 r2 r3 r4]
# Ax, Ay, Az - acceleration in G's
# r1, r2, r3, r4 - ground height under wheels in ft
#   1 = FL, 2 = FR, 3 = RL, 4 = RR

# BAND AID FIX OF INPUT ###
Accel = data_in[1:, 0:3]
Road = np.zeros((len(data_in[1:, 6]), 4))
sim_input = np.concatenate((Accel, Road), axis=1)


###


def t_round(t, input_step):
    return round(int(t / input_step) * input_step, 2)


setup_time = time.time()
print('Setup time:', setup_time - start_time)


# Simulate the system
def dX_dt(X, t):
    global A
    vect = A.dot(X)[:, np.newaxis] + B.dot(sim_input[np.where(data_in[1:, 6] == t_round(t, 0.05))].T)
    A = AUpdate(X[7], X[9], X[11], X[13])

    return vect.T.flatten()


x = integrate.odeint(dX_dt, x0.T, tspan, rtol=0.001, atol=0.001)
int_time = time.time()
print('Integration time:', int_time - setup_time)

'''
#Backwards Euler Solver ###################################################
A_bar = eye(length(x0))-A.*del_t
A_bar_inv = A_bar^-1

x = zeros(length(x0),length(tspan))
x(:,1) = x0 - E*input(:,1);

for i = 2:length(tspan)
    
    mean_input = input(:,i);
    
    x(:,i) = A_bar_inv*(x(:,i-1) + (del_t.*B_bar)*mean_input);
    
end
###########################################################################


# 4th Order Runge-Kutta ####################################################
# 
# x = zeros(length(x0),length(tspan))
# x(:,1) = x0 - E*input(:,1);
# 
# for i = 2:length(tspan)
#     
#     current_input = input(:,i-1);
#     k1 = A*x(:,i-1) + B_bar*current_input
#     
#     current_input = (input(:,i-1)+input(:,i))/2
#     x_temp = x(:,i-1)+del_t*k1/2
#     k2 = A*x_temp + B_bar*current_input
#     x_temp = x(:,i-1)+del_t*k2/2
#     k3 = A*x_temp + B_bar*current_input
#     
#     current_input = input(:,i);
#     x_temp = x(:,i-1)+del_t*k3
#     k4 = A*x_temp + B_bar*current_input
#         
#     x(:,i) = x(:,i-1) + del_t.*(k1+2*k2+2*k3+k4)/6
# 
# end
###########################################################################
'''

t = tspan

y1 = 12 * x[:, 0]  # CG height
y2 = (180 / math.pi) * x[:, 2]  # Pitch angle
y3 = - (180 / math.pi) * x[:, 4]  # Roll angle
u1 = x[:, 6]  # FL wheel height change
u2 = x[:, 8]  # FR wheel height change
u3 = x[:, 10]  # RL wheel height change
u4 = x[:, 12]  # RR wheel height change

# Damper Position
damper_u1 = u1 * IR_F
damper_u2 = u2 * IR_F
damper_u3 = u3 * IR_R
damper_u4 = u4 * IR_R

'''
# Input derivative

input_dot = zeros(size(input))

for i = 1:N-1
    
    if i < 3
        input_dot(:,i) = (input(:,i+1)-input(:,i))/del_t

elseif i > N-2
        input_dot(:,i) = (input(:,i+1)-input(:,i))/del_t

else
        input_dot(:,i) = (input(:,i-2)-8.*input(:,i-1)+8.*input(:,i+1)-input(:,i+2))/del_t/12

end

end
'''
# Wheel height velocities, Units: ft/s
ud1 = x[:, 7]
ud2 = x[:, 9]
ud3 = x[:, 11]
ud4 = x[:, 13]

# Define Chassis Points (z); Units: ft
z1 = x[:, 0] + (D.wheelbase * (1 - W.spr_dis)) * x[:, 3] + (D.track_F / 2) * x[:, 5]
z2 = x[:, 0] + (D.wheelbase * (1 - W.spr_dis)) * x[:, 3] - (D.track_F / 2) * x[:, 5]
z3 = x[:, 0] - (D.wheelbase * W.spr_dis) * x[:, 3] + (D.track_R / 2) * x[:, 5]
z4 = x[:, 0] - (D.wheelbase * W.spr_dis) * x[:, 3] - (D.track_R / 2) * x[:, 5]

# Define Chassis Point velocities (zd); Units: ft/s
zd1 = x[:, 2] + (D.wheelbase * (1 - W.spr_dis)) * x[:, 4] + (D.track_F / 2) * x[:, 6]
zd2 = x[:, 2] + (D.wheelbase * (1 - W.spr_dis)) * x[:, 4] - (D.track_F / 2) * x[:, 6]
zd3 = x[:, 2] - (D.wheelbase * W.spr_dis) * x[:, 4] + (D.track_R / 2) * x[:, 6]
zd4 = x[:, 2] - (D.wheelbase * W.spr_dis) * x[:, 4] - (D.track_R / 2) * x[:, 6]

print('Calc time:', time.time() - int_time)

# Figure 1
plt.plot(tspan, y1, label='CG Height [in]', linewidth=0.5)
plt.plot(tspan, y2, label='Pitch Angle [deg]')
plt.plot(tspan, y3, label='Roll Angle [deg]')
plt.xlabel("Time [s]")
plt.title("CG Height, Pitch, and Roll vs. Time")
plt.legend(loc='upper right')
plt.show()

'''
# Calculate wheel loads (add tire damping here if added later)
n_fl = W.vehicle_total * W.distribution/2 - K_Tire*(u1 - sim_input[4,:]) - Tire_Damping*(ud1-input_dot(4,:))
n_fr = W.vehicle_total*W.distribution/2 - K_Tire*(u2-input(5,:)) - Tire_Damping*(ud2-input_dot(5,:))
n_rl = W.vehicle_total*(1-W.distribution)/2 - K_Tire*(u3-input(6,:)) - Tire_Damping*(ud3-input_dot(6,:))
n_rr = W.vehicle_total*(1-W.distribution)/2 - K_Tire*(u4-input(7,:)) - Tire_Damping*(ud4-input_dot(7,:))

# Wheel forces from springs/dampers (not ARBs)
force_fl = wheelRateF * (u1 - z1) + wheelDampF * (ud1 - zd1) + W.sprung * W.spr_dis / 2
force_fr = wheelRateF * (u2 - z2) + wheelDampF * (ud2 - zd2) + W.sprung * W.spr_dis / 2
force_rl = wheelRateR * (u3 - z3) + wheelDampR * (ud3 - zd3) + W.sprung * (1 - W.spr_dis) / 2
force_rr = wheelRateR * (u4 - z4) + wheelDampR * (ud4 - zd4) + W.sprung * (1 - W.spr_dis) / 2

# total_n = n.fr+n.fl+n.rr+n.rl  # Check total vehicle weight

# Plot wheel loads over time
figure(2); hold on
plot(t,n.fr,'k',t,n.fl,'b',t,n.rr,'g',t,n.rl,'r')
legend('FR','FL','RR','RL')
xlabel('Time [s]')
ylabel('Wheel Loads [lbf]')
title('Wheel Loads vs. Time')
'''

# figure(3); plot(t,total_n) # Check total vehicle weight
# xlabel('Time [s]')

# Plot wheel forces from springs/dampers (not ARBs) over time
# figure(4); hold on
# plot(t,force.fr,'k',t,force.fl,'b',t,force.rr,'g',t,force.rl,'r')
# legend('FR','FL','RR','RL')
# xlabel('Time [s]')
# ylabel('Suspension Force at Wheels [lbf]')
# title('Force at Uprights from Springs/Dampers (no ARBs) vs. Time')

output = np.array([tspan, y1, y2, y3, damper_u1, damper_u2, damper_u3, damper_u4])
out_file = open("outNL.txt", 'w')
np.savetxt('outNL.txt', output.T)
