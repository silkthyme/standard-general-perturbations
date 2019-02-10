import sys
import math
from math import radians
import turtle

# Useful Conversions
REVOLUTIONS_TO_RADIANS = 2 * math.pi
DAYS_TO_MINUTES = 1440
MINUTES_TO_SECONDS = 60
ER_TO_KM = 6378.135

# Obtain input as stringc
raw_data = ''
tmp_line = ''
lines = []

while True:
	line = sys.stdin.readline()
	if line == '':
		break
	else:
		raw_data += line

i = 0
for char in raw_data:
	if char != '\n':
		tmp_line += char
	else:
		lines.append(tmp_line)
		i += 1
		tmp_line = ''

def initialize_constants(t, n0_dot, n0_double_dot, i0, omega0,e0, w0, M0, n0):
    # Initializing constants
    # print("\n------------INITIALIZING CONSTANTS------------\n")

    k_e = 0.0743669161 # er ^ (3/2) s ^ -1
    a1 = (k_e / n0) ** (2 / 3) 
    # print("a1 =", a1, "er")

    J2 = 0.00108262545
    J3 = -0.253881 * 10 ** -5 # paper
    # using the value from the paper
    JE = 5.413080 * 10 ** -4 * 2
    aE = 1 # radius of the earth
    cos_i0 = math.cos(i0)
    sin_i0 = math.sin(i0)
    delta1 = (3 / 4) * J2 * (aE * aE) / (a1 * a1) * (3 * cos_i0 * cos_i0 - 1) / (1 - e0 * e0) ** (3/2)
    # print("delta1 =", delta1, "[unitless]")
    a0 = a1 * (1 - 1 / 3 * delta1 - delta1 * delta1 - 134 / 81 * delta1 * delta1 * delta1)
    # print("a0 =", a0, "er")

    p0 = a0 * (1 - e0 * e0)
    # print("p0 =", p0, "er")

    q0 = a0 * (1 - e0)
    # print("q0 =", q0, "er")

    L0 = M0 + w0 + omega0
    # print("L0 =", L0, "radians")

    d_omega_dt = - (3 / 2) * J2 * (aE * aE) / (p0 * p0) * n0 * cos_i0
    # print("d_omega_dt =", d_omega_dt, "radians per min")

    dw_dt = 3 / 4 * J2 * (aE * aE) / (p0 * p0) * n0 * (5 * cos_i0 * cos_i0 - 1)
    # print("dw_dt =", dw_dt, "radians per min")

    return (k_e, J2, J3, aE, cos_i0, sin_i0, a0, q0, L0, d_omega_dt, dw_dt)

def secular_gravity_and_atmospheric_drag(t, n0, n0_dot, a0, q0, omega0, d_omega_dt, w0, dw_dt, L0):
    # Update for secular effects of atmospheric drag and gravitation
    # print("\n-----SECULAR GRAVITY AND ATMOSPHERIC DRAG-----\n")
    a = n0 + n0_dot * t + n0_double_dot / 2 * t * t
    a = a0 * (n0 / a) ** (2 / 3)
    # print("a =", a, "er")

    e = 10 ** -6
    if a > q0:
        e = 1 - q0 / a0
    # print("e =", e, "[unitless]")

    p = a * (1 - e * e)
    # print("p =", p, "er")

    omega_s0 = omega0 + d_omega_dt * t
    # print("omega_s0 =", omega_s0, "radians")

    w_s0 = w0 + dw_dt * t
    # print("w_s0 =", w_s0, "radians")

    Ls = L0 + (n0 + dw_dt + d_omega_dt) * t \
        + n0_dot * t * t / 2 \
        + n0_double_dot * t * t * t / 6 
    Ls = Ls % (2 * math.pi)
    # print("Ls =", Ls, "radians")

    return (a, e, p, omega_s0, w_s0, Ls)

def long_period_periodics(e, w_s0, Ls, J2, J3, aE, p, cos_i0, sin_i0):
    # print("\n------------LONG PERIOD PERIODICS------------\n")

    a_xNSL = e * math.cos(w_s0)
    # print("a_xNSL =", a_xNSL, "[unitless]")

    a_yNSL = e * math.sin(w_s0) - J3 / J2 / 2 * aE / p * sin_i0
    # print("a_yNSL =", a_yNSL, "[unitless]")

    L = Ls - J3 / J2 / 4 * aE / p * a_xNSL * sin_i0 * (3 + 5 * cos_i0) / (1 + cos_i0)
    L = L % (2 * math.pi)
    # print("L =", L, "radians")

    return (a_xNSL, a_yNSL, L)

def keplers_equation(L, omega_s0, a_xNSL, a_yNSL ):
    # print("\n-----------SOLVING KEPLERS EQUATION-----------\n")
    U = L - omega_s0
    U = U % (2 * math.pi)
    # print("U =", U, "radians")

    iteration = 0 # ITEM3 in FORTRAN
    E = U
    deltaE = 1 #TEM5 in FORTRAN
    while True:
        cosE = math.cos(E)
        sinE = math.sin(E)

        # print("Iteration", iteration, "    E =", E, "radians", "    deltaE =", deltaE)

        if abs(deltaE) < 10 ** -6 or iteration >= 10:
            break


        iteration += 1

        deltaE = -a_yNSL * sinE - a_xNSL * cosE + 1
        deltaE = (U - a_yNSL * cosE + a_xNSL * sinE - E) / deltaE
        absDeltaE = abs(deltaE) # TEM2 in FORTRAN
        if absDeltaE > 1:
            deltaE = absDeltaE / deltaE

        E = E + deltaE

    # print("E =", E, "radians")

    return (E, cosE, sinE)

def short_period_preliminary_quantities(a_xNSL, a_yNSL, cosE, sinE, k_e, a):
    # print("\n----SHORT PERIOD PRELIMINARY QUANTITIES----\n")
    
    eCosE = a_xNSL * cosE + a_yNSL * sinE
    # print("eCosE =", eCosE, "[unitless]")

    eSinE = a_xNSL * sinE - a_yNSL * cosE
    # print("eSinE =", eSinE, "[unitless]")

    eL_squared = a_xNSL * a_xNSL + a_yNSL * a_yNSL
    # print("eL_squared =", eL_squared, "[unitless]")

    pL = a * (1 - eL_squared)
    # print("pL =", pL, "er")

    r = a * (1 - eCosE)
    # print("r =", r, "er")

    r_dot = k_e * math.sqrt(a) / r * eSinE
    # print("r_dot =", r_dot, "er per min")

    r_v_dot = k_e * math.sqrt(pL) / r
    # print("r_v_dot=", r_v_dot, "er per min")

    temp = eSinE / (1 + math.sqrt(1 - eL_squared))
    sinU = a / r * (sinE - a_yNSL - a_xNSL * temp)
    cosU = a / r * (cosE - a_xNSL + a_yNSL * temp)
    # print("sinU =", sinU, "[unitless]")
    # print("cosU =", cosU, "[unitless]")

    u = math.atan2(sinU, cosU)
    if u < 0:                  # ensures that u is in range [0, 2pi)
        u += 2 * math.pi
    # print("u =", u, "radians")

    # print("    checking u...")
    # print("    sinU / cosU =", sinU/cosU)
    # print("    tan(u) =", math.tan(u))

    return (pL, r, r_dot, r_v_dot, sinU, cosU, u)

def update_for_short_periodics(cosU, sinU, pL, r , J2, aE, sin_i0, cos_i0, omega_s0, u):
    # print("\n-------UPDATE FOR SHORT PERIODICS-------\n")
    # more constants...
    sin2u = (cosU + cosU) * sinU
    cos2u = 1 - 2 * sinU * sinU

    r_k = r + J2 / 4 * (aE * aE) / pL * sin_i0 * sin_i0 * cos2u
    # print("r_k =", r_k, "er")

    u_k = u - J2 / 8 * (aE * aE) / (pL * pL) * (7 * cos_i0 * cos_i0 - 1) * sin2u
    # print("u_k =", u_k, "radians")

    omega_k = omega_s0 + 3 / 4 * J2 * (aE * aE) / (pL * pL) * cos_i0 * sin2u
    # print("omega_k =", omega_k, "radians")

    i_k = i0 + 3 / 4 * J2 * (aE * aE) / (pL * pL) * sin_i0 * cos_i0 * cos2u
    # print("i_k =", i_k, "radians")

    return (r_k, u_k, omega_k, i_k)

def orientation_vectors(u_k, omega_k, i_k):
    # print("\n------------ORIENTATION VECTORS------------\n")
    # constants ...
    sin_u_k = math.sin(u_k)
    cos_u_k = math.cos(u_k)
    sin_omega_k = math.sin(omega_k)
    cos_omega_k = math.cos(omega_k)
    sin_i_k = math.sin(i_k)
    cos_i_k = math.cos(i_k)

    Mx = -sin_omega_k * cos_i_k
    My = cos_omega_k * cos_i_k
    Mz = sin_i_k

    Nx = cos_omega_k
    Ny = sin_omega_k
    Nz = 0

    Ux = Mx * sin_u_k + Nx * cos_u_k
    Uy = My * sin_u_k + Ny * cos_u_k
    Uz = Mz * sin_u_k + Nz * cos_u_k

    Vx = Mx * cos_u_k - Nx * sin_u_k
    Vy = My * cos_u_k - Ny * sin_u_k
    Vz = Mz * cos_u_k - Nz * sin_u_k

    return ([Ux, Uy, Uz], [Vx, Vy, Vz])

def position_and_velocity(r_k, r_dot, r_v_dot, U, V):
    # print("\n------------POSITION AND VELOCITY------------\n")
    R = [r_k * x for x in U]

    R_dot = [r_dot * u + r_v_dot * v for (u,v) in zip(U,V)]

    # x_dot = r_dot * Ux + r_v_dot * Vx
    # y_dot = r_dot * Uy + r_v_dot * Vy
    # z_dot = r_dot * Uz + r_v_dot * Vz

    return (R, R_dot)

def convert_to_km_and_seconds(R, R_dot):
    # print("\n------------CONVERTING TO KM AND SECONDS-----------\n")
    R = [ER_TO_KM * x for x in R]

    R_dot = [ER_TO_KM / MINUTES_TO_SECONDS * x_dot for x_dot in R_dot]

    return (R, R_dot)

def sgp(t, n0_dot, n0_double_dot, i0, omega0,e0, w0, M0, n0):
    k_e, J2, J3, aE, cos_i0, sin_i0, a0, q0, L0, d_omega_dt, dw_dt \
        = initialize_constants(t, n0_dot, n0_double_dot, i0, omega0,e0, w0, M0, n0)

    a, e, p, omega_s0, w_s0, Ls \
        = secular_gravity_and_atmospheric_drag(t, n0, n0_dot, a0, q0, omega0, d_omega_dt, w0, dw_dt, L0)

    a_xNSL, a_yNSL, L \
        = long_period_periodics(e, w_s0, Ls, J2, J3, aE, p, cos_i0, sin_i0)

    E, cosE, sinE \
        = keplers_equation(L, omega_s0, a_xNSL, a_yNSL)

    pL, r, r_dot, r_v_dot, sinU, cosU, u \
        = short_period_preliminary_quantities(a_xNSL, a_yNSL, cosE, sinE, k_e, a)

    r_k, u_k, omega_k, i_k \
        = update_for_short_periodics(cosU, sinU, pL, r , J2, aE, sin_i0, cos_i0, omega_s0, u)

    U, V = \
        orientation_vectors(u_k, omega_k, i_k)
    
    R, R_dot = \
        position_and_velocity(r_k, r_dot, r_v_dot, U, V)

    R, R_dot = \
        convert_to_km_and_seconds(R, R_dot)
    
    return (R, R_dot)

def print_state_vectors(R, R_dot):
    print("r =")
    for x in R:
        print("    ", x)
    
    print("r_dot =")
    for x_dot in R_dot:
        print("    ", x_dot)

        
#print("{0}\n".format(lines))

satellite_name = format(lines[0].strip())
satellite_number = format(lines[1][2:7].strip())
classification = format(lines[1][7].strip())
launch_year = format(lines[1][8:11].strip())
launch_number = format(lines[1][11:14].strip())
peice = format(lines[1][14:18].strip())
epoch_year = format(lines[1][18:20].strip())
epoch_day = format(lines[1][20:23].strip())
epoch_day_fraction = format(lines[1][24:32].strip())



first_time_derivative_of_mean_motion  = str(lines[1][32:44].strip())
print("First Time Derivative of Mean Motion:            " + first_time_derivative_of_mean_motion)
n0_dot = lines[1][32:44].strip()
n0_dot = 2 * float(n0_dot) # rev per day squared and reverse 2 division
print("----n0_dot", n0_dot, "rev per day squared")
n0_dot = n0_dot * REVOLUTIONS_TO_RADIANS / (DAYS_TO_MINUTES * DAYS_TO_MINUTES) # radians per min squared

print("----n0_dot", n0_dot, "radians per min squared")

print("Second Time Derivative of Mean Motion:           "+str(lines[1][44:53].strip()))
temp = lines[1][44:53]
mantissa = temp[0:6]
power = temp[-3:]
mantissa = int(mantissa)
print("--------mantissa", mantissa)
power = int(power)
print("--------power", power)
n0_double_dot = 6 * (mantissa / 100000) * 10 ** power # reverse 6 division
print("----n0_double_dot", n0_double_dot, "rev per day cubed")
n0_double_dot *= REVOLUTIONS_TO_RADIANS / (DAYS_TO_MINUTES ** 3) # radians per min cubed
print("----n0_double_dot", n0_double_dot, "radians per min cubed")

bstar_drag_term = format(lines[1][53:62].strip())
print("BSTAR Drag Term:                                 " + str(bstar_drag_term))
ephemeris_type = format(lines[1][62:64].strip())
print("Ephemeris Type:                                  " + str(ephemeris_type))
element_set_number = format(lines[1][64:68].strip())
print("Element Set Number:                              " + str(element_set_number))
checksum = format(lines[1][68:70].strip())
print("Checksum:                                        " + str(checksum))

print("\n")
print("Satellite Number:                                " + str(satellite_number))
inclination = format(lines[2][8:17].strip())
print("Inclination:                                     " + str(inclination))
i0 = lines[2][8:17].strip()
i0 = float(i0)
print("----i0", i0, "degrees")
i0 = radians(i0)
print("----i0", i0, "radians")

print("Right Ascension of Ascending Node:               " + str(lines[2][17:26].strip()))
omega0 = lines[2][17:26].strip()
omega0 = float(omega0)
print("----omega0", omega0, "degrees")
omega0 = radians(omega0)
print("----omega0", omega0, "radians")

print("Eccentricity:                                    " + str(lines[2][26:34].strip()))
e0 = lines[2][26:34].strip()
e0 = float(e0) / (10 ** 7)
print("----e0", e0, "[unitless (decimal point assumed converted)]")

print("Argument of Perigee:                             " + str(lines[2][34:43].strip()))
w0 = lines[2][34:43].strip()
w0 = float(w0) 
print("----w0", w0, "degrees")
w0 = radians(w0)
print("----w0", w0, "radians")

print("Mean Anomaly:                                    " + str(lines[2][43:52].strip()))
M0 = lines[2][43:52].strip()
M0 = float(M0)
print("----M0", M0, "degrees")
M0 = radians(M0)
print("----M0", M0, "radians")

print("Mean Motion:                                     " + str(lines[2][52:63].strip()))
n0 = lines[2][52:63].strip()
n0 = float(n0)
print("----n0", n0, "rev per day")
n0 *= REVOLUTIONS_TO_RADIANS / DAYS_TO_MINUTES
print("----n0", n0, "radians per min")

print("Revolutions at Epoch:                            " + str(lines[2][63:68].strip()))
print("Checksum:                                        " + str(lines[2][68:70].strip()))




print("\n\n\n-------------SGP VALUES------------\n\n\n")
t = 0
while t <= 1440:
    R, R_dot = sgp(t, n0_dot, n0_double_dot, i0, omega0,e0, w0, M0, n0)
    print("Time since epoch =", t, "minutes")
    print_state_vectors(R, R_dot)

    #Print new version of tle after epoch t
    print("The TLE after time " + str(t) + " is:")
    print(satellite_name)
    print("1 " + str(satellite_number) + str(classification) + " " + str(launch_number) + str(launch_number) + str(peice) + " " + str(epoch_year) + str(epoch_day) + "." + str(int(epoch_day_fraction) + int(t)))
    
    t += 360