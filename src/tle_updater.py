import sys
import math
from math import radians

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
def sgp(t):
    # Initializing constants
    # print("\n------------INITIALIZING CONSTANTS------------\n")

    k_e = 0.0743669161 # er ^ (3/2) s ^ -1
    a1 = (k_e / n0) ** (2 / 3) 
    # print("a1 =", a1, "er")

    J2 = 0.00108262545
    # using the value from the paper
    JE = 5.413080 * 10 ** -4 * 2
    aE = 1 # radius of the earth
    cos_i0 = math.cos(i0)
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


    # more constants...
    J3 = -0.253881 * 10 ** -5 # paper
    sin_i0 = math.sin(i0)


    # print("\n------------LONG PERIOD PERIODICS------------\n")

    a_xNSL = e * math.cos(w_s0)
    # print("a_xNSL =", a_xNSL, "[unitless]")

    a_yNSL = e * math.sin(w_s0) - J3 / J2 / 2 * aE / p * sin_i0
    # print("a_yNSL =", a_yNSL, "[unitless]")

    L = Ls - J3 / J2 / 4 * aE / p * a_xNSL * sin_i0 * (3 + 5 * cos_i0) / (1 + cos_i0)
    L = L % (2 * math.pi)
    # print("L =", L, "radians")


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

    # print("\n-------UPDATE FOR SHORT PERIODICS-------\n")
    # more constants...
    sin2u = (cosU + cosU) * sinU
    cos2u = 1 - 2 * sinU * sinU
    pL_squared = pL * pL

    r_k = r + J2 / 4 * (aE * aE) / pL * sin_i0 * sin_i0 * cos2u
    # print("r_k =", r_k, "er")

    u_k = u - J2 / 8 * (aE * aE) / (pL * pL) * (7 * cos_i0 * cos_i0 - 1) * sin2u
    # print("u_k =", u_k, "radians")

    omega_k = omega_s0 + 3 / 4 * J2 * (aE * aE) / (pL * pL) * cos_i0 * sin2u
    # print("omega_k =", omega_k, "radians")

    i_k = i0 + 3 / 4 * J2 * (aE * aE) / (pL * pL) * sin_i0 * cos_i0 * cos2u
    # print("i_k =", i_k, "radians")

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
    # print("M =")
    # print("    ", Mx)
    # print("    ", My)
    # print("    ", Mz)

    Nx = cos_omega_k
    Ny = sin_omega_k
    Nz = 0
    # print("N =")
    # print("    ", Nx)
    # print("    ", Ny)
    # print("    ", Nz)

    Ux = Mx * sin_u_k + Nx * cos_u_k
    Uy = My * sin_u_k + Ny * cos_u_k
    Uz = Mz * sin_u_k + Nz * cos_u_k

    Vx = Mx * cos_u_k - Nx * sin_u_k
    Vy = My * cos_u_k - Ny * sin_u_k
    Vz = Mz * cos_u_k - Nz * sin_u_k

    # print("U =")
    # print("    ", Ux)
    # print("    ", Uy)
    # print("    ", Uz)

    # print("V =")
    # print("    ", Vx)
    # print("    ", Vy)
    # print("    ", Vz)

    # print("\n------------POSITION AND VELOCITY------------\n")

    x = r_k * Ux
    y = r_k * Uy
    z = r_k * Uz

    x_dot = r_dot * Ux + r_v_dot * Vx
    y_dot = r_dot * Uy + r_v_dot * Vy
    z_dot = r_dot * Uz + r_v_dot * Vz

    # print("r =")
    # print("    ", x)
    # print("    ", y)
    # print("    ", z)

    # print("r_dot =")
    # print("    ", x_dot)
    # print("    ", y_dot)
    # print("    ", z_dot)

    # print("\n------------CONVERTING TO KM AND SECONDS-----------\n")
    x *= ER_TO_KM
    y *= ER_TO_KM
    z *= ER_TO_KM

    x_dot *= ER_TO_KM / MINUTES_TO_SECONDS
    y_dot *= ER_TO_KM / MINUTES_TO_SECONDS
    z_dot *= ER_TO_KM / MINUTES_TO_SECONDS
    print("Time since epoch =", t, "minutes")
    print("r =")
    print("    ", x)   
    print("    ", y)   
    print("    ", z)   

    print("r_dot =")
    print("    ", x_dot)
    print("    ", y_dot)
    print("    ", z_dot)
    
#print("{0}\n".format(lines))

satellite_name = format(lines[0].strip())
satellite_number = format(lines[1][2:7].strip())
classification = format(lines[1][7].strip())
launch_year = format(lines[1][8:11].strip())
launch_number = format(lines[1][11:14].strip())
peice = format(lines[1][14:18].strip())
epoch_year = format(lines[1][18:20].strip())

epoch_time = float(lines[1][20:32].strip())
epoch_day = math.floor(epoch_time) 
print("Epoch Day:                                       " + str(epoch_day))
epoch_time-=epoch_day 
epoch_hour = 24 * epoch_time
print("Epoch Hour:                                      " + str(epoch_hour))
epoch_hour-=math.floor(epoch_hour)
epoch_minute = 60 * epoch_hour
print("Epoch Minute:                                    " + str(epoch_minute))
epoch_minute-=math.floor(epoch_minute)
epoch_second = 60 * epoch_minute
print("Epoch Second:                                    " + str(epoch_second))

print("First Time Derivative of Mean Motion:            " + str(lines[1][32:44].strip()))
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
print("Inclination:                                     " + inclination)
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
    sgp(t)
    t += 360
    #Print new version of tle after epoch t
    print("The TLE after time " + str(t) + " is:")
    print(satellite_name)
    print("1 " + str(satellite_number) + str(classification) + " " + str(launch_number) + str(launch_number) + str(peice))