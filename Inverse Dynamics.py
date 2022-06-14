import numpy as np
import math

# l1 = int(input('Length of 1st Link'))
# l2 = int(input('Length of 2nd Link'))
# l3 = int(input('Length of 3rd Link'))

Xe = int(input('X coordinate of End Effector: '))
Ye = int(input('Y coordinate of End Effector: '))

Y = int(input('Orientation of Robot: '))

l1 = l2 = l3 = 1

# m1 = int(input('Mass of 1st Link'))
# m2 = int(input('Mass of 2nd Link'))
# m3 = int(input('Mass of 3rd Link'))

m1 = m2 = m3 = 2

g = 9.81

x3 = Xe - l3*math.cos(Y*math.pi/180)
y3 = Ye - l3*math.cos(Y*math.pi/180)

alpha = math.acos((x3**2 + y3**2 - l1**2 - l2**2)/(2*l1*l2))
beta = math.asin((l2*math.sin(alpha))/math.sqrt(x3**2 + y3**2))

theta1a = math.atan(y3/x3) - beta
theta1b = math.atan(y3/x3) + beta

theta2a = (math.pi-alpha)
theta2b = -(math.pi-alpha)

theta3a = Y*math.pi/180 - (theta1a+theta2a)
theta3b = Y*math.pi/180 - (theta1b+theta2b)

Ja = np.array([[-Ye, l1*math.sin(theta1a)-Ye, -l3*math.sin(theta1a + theta2a + theta3a)], [Xe, Xe-l1*math.cos(theta1a), l3*math.cos(theta1a + theta2a + theta3a)]])
Jb = np.array([[-Ye, l1*math.sin(theta1b)-Ye, -l3*math.sin(theta1b + theta2b + theta3b)], [Xe, Xe-l1*math.cos(theta1b), l3*math.cos(theta1b + theta2b + theta3b)]])



theta1 = 90*math.pi/180
theta1_d = 10 #rad/sec
theta1_dd = 5 #rad/sec^2
theta2 = 90*math.pi/180
theta2_d = 10
theta2_dd = 5
theta3 = 90*math.pi/180
theta3_d = 10
theta3_dd = 5

theta = np.array([[theta1], [theta2], [theta3]])
theta_d = np.array([[theta1_d], [theta2_d], [theta3_d]])
theta_dd = np.array([[theta1_dd], [theta2_dd], [theta3_dd]])


def findJointTorque():
    # Mass Matrix
    M11 = ((l1**2)*m1)/3 + ((l2**2)*m2)/12 + (m3*(3*l1**2 + 6*l1*l2*math.cos(theta2) + 3*l1*l3*math.cos(theta2+theta3) + 3*l2**2 + 3*l2*l3*math.cos(theta3) + l3**2))/3

    M12 = ((l2**2)*m2)/12 + (m3*(6*l1*l2*math.cos(theta2) + 3*l1*l3*math.cos(theta2+theta3) + 6*l2**2 + 6*l2*l3*math.cos(theta3) + 2*l3**2))/6

    M13 = m3*l3*(3*l1*math.cos(theta2 + theta3) + 3*l2*math.cos(theta3) + 2*l3)/6

    M21 = l1*l2*m3*math.cos(theta2) + (l1*l3*m3*math.cos(theta2 + theta3))/2 + ((l2**2)*m2)/12 + (l2**2)*m3 + \
          l2*l3*m3*math.cos(theta3) + ((l3**2)*m3)/3

    M22 = ((l2**2)*m2*(math.sin(theta1+2*theta2) + math.sin(2*theta2))**2)/4 + ((l2**2)*m2*(math.cos(theta1 + 2*theta2) + math.cos(2*theta2))**2)/4 + \
        ((l2**2)*m2)/12 + ((l2**2)*m2)/12 + (l2**2)*m3 + l2*l3*m3*math.cos(theta3) + ((l3**2)*m3)/3

    M23 = (l3*m3*(3*l2*math.cos(theta3) + 2*l3))/6

    M31 = (l3*m3*(3*l1*math.cos(theta2 + theta3) + 3*l2*math.cos(theta3) + 2*l3))/6

    M32 = (l3*m3*(3*l2*math.cos(theta3)+ 2*l3))/6

    M33 = ((l3**2)*m3)/3

    M = np.array([[M11, M12, M13], [M21, M22, M23], [M31, M32, M33]])
    print('M = ', M)

    # Velocity Vector

    V11 = -2*l1*l3*m3*(math.sin(theta2))*theta1_d*theta2_d - l1*l2*m3*(theta2_d**2)*math.sin(theta2_d) - \
          l1*l3*m3*theta1_d*theta2_d*math.sin(theta2+theta3) - l1*l3*m3*theta1_d*theta3_d*math.sin(theta2 + theta3) -\
          (l1*l3*m3*(theta2_d**2)*math.sin(theta2 + theta3))/2 - l1*l3*m3*theta2_d*theta3_d*math.sin(theta2 + theta3) -\
            (l1*l3*m3*(theta3_d**2)*math.sin(theta2 + theta3))/2 + ((l2**2)*m2*(math.sin(theta1+2*theta2)*math.cos(2*theta2))*(theta2_d**2))/4 -\
        ((l2**2)*m2*math.sin(2*theta2)*math.cos(theta1+2*theta2)*theta2**2)/4 - l2*l3*m3*theta1_d*theta3_d*math.sin(theta3) - \
          l2*l3*m3*theta2_d*theta3_d*math.sin(theta3) - (l2*l3*m3*math.sin(theta3)*(theta3_d**2))/2

    V21 = l1*l2*m3*(theta1_d**2)*math.sin(theta2) - l1*l2*m3*math.sin(theta2)*theta1_d*theta2_d +\
          (l1*l3*m3*math.sin(theta2+theta3)*(theta1_d**2))/2 - (l1*l3*m3*math.sin(theta2+theta3)*(theta1_d*theta2_d))/2 - \
          ((l2**2)*m2*math.sin(theta1+2*theta2)*math.cos(2*theta2)*(theta1_d*theta2_d))/2 + \
          ((l2**2)*m2*math.cos(theta1+2*theta2)*math.sin(2*theta2)*(theta1_d*theta2_d))/2 - \
        l2*l3*m3*math.sin(theta3)*theta1_d*theta3_d - l2*l3*m3*math.sin(theta3)*theta2_d*theta3_d -\
          (l2*l3*m3*math.sin(theta3)*(theta3_d**2))/2

    V31 = l3*m3*(l1*math.sin(theta2+theta3)*(theta1_d**2) - l1*math.sin(theta2+theta3)*theta1_d*theta3_d + l2*math.sin(theta3)*(theta1_d**2) + \
        2*l2*math.sin(theta3)*theta1_d*theta2_d - l2*math.sin(theta3)*theta1_d*theta3_d + l2*math.sin(theta3)*(theta2_d**2) -l2*math.sin(theta3)*theta2_d*theta3_d)/2

    V = np.array([[V11], [V21], [V31]])
    print('V = ', V)

    # gravitational vector

    G11 = g*(l1*m1*math.cos(theta1) + 2*l1*m2*math.cos(theta1) + 2*l1*m3*math.cos(theta1) + l2*m2*math.cos(theta1 + theta2) + \
            2*l2*m3*math.cos(theta1+theta2) + l3*m3*math.cos(theta1+theta2+theta3))/2

    G21 = g*(l2*m2*math.cos(theta1+theta2) + 2*l2*m3*math.cos(theta1+theta2) + l3*m3*math.cos(theta1+theta2+theta3))/2

    G31 = g*(l3*m3*math.cos(theta1+theta2+theta3))/2

    G = np.array([[G11], [G21], [G31]])
    print('G = ', G)

    # Joint Torque T

    print('M*theta_dd = ', np.dot(M, theta_dd))
    T = np.dot(M, theta_dd) + V + G

    return T


print('T = ', findJointTorque())

