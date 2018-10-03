# import numpy as np
# test_list = [0.0 for i in range(5)]
# array = np.zeros(5)
#
# V = np.ones(5)
# W = 5*V
# print "W = ", W
#
# V = np.zeros((3,5))
# print "V = " , V
# np.shape(V)
#
# L = []
# for i in range(10):
#     L.append(i)
# L = np.array(L)
# print "L = ", L
# V = np.arange(0,6).reshape(2,3)
# print "V = ", V
# print "sum of V " + str(np.sum(V))
# print "mean of V " + str(np.mean(V))

#3 Molecular dymanics using numpy
import md_video as video
import md_header as md
import numpy as np
from matplotlib import pyplot as plt
n_atoms = 4
T_want = 3
dt = 0.0005
n_steps = 500000
rho  = 0.2
eps = np.ones((n_atoms, n_atoms))
eps = eps/10.0
eps[0,1] = 0.8
eps[1,0] = 0.8
symbols = n_atoms*["Ar "]
symbols[0] = "Ne "
symbols[1] = "Ne "
video.add_color(symbols)

def dister(num1, num2, R, box_width):
    # calculate the distance between 0 and 1
    X =R[0,num1]-R[0,num2]
    Y =R[1,num1]-R[1,num2]
    # Periodic boundary condition
    X -= box_width * np.rint(X/box_width)
    Y -= box_width * np.rint(Y/box_width)
    # calculate distance
    return np.sqrt(X**2 + Y**2)

def scale_temp(V, T_want, N_a):

    gamma = np.sqrt((T_want) / (( np.sum(V*V)) / (N_a)))
    return V * gamma
def simulation(n_atoms, n_steps, rho, T_want, eps, dt):
    dist = []
    t_together = 0.0
    t_apart = 0.0
    N_a = 2 * n_atoms
    R, V, F, box_width = md.initialize_particles(n_atoms, T_want, rho, eps)
    KinE = []
    PotE = []
    TotE = []
    Temp = []
    for i in range(n_steps):
        Rnew = R + V*dt + 0.5*dt**2*F

        if i%1000 ==0:
            video.add_frame(R[0],R[1], i)
        if i%100000 == 0:
            print i

        energy_potential, Fnew = md.lennard_jones(Rnew, box_width, eps)
        V = scale_temp(V, T_want, N_a)
        ke = (1.0/2.0) * np.sum(V*V)
        d = dister(0,1,R, box_width)
        dist.append(d)
        if d < 1.52:
            t_together += 1.0
        if d > 1.52:
            t_apart += 1.0
        KinE.append(ke)
        PotE.append(energy_potential)
        TotE.append(ke + energy_potential)
        Temp.append(2 * ke / (N_a))
        Vnew = V + 0.5*dt*(F + Fnew)

        V = Vnew.copy()
        R = Rnew.copy()
    video.save("colorparticle", periodic_boundary=True)
    return PotE, KinE, TotE, Temp, dist, t_together/t_apart


PotE, KinE, TotE, Temp, dist, k = simulation(n_atoms, n_steps, rho, T_want, eps, dt)

# plot energies
for i,label in zip([PotE, KinE, TotE],["Potential Energy", "Kinetic Energy", "Total Energy"]):
    plt.plot(range(n_steps), i, label=label)
plt.ylabel("Energy")
plt.xlabel("n_step")
plt.legend()
plt.show()

#plot temperature
plt.plot(range(n_steps), Temp, label="Temperature")
plt.legend()
plt.xlabel("n_step")
plt.ylabel("Temperature (K)")
plt.show()

#plot dist between 0 and 1 against time
plt.plot(range(n_steps), dist, label="distance")
plt.legend()
plt.xlabel("n_step")
plt.ylabel("distance")
plt.show()

print "k =", k
A = -T_want*np.log(k)
print "Free energy = ", A
print "T for dissociation = ", 0 = A/T
