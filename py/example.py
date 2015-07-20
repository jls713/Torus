import Torus_py
import numpy as np
import matplotlib.pyplot as plt

TorusPath = '../'

Gal = Torus_py.GalaxyPotential(TorusPath+"pot/PJM11_convenient.Tpot")
print "PJM11_convenient Phi(R=1,z=1)="+str(Gal(1., 1.))

# Initialise a logarithmic potential and plot
Log = Torus_py.LogarithmicPotential(240. / 977.775, 0.9, 1e-5, 0.)

x = np.linspace(0., 10., 1000)
plt.plot(x, map(lambda i: Log(i, 0.), x))
plt.xlabel(r'$x$')
plt.ylabel(r'$\Phi$')
plt.show()
plt.clf()


# Construct torus
J = np.array([0.1, 0.1, 1.8])

Donut = Torus_py.Torus()
Donut.make_torus(Log, J, 0.)

for Theta in np.random.uniform(0., 2. * np.pi, size=(1000, 3)):
    X = Donut.map3D(Theta)
    plt.plot([X[0]], [X[1]], 'k.')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.show()


# Define user potential


class logpot_python(Torus_py.WDPotential):

    def __init__(self, Vc, q):
        Torus_py.WDPotential.__init__(self)
        self.Vc2 = Vc * Vc
        self.q2 = q * q

    def Phi(self, Rz):
        return .5 * self.Vc2 * np.log(Rz[0] * Rz[0] + Rz[1] * Rz[1] / self.q2)

    def Forces(self, Rz):
        F = -self.Vc2 * Rz / (Rz[0] * Rz[0] + Rz[1] * Rz[1] / self.q2)
        F[1] = F[1] / self.q2
        return F

Log2 = logpot_python(240. / 977.775, 0.9)
Donut2 = Torus_py.Torus()
Donut2.make_torus(Log2, J, 0.)

for Theta in np.random.uniform(0., 2. * np.pi, size=(1000, 3)):
    X = Donut2.map3D(Theta)
    Y = Donut.map3D(Theta)
    plt.plot([X[0]], [X[1]], 'k.')
    plt.plot([Y[0]], [Y[1]], 'rx')
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')
plt.show()
