

class simbox:
  import numpy as np

  def __init__(self, n_atoms, r_min, vval, n_steps, box_width):
    import numpy as np
    import particle
    import imp
    imp.reload(particle)
    self.n_atoms = n_atoms
    self.r_min = r_min
    self.particles = []
    self.dt =  (2 * r_min) / (vval)
    self.n_steps = n_steps
    self.box_width = box_width
    self.partdist = []

    for i in range(n_atoms):
      x = np.random.uniform(box_width/2. + r_min, box_width - r_min)
      y = np.random.uniform(-box_width + r_min, box_width - r_min)

      vx = np.random.uniform(-vval, vval)
      vy = np.random.uniform(-vval, vval)

      v = np.array([vx,vy])

      self.particles.append(particle.particle(x, y, v, i))


  def disc(self):
    [i.disc() for i in self.particles]

  def run(self, pval):
    from matplotlib import pyplot as plt

    from matplotlib import animation, rc
    # animation function. This is called sequentially
    def animate(i):

      if i % pval == 0 and i != 0:
          print(i)
      #Check for collisions
      self.Colhandler()
      #Update positions
      for p in self.particles:
          p.Update(self.dt, self.box_width, self.r_min)

      line.set_data([i.getx() for i in self.particles],
                    [i.gety() for i in self.particles])

      return (line,)
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    #calculate radius in plot
    r = self.r_min
    rpix = (ax.transData.transform((0,r)) - ax.transData.transform((0,0)))[1]

    #generate particle plots
    line, = plt.plot([i.getx() for i in self.particles],
                    [i.gety() for i in self.particles], "ro",
                    markersize = rpix)
    plt.xlim(-1,1)
    plt.ylim(-1,1)

    # call the animator. blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, frames=self.n_steps, blit=True)

    anim.save('im2.mp4', writer=writer)
    print("done")

  def Colhandler(self):
      import numpy as np
      for p1 in self.particles:
          for p2 in self.particles:
              if p1 > p2:
                  dist = np.sqrt((p1.getx() - p2.getx())**2 \
                                 + (p1.gety() - p2.gety())**2)
                  if dist <= self.r_min * 2 and p1 != p2:
                      v1 = p1.getv().copy()
                      v2 = p2.getv().copy()

                      p1.v = v2
                      p2.v = v1

  def plot(self,xmin ,xmax, ymin, ymax):
    """plots box particles in plot with range
    xmin, xmax, ymin, ymax)"""
    from matplotlib import pyplot as plt
    plt.plot([i.getx() for i in self.particles],
             [i.gety() for i in self.particles], 'ro')


    if xmax != None or xmin != None:
      plt.xlim(xmin, xmax)

    if ymax != None or ymin != None:
      plt.ylim(ymin, ymax)
    plt.show()

class nosimbox(simbox):
    def partdistrecord(self):
        partdist = 0
        for p in self.particles:
            if p.getx() > 0:
                partdist += 1
        inserter = [partdist]
        self.partdist += inserter

    def run(self):
        from matplotlib import pyplot as plt
        from matplotlib import animation, rc
        # animation function. This is called sequentially
        for i in range(self.n_steps):
          if i % 1000 == 0 and i != 0:
              print(i)
          if partdist and i % 5 == 0:
              self.partdistrecord()

          #Check for collisions
          self.Colhandler()
          #Update positions
          for p in self.particles:
              p.Update(self.dt, self.box_width, self.r_min)
        if partdist:
            plt.close()

            plt.hist(self.partdist,bins=range(0, self.n_atoms), align="left")
            plt.title("partdist after %d t_steps" % self.n_steps)
            plt.savefig("dist_histogram.png")
            plt.show()

        print("done")

class Lennybox(simbox):
  def __init__(self, n_atoms, sigma, epsilon, m, n_steps):
    import numpy as np
    import particle
    import imp

    imp.reload(particle)
    self.n_atoms = n_atoms
    self.m = m
    self.epsilon = epsilon
    self.sigma = sigma
    self.tau = np.sqrt((self.m * self.sigma) / self.epsilon)
    self.particles = []
    self.dt =  0.001 * self.tau
    self.n_steps = n_steps
    box_width =  n_atoms * sigma*2
    self.box_width = box_width
    self.partdist = []
    vval =  self.sigma / self.tau

    for i in range(n_atoms):
      x = np.random.uniform(box_width/2. + sigma, box_width - sigma)
      y = np.random.uniform(-box_width + sigma, box_width - sigma)
      pos = np.array([x,y])

      vx = np.random.uniform(-1, 1) * vval
      vy = np.random.uniform(-1, 1) * vval
      v = np.array([vx,vy])


      pastpos = pos.copy() - v.copy()*self.dt

      self.particles.append(particle.ForceParticle(pastpos, pos, i, m, self.dt))

      #assign past forces to the particles
      forces = self.Forcecalc()
      [i.setpforce(forces[i.getnum()]) for i in self.particles]

    print("startT = " + str(self.Tcalc()))


  def Tcalc(self):
      sumvx = 0
      sumvy = 0
      for i in self.particles:
          sumvx += i.getv()[0]**2
          sumvy += i.getv()[1]**2
      KE = 0.5 * self.m *(sumvx/float(self.n_atoms) + sumvy/float(self.n_atoms))
      return(1.4e-23*KE)



  def run(self, pval):
    from matplotlib import pyplot as plt

    from matplotlib import animation, rc
    # animation function. This is called sequentially
    def animate(i):

      if i % pval == 0 and i != 0:
          print(i)

      forces = self.Forcecalc()
      #Update positions and velocities
      for p in self.particles:
          p.xUpdate(self.box_width, self.sigma, forces[p.getnum()])
          p.vUpdate(forces[p.getnum()])

      line.set_data([i.getx() for i in self.particles],
                    [i.gety() for i in self.particles])

      return (line,)
    # Set up formatting for the movie files
    Writer = animation.writers['ffmpeg']
    writer = Writer(fps=15, metadata=dict(artist='Me'), bitrate=1800)

    # First set up the figure, the axis, and the plot element we want to animate
    fig = plt.figure()
    ax = fig.add_subplot(111, aspect='equal')

    plt.xlim(-self.box_width,self.box_width)
    plt.ylim(-self.box_width,self.box_width)
    #calculate radius in plot
    r = self.sigma
    rpix = (ax.transData.transform((0,r)) - ax.transData.transform((0,0)))[1]

    #generate particle plots
    line, = plt.plot([i.getcpos()[0] for i in self.particles],
                    [i.getcpos()[1] for i in self.particles], "ro",
                    markersize = 3)


    # call the animator. blit=True means only re-draw the parts that have changed.
    anim = animation.FuncAnimation(fig, animate, frames=self.n_steps, blit=True)

    anim.save('im2.mp4', writer=writer)
    print("done" + "endT = " + str(self.Tcalc()))

  def plot(self,xmin ,xmax, ymin, ymax):
    """plots box particles in plot with range
    xmin, xmax, ymin, ymax)"""
    from matplotlib import pyplot as plt
    plt.plot([i.getcpos()[0] for i in self.particles],
             [i.getcpos()[1] for i in self.particles], 'ro')

    if xmax != None or xmin != None:
      plt.xlim(xmin, xmax)

    if ymax != None or ymin != None:
      plt.ylim(ymin, ymax)
    plt.show()

  def getsigma(self):
      return(self.sigma)

  def getepsilon(self):
      return(self.epsilon)

  def dist(self, pos1, pos2):
      import numpy as np
      return(np.sqrt((pos1[0]-pos2[0])**2 + (pos1[1]-pos2[1])**2))

  def Forcecalc(self):
      import sympy
      import numpy as np
      s = self.getsigma()
      e = self.getepsilon()
      x = sympy.symbols("x")
      forces = np.zeros((self.n_atoms,2))

      for i in self.particles:
          force = np.zeros(3)
          for j in self.particles:
              if i > j:
                  r = self.dist(i.getcpos(), j.getcpos())

                  force = np.float(sympy.diff(-4 * e * ((s/x)**12 - \
                          (s/x)**6),x).subs(x,r).evalf())

                  direction = (j.getcpos() - i.getcpos()) / r
                  forces[j.getnum()] += force * direction
                  forces[i.getnum()] -= force * direction
      return forces

class testbox(simbox):

  def __init__(self, n_atoms, r_min, vval, n_steps, box_width):
    from matplotlib import pyplot as plt
    import numpy as np
    from matplotlib import animation, rc
    import particle
    import imp
    self.n_atoms
    imp.reload(particle)
    self.r_min = r_min
    self.particles = []
    self.dt = (2 * r_min) / (vval)
    self.n_steps = n_steps
    self.box_width = box_width

    self.particles.append(particle.particle(-box_width / 2., 0., np.array([0.2,0.]), 0))
    self.particles.append(particle.particle(box_width / 2., 0., np.array([-0.2,0.]), 1))
