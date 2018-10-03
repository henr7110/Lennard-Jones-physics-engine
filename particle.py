
class particle:

  def __init__(self, xpos, ypos, v, num):
    import numpy as np
    self.pos = np.array([xpos, ypos])
    self.v = v
    self.num = num

  def disc(self):
    print(self.pos[0], self.pos[1], self.v)

  def getpos(self):
    return(self.pos)

  def getnum(self):
    return(self.num)

  def getx(self):
    return(self.pos[0])

  def gety(self):
    return(self.pos[1])

  def getv(self):
    return(self.v)


  def Update(self, dt, box_width, r_min):
    oldpos = self.getpos()
    oldv = self.getv()

    #if at boundary bounce
    #bounce of x-sides
    side = box_width
    if oldpos[0] < -side + r_min or oldpos[0] > side - r_min:
        self.v[0] = - oldv[0]
    if oldpos[1] < -side+ r_min or oldpos[1] > side - r_min:
        self.v[1] = -oldv[1]

    newpos = oldpos + self.v * dt

    self.pos = newpos

  def __eq__(self,other):
      if self.num == other.num:
          return True
      else:
          return False
  def __lt__(self,other):
      if self.num < other.num:
          return True
      else:
          return False

  def __gt__(self,other):
      if self.num > other.num:
          return True
      else:
          return False

  def __le__(self,other):
      if self.num <= other.num:
          return True
      else:
          return False
  def __ge__(self,other):
      if self.num >= other.num:
          return True
      else:
          return False

class ForceParticle (particle):
  def __init__(self, pastpos, currentpos, num, mass, dt):
    """all pos and v are 2D numpy arrays pastforce is the force at the pastpos
    from the initialization"""
    import numpy as np
    self.pastpos = pastpos
    self.currentpos = currentpos
    self.num = num
    self.mass = mass
    self.dt = dt

  def setpforce(self, pforce):
    self.pastforce = pforce

  def getv(self):
    return((self.getcpos() - self.getppos())/self.dt)

  def disc(self):
    pass

  def getcpos(self):
    return(self.currentpos)

  def getppos(self):
    return(self.pastpos)

  def getpforce(self):
      return(self.pastforce)

  def getmass(self):
      return(self.mass)

  def Update(self):
      pass

  def xUpdate(self, box_width, r_min, currentforce):
    """force contains the force on the particle used in the update step, pastpos
    is the position of the particle at t-dt"""

    currentpos = self.getcpos()
    oldpos = self.getppos()
    m = self.getmass()

    v = ((currentpos - oldpos) / self.dt)

    #if at boundary bounce
    #bounce of x-sides
    side = box_width
    if currentpos[0] < -side + r_min or currentpos[0] > side - r_min:
        self.v[0] = - currentv[0]
    if currentpos[1] < -side+ r_min or currentpos[1] > side - r_min:
        self.v[1] = -currentv[1]

    #update position
    newpos = currentpos + v*self.dt + (self.dt**2/(2*m) * currentforce)
    self.pos = newpos

  def vUpdate(self, currentforce):
    """updates velocity"""

    currentpos = self.getpos()
    oldpos = self.getpos()
    m = self.getmass()
    forcediff = currentforce - self.getpforce()

    pastv = ((currentpos - oldpos) / self.dt)
    self.v = pastv + (self.dt/(2*m)) * forcediff
    self.pastforce = currentforce
