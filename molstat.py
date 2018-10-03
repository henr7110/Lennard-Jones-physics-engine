# -*- coding: utf-8 -*-
from matplotlib import pyplot as plt
import numpy as np
from matplotlib import animation, rc
import boxes
from particle import particle
import imp
imp.reload(boxes)
n_atoms = 10
sigma = 0.34e-09 #m
t_steps = 1
epsilon = 1.65e-21 #J

m = 6.63e-25
box_width = 1.0
b = boxes.Lennybox(n_atoms, sigma, m, epsilon, t_steps)
b.plot(-b.box_width,b.box_width,-b.box_width,b.box_width)
b.sigma / b.dt

b.run(10)
print([i.getcpos() for i in b.particles])
print([i.getv() for i in b.particles])
b.plot(-b.box_width,b.box_width,-b.box_width,b.box_width)
