import copy
import numpy as np

symbols = []
trajectory = []


def add_color(symb):
	global symbols
	symbols = symb

def add_frame(R_x,R_y,n):
	global trajectory
	if n == 0:
		trajectory = []
	trajectory.append([copy.copy(R_x), copy.copy(R_y)])

def save(video_filename, box_width=10., periodic_boundary=False):
	global symbols
	n_atoms = len(trajectory[0][0])

	if len(symbols) == 0:
		symbols = n_atoms*["Ar "]

	with open(video_filename+'.xyz', 'w') as f:
		for positions_x,positions_y in trajectory:
			if periodic_boundary:
				positions_x -= box_width * np.floor(positions_x/box_width)
				positions_y -= box_width * np.floor(positions_y/box_width)
			f.write(str(n_atoms)+'\n')
			f.write(' \n')
			for symbol, x, y in zip(symbols,positions_x,positions_y):
				f.write(symbol+str(x)+" "+str(y)+' 0.0\n')