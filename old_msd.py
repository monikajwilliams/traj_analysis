#!/usr/bin/env python
import matplotlib
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA


def msd(
         # System Dimensionality
         nd=1,
	 axis=2,
        
         # Time Parameters

         # fs per step
	 timestep=50.0, 
	 interval_ps = 50.0,
         # ps to drop from beginning of simulation
	 drop_ps = 10.0,
         # length of msd segment
	 len_msd_ps = 100, 

         # Other
	 pDefect=True,
	 cell_len = 0.0,
	 mOrigins = True,
        
         # Filenames
	 fn_out = "msd_vals",
         fn_info = "msd_info.out",
	 fn_in = "Ostar",

         ):

	np.set_printoptions(threshold=10000)
	
	# => Unit Conversions <= #

	fs2ps = 0.001 

        # converting dropped part of trajectory into step index
	drop_ind = int(drop_ps/(timestep*fs2ps))
        # converting msd subtrajectory length into step length
	len_msd = int(len_msd_ps/(timestep*fs2ps))

	# => Loading Trajectory <= #

	Ostars = np.load(fn_in)
	print "Loaded Data"

        # Original length of simulation (steps)
	uncut_tmax = len(Ostars)
	# Length of simulation minus the dropped time
        times = np.arange(len(Ostars)-drop_ind)
	tmax = len(times)

        print "Total number of steps in simulation: %d" % (uncut_tmax)
        print "Total time of simulation: %2.2f ns" % (uncut_tmax*timestep*fs2ps/1000.0)

	if tmax > uncut_tmax:
		print "Length of simulation is shorter than length cut from beginning of simulation. Check for errors!"
		drop_ind = 0
		tmax = uncut_tmax
		times = np.arange(len(Ostars))

	if pDefect == True:
		Ostars = np.reshape(Ostars,(uncut_tmax,-1,3),order='F')

        # Eliminating position data in non-relevant dimensions
	if nd == 1:
	    for x in range(3):
	        if x != axis:
	            Ostars[:,:,x] *= 0.0
	if nd == 2:
	    for x in range(3):
	        if x != axis[0] and x != axis[1]:
	            Ostars[:,:,x] *= 0.0

	natoms = len(Ostars[0,:,0])

        # = > Multiple Origin Data <= #    
        	
	if mOrigins == True:
	    interval = int(interval_ps/(timestep*fs2ps)) 
	    print "interval = %d steps" % (interval)
	else:
	    interval = tmax 
	
	# => Calculating the Square Displacements <= # 

        # Looping over atoms
        # In case of a single proton defect natoms = 1

	atoms_origins = []
	b = 0
	for y in range(natoms):

            # Looping over trajectory
	    b += 1
	    origins = []
	    for x in range(drop_ind,tmax,interval): 
	        SDs = []
	        for Ostar in Ostars[x:x+len_msd,y]:
	            SDs.append((Ostar[0] - Ostars[x,y,0])**2 + 
	                       (Ostar[1] - Ostars[x,y,1])**2 + 
	                       (Ostar[2] - Ostars[x,y,2])**2) 
	        origins.append(SDs)
	    atoms_origins.append(origins)

	# => Averaging Square Displacements <= #

	print "Calculating MSD"
	atoms_origins = np.array(atoms_origins,ndmin=2)
	
	all_MSD = []
	for atom in atoms_origins:
	    n = 0
	    MSD = np.zeros_like(atoms_origins[0,0])
	    nsamples1 = 0
	    nsamples2 = 0
	    for origin in atom:
	        n += 1
	        nsamples = len(origin)
		if n == 0:
		   nsamples1 += nsamples
		if n == 1:
		   nsamples2 += nsamples
	        MSD[:nsamples] += origin
	    normalization = np.zeros_like(MSD)
	    nSDs = len(atom[0])
	    dSDs = nsamples1 - nsamples2

	    for i in range(n):
	        len(normalization)
	        normalization[0:(nSDs-(dSDs*(i)))] += 1
	    
	    nMSD = np.divide(MSD,normalization)
	    all_MSD.append(nMSD)

        # => Saving msd data <= # 

	all_MSD = np.array(all_MSD)
	all_MSD.dump(fn_out)

        # => Writing Analysis info <= #

	info = "Dropped Trajectory = %d\n " %(drop_ps)
	info1 = "Time interval between origins = %d\n " %(interval_ps)

	hs  = open(fn_info,"a")
	hs.write(info)
	hs.write(info1)
	hs.close() 
	
