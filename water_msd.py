#!/usr/bin/env python
import matplotlib
import os,sys
import numpy as np
import matplotlib.pyplot as plt
from traj_analysis import id_mols as im
from numpy import linalg as LA


def msd(
   # System Dimensionality
   nd=3,
   axis=2,
   
   # Time Parameters

   # fs per step
   timestep=5.0, 
   interval_ps = 1.00,
   # ps to drop from beginning of simulation
   drop_ps = 5.0,
   # length of msd segment
   len_msd_ps = 20.0, 

   # Other
   cell_len = 0.0,
   mOrigins = True,
   
   # Filenames
   fn_out = "O_msds",
   fn_info = "msd_O.out",
   fn_top = 'water.pdb',
   fn_traj='positions.xtc',

   cut = 1E23,
   len_chunk = 10,
   ):

   np.set_printoptions(threshold=10000)
   # => Unit Conversions <= #
   
   fs2ps = 0.001 
   ps2ns = 1.0/1000.0
   
   # converting dropped part of trajectory into step index
   drop_ind = int(drop_ps/(timestep*fs2ps))
   # converting msd subtrajectory length into step length
   len_msd = int(len_msd_ps/(timestep*fs2ps))

   atoms = np.array(im.get_coords(
       fn_top=fn_top,
       fn_traj = fn_traj,
       len_chunk = len_chunk,
       cut=cut,
       ))

   print "loaded trajectory"

   indsO = [ind for ind, val in enumerate(atoms[0,:,0]) if val == 8]
   Os = atoms[:,indsO,1:]
   
   # Original length of simulation (steps)
   tot_nsteps = len(Os)
   # Length of simulation minus the dropped time
   times = np.arange(len(Os)-drop_ind)
   tmax = len(times)
   
   print "Total number of steps in simulation: %d" % (tot_nsteps)
   print "Total time of simulation: %2.2f ns" % (tot_nsteps*timestep*fs2ps*ps2ns)
   
   if tmax > tot_nsteps:
   	print "Length of simulation is shorter than length cut from beginning of simulation. Check for errors!"
   	drop_ind = 0
   	tmax = tot_nsteps
   	times = np.arange(len(Os))
   
   # Eliminating position data in non-relevant dimensions
   if nd == 1:
       for x in range(3):
           if x != axis:
               Os[:,:,x] *= 0.0
   if nd == 2:
       for x in range(3):
           if x != axis[0] and x != axis[1]:
               Os[:,:,x] *= 0.0
   
   natoms = len(Os[0,:,0])
   
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
   print "looping over atoms"
   for y in range(natoms):
   
       # Looping over trajectory
       b += 1
       origins = []
       for x in range(drop_ind,tmax,interval): 
           SDs = []
           for O in Os[x:x+len_msd,y]:
               SDs.append((O[0] - Os[x,y,0])**2 + 
                          (O[1] - Os[x,y,1])**2 + 
                          (O[2] - Os[x,y,2])**2) 
           origins.append(SDs)
       atoms_origins.append(origins)
   
   # => Averaging Square Displacements <= #
   
   print "Calculating MSD"
   atoms_origins = np.array(atoms_origins,ndmin=2)
   print np.shape(atoms_origins)
   
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

