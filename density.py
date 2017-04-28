#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as mdt


def vol_dens(
            fn_traj = 'positions.xtc',
            fn_top = 'water2.pdb',
            len_chunk = 10,
            periodic = False,
            ):
            
	# => Load Info <= #
	
	top = mdt.load(fn_top).topology
	traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)

	# => Indices of O/H <= #
	
	atoms = [atom for atom in top.atoms]
	nwaters = (len(atoms) - 1)/3 #Assumes only one proton defect! 
	
	indsO = [ind for ind, val in enumerate(atoms) if val.name == 'O']
	n_O = len(indsO)
	
	pairs = []
	for ind,indO1 in enumerate(indsO):
	    for indO2 in indsO[ind+1:]:
	        pairs.append((indO1,indO2))


	# => Master Loop <= #
	
	times = []
        nchunk = 0.0
        step = 0

	# Loop over chunks
	for chunk in traj:  
            nchunk += 1
            if nchunk < 10.0:
                continue
            else:

	        # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
	        # Loop over timesteps in each chunk
	        for tind in range(chunk.xyz.shape[0]):
	            
	            times.append(chunk.time[tind])
	        
	            # Current frame  
	            # frame is np.array shape[nAtom,3]
	            frame = chunk.xyz[tind,:,:]
	 
	            # Find the Ostar
	            atomsO = frame[indsO,:]
                    nwaters = float(len(atomsO))

                    min_zO = np.min(atomsO[:,2])
                    max_zO = np.max(atomsO[:,2])
                    h = max_zO-min_zO

                    min_xO = np.min(atomsO[:,0])
                    max_xO = np.max(atomsO[:,0])
                    dist_x = max_xO-min_xO

                    r = max_xO/2.0
                    vol_dens = nwaters/(math.pi*(r**2.0)*h)
                    #jacobian = (4.0/3.0)*math.pi*(r**2)
                    jacobian = (4.0/3.0)*math.pi*r
                    lin_dens = nwaters/h

                    
                    break
                    
                break
	
	return lin_dens,1.0
