#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
from reference import quotes as q
import time
from datetime import date
import mdtraj as mdt

def rattle(
        fn_traj = 'positions.xtc',
        fn_top = 'water.pdb',
        fn_rattle = 'rattle_ostar',
        len_chunk = 10,
        timestep = 50,
        periodic = False,
        cutoff = 1.0, # ps
        cut=1E23,
        ):

	# => Load Info <= #
        fs2ps = 0.001
        cutoff /= timestep*fs2ps
	
	top = mdt.load(fn_top).topology
	traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)

	# => Indices of O/H <= #
	
	atoms = [atom for atom in top.atoms]
	nwaters = (len(atoms) - 1)/3 #Assumes only one proton defect! 
	
	indsO = [ind for ind, val in enumerate(atoms) if val.name == 'O']
	indsH = [ind for ind, val in enumerate(atoms) if val.name == 'H']
	n_O = len(indsO)
	n_H = len(indsH)
	
	pairs = []
	for indH in indsH:
	    for indO in indsO:
	        pairs.append((indH,indO))
	
        # => Initializing Counters <= #

        # Counts for tracking proton transfer events
        indsOstar = []
        Ostar = []

        # General Step Count
        step = 0
        nrattle = 0
        hop_count = 0
        prevOstar = 0

	# => Looping over Chunks <= #

	for chunk in traj:  

            # Establishing Periodic Boundary Conditions
            cell = chunk.unitcell_lengths

	    distances = mdt.compute_distances(chunk,pairs,opt=True,periodic=periodic).reshape((-1,n_H,n_O))
	    Hneighbors  = np.argmin(distances, axis=2) 
	    Oneighbors  = np.argmin(distances, axis=1) 
	
	    # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
	    # Loop over timesteps in each chunk
	    for tind in range(chunk.xyz.shape[0]):
	        
	        # Current frame  
	        # frame is np.array shape[nAtom,3]
	        frame = chunk.xyz[tind,:,:]
	 
	        # Find the Ostar
	        atomsO = frame[indsO,:]
	        atomsH = frame[indsH,:]
	
	        indOstar_local = [x for x in range(n_O) if list(Hneighbors[tind,:]).count(x) != 2][0]

	        if type(indOstar_local) != int:
	            print "Check trajectory! Wrong number of defects" 
	        indOstar = indsO[indOstar_local]
                if step < 2:
	            indsOstar.append(indOstar)
	            Ostar.append(frame[indOstar,:])
                    prevOstar = 1e23
                else: 
                    if periodic == True:
                        wrapped_Ostar[0] = np.mod((unwrapped_Ostar[0]/cell[tind,0]) , 1.0) * cell[tind,0]   
                        wrapped_Ostar[1] = np.mod((unwrapped_Ostar[1]/cell[tind,1]) , 1.0) * cell[tind,1] 
                        wrapped_Ostar[2] = np.mod((unwrapped_Ostar[2]/cell[tind,2]) , 1.0) * cell[tind,2]  
                    # No Hopping
                    if indOstar == indsOstar[-1]:
	                Ostar.append(frame[indOstar,:])
                        hop_count += 1
                    # Rattling
                    elif indOstar == prevOstar and hop_count < cutoff:
                        nrattle += 1
                        indsOstar[-1] = indOstar
                        Ostar[-hop_count:][:] = frame[indOstar,:]
	                Ostar.append(frame[indOstar,:])
                        hop_count = 0
                    # Transfer w/o Rattling
                    else:
	                Ostar.append(frame[indOstar,:])
                        prevOstar = indsOstar[-1]
                        hop_count += 1
                    
	            # Append to global target
	            indsOstar.append(indOstar)
                step += 1
                print "step: %d\r" % (step),
            if cut != 0.0:
                if step == cut:
                    break

        # => Dumping Data <= #        

        print "%d rattling events\n" % (nrattle)
	np.array(Ostar).dump(fn_rattle)

        q.h(9)

	return
