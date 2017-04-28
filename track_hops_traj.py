#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
from reference import quotes as q
import time
from datetime import date
import mdtraj as mdt

def drift_hops(
        fn_traj = 'positions.xtc',
        fn_top = 'water.pdb',
        fn_drift = 'drift_traj',
        fn_hops = 'hops_traj',
        fn_ostar = 'Ostar',
        fn_dhr = 'dhr_traj',
        fn_dr = 'dr_traj',
        len_chunk = 10,
        periodic = False,
        cut=1E23,
        cutoff = 0.3,
        wrapped=False,
        ):


	# => Load Info <= #
	
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

        # Count used for eliminating preliminary steps in trajectory
        nchunk = 0.0

        # General Step Count
        step = 0
        Ostars_drift = []
        Ostars_hops = []
        Ostars = []
        unwrapped_Ostars = []
        wrapped_Ostars = []
        indsOstars = []
        dhrs = []
        drs = []

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
                unwrapped_Ostar = frame[indOstar]
                wrapped_Ostar = [0.0,0.0,0.0]
                if step == 0:
                    Ostars_drift.append(unwrapped_Ostar)
                    Ostars_hops.append(unwrapped_Ostar)
                    Ostars.append(unwrapped_Ostar)
                    indsOstars.append(indOstar)
                else:    

                    if periodic == True:
                        wrapped_Ostar[0] = np.mod((unwrapped_Ostar[0]/cell[tind,0]) , 1.0) * cell[tind,0]   
                        wrapped_Ostar[1] = np.mod((unwrapped_Ostar[1]/cell[tind,1]) , 1.0) * cell[tind,1] 
                        wrapped_Ostar[2] = np.mod((unwrapped_Ostar[2]/cell[tind,2]) , 1.0) * cell[tind,2]  

                    # Drift movement of previous O*
                    if wrapped == False:
                        dx = frame[indsOstars[-1]][0]-unwrapped_Ostars[-1][0]
                        dy = frame[indsOstars[-1]][1]-unwrapped_Ostars[-1][1]
                        dz = frame[indsOstars[-1]][2]-unwrapped_Ostars[-1][2]                
                        dr = np.array([dx,dy,dz])
                    else:
                        dx = frame[indsOstars[-1]][0]-wrapped_Ostars[-1][0]
                        dy = frame[indsOstars[-1]][1]-wrapped_Ostars[-1][1]
                        dz = frame[indsOstars[-1]][2]-wrapped_Ostars[-1][2]                
                        dr = np.array([dx,dy,dz])

                        for ind,y in enumerate(dr):
                            if abs(y) > cell[tind,ind]/2.0:
                                if y < 0.0:
                                    dr[ind] += cell[tind,ind]
                                else:
                                    dr[ind] -= cell[tind,ind]

                    if indOstar != indsOstars[-1]:
                        # Hopping Movement without Drift Component
                        dhx = wrapped_Ostar[0]-wrapped_Ostars[-1][0]-dx
                        dhy = wrapped_Ostar[1]-wrapped_Ostars[-1][1]-dy
                        dhz = wrapped_Ostar[2]-wrapped_Ostars[-1][2]-dz
                       
                    else:
                        dhx = 0.0                  
                        dhy = 0.0
                        dhz = 0.0                 

                    dhr = np.array([dhx,dhy,dhz])

                    for ind,x in enumerate(dhr):
                        if abs(x) > cell[tind,ind]/2.0:
                            if x < 0.0:
                                dhr[ind] += cell[tind,ind]
                            else:
                                dhr[ind] -= cell[tind,ind]
        
                    # => Trouble Shooting Periodicity <= #

                    #if step > 10.0:
                    #    if abs(np.linalg.norm(np.array(dhr)+np.array(dr))) > 0.3:
                    #        print "Box Jump!\n"
                    #        ind = np.argmax(np.array(dhr))
                    #        if dhr < 0.0:
                    #            dhr[ind] = cell + dhr[ind]
                    #        else:
                    #            dhr[ind] = cell- dhr[ind]
                    #        print dhr
                    #        exit()

                        
                    dhrs.append(dhr)
                    drs.append(dr)

                    indsOstars.append(indOstar)
                    Ostars_drift.append(np.array(Ostars_drift[-1])+dr)
                    Ostars_hops.append(np.array(Ostars_hops[-1])+dhr)

                unwrapped_Ostars.append(np.array(unwrapped_Ostar))
                wrapped_Ostars.append(np.array(wrapped_Ostar))
                indsOstars.append(indOstar)

                step += 1
                print "step: %d\r" % (step),
            if step >= cut:
                break
            
        # => Dumping Data <= #        

	np.array(wrapped_Ostars).dump(fn_ostar)
	np.array(Ostars_hops).dump(fn_hops)
	np.array(Ostars_drift).dump(fn_drift)
	np.array(dhrs).dump(fn_dhr)
	np.array(drs).dump(fn_dr)

        q.h(9)

	return

    
