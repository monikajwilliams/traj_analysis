#!/usr/bin/env python
import matplotlib
import numpy as np
import time
import matplotlib.pyplot as plt
import mdtraj as mdt

def trkO(
        fn_top='water.pdb',
	fn_traj = 'positions.xtc',
	fn_info = "trkO_out.txt",
	fn_out = "Ostar",
	len_chunk = 10,
        periodic=None,
        cut=1E23,
        ):
	
	# => Load Info <= #
        old_ostar = np.load('old_Ostar')	
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

	# => Master Loop <= #
	
	indsOstar = []
	OstarZs = []
	Ostar = []
        step = 0

	# Loop over chunks

	for chunk in traj:
           
            # if the periodicity is left out of compute_distances
            # mdtraj determines if the unit cell information and periodicity 
            # of the system from the trajectory information
            # however, usually assumes periodic!!! 
 
            if periodic == None:
	        distances = mdt.compute_distances(chunk,pairs,opt=True).reshape((-1,n_H,n_O))
            else:
	        distances = mdt.compute_distances(chunk,pairs,opt=True,periodic=periodic).reshape((-1,n_H,n_O))
                
	    Hneighbors  = np.argmin(distances, axis=2)
	
	    # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
	    # Loop over timesteps in each chunk
	    for tind in range(chunk.xyz.shape[0]):
	        
#	        times.append(chunk.time[tind])
	    
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
	
	        # Append to global target
	        indsOstar.append(indOstar)
	        Ostar.append(frame[indOstar,:])

#                print "step: %d\r" % (step),
                print indOstar
                if old_ostar[step,0] != frame[indOstar,0]:
                    print old_ostar[step]
                    print frame[indOstar,:]
                if old_ostar[step,1] != frame[indOstar,1]:
                    print old_ostar[step]
                    print frame[indOstar,:]
                if old_ostar[step,2] != frame[indOstar,2]:
                    print old_ostar[step]
                    print frame[indOstar,:]
                step += 1
            if cut != 0.0:
                if step == cut:
                    break

        # => Saving Ostar positions <= #

	np.array(Ostar).dump(fn_out)
	info = "Number of water molecules = %d \n" % (nwaters)
	hs  = open(fn_info,"a")
	hs.write(info)
	hs.close() 

