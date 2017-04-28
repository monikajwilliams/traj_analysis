#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as mdt

# Calculating the radial distribution functions for various 
# water system pairs, including:
#           ostar-O
#           ostar-H
#           Hstar-O
#           O-O
#           O-H
#           H-H

def ostar_o(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	fn_out = "Ostar-O_dists",
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
        periodic=False,
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

	# => Master Loop <= #
	
	times = []
        OOdists = []
        nchunk = 0.0
        step = 0

	# Loop over chunks
	for chunk in traj:  
            nchunk += 1
            if nchunk < 10.0:
                continue
            else:
	        distances = mdt.compute_distances(chunk,pairs,opt=True,periodic=periodic).reshape((-1,n_H,n_O))
	        Hneighbors  = np.argmin(distances, axis=2) 
	        Oneighbors  = np.argmin(distances, axis=1) 
	
	        # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
	        # Loop over timesteps in each chunk
	        for tind in range(chunk.xyz.shape[0]):
	            
	            times.append(chunk.time[tind])
	        
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

                    # Calculating RDF

                    Ostar_O_inds = []
                    for indO in indsO:
                        Ostar_O_inds.append((indOstar,indO))
                        
                    OOdists.append(mdt.compute_distances(chunk[tind],Ostar_O_inds,periodic=periodic))
        
                    step += 1
                    print "step: %d\r" % (step),

	np.array(OOdists).dump(fn_out)

	return

def ostar_h(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	fn_out = "Ostar-H_dists",
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
        periodic = False,
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

	# => Master Loop <= #
	
	times = []
        OHdists = []
        nchunk = 0.0
        step = 0

	# Loop over chunks
	for chunk in traj:  
            nchunk += 1
            if nchunk < 10.0:
                continue
            else:
	        distances = mdt.compute_distances(chunk,pairs,opt=True,periodic=periodic).reshape((-1,n_H,n_O))
	        Hneighbors  = np.argmin(distances, axis=2) 
	        Oneighbors  = np.argmin(distances, axis=1) 
	
	        # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
	        # Loop over timesteps in each chunk
	        for tind in range(chunk.xyz.shape[0]):
	            
	            times.append(chunk.time[tind])
	        
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

                    # Calculating RDF

                    Ostar_H_inds = []
                    for indH in indsH:
                        Ostar_H_inds.append((indOstar,indH))
                        
                    OHdists.append(mdt.compute_distances(chunk[tind],Ostar_H_inds,periodic=periodic))
        
                    step += 1
                    print "step: %d\r" % (step),

	np.array(OHdists).dump(fn_out)

	return

def O_O(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	fn_out = "O-O_dists",
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
        periodic=False,
        cut = 0.0
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
	for ind,indO1 in enumerate(indsO):
	    for indO2 in indsO[ind+1:]:
	        pairs.append((indO1,indO2))

	# => Master Loop <= #
	
	times = []
        OOdists = []
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
	 
                    # Calculating RDF
                        
                    OOdists.append(mdt.compute_distances(chunk[tind],pairs,opt=True,periodic=periodic))
                    step += 1
                    print "step: %d\r" % (step),
            if cut != 0.0:
                if step == cut:
                    break

	np.array(OOdists).dump(fn_out)

	return

def O_H(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	fn_out = "O-H_dists",
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
        periodic=False,
        cut = 0.0,
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

	# => Master Loop <= #
	
	times = []
        OHdists = []
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
	 
                    # Calculating RDF
                        
                    OHdists.append(mdt.compute_distances(chunk[tind],pairs,opt=True,periodic=periodic))
                    step += 1
                    print "step: %d\r" % (step),
            if cut != 0.0:
                if step == cut:
                    break

	np.array(OHdists).dump(fn_out)

	return

def H_H(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	fn_out = "H-H_dists",
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
        periodic=False,
        cut = 0.0,
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
	for indH1 in indsH:
	    for indH2 in indsH:
                if indH1 != indH2: 
	            pairs.append((indH1,indH2))
                else:
                    continue

	# => Master Loop <= #
	
	times = []
        HHdists = []
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
	 
                    # Calculating RDF
                        
                    HHdists.append(mdt.compute_distances(chunk[tind],pairs,opt=True,periodic=periodic))
                    step += 1
                    print "step: %d\r" % (step),
            if cut != 0.0:
                if step == cut:
                    break

	np.array(HHdists).dump(fn_out)

	return

def hstar_o(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	fn_out = "Hstar-O_dists",
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
        periodic=False,
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

	# => Master Loop <= #
	
	times = []
        HstarOdists = []
        nchunk = 0.0
        step = 0

	# Loop over chunks
	for chunk in traj:  
            nchunk += 1
            if nchunk < 10.0:
                continue
            else:
	        distances = mdt.compute_distances(chunk,pairs,opt=True,periodic=periodic).reshape((-1,n_H,n_O))
	        Hneighbors  = np.argmin(distances, axis=2) 
	        Oneighbors  = np.argmin(distances, axis=1) 
	
	        # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
	        # Loop over timesteps in each chunk
	        for tind in range(chunk.xyz.shape[0]):
	            
	            times.append(chunk.time[tind])
	        
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
                    
                    # => Determing H-star based on minimum delta-coordinate with 2nd nearest oxygen <= #
                    all_protons =  distances[tind,:,indOstar_local]

                    # local indices of all three hydrogens attached to O-star
                    ind3Ps = np.argsort(all_protons)[:3]

                    deltas = []
                    for P in ind3Ps:
                        # distance between potential  H-star and O-star
                        d1 = distances[tind,P,indOstar_local]

                        # all distances between potential H-star and all oxygens
                        HstarOs = distances[tind,P,:]
                    
                        #local index of 2nd nearest oxygen to potential H-star
                        indO2_local = np.argsort(HstarOs)[1]

                        # distance between potential H-star and 2nd nearest  oxygen
                        d2 = HstarOs[indO2_local]
                        delta_test = d2-d1
                        if delta_test < 0.0:
                            print "delta is negative!"
                            exit()
                        deltas.append(delta_test)

                    delta_local = min(deltas)

                    # Local index of H-star
                    indPstar_local = ind3Ps[np.argmin(deltas)]

                    # Distances between H-star and all oxygens
                    HstarOs = distances[tind,indPstar_local,:]
                    HstarOdists.append(HstarOs)

                    step += 1
                    print "step: %d\r" % (step),

	np.array(HstarOdists).dump(fn_out)

	return

