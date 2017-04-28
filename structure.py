#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
import quotes as q
import mdtraj as mdt

def structure(nwaters=100,
	fn_traj = 'positions.xtc',
        fn_top='water2.pdb',
	infoname = "info.txt",
	filename1 = "Ostar",
	filename2 = "Pstar",
	filename3 = "delta",
	filename4 = "tau",
	filename5 = "dOO",
	filename6 = "angles",
        dump_Ostar   = False,  
        dump_Pstar   = False,
        dump_dOO     = False,
        dump_angles  = False,
        dump_tau     = False, 
        dump_delta   = False,
	len_chunk = 10,
        nchunk = 10, #number of chunks to skip
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

        pairsOO = []
	for indO in indsO:
	    for indO2 in indsO:
	        pairsOO.append((indO,indO2))
	# => Master Loop <= #
	
	times = []
	indsOstar = []
	indsPstar = []
	Ostar = []
	Pstar = []
        delta = []
        tau = []
        dOO = []
        angle = []
	# Loop over chunks
        nchunk = 0.0
	for chunk in traj:  
            nchunk += 1
            if nchunk < 10.0:
                continue
            else:
	        distances = mdt.compute_distances(chunk,pairs,opt=True).reshape((-1,n_H,n_O))
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
                    
                    #Determing H-star based on minimum delta-coordinate with 2nd nearest oxygen
                    
                    all_protons =  distances[tind,:,indOstar_local]

                    #local indices of all three hydrogens attached to O-star
                    ind3Ps = np.argsort(all_protons)[:3]

                    deltas = []
                    for P in ind3Ps:
                        #distance between potential  H-star and O-star
                        d1 = distances[tind,P,indOstar_local]

                        #all distances between potential H-star and all oxygens
                        HstarOs = distances[tind,P,:]
                    
                        #local index of 2nd nearest oxygen to potential H-star
                        indO2_local = np.argsort(HstarOs)[1]

                        #distance between potential H-star and 2nd nearest  oxygen
                        d2 = HstarOs[indO2_local]
                        delta_test = d2-d1
                        if delta_test < 0.0:
                            print "delta is negative!"
                            exit()
                        deltas.append(delta_test)

                    delta_local = min(deltas)

                    #local index of H-star
                    indPstar_local = ind3Ps[np.argmin(deltas)]

                    #all distances between H-star and all oxygens
                    HstarOs = distances[tind,indPstar_local,:]
                    
                    #local index of 2nd nearest oxygen to potential H-star
                    indO2_local = np.argsort(HstarOs)[1]
                    indO2 = indsO[indO2_local] 
                    pairsOO = [[indO2,indOstar]]
                    distOO2 = mdt.compute_distances(chunk[tind],pairsOO)
                    tau_local = float(delta_local/distOO2)

                    indPstar = indsH[indPstar_local]	
                    triple = [[indOstar,indPstar,indO2]]
                    angle_local = math.degrees(mdt.compute_angles(chunk[tind],triple))

	            # Append to global target
	            indsOstar.append(indOstar)
                    angle.append(angle_local)
	            Ostar.append(frame[indOstar,:])
	            Pstar.append(frame[indPstar,:])

                    np.append(delta,delta_local)
                    tau.append(tau_local)
                    dOO.append(distOO2)
                    
        if dump_Ostar == True:                
	    np.array(Ostar).dump(filename1)

        if dump_Pstar == True:                
	    np.array(Pstar).dump(filename2)

        if dump_delta == True:                
	    np.array(delta).dump(filename3)

        if dump_tau == True:                
	    np.array(tau).dump(filename4)

        if dump_dOO == True:                
	    np.array(dOO).dump(filename5)

        if dump_angles == True:                
	    np.array(angle).dump(filename6)

	info = "Number of water molecules = %d \n" % (nwaters)
	hs  = open(infoname,"a")
	hs.write(info)
	hs.close() 
        q.h(6)
	return

