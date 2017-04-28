#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
import quotes as q
import time
from datetime import date
import mdtraj as mdt

def structure(
        in_opts = {},
        ):


        # => Default Options <= #

        options = { 

        # > Filenames < #    

        # Trajectory
        'fn_traj' : 'positions.xtc',
        # Topology
        'fn_top' : 'water2.pdb',
        # Output file 
        'fn_out' : 'EZ_structure.out',
        # Max delta coordinate
        'fn_max_delta' : 'max_delta',
        # Min delta coordinate
        'fn_min_delta' : 'min_delta',
        # Max O*-O distance for nearest neighbor
        'fn_max_dOO' : 'max_dOO',
        # Min O*-O dostamce for nearest neighbor
        'fn_min_dOO' : 'min_dOO',
        # O1-O*-O2 angle
        'fn_EZ_angle' : 'EZ_angle',

        # > MSD sampling options < #

        # Length of chunks in loading trajectory (steps)
        'len_chunk' : 10,
        # Timestep of trajectory (fs)
        'timestep' : 50,
        # Number of chunks to skip from the beginning of the trajectory
        'nchunk' : 10,
        }

        # => Override Default Options <= #

        for key,val in in_opts.items():
           if key not in options.keys():
              raise ValueError('%s option unavailable, RTFM' % (key))
           if options[key] != val:
              print "Default Option Overridden: %s = %s" % (key,str(val))
           options[key] = val
        print '\n'
       
    
        # => Declaring Variables from Options <= #
        
       # len_MSDps = options['len_MSDps']
        len_chunk = options['len_chunk']
        nchunk = options['nchunk']
        timestep = options['timestep']
        #interval_ps = options['interval_ps']

        fn_traj  = options['fn_traj' ]  
        fn_top   = options['fn_top'  ]
        fn_out   = options['fn_out'  ]
        fn_max_delta = options['fn_max_delta']
        fn_min_delta = options['fn_min_delta']
        fn_max_dOO = options['fn_max_dOO']
        fn_min_dOO = options['fn_min_dOO']
        fn_EZ_angle  = options['fn_EZ_angle']

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
        
        # => Declaring arrays <= #

        times = []	
        max_delta = [] 
        min_delta = [] 
        max_dOO = [] 
        min_dOO = [] 
        EZ_angle  = [] 
        
        step = 0

	# => Master Loop <= #

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

                    # The three delta coordinates
                    deltas = []
                    # Global indices for the three interacting Oxygens
                    O_inds = []

                    for P in ind3Ps:
                        #distance between potential  H-star and O-star
                        d1 = distances[tind,P,indOstar_local]

                        #all distances between potential H-star and all oxygens
                        HstarOs = distances[tind,P,:]
                    
                        #local index of 2nd nearest oxygen to potential H-star
                        indO2_local = np.argsort(HstarOs)[1]
                        O_inds.append(indsO[indO2_local])

                        #distance between potential H-star and 2nd nearest  oxygen
                        d2 = HstarOs[indO2_local]
                        delta_test = d2-d1
                        if delta_test < 0.0:
                            print "delta is negative!"
                            exit()
                        deltas.append(delta_test)

                    min_delta.append(min(deltas))
                    max_delta.append(max(deltas))
                    indO1 = O_inds[np.argsort(deltas)[0]]
                    indO2 = O_inds[np.argsort(deltas)[1]]
                    indO3 = O_inds[np.argsort(deltas)[2]]
                    O_neighbors = [indO1,indO2,indO3]

                    triple = [[indO1,indOstar,indO2]]
                    EZ_angle.append(math.degrees(mdt.compute_angles(chunk[tind],triple)))

                    dOOs = []
                    for O in O_neighbors:
                        pair = [[indOstar,O]]
                        dOOs.append(mdt.compute_distances(chunk[tind],pair))

                    min_dOO.append(min(dOOs))
                    max_dOO.append(max(dOOs))

                    step += 1
                    print "Step: %d" % (step)

        np.array(EZ_angle).dump(fn_EZ_angle)
        np.array(min_delta).dump(fn_min_delta)
        np.array(max_delta).dump(fn_max_delta)

        np.array(min_delta).dump(fn_min_dOO)
        np.array(max_delta).dump(fn_max_dOO)

        timestamp = "Date: %s\n" % (str(date.today()))
        info0 = "Options selected for simulation: \n"

	hs  = open(fn_out,"a")
        hs.write(info0)
        for key,val in options.items():
           opt = "%s = %s\n" % (key,str(val))
           hs.write(opt)
	hs.close() 
        q.h(6)

	return

