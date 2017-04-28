#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from numpy import linalg as LA
import mdtraj as mdt

def trk_LD(in_opts
          ):

        # => Default Options <= #

        options = { 

        # > Filenames < #    

        # Trajectory
        'fn_traj' : 'positions.xtc',
        # Topology
        'fn_top' : 'water2.pdb',
        # Output file 
        'fn_out' : 'trk_L.out',
        # Position data of L defect
        'fn_L' : 'L_defect',

        # > MSD sampling options < #
        
        # Number of chunks to skip from the beginning of the trajectory
        'nchunk' : 10,
        # Periodic Boundary Conditions
        'periodic' : False,
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
        
        len_chunk = options['len_chunk']
        periodic = options['periodic']

        fn_traj  = options['fn_traj' ]  
        fn_top   = options['fn_top'  ]
        fn_out   = options['fn_out'  ]
        fn_L = options['fn_L']
        fn_D = options['fn_D']

	# => Load Info <= #
	
	top = mdt.load(fn_top).topology
	traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)

	# => Indices of O/H <= #
	
	atoms = [atom for atom in top.atoms]
	
	indsO = [ind for ind, val in enumerate(atoms) if val.name == 'O']
	indsH = [ind for ind, val in enumerate(atoms) if val.name == 'H']
        inds_all = [ind for ind, val in enumerate(atoms)]
	n_O = len(indsO)
	n_H = len(indsH)
	n_atoms = n_H + n_O
	
	pairs = []
	for ind_all in inds_all:
	    for indO in indsO:
                if ind_all != indO:
	            pairs.append((ind_all,indO))

	
        # => Initializing Counters <= #

        L_defect = []
        D_defect = []

	# => Loop over chunks <= #

	for chunk in traj:

	    distances = mdt.compute_distances(chunk,pairs,opt=True).reshape((-1,n_atoms-1,n_O))
            print n_O
            print np.shape(distances)
            exit()
	    Hneighbors  = np.argmin(distances, axis=2)
	
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


def test():

    trk_LD()    

