#!/usr/bin/env python
import matplotlib
import math
import numpy as np
import matplotlib.pyplot as plt
import quotes as q
import time
from datetime import date
import mdtraj as mdt

def msd_hops(
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
        'fn_out' : 'msd_hops.out',
        # MSD data of defecting hopping
        'fn_hops' : 'hops',
        # MSD data of water drift
        'fn_drift' : 'drift',
        # MSD data of transfers
        'fn_transfers' : 'transfers',

        # > MSD sampling options < #
        
        # Length of MSD subtrajectory (ps)
        'len_MSDps' :  100,
        # Interval between MSD subtrajectories (ps) 
        'interval_ps' : 100,
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
        
        len_MSDps = options['len_MSDps']
        len_chunk = options['len_chunk']
        timestep = options['timestep']
        interval_ps = options['interval_ps']

        fn_traj  = options['fn_traj' ]  
        fn_top   = options['fn_top'  ]
        fn_out   = options['fn_out'  ]
        fn_hops  = options['fn_hops' ]
        fn_drift = options['fn_drift']

        if interval_ps < len_MSDps:
            print "WARNING: Length of subtraj is longer than interval, samples will be correlated\n"

        # => Unit Conversions from ps to steps <= #
	fs2ps = 0.001 
	len_MSD = int(len_MSDps/(timestep*fs2ps))
        interval = int((interval_ps/(fs2ps*timestep)))

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
	
	times = []
        transfer_steps = []
        all_msd_hops = np.zeros(len_MSD)
        all_msd_dOs = np.zeros(len_MSD)
        originO = []
        origindO = []
        dhops = []
        dOs = []

        # => Initializing Counters <= #

        # Count used for eliminating preliminary steps in trajectory
        nchunk = 0.0

        # Counts for tracking proton transfer events
        Ochange_count = 0
        prevOstar = 0

        # Counts for MSD averaging
        interval_step_dO = 0
        interval_step_hops = 0
        seg_count_dOs = 0
        seg_count_hops = 0

        # General Step Count
        step = 0

	# => Looping over Chunks <= #

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
                    Ostar = frame[indOstar]
    
                    # Update Counts
                    if prevOstar == 0:
                        prevOstar = indOstar
                        originO = Ostar
                        oldOstar = Ostar 
                        origindO = Ostar
    
                    elif prevOstar != indOstar:
                        if interval_step_hops < len_MSD: 
                            dH = (Ostar[0]-originO[0])**2 + (Ostar[1]-originO[1])**2 + (Ostar[2]-originO[2])**2  
                            print "<=====(TRANSFER)=====>"

                            dhops.append(dH)
    
                            prevOstar = indOstar
                            Ochange_count += 1
                            interval_step_hops += 1
                            origindO[0] += (Ostar[0] - oldOstar[0])
                            origindO[1] += (Ostar[1] - oldOstar[1])
                            origindO[2] += (Ostar[2] - oldOstar[2])
                           
                            transfer_steps.append(1) 

                        else:
                            originO = Ostar
                            all_msd_hops += dhops

                            dhops = []
                            seg_count_hops += 1
                            interval_step_hops = 0
    
                    else:
                        if interval_step_dO < len_MSD: 
                            dO = (Ostar[0]-origindO[0])**2 + (Ostar[1]-origindO[1])**2 + (Ostar[2]-origindO[2])**2  
                            print dO
                            dOs.append(dO)

                            interval_step_dO += 1

                        else:
                            origindO = Ostar
                            all_msd_dOs += dOs

                            dOs = []
                            seg_count_dOs += 1
                            interval_step_dO = 0

                        transfer_steps.append(0) 

                    step += 1
                    oldOstar = Ostar
                    print "step: %d" % (step) 

        # => Normalizing MSDs <= #

        all_msd_hops /= seg_count_hops
        all_msd_dOs /= seg_count_dOs
        
        # => Dumping Data <= #        

	np.array(all_msd_hops).dump(fn_hops)
	np.array(all_msd_dOs).dump(fn_drift)
	np.array(transfer_steps).dump(fn_transfers)

        # => Writing Output <=# 

        timestamp = "Date: %s\n" % (str(date.today()))
        info0 = "Calculation of MSD for proton hops and water drift\n"
	info1 = "Number of Ostar transfers = %d \n" % (Ochange_count)
	info2 = "Number of total steps = %d \n \n" %     (step)
        info3 = "Options selected for simulation: \n"

	hs  = open(fn_out,"a")
        hs.write("<=====(msd_hops)=====>")
	hs.write(timestamp)
	hs.write(info0)
	hs.write(info1)
	hs.write(info2)
	hs.write(info3)
        for key,val in options.items():
           opt = "%s = %s\n" % (key,str(val))
           hs.write(opt)
        hs.write("\n")
	hs.close() 
        q.h(9)
	return

def test():
    
    options = {
        'len_MSDps' :  15,
        'interval_ps' : 15,
        }
    msd_hops(options)

#test()
    
