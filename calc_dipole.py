#!/usr/bin/env python
import matplotlib
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as mdt

def dipole(
        fn_top='water.pdb',
        fn_charges='charges',
	fn_traj = 'positions.xtc',
	fn_info = "trkO_out.txt",
	fn_dipoles = "dipoles",
	fn_fields = "fields",
	len_chunk = 10,
        periodic=None,
        ):
	
	# => Load Info <= #
	
	charges = np.load(fn_charges)
	top = mdt.load(fn_top).topology
	traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)

	# => Indices of O/H <= #
	
	atoms = [atom for atom in top.atoms]
	nwaters = (len(atoms) - 1)/3 #Assumes only one proton defect! 
        natoms = len(atoms)
	
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
	indsOstar = []
	OstarZs = []
	Ostar = []
        dipoles = []
        fields = []
        step = 0

	# Loop over chunks

	for chunk in traj:
           
            # if the periodicity is left out of compute_distances
            # mdtraj determines if the unit cell information and periodicity 
            # of the system from the trajectory information
 
            if periodic == None:
	        distances = mdt.compute_distances(chunk,pairs,opt=True).reshape((-1,n_H,n_O))
            else:
	        distances = mdt.compute_distances(chunk,pairs,opt=True,periodic=periodic).reshape((-1,n_H,n_O))
                
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


                #local indices of all three hydrogens attached to O-star
                all_protons =  distances[tind,:,indOstar_local]
                ind3Ps = np.argsort(all_protons)[:3]
                indH1 = indsH[ind3Ps[0]]
                indH2 = indsH[ind3Ps[1]]
                indH3 = indsH[ind3Ps[2]]
	
	        # Append to global target
	        indsOstar.append(indOstar)
                Ostar = frame[indOstar,:]

                pairs_Ostar = []

                all_inds = [ind for ind,val in enumerate(atoms) if ind != indOstar and ind != indH1 and ind != indH2 and ind != indH3] 
                all_inds = np.array(all_inds)
               
                Ostar = np.array(Ostar)
                Ostar_distances = [atom_z-Ostar[2] for ind,atom_z in enumerate(frame[all_inds,2])]
                dipole = np.sum([val*charges[tind,all_inds[ind]] for ind,val in enumerate(Ostar_distances)])
                dipoles.append(dipole)

                field = np.sum([charges[tind,all_inds[ind]]*(func_field(Ostar[2],atom_z)) for ind,atom_z in enumerate(frame[all_inds,2])])
                fields.append(field)
                step += 1
                print "step: %d\r" % (step),

        # => Saving Ostar positions <= #

	np.array(dipoles).dump(fn_dipoles)
	np.array(fields).dump(fn_fields)

	return

def func_field(ostarz,atomz):

    value = (ostarz-atomz)/(abs(ostarz-atomz)**3)
        

    return value
