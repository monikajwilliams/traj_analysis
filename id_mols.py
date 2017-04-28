#!/usr/bin/env python -u
import time
import math
import numpy as np
import matplotlib.pyplot as plt
import mdtraj as mdt

def find_bonds(
    Zxyz,
    Rscale = 1.0,
    ):

    # true vdw radii in A
    # Z : radii
    R_map = {
        1 : 1.09,
        6 : 1.70, 
        7 : 1.55,
        8 : 1.52,
        9 : 1.47,
        11 : 2.27,
        15 : 1.80,
        16 : 1.80,
        17 : 1.75,
        }
    
    # scaled covalent radius
    R_map = { key : 0.5*Rscale*val for key,val in R_map.items()}
    
    # number of atoms 
    nP = Zxyz.shape[0]
    # array of vdw radii
    R = np.array([R_map[int(Zxyz[P,0])] for P in range(nP)])

    # the grid-spacing for hash
    H = 2.0*np.max(R)

    # the hash
    hash = {}
    for P in range(nP):
        key = (
            int(math.floor(Zxyz[P,1]/H)),
            int(math.floor(Zxyz[P,2]/H)),
            int(math.floor(Zxyz[P,3]/H)),
            )

        hash.setdefault(key,[]).append(P)
        
    # find the bonds    
    bonds = []
    for P in range(nP):
        RP = R[P]
        key = (
            int(math.floor(Zxyz[P,1]/H)),
            int(math.floor(Zxyz[P,2]/H)),
            int(math.floor(Zxyz[P,3]/H)),
            )

        D = 1
        for nx in range(-D,+D+1):
            for ny in range(-D,+D+1):
                for nz in range(-D,+D+1):

                        key2 = (
                            key[0] + nx, 
                            key[1] + ny, 
                            key[2] + nz, 
                            )

                        blist = hash.get(key2,None)
                        if blist == None: 
                            continue

                        for Q in blist:
                            if P == Q:
                                continue
                            RQ = R[Q]
                            RPQ = math.sqrt(
                                (Zxyz[P,1]-Zxyz[Q,1])**2 + 
                                (Zxyz[P,2]-Zxyz[Q,2])**2 + 
                                (Zxyz[P,3]-Zxyz[Q,3])**2 
                                )
                            if RPQ < RP + RQ:
                                bonds.append((P,Q))

    return bonds

def find_cycles(
    inds,
    bonds,
    ):

    bonds2 = {}
    for bond in bonds:
        bonds2.setdefault(bond[0],[]).append(bond[1])

    clusters = [] 
    clustered = set()
    for P in inds:
        if P in clustered:
            continue
        queue = [P]
        cluster = []
        while len(queue):
            Q = queue.pop(0)
            cluster.append(Q)
            if not Q in bonds2:
                continue
            for R in bonds2[Q]:
                if not R in cluster:
                    queue.append(R)
        clusters.append(cluster)
        for Q in cluster: 
            clustered.add(Q)

    return clusters
   
        
def find_frags(
    Zxyz,
    Rscale = 1.0,
    ):

    bonds = find_bonds(Zxyz,Rscale)
    inds = list(range(Zxyz.shape[0]))
    frags = find_cycles(inds,bonds)

    return frags

def get_unitcell(
    fn_top='water.pdb',
    fn_traj = 'positions.xtc',
    len_chunk = 1,
    cut=100,
    ):

    # => Load Info <= #
    
    top = mdt.load(fn_top).topology
    traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)

    cells = []
    step = 0
    for chunk in traj:
       
        # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
        # Loop over timesteps in each chunk
        for tind in range(chunk.xyz.shape[0]):
            #cells.append(chunk.unitcell_lengths)
            if step == 0:
                cells = chunk.unitcell_lengths
            else:
                cells = np.vstack((cells,chunk.unitcell_lengths))
            step += 1
        if step == cut:
            break

    return cells*10.0

def get_coords(
    fn_top='water.pdb',
    fn_traj = 'positions.xtc',
    len_chunk = 1,
    cut=100,
    ):
    
    # => Load Info <= #
    
    top = mdt.load(fn_top).topology
    traj = mdt.iterload(fn_traj, top=top, chunk=len_chunk)
    
    # => Indices of O/H <= #
    
    atoms = [atom for atom in top.atoms]
    ids = [val.name for val in atoms]
    atom_ids = {
    	'H': 1,
    	'O': 8,
    	'C': 6,
    	'S': 16,
    	'Na': 11,
        }
    
    Z = []
    for i in ids:    
        Z.append(atom_ids[i])
    Z =  np.array(Z).reshape(len(ids))
    
    # => Master Loop <= #
    step = 0
    
    # Loop over chunks
    Zxyz = []
    for chunk in traj:
       
        # chunk.xyz is np.array shape[nTinChunk,nAtom,3]
        # Loop over timesteps in each chunk
        for tind in range(chunk.xyz.shape[0]):
            Zxyz.append(np.column_stack((Z,chunk.xyz[tind,:,:]*10.0)))
            step += 1
        if step == cut:
            break
    
    return Zxyz

def make_whole(
    Zxyzs, # Takes only 1 frame
    cell,
    Rscale=1.0,
    ):

    nP = len(Zxyzs)
    cut = nP

    # cell dimensions
    #dims[0,0] = np.min(Zxyzs[:,1])
    #dims[1,0] = np.min(Zxyzs[:,2])
    #dims[2,0] = np.min(Zxyzs[:,3])

    dims = np.zeros((3,2))
    dims[0,0] = 0.0 
    dims[1,0] = 0.0 
    dims[2,0] = 0.0 
    dims[0,1] = dims[0,0] + cell[0]
    dims[1,1] = dims[1,0] + cell[1]
    dims[2,1] = dims[2,0] + cell[2]

    # replicate system
    px_Zxyzs = np.copy(Zxyzs[np.argsort(Zxyzs[:,1])[:cut]])   
    py_Zxyzs = np.copy(Zxyzs[np.argsort(Zxyzs[:,2])[:cut]])   
    pz_Zxyzs = np.copy(Zxyzs[np.argsort(Zxyzs[:,3])[:cut]])   
    nx_Zxyzs = np.copy(Zxyzs[np.argsort(Zxyzs[:,1])[-cut:]])  
    ny_Zxyzs = np.copy(Zxyzs[np.argsort(Zxyzs[:,2])[-cut:]])  
    nz_Zxyzs = np.copy(Zxyzs[np.argsort(Zxyzs[:,3])[-cut:]])  
    px_Zxyzs[:,1] += cell[0]
    py_Zxyzs[:,2] += cell[1]
    pz_Zxyzs[:,3] += cell[2]
    nx_Zxyzs[:,1] -= cell[0]
    ny_Zxyzs[:,2] -= cell[1]
    nz_Zxyzs[:,3] -= cell[2]

    stack_Zxyzs = np.vstack((
        Zxyzs,
        px_Zxyzs,
        py_Zxyzs,
        pz_Zxyzs,
        nx_Zxyzs,
        ny_Zxyzs,
        nz_Zxyzs,
        ))

    frags = find_frags(stack_Zxyzs,Rscale=Rscale)
    new_Zxyzs = np.zeros_like(Zxyzs)
    count = 0
    for ind,frag in enumerate(frags):
        zxyzs = stack_Zxyzs[frag]
        cxyz = com(zxyzs)
        if np.all(np.greater_equal(cxyz,dims[:,0])) and np.all(np.less_equal(cxyz,dims[:,1])):
            size = np.shape(zxyzs)[0]
            new_Zxyzs[count:count+size,:] += zxyzs
            count += size
    return new_Zxyzs

def com(
    Zxyzs,
    ):

    # Calculates the center of mass
    # for a system of atoms
    
    mass_data_ = {
        1  :  1.00782503207,
        2  :  4.00260325415,
        3  :  7.016004548,
        4  :  9.012182201,
        5  :  11.009305406,
        6  :  12.0000000,
        7  :  14.00307400478,
        8  :  15.99491461956,
        9  :  18.998403224,
        10 :  19.99244017542,
        11 :  22.98976928087,
        12 :  23.985041699,
        13 :  26.981538627,
        14 :  27.97692653246,
        15 :  30.973761629,
        16 :  31.972070999,
        17 :  34.968852682,
        18 :  39.96238312251,
        } 

    natoms = len(Zxyzs)
    tot_mass = 0.0
    x_mass = 0.0
    y_mass = 0.0
    z_mass = 0.0
   
    for ind,Zxyz in enumerate(Zxyzs):
        m = mass_data_[Zxyz[0]]
        tot_mass += m
        x_mass += m*Zxyz[1] 
        y_mass += m*Zxyz[2]        
        z_mass += m*Zxyz[3]
    
    x_mass /= tot_mass
    y_mass /= tot_mass
    z_mass /= tot_mass

    com = [x_mass,y_mass,z_mass]

    return com
    

def distances(atoms1,atoms2):

    xi, xj = np.meshgrid(atoms1[:,0],atoms2[:,0],indexing='ij')
    yi, yj = np.meshgrid(atoms1[:,1],atoms2[:,1],indexing='ij')
    zi, zj = np.meshgrid(atoms1[:,2],atoms2[:,2],indexing='ij')
    dx = xi - xj
    dy = yi - yj
    dz = zi - zj
    dr = np.sqrt(np.square(dx) + np.square(dy) + np.square(dz))
    
    return dr
                
def get_ostar(
    start = 2000,
    stop = 2010,
    fn_out = 'test.npy',        
    fn_top = 'water.pdb',
    fn_traj='positions.xtc',
    Rscale = 1.0,
    unwrapped=False,
    periodic=False,
    ):

    Ostars = []
    print "Loading Unit Cell Vectors"
    cells = get_unitcell(
        fn_top=fn_top,
        fn_traj=fn_traj,
        cut=stop,
        )
    print "Loading Coordinates"
    Zxyzs = get_coords(
        fn_top=fn_top,
        fn_traj=fn_traj,
        cut=stop,
        )
    indsO = [ind for ind,val in enumerate(Zxyzs[0][:,0]) if val == 8] 
    nO = len(indsO)
    dstep = stop-start
    clusters = []
    
    for step in range(start,stop,1):
        Zxyz = Zxyzs[step]
        cell = cells[step]
        if unwrapped == True:
           Zxyz[:,1:] = np.mod(Zxyz[:,1:] / cell,1) * cell 
        if periodic == True:
            Zxyz = make_whole(Zxyz,cell,Rscale=Rscale)
        mols = find_frags(Zxyz,Rscale=Rscale)
        hydronium = np.array([val for ind,val in enumerate(mols) if len(val) == 4])
        if np.shape(hydronium)[0] != 0:
            hydronium = hydronium[0]
            Ostar = np.array(Zxyz[[hydronium[ind] for ind,val in enumerate(hydronium) if Zxyz[hydronium[ind],0] == 8][0],:])
            Ostars.append(Ostar)

        # in the case that the excess proton is not found in standard
        # hydronium format
        else:
            Hind = np.array([val for ind,val in enumerate(mols) if len(val) != 3])

            # if the excess proton is a single atom
            if np.shape(Hind)[1] == 1:
                atoms2 = Zxyz[indsO,1:]
                atoms1 = Zxyz[Hind,1:][0]
                dr = distances(atoms1,atoms2)
                Ostar = np.array(Zxyz[indsO[np.argmin(dr)],:])
                Ostars.append(Ostar)

            # if the excess proton resides in a larger structure
            else:
                hs = Zxyz[[Hind[0,ind] for ind,val in enumerate(Hind[0,:]) if Zxyz[Hind[0,ind],0] == 1],:]
                os = Zxyz[[Hind[0,ind] for ind,val in enumerate(Hind[0,:]) if Zxyz[Hind[0,ind],0] == 8],:]
                dr = distances(hs[:,1:],os[:,1:]) 
                hstar = np.argmin(dr[:,0]+dr[:,1])
                ostar = np.argmin(dr[hstar,:])
                Ostar = np.array(os[ostar,:])
                Ostars.append(Ostar)
                clusters.append(Zxyz[Hind][0])
        print '%d/%d\r' % (step-start,dstep),
        

    np.array(Ostars)[:,1:].dump(fn_out)
    np.array(clusters).dump('clusters.npy')

    return 

