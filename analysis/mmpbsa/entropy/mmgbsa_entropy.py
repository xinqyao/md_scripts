#! /usr/bin/env python

import os, sys, glob, pickle, gzip
import numpy as np
from numpy import arange, dtype
import scipy.stats
from scipy.stats import spearmanr, kendalltau, gaussian_kde
import math
import itertools
import pymbar

BETA = 1./(8.3144621E-3/4.184*300) # mol/kcal
#rmsd_filters = [3.0, np.inf]

if not 'system_dir' in locals().keys():
    system_dir = '/home/cli78/MM_PBSA/bromodomains/'

#########################
# Parse input arguments #
#########################

import argparse
parser = argparse.ArgumentParser()
parser.add_argument('--reference_data', \
    default=system_dir + 'experimental.dat', \
    help='File with reference free energies')

parser.add_argument('--sampler_dir', \
    default=system_dir, \
    help='Directory with sampler energies')
parser.add_argument('--sampler_total_E_FN', \
    default='trajectory.log', \
    help='File name for sampler complex energies')

parser.add_argument('--nm_dir', \
    default=system_dir, \
    help='Directory with normal modes analysis')
parser.add_argument('--nm_RL_FN', \
    default='complex_nm.out', \
    help='File name for normal modes of the complex')
parser.add_argument('--nm_R_FN', \
    default='receptor_nm.out', \
    help='File name for normal modes of the receptor')
parser.add_argument('--nm_L_FN', \
    default='ligand_nm.out', \
    help='File name for normal modes of the ligand')

parser.add_argument('--GBSA_dir', \
    default=system_dir, \
    help='Directory with GBSA energies')
parser.add_argument('--GBSA_RL_FN', \
    default='trajectory.log', \
    help='File name for GBSA complex energies')
parser.add_argument('--GBSA_R_FN', \
    default='protein.dat', \
    help='File name for GBSA receptor energies')
parser.add_argument('--GBSA_L_FN', \
    default='ligand.dat', \
    help='File name for GBSA ligand energies')

parser.add_argument('--ligand_dir',
    default=system_dir, \
    help='Directory with ligand RMSDs, prmtop, and dcd files')
parser.add_argument('--ligand_RMSD_FN', \
    default='ligand_rmsd.dat', help='File name for ligand rmsds')
parser.add_argument('--ligand_prmtop_FN', \
    default='vacuum.prmtop', help='File name for ligand prmtop')
parser.add_argument('--ligand_dcd_FN', \
    default='ligand_aligned.dcd', \
    help='Trajectory of the ligand from the simulation of the complex')

parser.add_argument('--skip_analysis', \
    action='store_true', \
    help='Does not do analysis')
parser.add_argument('-f', help='Dummy argument for ipython')

args = parser.parse_args()
del argparse

#######################
# Read reference data #
#######################

#F = open(args.reference_data)
#dat = F.readlines()
#F.close()
#
#ref_dG = {}
#for line in dat:
#    (id, dG) = line.split()
#    ref_dG[id] = float(dG)
#
#ids ref_dG.keys()
#del dat, F

ids = ['sys1']

#############
# Load data #
#############

def load_potential_energy(FN):
    """
    Loads a file with potential energies, returning them in kcal/mol
    """
    if os.path.isfile(FN):
        F = open(FN, 'r')
        dat = F.read().strip()
        F.close()
        if FN.endswith('.log'):
            lines = dat.split('\n')
            return np.array([float(line.split('\t')[3]) for line in lines[1:]])/4.184
        elif FN.endswith('.dat'):
            lines = dat.split('\n')
#            return np.array([float(line) for line in lines])/4.184
            return np.array([float(line) for line in lines])
        elif FN.endswith('.0'):
            frames = dat.split('Processing frame')[1:]
            energies = np.array(\
                [f.replace('1-4 ','1-4_')[f.find('BOND'):].split()[2:35:3] \
                for f in frames],dtype=float)
            # return np.sum(energies,1)/4.184
            return np.sum(energies[:,[3,4,5,9,10]],1)/4.184
        elif FN.endswith('.out'):
            # The normal modes output is in cal/mol/K
            # The conversion is to kcal/mol at 298.15 K
            if dat=='':
                print 'Empty file '+FN
                return np.array(np.nan)
            entropies = np.array([d.split()[-1] for d in dat.split('\n')], \
                dtype=float)*298.15/1000.
            return entropies
            pass
        else:
            raise Exception('Extension not recognized')
    else:
        print 'Missing file '+FN
        return np.array(np.nan)

#def load_RMSD(FN, nvals):
#    """
#    Loads a file with ligand RMSDs, in Angstroms
#    """
#    if os.path.isfile(FN):
#        F = open(FN, 'r')
#        dat = F.read().strip()
#        F.close()
#        return np.array([np.float(d.split()[1]) for d in dat.split('\n')[1:]])
#    else:
#        print 'Missing file ' + FN
#        return np.zeros(nvals)

def get_external_dof(traj):
    '''
    Calculate external degrees of freedom:
    (1) the center of mass and
    (2) the euler angle based on the principle axes of rotation
    '''
    import mdtraj
    com = mdtraj.compute_center_of_mass(traj)
    inertial_tensors = mdtraj.geometry.compute_inertia_tensor(traj)
    
    euler_angles = []
    # Principal axes are eigenvectors of the interial tensors
    # This makes the first principal axis Z, second Y, and third X
    for inertial_tensor in inertial_tensors:
        [w,v] = np.linalg.eig(inertial_tensor)
        (X,Y,Z) = [v[:,i] for i in np.argsort(w)]
        # Euler angles are determined from principal axes
        # https://en.wikipedia.org/wiki/Euler_angles
        euler_angles.append([np.arctan2(Z[0],Z[1]), \
                             np.arccos(Z[2]), \
                             np.arctan2(X[2],Y[2])])
    euler_angles = np.array(euler_angles)
    return (com, euler_angles)

Es = dict([(id,{}) for id in ids])

for id in ids:
#    Es[id]['sampler'] = load_potential_energy(\
#        os.path.join(args.sampler_dir,id,args.sampler_total_E_FN))
#    for solvent in ['GBSA','nm']:
    for solvent in ['GBSA']:
        dir = os.path.join(getattr(args,solvent+'_dir'),id)
        for moeity in ['L','R','RL']:
            Es[id]['%s_%s'%(solvent,moeity)] = load_potential_energy(\
                    os.path.join(dir, getattr(args,'%s_%s_FN'%(solvent,moeity))))
#    Es[id]['RMSD'] = load_RMSD(\
#        os.path.join(args.ligand_dir,id,args.ligand_RMSD_FN), \
#        len(Es[id]['sampler']))
    if np.any([np.any(np.isnan(Es[id][key])) for key in Es[id].keys()]):
        print 'Removing '+id
        del Es[id]
        ids = [item for item in ids if not item==id]

import mdtraj
from netCDF4 import Dataset
import sys
sys.stdout.write('Getting external degrees of freedom for')

if not os.path.isdir('external_dof'):
    os.mkdir('external_dof')
for id in ids:
    sys.stdout.write(' ' + id)

    loaded = False
    FN = os.path.join('external_dof',id + '.nc')
    if os.path.isfile(FN):
        try:
            nc = Dataset(FN,'r')
            Es[id]['com'] = nc.variables['com'][:]
            Es[id]['euler_angles'] = nc.variables['euler_angles'][:]
            nc.close()
            loaded = True
        except:
            loaded = False
    if not loaded:
        traj_FN = os.path.join(args.ligand_dir,id,args.ligand_dcd_FN)
        if not os.path.isfile(traj_FN):
            print 'Ligand trajectory missing in '+traj_FN
            (Es[id]['com'], Es[id]['euler_angles']) = (np.zeros(3), np.zeros(3))
            continue
        prmtop_FN = os.path.join(args.ligand_dir,id,args.ligand_prmtop_FN)
        traj = mdtraj.load(traj_FN, top=prmtop_FN)
        (Es[id]['com'], Es[id]['euler_angles']) = get_external_dof(traj)
        nc = Dataset(FN,'w',format='NETCDF4')
        nc.createDimension('n_frames', Es[id]['com'].shape[0])
        nc.createDimension('n_coords', 3)
        nc.createVariable('com','f8',('n_frames','n_coords'))
        nc.createVariable('euler_angles','f8',('n_frames','n_coords'))
        nc.variables['com'][:] = Es[id]['com']
        nc.variables['euler_angles'][:] = Es[id]['euler_angles']
        nc.close()
sys.stdout.write('\n\n')

print 'Number of energies:'
#print '\tSampler\t|\tGBSA\t\t\t|\tnm\t\t\t|\tCOM\t|\tRMSD'
print '\t|\tGBSA\t\t'
for id in ids:
#    count_str = id + '\t%d'%(len(Es[id]['sampler'])) + '\t|\t'
    count_str = '\t|\t'.join(['\t'.join([\
        'nan' if np.isnan(Es[id]['%s_%s'%(solvent,moeity)]).any() \
        else '%d'%len(Es[id]['%s_%s'%(solvent,moeity)]) \
#            for moeity in ['L','R','RL']]) for solvent in ['GBSA','nm']])
            for moeity in ['L','R','RL']]) for solvent in ['GBSA']])
#    count_str += '\t|\t' + '%d'%len(Es[id]['com'])
#    count_str += '\t|\t' + '%d'%len(Es[id]['RMSD'])
    print count_str
print

#################################
# Binding free energy estimates #
#################################

def estimate_dG(Psi, nm_TdS, DeltaU=None):
    """
    Returns an array of free energy estimates:
    based on cumulants 1-4, EXP, and then NM.
    Psi are interaction energies.
    nm_TdS are entropy terms from normal modes analysis.
    DeltaU are reweighting factors for individual configurations.
    """
    if DeltaU is None:
        cumulants = [np.nan] + [scipy.stats.kstat(Psi,n) for n in range(1,5)]
        dG = np.cumsum([BETA**(n-1)*cumulants[n]/math.factorial(n) \
            for n in range(1,5)])
        exp = (np.log(np.mean(np.exp(BETA*Psi-max(BETA*Psi))))+max(BETA*Psi))/BETA
    else:
        num = Psi - DeltaU
        den = -DeltaU
        cumulants_num = [np.nan] + [scipy.stats.kstat(num,n) for n in range(1,5)]
        cumulants_den = [np.nan] + [scipy.stats.kstat(den,n) for n in range(1,5)]
        dG = np.cumsum([BETA**(n-1)*(cumulants_num[n] - cumulants_den[n])/\
            math.factorial(n) for n in range(1,5)])
        exp = (np.log(np.average(np.exp(BETA*num-max(BETA*num)))) + max(BETA*num) \
            - np.log(np.average(np.exp(BETA*den-max(BETA*den)))) - max(BETA*den))/BETA
    return np.array(list(dG) + [exp, dG[0] - nm_TdS])

def calc_confinement_fe(com, euler, nbins=100):
    """
    Obtains the free energy of
    (1) putting the unbound ligand into the binding site, dG_xi
    (2) adding constraints to the unbound ligand and
        removing them from the bound complex, dG_c
    """
    # Scott's rule for Gaussian kernel density estimate bandwidth
    npts = com.shape[0]
    scotts_factor = np.power(npts,(-1./5))
    dG_xi = []
    dG_c = []
    U_c = []

    # Translational degrees of freedom

    # The standard state volume for a single molecule
    # in a box of size 1 L is 1.66053928 nanometers**3
    L_o = 1.66053928**(1/3.)

    def calc_translational(x):
        tail = (max(x)-min(x))*0.1
        edges = np.linspace(min(x)-tail, max(x)+tail, nbins)
        delta = (edges[1]-edges[0])
        centers = edges + delta/2.

        # Obtain confinement potential
        kernel = gaussian_kde(x)
        rho = kernel(centers)
        U = -np.log(rho)/BETA

        # Obtain confinement energies for bound state
        inds = np.array(np.floor(\
            (x-edges[0])/delta),dtype=np.int)
        inds[inds==nbins] = nbins-1
        DeltaU = U[inds]

        # Estimate confinement free energies
        # For dG_cL, use numerical integral for expectation
        dG_cL = -1/BETA*np.log(np.mean(np.exp(-BETA*U)))
        # For dG_cRL, use sample mean estimator
        dG_cRL = -1/BETA*np.log(np.mean(\
            np.exp(-BETA*DeltaU+min(BETA*DeltaU)))) + min(DeltaU)

        # Standard state correction
        L = edges[-1]-edges[0]+delta
        dG_xi = -1/BETA * np.log(L / L_o)

        return dG_xi, dG_cL - dG_cRL, DeltaU

    # Performs principal components analysis
    com_c = com - np.mean(com, 0)
    [w,v] = np.linalg.eig(np.dot(com_c.T,com_c)/com_c.shape[0])
    com_pca = np.dot(com_c,v)

    for dim in range(3):
        dG_xi_d, dG_c_d, U_c_d = calc_translational(com_pca[:,dim])
        dG_xi.append(dG_xi_d)
        dG_c.extend([dG_xi_d + dG_c_d])
        U_c.append(U_c_d)

    # Rotational degrees of freedom

    edges_2pi = np.linspace(-np.pi,np.pi,nbins)
    edges_pi = np.linspace(0,np.pi,nbins)

    def calc_rotational(x, edges, sine=False):
        delta = edges[1]-edges[0]
        centers = edges[:-1] + delta/2.

        # Obtain confinement potential
        # See https://stackoverflow.com/questions/18921419/\
        # implementing-a-2d-fft-based-kernel-density-estimator-in-python-and-comparing-i
        H = np.histogram(x,edges,density=True)[0]
        sigma = scotts_factor*10.*delta
        ker = np.exp(-(centers-centers[int(len(centers)/2)])**2/\
                     2/sigma**2)/sigma/np.sqrt(2*np.pi)
        rho = np.abs(np.fft.fftshift(np.fft.ifft(\
            np.fft.fft(H)*np.fft.fft(ker))).real)
        rho /= np.sum(rho)*delta
        U = -np.log(np.array(rho))/BETA

        # Obtain confinement energies for bound state
        inds = np.array(np.floor(\
            (x-edges[0])/delta),dtype=np.int)
        inds[inds==nbins] = nbins-1
        DeltaU = U[inds]

        # Estimate confinement free energies
        # For dG_cL, use numerical integral for expectation
        if sine:
            dG_cL = -1/BETA*np.log(np.sum(np.sin(centers)*np.exp(-BETA*U))/\
                                   np.sum(np.sin(centers)))
        else:
            dG_cL = -1/BETA*np.log(np.mean(np.exp(-BETA*U)))
        # For dG_cRL, use sample mean estimator
        dG_cRL = -1/BETA*np.log(np.mean(\
            np.exp(-BETA*DeltaU+min(BETA*DeltaU)))) + min(DeltaU)

        return dG_cL - dG_cRL, DeltaU

    dG_c_d, U_c_d = calc_rotational(euler[:,0], edges_2pi)
    dG_c.append(dG_c_d)
    U_c.append(U_c_d)
    dG_c_d, U_c_d = calc_rotational(euler[:,1], edges_pi, sine=True)
    dG_c.append(dG_c_d)
    U_c.append(U_c_d)
    dG_c_d, U_c_d = calc_rotational(euler[:,2], edges_2pi)
    dG_c.append(dG_c_d)
    U_c.append(U_c_d)

    return np.sum(dG_xi), np.sum(dG_c), np.sum(U_c,0)

dG = dict([(id,{}) for id in ids])

# Actually does free energy estimates
for id in ids:
    for solvent in ['GBSA']:
#        for rmsd_filter in rmsd_filters:
#            stride = len(Es[id]['RMSD'])/len(Es[id][solvent+'_RL'])
#            ok_rmsd = Es[id]['RMSD'][::stride]<rmsd_filter
#            Psi = (Es[id][solvent+'_RL'] - Es[id][solvent+'_R'] \
#                    - Es[id][solvent+'_L'])[ok_rmsd]
            Psi = (Es[id][solvent+'_RL'] - Es[id][solvent+'_R'] \
                    - Es[id][solvent+'_L'])

#            # Entropy from normal modes analysis
#            if not np.isnan(Es[id]['nm_RL']).any():
#                nm = Es[id]['nm_RL'] - Es[id]['nm_R'] - Es[id]['nm_L']
#                stride = len(ok_rmsd)/len(Es[id]['nm_RL'])
#                nm_TdS = np.mean(nm[ok_rmsd[::stride]])
#            else:
#                nm_TdS = np.nan
            nm_TdS = 0.0

            # Free energy estimate without binding site correction
            # or ligand confinement free energy
            # and without reweighting
#            dG[id]['%s_no_rw_%.2f_o'%(solvent,rmsd_filter)] = \
#                estimate_dG(Psi, nm_TdS)
            dG[id]['%s_no_rw_o'%solvent] = \
                estimate_dG(Psi, nm_TdS)
            # or with reweighting
#            stride = len(Es[id]['sampler'])/len(Es[id][solvent+'_RL'])
#            DeltaU = (Es[id][solvent+'_RL'] - Es[id]['sampler'][::stride])[ok_rmsd]
#            dG[id]['%s_rw_%.2f_o'%(solvent,rmsd_filter)] = \
#                estimate_dG(Psi, nm_TdS, DeltaU)

            # External dof corrections
            (dG_xi, dG_c, U_c) = calc_confinement_fe(\
                Es[id]['com'], \
                Es[id]['euler_angles'])
#                Es[id]['com'][ok_rmsd[:Es[id]['com'].shape[0]]], \
#                Es[id]['euler_angles'][ok_rmsd[:Es[id]['euler_angles'].shape[0]]])

            # With binding site correction
            dG[id]['%s_no_rw_xi'%solvent] = \
                estimate_dG(Psi, nm_TdS) + dG_xi
            # With ligand confinement free energy
            dG[id]['%s_no_rw_c'%solvent] = \
                estimate_dG(Psi[(Psi.shape[0] - U_c.shape[0]):], nm_TdS, \
                U_c) + dG_c
#            dG[id]['%s_no_rw_%.2f_xi'%(solvent, rmsd_filter)] = \
#                estimate_dG(Psi, nm_TdS) + dG_xi
#            dG[id]['%s_rw_%.2f_xi'%(solvent, rmsd_filter)] = \
#                estimate_dG(Psi, nm_TdS, DeltaU) + dG_xi
#            # With ligand confinement free energy
#            dG[id]['%s_no_rw_%.2f_c'%(solvent, rmsd_filter)] = \
#                estimate_dG(Psi[(Psi.shape[0] - U_c.shape[0]):], nm_TdS, \
#                U_c) + dG_c
#            dG[id]['%s_rw_%.2f_c'%(solvent, rmsd_filter)] = \
#                estimate_dG(Psi[(Psi.shape[0] - U_c.shape[0]):], nm_TdS, \
#                DeltaU[(Psi.shape[0] - U_c.shape[0]):] + U_c) + dG_c

for confinement in ['o','xi','c']:
#    key = 'GBSA_no_rw_%.2f_%s'%(rmsd_filters[0], confinement)
    key = 'GBSA_no_rw_%s'%confinement
    print 'Free energies for '+key
    for id in ids:
        print id + '\t' + ' '.join(['%6.2f'%fe \
            for fe in dG[id][key]])
    print

########################
# Performance analysis #
########################

from scipy.stats import spearmanr, kendalltau

def PearsonR(x, y):
    return np.corrcoef(x,y)[0,1]

def SpearmanR(x, y):
    return spearmanr(x,y)[0]

def KendallTau(x, y):
    return kendalltau(x,y)[0]

def RMSE(x,y):
    return np.sqrt(np.mean((x-y)**2))

def aRMSE(x,y):
    mean_diff = np.mean(x)-np.mean(y)
    return np.sqrt(np.mean((x-y-mean_diff)**2))

def bootstrap_std(x,y,func):
    vals = []
    rep = 0
    while rep<1000:
        inds = np.random.randint(len(x), size=len(x))
        val = func(x[inds],y[inds])
        if not np.isnan(val):
            vals.append(func(x[inds],y[inds]))
            rep += 1
    return np.std(vals)
metrics = ['PearsonR', 'SpearmanR', 'KendallTau', 'RMSE', 'aRMSE']

if not args.skip_analysis:
    dG_r = np.array([ref_dG[id] for id in ids])
    print 'Performance analysis:'
    print
    print '&&%17s'%'' + ''.join(['&%13s'%t \
        for t in ['1','2','3','4','EXP','1+nm']]) + '\\\\'

    solvent = 'GBSA' # also could be 'PBSA'
    rw = 'no_rw' # also could be 'rw'
    for rmsd_filter in rmsd_filters:
        for confinement in ['_o','_xi','_c']:
            prefix = repr(rmsd_filter) + '&' + \
                {'_o':' No', '_xi':'Site', '_c':'Yes'}[confinement] + '&'
            tables = {}
            for metric in metrics:
                tables[metric] = []
            rows = {}
            for metric in metrics:
                rows[metric] = [prefix + '%11s'%metric]
            key = '%s_%s_%.2f%s'%(solvent,rw,rmsd_filter,confinement)
            for order_ind in range(6):
                dG_t = np.array([dG[id][key][order_ind] for id in ids])
                rows['PearsonR'].append('&%6.2f (%.2f)'%(\
                    PearsonR(dG_r,dG_t),bootstrap_std(dG_r,dG_t,locals()['PearsonR'])))
                rows['SpearmanR'].append('&%6.2f (%.2f)'%(\
                    SpearmanR(dG_r,dG_t),bootstrap_std(dG_r,dG_t,locals()['SpearmanR'])))
                rows['KendallTau'].append('&%6.2f (%.2f)'%(\
                    KendallTau(dG_r,dG_t),bootstrap_std(dG_r,dG_t,locals()['KendallTau'])))
                rows['RMSE'].append('&%6.2f (%.2f)'%(\
                    RMSE(dG_r,dG_t),bootstrap_std(dG_r,dG_t,locals()['RMSE'])))
                rows['aRMSE'].append('&%6.2f (%.2f)'%(\
                    aRMSE(dG_r,dG_t),bootstrap_std(dG_r,dG_t,locals()['aRMSE'])))
            for metric in metrics:
                tables[metric].append(rows[metric][0] + \
                    ''.join(['%14s'%r for r in rows[metric][1:]]) + '\\\\')

            # Print out a table
            for metric in metrics:
                print '\n'.join(tables[metric])
            print
