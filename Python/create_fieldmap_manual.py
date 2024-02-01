#! /usr/bin/python
"""
create fieldmap for randomly oriented vessels with a given oxygenation level

@author: Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>

Reference for the fieldmap equation:
Báez-Yánez et al. "The impact of vessel size, orientation and intravascular contribution on the neurovascular fingerprint of BOLD bSSFP fMRI". In NeuroImage (2017)(Vol. 163, pp. 13-23).
"""

import os
import sys
import numpy as np
import scipy.io
from math import sin, cos
from tqdm import tqdm


########  define parameters ########
FoV = np.asarray([3000, 3000, 3000])  # Field of view in um
Res = np.asarray([1500, 1500, 1500])  # Resolution in voxels
# update Field of view to account for the voxel size and center the field map
FoV2 = (FoV - FoV/Res) / 2

Y    = [0.77, 0.85] # oxygenation level of blood, 0.77, 0.85
Htc  = 0.4              # hematocrit level
dChi = 0.273e-6         # susceptibility difference between fully deoxygeneted blood and tissue
vessel_radius = 150      # vessel diameter in um

Theta = [np.pi/2] # angle between axis of vessel and B0 (in radian), B0 is aligned with z axis
Alpha = [0]  # inplace angle of the vessel axis, with respect to x axis (in radian)

# vessel, is a strait line defined by two points, first point is randomly generated in cartesian coordinate, second point is generated in spherical coordinate to address orientation respect to B0
p1z = int(sys.argv[1]) # vessel center in z direction
P1  = [np.asarray([0, 0, p1z], dtype=np.float32)]  # first point

########  prepare ########
print('create a grid of points...', end='')
xv, yv, zv = np.meshgrid(np.linspace(-FoV2[0], FoV2[0], Res[0], endpoint=True, dtype=np.float32),
                         np.linspace(-FoV2[1], FoV2[1], Res[1], endpoint=True, dtype=np.float32),
                         np.linspace(-FoV2[2], FoV2[2], Res[2], endpoint=True, dtype=np.float32), indexing='ij')
print('step1...', end='')
p3 = np.vstack((xv.flatten(order='F'), yv.flatten(order='F'), zv.flatten(order='F'))).T
print('step2...done.')
xv = None  # save memory
yv = None
zv = None

for Y in Y:
    field_map   = np.zeros(p3.shape[0], dtype=np.float32)
    vessel_mask = np.zeros_like(field_map, dtype=bool)

    for (theta, alpha, p1) in zip(Theta, Alpha, P1):
        p21 = np.asarray([sin(theta) * cos(alpha), sin(theta) * sin(alpha), cos(theta)], dtype=np.float32)  # p2 - p1
        dW = np.array([], dtype=np.float32)
        D  = np.array([], dtype=np.float32)
        # break p3 to smaller chunks to save memory
        p3chunks = np.array_split(p3, p3.shape[0]//1e7, axis=0)
        for p3c in tqdm(p3chunks):
            # vector from projection of p3 in p21 to p3
            t = np.matmul((p3c - p1), p21) / np.linalg.norm(p21) # A.B = |A||B|cos(theta)  --> t = |B|cos(theta) = A.B / |A|
            t = t[..., np.newaxis] # add a new axis to t to be able to multiply it with p21    
            v1_u = p3c - (p1 + t * p21)
            t = None  # save memory
            d = np.linalg.norm(v1_u, axis=1) # distance from p3 (all spatial points) to p21 (cylinder axis)
            
            # reference vector of current cylinder (vessel), a vector prependicular to the cylinder axis and its Y element is zero
            # suppose this is a vector of form [1,0,x], connects point P3 to p1
            v2_u = np.asarray([1., 0., 0.], dtype=np.float32)
            # update the x (3rd) element of the v2_u. We did a dot product between v and p21 to get the x element of v2_u. Dot product should be zero since two vectors are perpendicular. We can solve for x element of v2_u.
            v2_u[2] = -v2_u[0] * p21[0]/p21[2] # v2u.p21 = 0  --> v2u[0]*p21[0] + v2u[1]*p21[1] + v2u[2]*p21[2] = 0 --> setting v2u[1]=0, then v2u[2] = -v2u[0]*p21[0]/p21[2]
            cos2phi = np.divide(np.dot(v1_u, v2_u), d) / np.linalg.norm(v2_u) # cos(theta) = A.B / |A||B|
            cos2phi = cos2phi * cos2phi
            v1_u = None  # save memory

            dW_temp = 2*np.pi*(1-Y)*dChi*Htc * (vessel_radius**2) / (d**2) * (sin(theta)**2) * (2*cos2phi - 1)  
            cos2phi = None  # save memory

            dW = np.hstack([dW, dW_temp]) if dW.size else dW_temp   
            D  = np.hstack([D, d]) if D.size else d
            

        f_in = 2*np.pi*(1-Y)*dChi*Htc * (cos(theta)**2 - 1/3)
        # print('RAM Used7 (GB):', psutil.virtual_memory()[3]/1024/1024/1024)

        # maks for the current vessel
        mask_local = np.zeros_like(vessel_mask, dtype=bool)
        np.putmask(mask_local, (D < vessel_radius), True)

        vessel_mask = np.logical_or(vessel_mask, mask_local)  # expand vessel mask
        np.putmask(dW, vessel_mask, 0)  # exclude fields inside vessels

        field_map += dW
        # set the values of fieldmap inside current vessel
        np.putmask(field_map, mask_local, f_in)

    BV  = np.sum(vessel_mask) / np.prod(Res)
    print(f'Blood Volume = {BV*100:.1f}%')

    field_map = field_map.reshape(Res, order='F')
    vessel_mask = vessel_mask.reshape(Res, order='F')
    # filename = f'./Y{Y}.mat'
    # print(f'Saving to {os.path.abspath(filename)}')
    # scipy.io.savemat(filename, mdict={'fieldmap': field_map, 'mask': vessel_mask, 'FoV': FoV * 1e-6, 'Res': Res, 'vessel_radius': vessel_radius*1e-6, 'BV': BV, 'Y': Y})

    filename = f'./Y{Y}_d{p1z}.bin'
    print(f'Saving to {os.path.abspath(filename)}')
    with open(filename, 'wb') as f:
        f.write(Res.astype('uint32').tobytes(order='F'))
        f.write((FoV*1e-6).astype('float32').tobytes(order='F'))
        f.write(field_map.astype('float32').tobytes(order='F'))
        f.write(vessel_mask.astype('uint8').tobytes(order='F'))