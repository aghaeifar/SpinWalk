"""
create fieldmap for randomly oriented vessels with a given oxygenation level

@author: Ali Aghaeifar <ali.aghaeifar@tuebingen.mpg.de>

Reference for the fieldmap equation:
Báez-Yánez et al. "The impact of vessel size, orientation and intravascular contribution on the neurovascular fingerprint of BOLD bSSFP fMRI". In NeuroImage (2017)(Vol. 163, pp. 13-23).
"""


import numpy as np
import scipy as sp
from math import sin, cos

# define parameters
FoV  = np.asarray([1e3, 1e3, 1e3]) # Field of view in um
Res  = np.asarray([256, 256, 256]) # Resolution in voxels
FoV2 = (FoV - FoV/Res) / 2 # update Field of view to account for the voxel size and center the field map

is_orientation_random = False 
theta_B0 = np.pi/2     # angle between axis of vessel and B0 if is_orientation_random == True
Y     = 0.75  # oxygenation level of blood
Htc   = 0.4   # hematocrit level
dChi  = 0.273e-6 # susceptibility difference between fully deoxygeneted blood and tissue
BV    = 0.03  # blood volume fraction
vessel_radius = 50 # vessel diameter in um


# create a grid of points
xv, yv, zv = np.meshgrid(np.linspace(-FoV2[0], FoV2[0], Res[0], endpoint=True), 
                         np.linspace(-FoV2[1], FoV2[1], Res[1], endpoint=True),
                         np.linspace(-FoV2[2], FoV2[2], Res[2], endpoint=True), indexing='ij')
p3 = np.vstack((xv.flatten(order='F'), yv.flatten(order='F'), zv.flatten(order='F')), dtype=np.float32).T
xv = None  # save memory
yv = None
zv = None


# create a fieldmap
dW = np.zeros(p3.shape[0], dtype=np.float32)
vessel_mask = np.zeros(Res, dtype=bool).flatten(order='F')

while True:
    # vessel, is a strait line defined by two points, first point is randomly generated in cartesian coordinate, second point is generated in spherical coordinate to address orientation respect to B0
    p1 = np.random.uniform(-FoV2, FoV2, 3) # first point
    if is_orientation_random == True:
        phi   = np.random.uniform(0, 2*np.pi, 1)
        theta_B0 = np.random.uniform(0, np.pi/2, 1)
    else:
        phi = 0
    
    p21 = np.asarray([sin(theta_B0) * cos(phi), sin(theta_B0) * sin(phi), cos(theta_B0)]) # p2 - p1
    p2  = p21 + p1 # second point    

    # project p3 (all spatial points) to the line defined by p1 and p2
    t = np.sum((p3 - p1) * (p2 - p1), axis=1) / np.sum((p1-p2)**2)
    t = t[..., np.newaxis]
    projection = p1 + t * (p2 - p1)
    v1_u = p3 - projection  # vector from projection point to p3
    d = np.linalg.norm(v1_u, axis=1)  # distance from projection point (cylinder axis) to p3 (all spatial points)
    # d = np.linalg.norm(np.cross(p2-p1,p3-p1), axis=1)/np.linalg.norm(p2-p1) # distance to vessel

    # reference vector of current cylinder (vessel), a vector prependicular to the cylinder axis and its Y element is zero
    v2_u = np.asarray([1., 0., 0.], dtype = np.float32) # suppose this is a vector of form [1,0,x], connects point P3 to p1
    v2_u[2] = -v2_u[0] * p21[0]/p21[2] # update the x (3rd) element of the v2_u. We did a dot product between v and p21 to get the x element of v2_u. Dot product should be zero since two vectors are perpendicular. We can solve for x element of v2_u.
    cos2theta = np.divide(np.dot(v1_u, v2_u), np.linalg.norm(v1_u, axis=1)) / np.linalg.norm(v2_u)
    cos2theta = cos2theta * cos2theta
    
    dW_temp = 2*np.pi*(1-Y)*dChi*Htc * (vessel_radius**2) / (d**2) * (sin(theta_B0)**2) * (2*cos2theta - 1)
    f_in = 2*np.pi*(1-Y)*dChi*Htc * (cos(theta_B0)**2 - 1/3)
        
    mask_local = np.zeros_like(vessel_mask, dtype = bool) # maks for the current vessel
    np.putmask(mask_local, (d < vessel_radius), True)  
    vessel_mask = np.logical_or(vessel_mask, mask_local) # expand vessel mask
    
    np.putmask(dW_temp, vessel_mask, 0) # exclude fields inside vessels
    dW += dW_temp
    np.putmask(dW, mask_local, f_in)  # set the values of fieldmap inside current vessel
        
    BV_t = np.sum(vessel_mask) / np.prod(Res)
    if BV_t > BV:
        break
   
dW = dW.reshape(Res)
vessel_mask = vessel_mask.reshape(Res)
filename = 'fieldmap.mat'
if is_orientation_random == False:
    filename = 'fieldmap_ori=' + "{:.1f}".format(theta_B0*180/np.pi) + '.mat'
sp.io.savemat(filename, mdict={'fieldmap': dW, 'mask': vessel_mask, 'FoV': FoV*1e-6, 'Res': Res, 'vessel_radius': vessel_radius*1-6, 'BV': BV, 'Y': Y})
