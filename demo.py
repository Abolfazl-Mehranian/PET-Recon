
import numpy as np
import matplotlib.pyplot as plt
from geometry.BuildGeometry import BuildGeometry

#%% Siemense mMR PET
mmr = BuildGeometry('mmr')

save_dir='J:\MyPyWorkSapce\mmr2008\\'
isPrecomputed = False
if not isPrecomputed:
    mmr.buildSystemMatrixUsingSymmetries(save_dir)

img = mmr.buildPhantom(model=0)    
    
#2D example
mmr.loadSystemMatrix(save_dir,is3d=False)

img2d = img[:,:,50]
y2d = mmr.forwardProject2D(img2d,psf=0.2)
imgx2d = mmr.MLEM2D(y2d,psf=0.2,niter = 50)
plt.figure()
plt.subplot(1,2,1),plt.imshow(img2d),plt.title('Original',fontsize=15)
plt.subplot(1,2,2),plt.imshow(imgx2d),plt.title('Reconstructed',fontsize=15)

# 3D example 
mmr.loadSystemMatrix(save_dir)
y3d = mmr.forwardProject3D(img2d,psf=0.2)
imgx3d = mmr.MLEM3D(y3d,psf=0.2,niter = 5)
plt.figure()
plt.subplot(1,2,1),plt.imshow(img[:,:,50]),plt.title('Original',fontsize=15)
plt.subplot(1,2,2),plt.imshow(imgx3d[:,:,50]),plt.title('Reconstructed',fontsize=15)

#%% Siemense mCT TOF PET

mct = BuildGeometry('mct')

save_dir='J:\MyPyWorkSapce\mct1104\\'
isPrecomputed = False
if not isPrecomputed:
    mct.buildSystemMatrixUsingSymmetries(save_dir)

img = mct.buildPhantom(model=0)    
    
#2D example
mct.loadSystemMatrix(save_dir,is3d=False)

img2d = img[:,:,50]


y2d_tof = mct.forwardProject2D(img2d,psf=0.2,tof=True)
y2d = mct.forwardProject2D(img2d,psf=0.2) # y2d = np.sum(y2d_tof,axis=2)

imgx2d = mct.MLEM2D(y2d,psf=0.2,niter=20)
imgx2d_tof = mct.MLEM2D(y2d_tof,psf=0.2,niter=20,tof=True)

plt.figure()
plt.subplot(1,3,1),plt.imshow(img2d),plt.title('Original',fontsize=15)
plt.subplot(1,3,2),plt.imshow(imgx2d),plt.title('Reconstructed',fontsize=15)
plt.subplot(1,3,3),plt.imshow(imgx2d_tof),plt.title('TOF reconstructed',fontsize=15)

# 3D example 
mct.loadSystemMatrix(save_dir)

y3d_tof = mct.forwardProject2D(img2d,psf=0.2,tof=True)
y3d = np.sum(y3d_tof,axis=3)

imgx3d = mct.MLEM3D(y3d,psf=0.2,niter=20)
imgx3d_tof = mct.MLEM3D(y3d_tof,psf=0.2,niter=20,tof=True)

plt.figure()
plt.subplot(1,3,1),plt.imshow(img),plt.title('Original',fontsize=15)
plt.subplot(1,3,2),plt.imshow(imgx3d),plt.title('Reconstructed',fontsize=15)
plt.subplot(1,3,3),plt.imshow(imgx3d_tof),plt.title('TOF reconstructed',fontsize=15)
