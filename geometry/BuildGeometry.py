# -*- coding: utf-8 -*-
"""
Created on Sun Apr 29 23:45:29 2018

@author: abm15
"""

mct = {'model_number':          1104,
       'circularGantry':        1,
       'nBuckets':              48,
       'nBlockRings':           4,
       'nBlockPerRing':         48,
       'nPhysCrystalsPerBlock': 13,
       'useVirtualCrystal':     1,
       'detectorRadiusCm':      42.76,
       'sinogramDOIcm':         0.67,
       'LORDOIcm':              0.96,
       'nRadialBins':           400,
       'nMash':                 2,
       'rCrystalDimCm':         2.0,
       'xCrystalDimCm':         0.40728,
       'zCrystalDimCm':         0.4050,
       'transaxialFovCm':       69.7266,
       'span':                  11,
       'nSegments':             9,
       'maxRingDiff':           49,
       'nTofBins':              13,
       'coinciWindowWidthNsec': 4.0625,
       'tofResolutionNsec':     0.580,
       'tofOffsetNsec':         0.039            
       }

mmr = {'model_number':          2008,
       'circularGantry':        1,
       'nBuckets':              224,
       'nBlockRings':           8,
       'nBlockPerRing':         56,
       'nPhysCrystalsPerBlock': 8,
       'useVirtualCrystal':     1,
       'detectorRadiusCm':      32.8,
       'sinogramDOIcm':         0.67,
       'LORDOIcm':              0.96,
       'nRadialBins':           344,
       'nMash':                 1,
       'rCrystalDimCm':         2.0,
       'xCrystalDimCm':         0.41725,
       'zCrystalDimCm':         0.4063,
       'transaxialFovCm':       60.0,
       'span':                  11,
       'nSegments':             11,
       'maxRingDiff':           60,
       'nTofBins':              1,
       'coinciWindowWidthNsec': 5.85938,
       'tofResolutionNsec':     5.85938,
       'tofOffsetNsec':         0            
       }


import numpy as np


class Attr():
    def __setattr__(self, name, value):
         self.__dict__[name] = value
    def __getitem__(self, name):
        return self[name]

class BuildGeometry:
    def __init__(self,scannerModel):
        
        self.scanner = Attr()
        self.sinogram = Attr()
        self.image = Attr()

        if scannerModel.lower() == 'mct':
            self.__computeGantryInfo(mct)
        elif scannerModel.lower() == 'mmr':
            self.__computeGantryInfo(mmr)  


    def __computeGantryInfo(self,g):
        self.__load_gantry_dict(g)
        self.scanner.nCrystalsPerBlock = self.scanner.nPhysCrystalsPerBlock + self.scanner.useVirtualCrystal
        self.scanner.nCrystalsPerRing = self.scanner.nBlockPerRing * self.scanner.nCrystalsPerBlock
        if self.scanner.model_number == 1104:
            self.scanner.nCrystalRings = self.scanner.nBlockRings * self.scanner.nPhysCrystalsPerBlock + (self.scanner.nBlockRings-1)*self.scanner.useVirtualCrystal
        elif self.scanner.model_number ==2008:
            self.scanner.nCrystalRings = self.scanner.nBlockRings * self.scanner.nPhysCrystalsPerBlock
        self.scanner.effDetectorRadiusCm = self.scanner.detectorRadiusCm + self.scanner.LORDOIcm
        self.scanner.isTof = self.sinogram.nTofBins>1
        self.scanner.TofBinWidthNsec = self.scanner.coinciWindowWidthNsec/self.sinogram.nTofBins
        self.scanner.planeSepCm = self.scanner.zCrystalDimCm/2.0        
        self.sinogram.nAngularBins = self.scanner.nCrystalsPerRing//2//self.sinogram.nMash
        self.image.matrixSize = [self.sinogram.nRadialBins,self.sinogram.nRadialBins,2*self.scanner.nCrystalRings-1]
        self.image.voxelSizeCm = [self.scanner.xCrystalDimCm/2.0,self.scanner.xCrystalDimCm/2.0,self.scanner.planeSepCm]
        
    def __load_gantry_dict(self,g):
        self.scanner.model_number = g['model_number']
        self.scanner.circularGantry = g['circularGantry']
        self.scanner.nBuckets = g['nBuckets']
        self.scanner.nBlockRings = g['nBlockRings']
        self.scanner.nBlockPerRing = g['nBlockPerRing']
        self.scanner.nPhysCrystalsPerBlock = g['nPhysCrystalsPerBlock']
        self.scanner.useVirtualCrystal = g['useVirtualCrystal']
        self.scanner.detectorRadiusCm = g['detectorRadiusCm']
        self.scanner.sinogramDOIcm = g['sinogramDOIcm']
        self.scanner.LORDOIcm = g['LORDOIcm']
        self.scanner.rCrystalDimCm = g['rCrystalDimCm']
        self.scanner.xCrystalDimCm = g['xCrystalDimCm']
        self.scanner.zCrystalDimCm = g['zCrystalDimCm']
        self.scanner.transaxialFovCm = g['transaxialFovCm']
        self.scanner.maxRingDiff = g['maxRingDiff']
        self.scanner.coinciWindowWidthNsec = g['coinciWindowWidthNsec']
        self.scanner.tofResolutionNsec = g['tofResolutionNsec']
        self.scanner.tofOffsetNsec = g['tofOffsetNsec']        
        self.sinogram.nRadialBins = g['nRadialBins']
        self.sinogram.nMash = g['nMash']
        self.sinogram.span = g['span']
        self.sinogram.nSegments = g['nSegments']
        self.sinogram.nTofBins = g['nTofBins']

    def description(self):
        print('to do')
        
    def buildMichelogram(self):
        a = np.transpose(np.arange(1,self.scanner.nCrystalRings**2+1).reshape(self.scanner.nCrystalRings,self.scanner.nCrystalRings))
        b = np.arange(-1*self.scanner.maxRingDiff,self.scanner.maxRingDiff + 1).reshape(self.sinogram.nSegments,self.sinogram.span)
        direction = self.sinogram.nSegments//2
        isodd = np.remainder(b[direction,0],2)
        Segments = []
        maxNumberOfPlanesPerSeg = np.zeros([self.sinogram.nSegments, 2],dtype='int16')
        
        for j in range(self.sinogram.nSegments):
            diagonalsPerSegment = []
            for i in range(self.sinogram.span):
                diagonalsPerSegment.append(np.diag(a,k=b[j,i]))
            if j == direction and isodd:
                c=0; k=1
            else:
                c=1; k=0
            oddPlanes,maxNumberOfPlanesPerSeg[j,0] = self.__zero_pad(diagonalsPerSegment[k::2]) # Odd planes
            evenPlanes,maxNumberOfPlanesPerSeg[j,1] = self.__zero_pad(diagonalsPerSegment[c::2]) # Even planes
            OddEvenPlanesPerSegment = np.empty((np.sum(maxNumberOfPlanesPerSeg[j,:]), ), dtype=object)
            OddEvenPlanesPerSegment[0::2] = self.__zero_trim(oddPlanes)
            OddEvenPlanesPerSegment[1::2] = self.__zero_trim(evenPlanes)
            Segments.append(OddEvenPlanesPerSegment)
  
        self.sinogram.numberOfPlanesPerSeg = np.sum(maxNumberOfPlanesPerSeg,axis=1)
        self.sinogram.totalNumberOfSinogramPlanes = np.sum(self.sinogram.numberOfPlanesPerSeg)
        return Segments

    def plotMichelogram(self,showRingNumber=0):
        import matplotlib.pyplot as plt
        Segments = self.buildMichelogram()
        grid = np.zeros([self.scanner.nCrystalRings**2,1],dtype='int16')
        nS = self.sinogram.nSegments
        colourPerSeg = np.concatenate([np.arange(0,(nS-1)/2), [(nS-1)/2 +1], np.arange((nS-1)/2 -1,-1,-1)]) + 1
        
        for i in range(nS):
            idx = np.concatenate(Segments[i][:])-1
            grid[idx] = colourPerSeg[i]
        grid = grid.reshape([self.scanner.nCrystalRings,self.scanner.nCrystalRings])    
        ringNumber = np.arange(1,grid.size + 1)
        plt.imshow(grid, aspect='equal')
        if showRingNumber == 1:
            k = 0
            for (j, i), _ in np.ndenumerate(grid):
                label = '{}'.format(ringNumber[k])
                plt.text(i,j,label,ha='center',va='center',fontsize=12)
                k+=1 
        ax = plt.gca();
        ax = plt.gca();
        ax.set_xticks(np.arange(0, self.scanner.nCrystalRings, 1));
        ax.set_yticks(np.arange(0, self.scanner.nCrystalRings, 1));
        ax.set_xticklabels(np.arange(1, self.scanner.nCrystalRings+1, 1));
        ax.set_yticklabels(np.arange(1, self.scanner.nCrystalRings+1, 1));
        ax.set_xticks(np.arange(-.5, self.scanner.nCrystalRings, 1), minor=True);
        ax.set_yticks(np.arange(-.5, self.scanner.nCrystalRings, 1), minor=True);
        ax.grid(which='minor', color='k', linestyle='-', linewidth=2)
        plt.title('Michelogram: Span = {}, nSegments = {}'.format(self.sinogram.span,self.sinogram.nSegments),fontsize=18)
        plt.tight_layout()
        plt.show()

    def LorsAxialCoor(self):
        #Axial coordinate of LORs in each segment
        Segments = self.buildMichelogram()
        z_axis = self.scanner.zCrystalDimCm* np.arange(-(self.scanner.nCrystalRings -1)/2,(self.scanner.nCrystalRings -1)/2+1)
        axialCoorPerSeg = []
        
        for i in range(self.sinogram.nSegments):
            zy = np.zeros([self.sinogram.numberOfPlanesPerSeg[i],4])
            for j in range(self.sinogram.numberOfPlanesPerSeg[i]):
                ii,jj = self.__col2ij(Segments[i][j],self.scanner.nCrystalRings)
                zy[j,0] = np.mean(z_axis[ii])
                zy[j,1] = np.mean(z_axis[jj])
                zy[j,2] = self.scanner.effDetectorRadiusCm
                zy[j,3] = -self.scanner.effDetectorRadiusCm
            axialCoorPerSeg.append(zy)
        return axialCoorPerSeg,z_axis
    
    def plotLorsAxialCoor(self,plotSeparateSegmentsToo=0):
        import matplotlib.pyplot as plt
        axialCoorPerSeg,z_axis = self.LorsAxialCoor()
        
        plt.figure()
        for j in range(self.sinogram.nSegments):
            for i in range(len(axialCoorPerSeg[j])):
                plt.plot(axialCoorPerSeg[j][i,0:2],axialCoorPerSeg[j][i,2:4],color='green', linestyle='solid')  
        plt.plot(z_axis,self.scanner.effDetectorRadiusCm*np.ones([len(z_axis),]),'bs',fillstyle='none',markeredgewidth=2)
        plt.plot(z_axis,-self.scanner.effDetectorRadiusCm*np.ones([len(z_axis),]),'bs',fillstyle='none',markeredgewidth=2)
        plt.xlabel('Axial Distance (cm)',fontsize=18)
        plt.ylabel('Radial Distance (cm)',fontsize=18)
        plt.title('All Segments',fontsize=18)

        if plotSeparateSegmentsToo==1:
            ii = (self.sinogram.nSegments)//2
            order = np.zeros([self.sinogram.nSegments,],dtype='int16')
            order[0] = ii
            order[1::2] = np.arange(ii+1,self.sinogram.nSegments,dtype='int16')
            idx = order!=0
            order[2::2] = np.arange(ii-1,-1,-1,dtype='int16')
            q = 0       
            for j in range(self.sinogram.nSegments):
                if idx[j]:
                    if q==0:
                        seg_title = "Segment: {}"
                    else:
                        seg_title = "Segment: $\pm$ {}"
                    plt.figure()
                    plt.plot(z_axis,self.scanner.effDetectorRadiusCm*np.ones([len(z_axis),]),'bs',fillstyle='none',markeredgewidth=2)
                    plt.plot(z_axis,-self.scanner.effDetectorRadiusCm*np.ones([len(z_axis),]),'bs',fillstyle='none',markeredgewidth=2)
                    plt.xlabel('Axial Distance (cm)',fontsize=18)
                    plt.ylabel('Radial Distance (cm)',fontsize=18)
                    plt.title(seg_title.format(q),fontsize=18)
                    q+=1
                ii = order[j]
                for i in range(len(axialCoorPerSeg[ii])):
                    plt.plot(axialCoorPerSeg[ii][i,0:2],axialCoorPerSeg[ii][i,2:4],color='green', linestyle='solid')

    def LorsTransaxialCoor(self,startXtal=None):
        #Transaxial coordinate of LORs in each segment
        if startXtal is None:
            if self.scanner.model_number ==1104:
                startXtal = 68
            elif self.scanner.model_number ==2008:
                startXtal = 40     
        else:  
            self.sinogram.startXtal = startXtal
  
        p = np.linspace(2*np.pi,0,self.scanner.nCrystalsPerRing+1)
        centerCm = np.zeros([self.scanner.nCrystalsPerRing,2])
        centerCm[:,0]= self.scanner.effDetectorRadiusCm*np.cos(p[1::])
        centerCm[:,1]= self.scanner.effDetectorRadiusCm*np.sin(p[1::])
        isVirtualCrystal = np.zeros([self.scanner.nCrystalsPerRing,],dtype='bool');
        idx = np.arange(self.scanner.nPhysCrystalsPerBlock+1,self.scanner.nCrystalsPerRing+self.scanner.nPhysCrystalsPerBlock+1,
                        self.scanner.nPhysCrystalsPerBlock+1)
        isVirtualCrystal[idx-1]= 1
        
        increment = np.zeros([self.scanner.nCrystalsPerRing,2],dtype='int16')
        increment[0::2,0] = np.arange(1,self.scanner.nCrystalsPerRing/2+1)
        increment[1::2,0] = increment[0::2,0]+1
        increment[0::2,1] = np.arange(self.scanner.nCrystalsPerRing/2+1,self.scanner.nCrystalsPerRing+1)
        increment[1::2,1] = increment[0::2,1]

        halfNumberOfRadialBins = self.sinogram.nRadialBins//2+1 # Before interleaving
        R = np.empty((self.scanner.nCrystalsPerRing,3), dtype=object)
        V = np.zeros([halfNumberOfRadialBins,2],dtype='int16')
        
        for ii in range(self.scanner.nCrystalsPerRing):  
            s1 = (startXtal + np.arange(0, halfNumberOfRadialBins,dtype='int16')) - increment[ii,0]
            s2 = (startXtal + np.arange(0, halfNumberOfRadialBins,dtype='int16')) - increment[ii,1]
            s1 = self.__rem_p(s1, self.scanner.nCrystalsPerRing)-1
            s2 = self.__rem_p(s2, self.scanner.nCrystalsPerRing)-1
            s2 = s2[::-1]
            P1 = centerCm[s1,:]
            P2 = centerCm[s2,:]
            V = 0*V
            V[:,0] = isVirtualCrystal[s1]
            V[:,1] = isVirtualCrystal[s2]
            R[ii,0] = P1
            R[ii,1] = P2
            R[ii,2] = V    
            
        if 0: # test rotation
            import matplotlib.pyplot as plt
            plt.figure()
            for ii in range(0,self.scanner.nCrystalsPerRing,1):
                plt.clf()
                P1 = R[ii,0]
                P2 = R[ii,1]
                V = R[ii,2]
                for j in range(halfNumberOfRadialBins):
                    if V[j,0]:
                        s = 'og'
                    else:
                        s = '.b'
                    plt.plot(P1[j,0],P1[j,1],s)
                    if V[j,1]:
                        f = 'og'
                    else:
                        f = '.r'
                    plt.plot(P2[j,0],P2[j,1],f)
                    plt.axis('equal')
                    plt.xlim((-self.scanner.transaxialFovCm,self.scanner.transaxialFovCm))
                    plt.ylim((-self.scanner.transaxialFovCm,self.scanner.transaxialFovCm))
                    plt.plot([0,0],[-60,60],'--k')
                    plt.plot([-60,60],[0,0],'--k')
                    plt.plot(centerCm[:,0],centerCm[:,1],c='k',ls='--',lw=0.5)
                    plt.show()
                    plt.pause(0.1)
        # Interleaving
        xy1 = np.zeros([self.scanner.nCrystalsPerRing//2,self.sinogram.nRadialBins,2])
        xy2 = np.zeros([self.scanner.nCrystalsPerRing//2,self.sinogram.nRadialBins,2])
        gaps = np.zeros([self.scanner.nCrystalsPerRing//2,self.sinogram.nRadialBins],dtype='int16')#
        
        for i in range(self.scanner.nCrystalsPerRing//2):
            idx = self.scanner.nCrystalsPerRing-(2*i+2)
            P1=R[idx,0]
            P2=R[idx,1]
            xy1[i,0:self.sinogram.nRadialBins:2,:] = P1[0:-1,:]
            xy1[i,1:self.sinogram.nRadialBins:2,:] = (P1[0:-1,:]+P1[1::,:])/2
            xy2[i,0:self.sinogram.nRadialBins:2,:] = P2[0:-1,:]
            xy2[i,1:self.sinogram.nRadialBins:2,:] = (P2[0:-1,:]+P2[1::,:])/2
            a = np.sum(R[idx+1,2],axis=1).reshape(-1,1)>0
            b = np.sum(R[idx ,2],axis=1).reshape(-1,1)>0
            c = np.concatenate([a,b],axis=1).flatten()
            gaps[i,:]= c[1:-1]

        if self.sinogram.nMash==2:
            xy1 = (xy1[0::2,:,:] + xy1[1::2,:,:])/2
            xy2 = (xy2[0::2,:,:] + xy2[1::2,:,:])/2
            gap = np.zeros([self.sinogram.nAngularBins,self.sinogram.nRadialBins],dtype='int16')
            for i in range(self.sinogram.nAngularBins):
                gap[i,:] = np.sum(gaps[2*i:2*i+2,:],axis=0)
            gaps = gap
        gaps = np.transpose(gaps)
        
        # Calculate angular sampling
        centalBin = self.sinogram.nRadialBins//2-1
        p1 = xy1[0,centalBin,:]
        p2 = xy2[0,centalBin,:]
        lor1 = p2-p1
        p1 = xy1[1,centalBin,:]
        p2 = xy2[1,centalBin,:]
        lor2 = p2-p1
        CosTheta = np.dot(lor1,lor2)/(np.linalg.norm(lor1)*np.linalg.norm(lor2))
        self.sinogram.angSamplingDegrees = np.arccos(CosTheta)*180/np.pi                
        return xy1, xy2, gaps 

    def plotLorsTransaxialCoor(self):
        import matplotlib.pyplot as plt
        xy1, xy2, gaps = self.LorsTransaxialCoor()
        
        plt.figure()
        for i in range(self.sinogram.nAngularBins//4):
            plt.clf()
            plt.plot(np.array([xy1[i,:,0],xy2[i,:,0]]), np.array([xy1[i,:,1],xy2[i,:,1]]),c='green',ls='-',lw=0.5)
            if self.sinogram.nMash==2:
                idx = gaps[:,i]>1
            else:
                idx = gaps[:,i]>0
            
            plt.plot(np.array([xy1[i,idx,0],xy2[i,idx,0]]), np.array([xy1[i,idx,1],xy2[i,idx,1]]),c='blue',ls='-',lw=0.75)
            plt.axis('square')
            lim = self.scanner.transaxialFovCm*3/4
            plt.xlim((-lim,lim))
            plt.ylim((-lim,lim))
            plt.plot([0,0],[-lim,lim],'--k')
            plt.plot([-lim,lim],[0,0],'--k')
            plt.title('Angle: {}'.format(str(i+1)),fontsize=15)
            plt.show()
            plt.pause(0.1)  

    def Lors3DEndPointCoor(self,reduce4symmetries = 0):
        axialCoorPerSeg,_= self.LorsAxialCoor()
        xy1, xy2, gaps = self.LorsTransaxialCoor()
        xyz01 = np.zeros([self.sinogram.nAngularBins,self.sinogram.nRadialBins,3,self.sinogram.totalNumberOfSinogramPlanes],dtype='float32')
        xyz02 = np.zeros([self.sinogram.nAngularBins,self.sinogram.nRadialBins,3,self.sinogram.totalNumberOfSinogramPlanes],dtype='float32')
        k = 0
        flag = 1
        centralSegment = (self.sinogram.nSegments)//2
        for j in range(self.sinogram.nSegments):
            for i in range(self.sinogram.numberOfPlanesPerSeg[j]):
                z1 = axialCoorPerSeg[j][i,0]*np.ones([self.sinogram.nAngularBins,self.sinogram.nRadialBins])
                z2 = axialCoorPerSeg[j][i,1]*np.ones([self.sinogram.nAngularBins,self.sinogram.nRadialBins])
                if j>centralSegment: # for angles greater than 180
                    tmp = z1
                    z1 = z2
                    z2 = tmp
                    if flag:
                        xy1 = -xy1
                        xy2 = -xy2
                        flag = 0
                xyz01[:,:,0:2,k] = xy1
                xyz01[:,:,2,k] = z1
                xyz02[:,:,0:2,k] = xy2
                xyz02[:,:,2,k] = z2
                k+=1
        # sort the coordinate for segments 0 +1 -1 +2 -2,...
        cumulativePlaneNumber = np.cumsum(self.sinogram.numberOfPlanesPerSeg)
        planeRange = np.zeros([len(cumulativePlaneNumber),2],dtype='int16')
        planeRange[0,0] = 0
        planeRange[1:,0] = cumulativePlaneNumber[0:-1]
        planeRange[:,1] = cumulativePlaneNumber
        
        o = np.zeros([self.sinogram.nSegments,],dtype='int16')
        o[0::2] = np.arange(centralSegment,self.sinogram.nSegments)
        o[1::2] = np.arange(centralSegment-1,-1,-1)
        
        newCumulativePlaneNumber = np.cumsum(self.sinogram.numberOfPlanesPerSeg[o])
        newPlaneRange = np.zeros([len(newCumulativePlaneNumber),2],dtype='int16')
        newPlaneRange[0,0] = 0
        newPlaneRange[1:,0] = newCumulativePlaneNumber[0:-1]
        newPlaneRange[:,1] = newCumulativePlaneNumber
        
        self.sinogram.numberOfPlanesPerSeg = self.sinogram.numberOfPlanesPerSeg[o]
        self.sinogram.originalSegmentOrder = o
        
        if self.scanner.model_number ==2008:
            S = 0*newPlaneRange
            S[0,:] = newPlaneRange[0,:]
            for i in range(centralSegment):
                S[2*i+1,:] = newPlaneRange[2*i+2,:]
                S[2*i+2,:] = newPlaneRange[2*i+1,:]                      
            newPlaneRange = S
        self.sinogram.planeRange = newPlaneRange
        xyz1 = 0*xyz01
        xyz2 = 0*xyz02
        for i in range(self.sinogram.nSegments):
            xyz1[:,:,:,newPlaneRange[i,0]:newPlaneRange[i,1]]= xyz01[:,:,:,planeRange[o[i],0]:planeRange[o[i],1]]
            xyz2[:,:,:,newPlaneRange[i,0]:newPlaneRange[i,1]]= xyz02[:,:,:,planeRange[o[i],0]:planeRange[o[i],1]]
        if reduce4symmetries==1:
            self.calculateAxialSymmetries()
            xyz1 = xyz1[0:self.sinogram.nAngularBins//2,:,:,self.sinogram.uniqueAxialPlanes-1]
            xyz2 = xyz2[0:self.sinogram.nAngularBins//2,:,:,self.sinogram.uniqueAxialPlanes-1]

        return xyz1, xyz2, newPlaneRange

    def plotLors3DEndPointCoor(self,planeNumber=151):
        import matplotlib.pyplot as plt
        from mpl_toolkits.mplot3d import axes3d
        
        xyz1, xyz2, _ = self.Lors3DEndPointCoor()
        p = np.linspace(2*np.pi,0,self.scanner.nCrystalsPerRing+1)
        centerCm = np.zeros([self.scanner.nCrystalsPerRing,2])
        centerCm[:,0]= self.scanner.effDetectorRadiusCm*np.cos(p[1::])
        centerCm[:,1]= self.scanner.effDetectorRadiusCm*np.sin(p[1::])
                
        fig = plt.figure()
        ax = fig.add_subplot(111,projection='3d')
        i = 0
        x1 = xyz1[i,:,0,planeNumber]
        y1 = xyz1[i,:,1,planeNumber]
        z1 = xyz1[i,:,2,planeNumber]
        x2 = xyz2[i,:,0,planeNumber]
        y2 = xyz2[i,:,1,planeNumber]
        z2 = xyz2[i,:,2,planeNumber]

        for j in range(0,self.sinogram.nRadialBins,3):
            ax.plot(np.array([x1[j],x2[j]]),np.array([y1[j],y2[j]]),np.array([z1[j],z2[j]]),c='green',ls='-',lw=0.75)
        ax.plot(centerCm[0::3,0],centerCm[0::3,1],zs = z1[0],c='blue', markersize = 5, marker = 's', ls='None')
        ax.plot(centerCm[0::3,0],centerCm[0::3,1],zs = z2[0],c='blue', markersize = 5, marker = 's', ls='None')
        plt.show()        

    def calculateAxialSymmetries(self):                   
        newPlaneRange = self.sinogram.planeRange
        newPlaneRange[:,0] +=1                  
        l = self.sinogram.span//2 + 1   
        c = self.sinogram.nSegments//2           
        K = np.zeros([c,l],dtype='int16')
        K[0:c,0] = newPlaneRange[1::2,0] 
        K[0:c,1] = K[0:c,0] + 1   
        for i in range(2,l):
            K[0:,i] = K[0:,i-1]+2             
        self.sinogram.uniqueAxialPlanes = np.concatenate([[1],K.flatten()])

        ## calculate the translational and mirror symmetries
        b = newPlaneRange.flatten()
        b = np.reshape(b[2::],(4,c),order='F').transpose()
        n = self.sinogram.span-1 
        I = np.zeros([n,4],dtype='int16')
        x = np.arange(n)
        P = []
        for i in range(c):
            a = b[i,:]
            I = 0*I
            I[:,0] = a[0] + x
            I[:,1] = a[1] - x
            I[:,2] = a[2] + x
            I[:,3] = a[3] - x
            P.append(I)
        P = np.concatenate(P[:],axis=0)
        
        symID = np.zeros([P.shape[0],],dtype='int16')
        for i in range(1,len(self.sinogram.uniqueAxialPlanes)):
            j = P[:,0] == self.sinogram.uniqueAxialPlanes[i]
            symID[j] = i+1
        i = np.array(np.nonzero(symID==0),dtype='int16')
        symID[i] = symID[i-1]
        
        Ax = np.zeros([self.sinogram.totalNumberOfSinogramPlanes,],dtype='int16')
        mirror = np.ones([self.sinogram.totalNumberOfSinogramPlanes,],dtype='int16')
        Ax[0:self.sinogram.numberOfPlanesPerSeg[0]] = 1
        idx = n*np.arange(1,c+1)
        
        for i in range(len(symID)):
            if np.any((i+1)==idx):
                ii = np.concatenate([np.arange(P[i,0],P[i,1]+1),np.arange(P[i,2],P[i,3]+1)])
                Ax[ii-1] = symID[i]
                ii = np.arange(P[i,2],P[i,3])
                mirror[ii-1] = -1
            Ax[P[i,:]-1] = symID[i]
            mirror[P[i,2:4]-1] = -1

        offset = np.zeros([self.sinogram.totalNumberOfSinogramPlanes,],dtype='int16')
        for i in range(self.sinogram.totalNumberOfSinogramPlanes):
            if mirror[i]==1:
                offset[i] = (i+1) - self.sinogram.uniqueAxialPlanes[Ax[i]-1]
            else:
                j = np.array(np.nonzero(symID==(Ax[i]))[0])[0]
                x = P[j,0:3]
                offset[i] = (self.image.matrixSize[2]-1) -(x[1]-x[0]) + ((i+1) - x[2])
        planeMirrorTranslation = np.zeros([self.sinogram.totalNumberOfSinogramPlanes,3],dtype='int16')
        planeMirrorTranslation[:,0] = Ax
        planeMirrorTranslation[:,1] = mirror
        planeMirrorTranslation[:,2] = offset
        self.sinogram.planeMirrorTranslation = planeMirrorTranslation

              
    def calculateSystemMatrixPerPlane(self,xyz1,xyz2,I,reconFovRadious=None):
        if reconFovRadious is None:
            reconFovRadious = self.scanner.transaxialFovCm/2.5
        def paramInterSectPoint(p1x,p2x, amin,axmin,amax,axmax, tx,bx,dx,Nx):
            if  p1x < p2x:
                if amin == axmin:
                    imin = 1
                else:
                    imin = np.ceil(((p1x + amin*tx) - bx)/dx)
                if amax == axmax:
                    imax = Nx-1
                else:
                    imax = np.floor(((p1x + amax*tx) - bx)/dx)
                ax = (bx + np.arange(imin,imax+1)*dx - p1x ) / tx
            else:
                if amin == axmin:
                    imax = Nx-2
                else:
                    imax = np.floor(((p1x + amin*tx) - bx)/dx)
                if amax == axmax:
                    imin = 0
                else:
                    imin = np.ceil(((p1x + amax*tx) - bx)/dx)
                ax = (bx + np.arange(imax,imin-1,-1)*dx - p1x ) / tx
            return ax
        # Start Siddon tracing
        Nx = self.image.matrixSize[0] + 1
        Ny = self.image.matrixSize[1] + 1
        Nz = self.image.matrixSize[2] + 1
        dx = self.image.voxelSizeCm[0]
        dy = self.image.voxelSizeCm[1]
        dz = self.image.voxelSizeCm[2]
        bx = -(Nx-1)*dx/2
        by = -(Ny-1)*dy/2
        bz = -(Nz-1)*dz/2
        vCenter_y = dx*np.arange(-(Nx-2)/2,(Nx-2)/2+1)
        vCenter_x = dy*np.arange(-(Ny-2)/2,(Ny-2)/2+1)
        vCenter_z = dz*np.arange(-(Nz-2)/2,(Nz-2)/2+1)
        thresholdOfWeakLors = 50
        sMatrix = np.zeros((self.sinogram.nAngularBins//2,self.sinogram.nRadialBins), dtype=object)
        if self.scanner.isTof:
            from scipy.stats import norm
            tofBinBoundaries = np.linspace(-self.scanner.coinciWindowWidthNsec/2,self.scanner.coinciWindowWidthNsec/2,self.sinogram.nTofBins+1)
            sigma = self.scanner.tofResolutionNsec/np.sqrt(np.log(256))
            tofMatrix = np.zeros((self.sinogram.nAngularBins//2,self.sinogram.nRadialBins), dtype=object)
    
        for ang in range(self.sinogram.nAngularBins//2):
            for rad in range(self.sinogram.nRadialBins): 
    
                p1x = xyz1[ang,rad,0,I]
                p1y = xyz1[ang,rad,1,I]
                p1z = xyz1[ang,rad,2,I]
                p2x = xyz2[ang,rad,0,I]
                p2y = xyz2[ang,rad,1,I]
                p2z = xyz2[ang,rad,2,I]
                tx = p2x - p1x
                if tx == 0: 
                    p2x += 1e-2
                    tx = p2x - p1x
                ax = (bx + np.array([0 , Nx-1])*dx - p1x)/ tx
                axmin = ax.min()
                axmax = ax.max()
                ty = p2y - p1y
                if ty == 0: 
                    p2y += 1e-2
                    ty = p2y - p1y
                ay = (by + np.array([0 , Ny-1])*dy - p1y)/ ty
                aymin = ay.min()
                aymax = ay.max()
                tz = p2z - p1z
                if tz == 0: 
                    p1z += 1e-2
                    tz = p2z - p1z
                az = (bz + np.array([0 , Nz-1])*dz - p1z)/ tz
                azmin = az.min()
                azmax = az.max()
    
                amin = np.array([0,axmin, aymin, azmin]).max()
                amax = np.array([1,axmax, aymax, azmax]).min()
                
                if amin < amax:
                    ax = paramInterSectPoint(p1x,p2x, amin,axmin,amax,axmax, tx,bx,dx,Nx)
                    ay = paramInterSectPoint(p1y,p2y, amin,aymin,amax,aymax, ty,by,dy,Ny)
                    az = paramInterSectPoint(p1z,p2z, amin,azmin,amax,azmax, tz,bz,dz,Nz)
                    a = np.unique(np.concatenate(([[amin],ax,ay,az,[amax]])))
                    k = np.arange(len(a)-1)
                    im = np.floor(((p1x + ((a[k+1] + a[k])/2)*tx) - bx)/dx)
                    jm = np.floor(((p1y + ((a[k+1] + a[k])/2)*ty) - by)/dy)
                    km = np.floor(((p1z + ((a[k+1] + a[k])/2)*tz) - bz)/dz)
                    LorIntersectionLength = (a[k+1]-a[k])*np.sqrt(tx**2 + ty**2 + tz**2)*1e4/dx #normaized by pixel size
                    
                    M = np.stack([im, jm, km, LorIntersectionLength]).transpose().astype('int16')
                    #remove weak LOR interactions
                    weaksID = M[:,3]>thresholdOfWeakLors
                    M = M[weaksID,:]
                    if reconFovRadious!=0:# remove LOR interactions out of reduced reconstruction FOV 
                        IdsToKeep = (vCenter_y[M[:,0]]**2 + vCenter_x[M[:,1]]**2)<(reconFovRadious)**2
                        M = M[IdsToKeep,:]

                    if M.size!=0:
                        sMatrix[ang,rad] = M
                        if self.scanner.isTof:
                            VoxelCenters = np.stack([vCenter_y[M[:,0]],vCenter_x[M[:,1]],vCenter_z[M[:,2]]]).transpose()
                            nEmiPoint = VoxelCenters.shape[0]
                            endPoint1 = np.tile(np.array([p1x, p1y, p1z]),(nEmiPoint,1))
                            endPoint2 = np.tile(np.array([p2x, p2y, p2z]),(nEmiPoint,1))
                            dL = np.sqrt(np.sum((endPoint2 - VoxelCenters)**2 ,axis=1)) - np.sqrt(np.sum((endPoint1 - VoxelCenters)**2 ,axis=1))
                            dT = dL/30 - self.scanner.tofOffsetNsec
                            tofWeights = np.zeros([nEmiPoint,self.sinogram.nTofBins])
                            
                            for q in range(nEmiPoint):
                                tmp = norm.cdf(tofBinBoundaries,dT[q],sigma)
                                tofWeights[q,:] = tmp[1:] - tmp[:-1]
                            tofMatrix[ang,rad] = (tofWeights*1e4).astype('int16')
                        
        if not self.scanner.isTof:
            tofMatrix = 0
        return sMatrix,tofMatrix    
    
    def buildSystemMatrixUsingSymmetries(self,save_dir=None, reconFovRadious=None, ncores = 1):
        import os
        if save_dir is None:
            save_dir = os.getcwd()
        if not os.path.exists(save_dir):
            os.makedirs(save_dir)
        print("save path: {}".format(save_dir))
        if reconFovRadious is None:
            reconFovRadious = self.scanner.transaxialFovCm/2.5
        xyz1, xyz2, _ = self.Lors3DEndPointCoor(1)
        if ncores==1:
            for i in range(len(self.sinogram.uniqueAxialPlanes)):
                print(i)
                geoMatrix,tofMatrix = self.calculateSystemMatrixPerPlane(xyz1,xyz2,i,reconFovRadious)
                np.save(save_dir+'geoMatrix-'+str(i)+'.npy', geoMatrix)
                if self.scanner.isTof:
                    np.save(save_dir+'tofMatrix-'+str(i)+'.npy', tofMatrix)
        else:
            print('to do: multiprocessing')
    
    def loadSystemMatrix(self,save_dir,is3d = True, tof=True):
        import time
        """
        if case usage of: 
            PET = BuildGeometry('mct') 
            PET.loadSystemMatrix(save_dir)
            call "self.Lors3DEndPointCoor(1)" internally to set sinogram attributes 
        """
        tic = time.time()
        if not hasattr(self.sinogram,'uniqueAxialPlanes'):
            self.Lors3DEndPointCoor(1);
        self.geoMatrix = []
        if tof and self.scanner.isTof:
            self.tofMatrix = []
        N = len(self.sinogram.uniqueAxialPlanes) if is3d else 1
        for i in range(N):
            self.geoMatrix.append(np.load(save_dir+'geoMatrix-'+str(i)+'.npy'))
            if tof and self.scanner.isTof:
                self.tofMatrix.append(np.load(save_dir+'tofMatrix-'+str(i)+'.npy'))
        print('loaded in: {} sec.'.format(time.time()-tic))
    
    def forwardProject3D(self,img3d,tof=False, psf=0):
        import time
        if tof and not self.scanner.isTof:
           raise ValueError("The scanner is not TOF") 
        nUniqueAxialPlanes = len(self.sinogram.uniqueAxialPlanes)
        allPlanes = []
        for i in range(len(self.sinogram.uniqueAxialPlanes)):
            allPlanes.append(np.nonzero(self.sinogram.planeMirrorTranslation[:,0] == i+1)[0])
            
        img3d = self.gaussFilter(img3d.flatten('F'),psf)
        dims = [self.sinogram.nRadialBins,self.sinogram.nAngularBins,self.sinogram.totalNumberOfSinogramPlanes]
        if tof: dims.append(self.sinogram.nTofBins)
        y = np.zeros(dims,dtype='float')
        matrixSize = self.image.matrixSize
        q = self.sinogram.nAngularBins//2
        planeMirrorTranslation = self.sinogram.planeMirrorTranslation
        prod = lambda x, y: x.reshape(-1,1).dot(y.reshape(1,-1)).T.astype('int32')
        
        tic = time.time()
        for i in range(self.sinogram.nAngularBins//2):
            for j in range(self.sinogram.nRadialBins):
                for p in range(nUniqueAxialPlanes):
                    M0 = self.geoMatrix[p][i,j]
                    if not np.isscalar(M0):
                        M = M0[:,0:3].astype('int32')
                        G = M0[:,3]/1e4
                        H = planeMirrorTranslation[allPlanes[p],:]
                        idxAxial = matrixSize[0]*matrixSize[1]*(prod(H[:,1],M[:,2]) + H[:,2])
                        idx1 = (M[:,0] + M[:,1]*matrixSize[0]).reshape(-1,1) + idxAxial
                        idx2 = (M[:,1] + matrixSize[0]*(matrixSize[0]-1-M[:,0])).reshape(-1,1) + idxAxial
                        if tof:
                            W = self.tofMatrix[p][i,j]/1e4
                            y[j,i,allPlanes[p],:] = (G.reshape(-1,1)*img3d[idx1]).T.dot(W)
                            y[j,i+q,allPlanes[p],:] = (G.reshape(-1,1)*img3d[idx2]).T.dot(W)
                        else:
                            y[j,i,allPlanes[p]] = G.dot(img3d[idx1])
                            y[j,i+q,allPlanes[p]] = G.dot(img3d[idx2]) 
        print('forward-projected in: {} sec.'.format(time.time()-tic))                    
        return y 

    def MLEM3D(self, Prompts,img=None,RS=None,niter=100, AN=None, tof=False, psf=0):
        import time
        if tof and not self.scanner.isTof:
               raise ValueError("The scanner is not TOF") 
        if img is None:
            img = np.ones(self.image.matrixSize,dtype='float')
        img = img.flatten('F')  
        sensImage = np.zeros_like(img)
        if np.ndim(Prompts)!=4:
            tof = False
        if RS is None:
            RS = 0*Prompts
        if AN is None:
            AN = np.ones([self.sinogram.nRadialBins,self.sinogram.nAngularBins,self.sinogram.totalNumberOfSinogramPlanes],dtype='float')
    
        matrixSize = self.image.matrixSize
        q = self.sinogram.nAngularBins//2
        Flag = True    
    
        nUniqueAxialPlanes = len(self.sinogram.uniqueAxialPlanes)
        allPlanes = []
        for i in range(len(self.sinogram.uniqueAxialPlanes)):
            allPlanes.append(np.nonzero(self.sinogram.planeMirrorTranslation[:,0] == i+1)[0])
        planeMirrorTranslation = self.sinogram.planeMirrorTranslation
        prod = lambda x, y: x.reshape(-1,1).dot(y.reshape(1,-1)).T.astype('int32')
        
        tic = time.time()
        for n in range(niter):    
            if np.any(psf!=0):
                imgOld = self.gaussFilter(img,psf,True)
            else:
                imgOld = img
            backProjImage = 0*img 
            for i in range(self.sinogram.nAngularBins//2):
                for j in range(self.sinogram.nRadialBins):
                    for p in range(nUniqueAxialPlanes):
                        M0 = self.geoMatrix[p][i,j]
                        if not np.isscalar(M0):
                            M = M0[:,0:3].astype('int32')
                            G = M0[:,3]/1e4
                            H = planeMirrorTranslation[allPlanes[p],:]
                            idxAxial = matrixSize[0]*matrixSize[1]*(prod(H[:,1],M[:,2]) + H[:,2])
                            idx1 = (M[:,0] + M[:,1]*matrixSize[0]).reshape(-1,1) + idxAxial
                            idx2 = (M[:,1] + matrixSize[0]*(matrixSize[0]-1-M[:,0])).reshape(-1,1) + idxAxial
                            if tof:
                                W = self.tofMatrix[p][i,j]/1e4
                                f1 = AN[j,i,allPlanes[p]].reshape(-1,1)*(Prompts[j,i,allPlanes[p],:]/ ((G.reshape(-1,1)*img[idx1]).T.dot(W) + RS[j,i,allPlanes[p],:]+1e-5))
                                f2 = AN[j,i+q,allPlanes[p]].reshape(-1,1)*(Prompts[j,i+q,allPlanes[p],:]/ ((G.reshape(-1,1)*img[idx2]).T.dot(W) + RS[j,i+q,allPlanes[p],:]+1e-5))                           
                                backProjImage[idx1] += G.reshape(-1,1)*W.dot(f1.T)
                                backProjImage[idx2] += G.reshape(-1,1)*W.dot(f2.T)
                            else:
                                f1 = AN[j,i,allPlanes[p]]*(Prompts[j,i,allPlanes[p]]/(G.dot(imgOld[idx1])+RS[j,i,allPlanes[p]]+1e-5))
                                f2 = AN[j,i+q,allPlanes[p]]*(Prompts[j,i+q,allPlanes[p]]/(G.dot(imgOld[idx2])+RS[j,i+q,allPlanes[p]]+1e-5))
                                backProjImage[idx1] += G.reshape(-1,1).dot(f1.reshape(1,-1))
                                backProjImage[idx2] += G.reshape(-1,1).dot(f2.reshape(1,-1)) 
                            if Flag:
                                if tof:
                                    GW = G*np.sum(W,axis = 1)
                                    sensImage[idx1] += GW.reshape(-1,1).dot(AN[j,i,allPlanes[p]].reshape(1,-1)) 
                                    sensImage[idx2] += GW.reshape(-1,1).dot(AN[j,i+q,allPlanes[p]].reshape(1,-1)) 
                                else:
                                    sensImage[idx1] += G.reshape(-1,1).dot(AN[j,i,allPlanes[p]].reshape(1,-1)) 
                                    sensImage[idx2] += G.reshape(-1,1).dot(AN[j,i+q,allPlanes[p]].reshape(1,-1)) 
            if np.any(psf!=0) and Flag:
                sensImage = self.gaussFilter(sensImage,psf,True)
            Flag = False
            img = imgOld*backProjImage/(sensImage+1e-5)
                          
        print('forward-projected in: {} sec.'.format((time.time()-tic)/60))                    
        return img.reshape(matrixSize,order='F') 
    
    def forwardProject2D(self,img2d,tof=False, psf=0):
        if tof and not self.scanner.isTof:
           raise ValueError("The scanner is not TOF") 
        import time
        img2d = self.gaussFilter(img2d.flatten('F'),psf)
        dims = [self.sinogram.nRadialBins,self.sinogram.nAngularBins]
        if tof: dims.append(self.sinogram.nTofBins)
        y = np.zeros(dims,dtype='float')
        matrixSize = self.image.matrixSize
        q = self.sinogram.nAngularBins//2
        tic = time.time()
        for i in range(self.sinogram.nAngularBins//2):
            for j in range(self.sinogram.nRadialBins):
                M0 = self.geoMatrix[0][i,j]
                if not np.isscalar(M0):
                    M = M0[:,0:3].astype('int32')
                    G = M0[:,3]/1e4
                    idx1 = M[:,0] + M[:,1]*matrixSize[0]
                    idx2 = M[:,1] + matrixSize[0]*(matrixSize[0]-1-M[:,0])
                    if tof:
                        W = self.tofMatrix[0][i,j]/1e4
                        y[j,i,:] = (G*img2d[idx1]).dot(W)
                        y[j,i+q,:] =(G*img2d[idx2]).dot(W)
                    else:
                        y[j,i] = G.dot(img2d[idx1])
                        y[j,i+q] = G.dot(img2d[idx2]) 
        print('forward-projected in: {} sec.'.format(time.time()-tic))                    
        return y
#        
    def MLEM2D(self, Prompts,img=None,RS=None,niter=100, AN=None, tof=False, psf=0):
        import time
        if tof and not self.scanner.isTof:
               raise ValueError("The scanner is not TOF") 
        if img is None:
            img = np.ones(self.image.matrixSize[:2],dtype='float')
        img = img.flatten('F')  
        sensImage = np.zeros_like(img)
        if np.ndim(Prompts)!=3:
            tof = False
        if RS is None:
            RS = 0*Prompts
        if AN is None:
            AN = np.ones([self.sinogram.nRadialBins,self.sinogram.nAngularBins],dtype='float')
    
        matrixSize = self.image.matrixSize
        q = self.sinogram.nAngularBins//2
        Flag = True
        tic = time.time()
        for n in range(niter):
            if np.any(psf!=0):
                imgOld = self.gaussFilter(img,psf)
            else:
                imgOld = img
            backProjImage = 0*img        
            for i in range(self.sinogram.nAngularBins//2):
                for j in range(self.sinogram.nRadialBins):
                    M0 = self.geoMatrix[0][i,j]
                    if not np.isscalar(M0):
                        M = M0[:,0:3].astype('int32')
                        G = M0[:,3]/1e4
                        idx1 = M[:,0] + M[:,1]*matrixSize[0]
                        idx2 = M[:,1] + matrixSize[0]*(matrixSize[0]-1-M[:,0])
                        if tof:
                            W = self.tofMatrix[0][i,j]/1e4
                            backProjImage[idx1] += G*AN[j,i]*W.dot(Prompts[j,i,:]/((G*img[idx1]).dot(W)+RS[j,i,:]+1e-5))
                            backProjImage[idx2] += G*AN[j,i+q]*W.dot(Prompts[j,i+q,:]/((G*img[idx2]).dot(W)+RS[j,i+q,:]+1e-5))
                        else:
                            backProjImage[idx1] += G*AN[j,i]*(Prompts[j,i]/(G.dot(imgOld[idx1])+RS[j,i]+1e-5))
                            backProjImage[idx2] += G*AN[j,i+q]*(Prompts[j,i+q]/(G.dot(imgOld[idx2])+RS[j,i+q]+1e-5)) 
                        if Flag:
                            if tof:
                                GW = G*np.sum(W,axis = 1)
                                sensImage[idx1] += GW*AN[j,i]
                                sensImage[idx2] += GW*AN[j,i+q]
                            else:
                                sensImage[idx1] += G*AN[j,i]
                                sensImage[idx2] += G*AN[j,i+q]
            if np.any(psf!=0) and Flag:
                sensImage = self.gaussFilter(sensImage,psf)
            Flag = False
            img = imgOld*backProjImage/(sensImage+1e-5)
                        
        print('reconstructed in: {} min.'.format((time.time()-tic)/60))                    
        return img.reshape(matrixSize[:2],order='F')


        
    ''' PRIVATE HELPERS '''          
    def __zero_pad(self,y):
        maxNumberOfPlanes = [len(y[i]) for i in range(len(y))]
        PlaneNumbers = np.zeros([len(y),np.max(maxNumberOfPlanes)],dtype='int16')
        for i in range(len(y)):
            ii = (np.max(maxNumberOfPlanes) - len(y[i]))//2
            if ii==0:
                PlaneNumbers[i,:] = y[i]
            else:
                PlaneNumbers[i,ii:-ii] = y[i]
        return PlaneNumbers, np.max(maxNumberOfPlanes)
    
    def __zero_trim(self,y):
        out = []
        for i in range(len(y[0])):
            tmp = y[:,i]
            out.append(tmp[np.nonzero(tmp)])
        return out  

    def __col2ij(self,m,n):
        if np.max(m) > n**2:
            raise ValueError("m is greater than the max number of elements") 
        j = np.ceil(m/n)-1
        i = m-j*n-1
        return i.astype(int),j.astype(int) 
    
    def __rem_p(self,x,nx):
        for i in range(len(x)):
            while x[i] < nx:
                x[i] += nx
            while x[i] > nx:
                x[i] -= nx
        return x

    def gaussFilter(self,img,fwhm,is3d=False):
        # 3D aniso/isotropic Gaussian filtering
        from scipy import ndimage
        fwhm = np.array(fwhm)
        if np.all(fwhm==0):
            return img
        if is3d:
            img = img.reshape(self.image.matrixSize,order='F')
            voxelSizeCm = self.image.voxelSizeCm
        else:
            img = img.reshape(self.image.matrixSize[:2],order='F')
            voxelSizeCm = self.image.voxelSizeCm[:2]
        if fwhm.shape==1:
            if is3d:
                fwhm=fwhm*np.ones([3,])
            else:
                fwhm=fwhm*np.ones([2,])
        sigma=fwhm/voxelSizeCm/np.sqrt(2**3*np.log(2))
        imOut = ndimage.filters.gaussian_filter(img,sigma)
        return imOut.flatten('F')        
        
    def buildPhantom(self,model = 0,display = True):
        if model==0: # shepp-logan like phantom
            x = np.arange(-self.image.matrixSize[0]//2, self.image.matrixSize[0]//2, 1)
            y = np.arange(-self.image.matrixSize[1]//2, self.image.matrixSize[1]//2, 1)
            z = np.arange(-self.image.matrixSize[2]//2, self.image.matrixSize[2]//2, 1)
            xx, yy,zz = np.meshgrid(x, y,z)
            
            z1 = (xx**2+(yy/1.3)**2+((zz)/0.8)**2)<80**2
            z2 = 2.5*(((xx+20)**2+(yy+20)**2+(zz/3)**2)<10**2)
            z3 = 3*(((xx/0.45)-60)**2+(yy+20)**2+(zz/3)**2)<60**2
            z4 = 2*(((xx)**2+(yy/0.8-30)**2+(zz/3)**2)<15**2)
            img = z1.astype('float') + z2.astype('float')+ z3.astype('float')+ z4.astype('float')
        else:
            raise ValueError("unknown phantom")
        
        if display:
            import matplotlib.pyplot as plt
            slices = np.sort(np.random.randint(self.image.matrixSize[2],size=4))
            plt.figure()
            plt.subplot(2,2,1),plt.imshow(img[:,:,slices[0]]),plt.title('Slice: {}'.format(str(slices[0])),fontsize=15)
            plt.subplot(2,2,2),plt.imshow(img[:,:,slices[1]]),plt.title('Slice: {}'.format(str(slices[1])),fontsize=15)
            plt.subplot(2,2,3),plt.imshow(img[:,:,slices[2]]),plt.title('Slice: {}'.format(str(slices[2])),fontsize=15)
            plt.subplot(2,2,4),plt.imshow(img[:,:,slices[3]]),plt.title('Slice: {}'.format(str(slices[3])),fontsize=15)
            plt.show()
        return img