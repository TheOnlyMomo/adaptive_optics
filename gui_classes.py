import numpy as np

class sci_cameras:
    def __init__(self,name,n_pixels,qe,bits,fw,dc,rn,sens):
        self.name = name
        self.n_pixels = n_pixels
        self.qe = qe
        self.bits = bits
        self.fw = fw
        self.dc = dc
        self.rn = rn
        self.sens = sens
    
    def saturation(self):
        return (self.n_pixels)*((2**self.bits)*self.sens)/self.qe
    
class wfs_cameras:
    def __init__(self,name,n_pixels,qe,bits,fw,dc,rn,sens):
        self.name = name
        self.n_pixels = n_pixels
        self.qe = qe
        self.bits = bits
        self.fw = fw
        self.dc = dc
        self.rn = rn
        self.sens = sens
    
    def saturation(self):
        return (self.n_pixels)*((2**self.bits)*self.sens)/self.qe
    
class telescopes:
    def __init__(self, name, diameter, fov):
        self.name = name
        self.diameter = diameter
        self.fov = fov # In degrees


    def area(self):
        return (np.pi*(self.diameter/2)**2)*10000 # Converting to cm2 for CGS units
    
class shwf_sensor:
    def __init__(self,name,n_subaps, subap_fov):
        self.name = name
        self.n_subaps = n_subaps
        self.subap_fov = subap_fov

    def pxl_per_subap(self,camera_pixels):
        return camera_pixels/self.n_subaps
        
    
ASI6200 = sci_cameras('ASI6200',4096*2160, 0.9,16,352618,0.02,5, 0.8)
ASI6200_128 = wfs_cameras('ASI6200',128*128, 0.9,16,1504500,0.02,5, 0.8)
Marana42B6 = wfs_cameras('Marana42B6',128*128, 0.95,16,55000,0.15,1.6, 0.8)
telescope_1_2 = telescopes("1.2 M, FOV: 0.4",1.2,0.4)
telescope_0_75 = telescopes("0.75M, FOV: 0.4",0.75,0.4)
telescope_4 = telescopes("0.75M, FOV: 0.4",4,0.4)
default_shwf = shwf_sensor("SHWF 64 Subap",64, 10)
