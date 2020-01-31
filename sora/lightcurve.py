import numpy as np
import pylab as pl
import astropy.units as u


class LightCurve():
    '''
    Docstring
    Define a Light Curve
    '''
    def __init__(self,time,flux,**kwargs):
        if len(flux) != len(time):
            print('The time and the flux should have the same length')
            return
        self.flux = flux
        self.time = time
        self.exptime = None
        self.model = np.ones(len(self.time))
        self.lambda_0 = 0.70 #microns #### *u.micrometer.to('km')
        self.delta_lambda = 0.30 #microns #### *u.micrometer.to('km')
        self.sigma = None
        self.tref = None
        return
    
    def set_lightcurve(self,time,flux):
        '''
        Set the light curve
        Inputs:
        time = numpy array, in seconds relative to reference time (self.tref)
        flux = numpy array, in flux
        '''
        self.time = time
        self.flux = flux
        return
    
    def set_filter(self,lambda_0,delta_lambda):
        '''
        Set the filter utilized in the observation centred at lambda0 with a width of delta_lambda. 
        Inputs:
        lambda_0 = float, in microns
        delta_lambda = float, in microns
        '''
        self.lambda_0 = lambda_0 #microns
        self.delta_lambda = delta_lambda #microns
        return
    
    def set_exposure(self,exp):
        '''
        Set the exposure time utilized
        Inputs:
        exp = float, in seconds
        '''
        self.exptime = np.float(exp)
        return
    
    def set_sigma_noise(self,sigma):
        '''
        Set the standard deviation of the light-curve (sigma)
        Inputs:
        sigma = float, in flux
        '''
        self.sigma = sigma
        return
    
    def calc_sigma_noise(self):
        '''
        Determine the standard deviation of the light-curve (sigma)
        '''
        if (np.any(self.flux) != None):
            cut = np.absolute(self.flux-1) < np.std(self.flux)
            for loop in range(10):
                cut = np.absolute(self.flux-1) < 3*np.std(self.flux[cut])
            self.sigma = np.std(self.flux[cut])
        return
    
    
    def calc_magnitude_drop(self,star,obj):
        '''
        Determine the magnitude drop
        Inputs:
        star = Star()
        obj = Ephemerides()
        '''
        ###DEVE SER MODIFICADO PARA ENTRAR NO OCCULTATION [????]
        self.contrast= 1/(1+(10**((star.magV-obj.magV)*0.4)))
        self.mag_combined = star.magV-(2.5*(np.log10(1/self.contrast)))
        self.mag_drop = obj.magV - self.mag_combined
        self.bottom_flux = 10**((self.mag_combined - obj.magV)*0.4)
        return
    
    
    def calc_occmodel(self,t_ingress,t_egress,opa_ampli=1.0,baseflux=1.0,depht=0.0):
        #INCLUIR MODELO CORRETO!!!!!!!
        if (np.any(self.time) != None):
            self.model = baseflux*np.ones(len(self.time))
            event_id = (self.time > t_ingress) & (self.time < t_egress)
            self.model[event_id] = depht
        return

    def plot_lc(self,fig_name=None):
        '''
        Plot the light curve if you want to save the file, the fig_name should be different than None
        '''
        if (np.any(self.time) != None):
            pl.close()
            pl.figure(figsize=[8.4, 3.4])
            pl.plot(self.time,self.flux,'k.-',label='Obs.',zorder=1)
            if (np.any(self.model) != None):
                pl.plot(self.time,self.model,'r-',label='Model',zorder=2)
                pl.scatter(self.time,self.model, s=50, facecolors='none', edgecolors='r',zorder=3)
            pl.tight_layout()
            pl.xlabel('Time [seconds]',fontsize=20)
            pl.ylabel('Relative Flux',fontsize=20)
            pl.legend(fontsize=20)
            pl.xticks(fontsize=20)
            pl.yticks(fontsize=20)
            if (fig_name == None):
                pl.show()
            else:
                pl.savefig(fig_name,format='png', dpi=200)
        else:
            print('There is no lightcurve to plot')
        return

