import numpy as np
import matplotlib.pylab as pl
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
from astropy.time import Time
import scipy.special as scsp
import scipy.stats as scst
from datetime import datetime
from scipy.odr import odrpack as odr
from scipy.odr import models
from .extra import ChiSquare
import os
import warnings

def calc_fresnel(distance,lambida):
    """ Returns the fresnel scale.
    ----------
    Parameters
    ----------
    distance  (int, float, array): distances, in km.
    lambida   (int, float, array): wavelength, in km.
    ----------
    Return
    fresnel_scale  (float, array): fresnel scale, in km
    """
    return np.sqrt(lambida*distance/2)


def fit_pol(x,y,deg):
    """ Fit a polynom to the data.
    ----------
    Parameters
    ----------
    x  (array): 
    y  (array):
    deg  (int): order of the polynom 
    ----------
    Return
    ----------
    param (array): Array with the fitted values
    param_err (array): Array with the errors of the fitted values    
    """
    func   = models.polynomial(deg)
    mydata = odr.Data(x, y)
    myodr  = odr.ODR(mydata, func,maxit=200)
    myodr.set_job(fit_type=2)
    fit = myodr.run()
    param = fit.beta[::-1]                      
    param_err = fit.sd_beta[::-1]
    return param, param_err


class LightCurve():
    __names = []
    def __init__(self, name, **kwargs):
        """ Define a Light Curve

        Parameters:
            name (str): The name of the LightCurve. (required)
                Each time an LightCurve object is defined,
                the name must be different.
            tref (Time,str,float): Instant of reference. It can be
                in Julian Date, string in ISO format or Time object.
                Required only if LightCurve have input fluxes and
                given time is not in Julian Date.
            lambda (int,float): The band pass of observation (not required)
                in microns
            delta_lambda (int,float): The band width (not required.)

            Input data must be one of the 4 options below:
            file (str): a file with the time and flux in the first
                and second columns, respectively. A third columns
                with error in flux can also be given.

            IF file is not given:
            time: time must be a list of times, in seconds from tref,
                or Julian Date, or a Time object.
            flux: flux must be a list of fluxes. It must have the
                same lenght as time.
            dflux: if file not given, dflux must be a list of fluxes errors.
                It must have the same lenght as time. (not required)

            IF time and flux are not given.
            For a positive occultation
            immersion: The instant of immersion.
            emersion: The instant of emersion
            immersion_err: Immersion time uncertainty
            emersion_err: Emersion time uncertainty

            For a negative occultation
            initial_time: The initial time of observation
            end_time: The end time of observation.
        """
        input_done = False
        self.dflux = None
        self.__name = name
        self.flux = None
        self.time_model = None
        if self.__name in self.__names:
            raise ValueError('name {} already defined for another LightCurve object. Please choose a different one.'.format(self.__name))
        if 'immersion' in kwargs and 'emersion' in kwargs:
            if type(kwargs['immersion']) == str:
                self._immersion = Time(kwargs['immersion'])
            else:
                self._immersion = kwargs['immersion']
            if type(kwargs['emersion']) == str:
                self._emersion = Time(kwargs['emersion'])
            else:
                self._emersion = kwargs['emersion']
            self.immersion_err = 0.0
            if 'immersion_err' in kwargs:
                self.immersion_err = kwargs['immersion_err']
            self.emersion_err = 0.0
            if 'emersion_err' in kwargs:
                self.emersion_err = kwargs['emersion_err']
            input_done = True
        if 'initial_time' in kwargs and 'end_time' in kwargs:
            if type(kwargs['initial_time']) == str:
                self.initial_time = Time(kwargs['initial_time'])
            else:
                self.initial_time = kwargs['initial_time']
            if type(kwargs['end_time']) == str:
                self.end_time = Time(kwargs['end_time'])
            else:
                self.end_time = kwargs['end_time']
            input_done = True
        if not input_done:
            try:
                self.set_flux(**kwargs)
            except:
                raise ValueError('No allowed input conditions satisfied. Please refer to the tutorial.')
        self.lambda_0 = kwargs.get('lambda', 0.70)
        self.delta_lambda = kwargs.get('delta_lambda', 0.30)
        self.dt = 0.0
        self.__names.append(self.__name)

    @property
    def fresnel_scale(self):
        lamb  = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist  = self.dist*u.au.to('km')
        fresnel_scale_1 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist,lamb+dlamb/2.0)
        fresnel_scale   = (fresnel_scale_1 + fresnel_scale_2)/2.0
        return fresnel_scale

    @property
    def name(self):
        return self.__name

    @property
    def immersion(self):
        return self._immersion + self.dt*u.s

    @property
    def emersion(self):
        return self._emersion + self.dt*u.s

    def check_names(self):
        return self.__names
    
    def set_flux(self, **kwargs):
        """ Set the flux for the LightCurve

        Parameters:
            exptime (int,float): The exposure time of the observation. (required)
            file (str): a file with the time and flux in the first
                and second columns, respectively. A third columns
                with error in flux can also be given.
            time: if file not given, time must be a list of times,
                in seconds from tref, or Julian Date, or a Time object.
            flux: if file not given, flux must be a list of fluxes.
                It must have the same lenght as time.
            dflux: if file not given, dflux must be a list of fluxes errors.
                It must have the same lenght as time.
            tref (Time,str,float): Instant of reference. It can be
                in Julian Date, string in ISO format or Time object.
        """
        input_done = False
        if 'file' in kwargs:
            if not os.path.isfile(kwargs['file']):
                raise ValueError('{} not found'.format(kwargs['file']))
            try:
                time, self.flux, self.dflux = np.loadtxt(kwargs['file'], usecols=[0,1,2], unpack=True)
                self.flux_obs = self.flux
            except:
                time, self.flux = np.loadtxt(kwargs['file'], usecols=[0,1], unpack=True)
                self.flux_obs = self.flux
            else:
                raise ValueError('Input file must have 2 or 3 columns')
            input_done = True
        if 'time' in kwargs and 'flux' in kwargs:
            if input_done:
                raise ValueError('Only one type of input can be given. Please refer to the tutorial.')
            self.flux = kwargs['flux']
            time = kwargs['time']
            if len(self.flux) != len(time):
                raise ValueError('time and flux must have the same length')
            if 'dflux' in kwargs:
                self.dflux = kwargs['dflux']
                if len(self.flux) != len(self.dflux):
                    raise ValueError('dflux must have the same length as flux and time')
            input_done = True
        if 'exptime' not in kwargs:
            raise ValueError('exptime not defined')
        self.exptime = kwargs['exptime']
        if 'tref' in kwargs:
            try:
                if type(kwargs['tref']) in [Time,str]:
                    self.tref = Time(kwargs['tref'])
                elif type(kwargs['tref']) in [int,float]:
                    self.tref = Time(kwargs['tref'], format='jd')
            except:
                raise ValueError('tref must be an Julian Date or ISO format Date')
        if 'time' in locals():
            if type(time) == Time:
                if not 'tref' in kwargs:
                    self.tref = Time(time[0].iso.split(' ')[0] + ' 00:00:00.000')
            elif all(time > 2400000):
                time = Time(time, format='jd')
                if not 'tref' in kwargs:
                    self.tref = Time(time[0].iso.split(' ')[0] + ' 00:00:00.000')
            elif not hasattr(self,'tref'):
                raise ValueError('tref must be given')
            else:
                time = self.tref + time*u.s
            self.time = (time - self.tref).sec
            order = np.argsort(self.time)
            self.model = np.ones(len(self.time))
            self.flux = self.flux[order]
            self.flux_obs = self.flux
            self.time = self.time[order]
            if self.dflux != None:
                self.dflux = self.dflux[order]
            self.initial_time = np.min(time)
            self.end_time = np.max(time)
            self.cycle = np.median(self.time[1:] - self.time[:-1])
            if self.cycle < self.exptime:
                warnings.warn('Exposure time ({:0.4f} seconds) higher than Cycle time ({:0.4f} seconds)'.format(self.exptime,self.cycle))

    def set_vel(self,vel):
        '''
        Set the occultation velocity
        Inputs:
        vel = float, in km/s
        '''
        if type(vel) == u.quantity.Quantity:
            vel = vel.to(u.km/u.s).value
        elif type(vel) in [float,int]:
            pass
        else:
            raise TypeError('vel must be an integer, a float or an Astropy Unit object')
        self.vel = np.absolute(vel)

    def set_dist(self,dist):
        '''
        Set the object distance
        Inputs:
        dist = float, in km
        '''
        if type(dist) == u.quantity.Quantity:
            dist = dist.to(u.AU).value
        elif type(dist) in [float,int]:
            pass
        else:
            raise TypeError('dist must be an integer, a float or an Astropy Unit object')
        self.dist = dist

    def set_star_diam(self,d_star):
        '''
        Set the star diameter
        Inputs:
        diam = float, in km
        '''
        if type(d_star) == u.quantity.Quantity:
            d_star = d_star.to(u.km).value
        elif type(d_star) in [float,int]:
            pass
        else:
            raise TypeError('d_star must be an integer, a float or an Astropy Unit object')
        self.d_star = d_star

    def set_filter(self,lambda_0,delta_lambda):
        '''
        Set the filter bandwidth in microns
        Inputs:
        lambda_0 = float, in microns
        delta_lambda = float, in microns
        '''
        if type(lambda_0) == u.quantity.Quantity:
            lambda_0 = lambda_0.to(u.micrometer).value
        elif type(lambda_0) in [float]:
            pass
        else:
            raise TypeError('lambda_0 must be a float or an Astropy Unit object')
        self.lambda_0 = lambda_0
        if type(delta_lambda) == u.quantity.Quantity:
            delta_lambda = delta_lambda.to(u.micrometer).value
        elif type(delta_lambda) in [float]:
            pass
        else:
            raise TypeError('delta_lambda must be a float or an Astropy Unit object')
        self.delta_lambda = delta_lambda

        
    def calc_magnitude_drop(self,mag_star,mag_obj):
        """
        Determine the magnitude drop
        ----------
        Parameters
        ----------
        mag_star (float): Stellar magnitude.
        mag_obj  (float): Object apparent magnitude to the date.
        ----------
        Returns
        ----------
        mag_drop (float): Magnitude drop for the given magnitudes     
        bottom_flux (float): Normalized bottom flux for the given magnitudes     
        """
        contrast= 1/(1+(10**((mag_star-mag_obj)*0.4)))
        mag_combined = mag_star-(2.5*(np.log10(1/contrast)))
        mag_drop = mag_obj - mag_combined
        bottom_flux = 10**((mag_combined - mag_obj)*0.4)
        self.mag_drop = mag_drop
        self.bottom_flux = bottom_flux
        return 
 
    def normalize(self,poly_deg=None,mask=None,**kwargs):
        """ Returns the fresnel scale.
        ----------
        Parameters
        ----------
        poly_deg  (int): degree of the polynom to be fitted
        mask   (array of Bolleans): which values to be fitted 
        ----------
        kwargs
        ----------
        flux_min (int,float): event flux to be setted as 0.0
        flux_max (int,float): baseline flux to be setted as 1.0        
        plot (Bollean): If True plot the steps for visual aid  
        """
        #Create a mask where the polynomial fit will be done
        if np.all(self.flux) == None:
            raise ValueError('Normalization is only possible when a LightCurve is instatiated with time and flux.')
        self.reset_flux()
        flux_min = 0.0
        flux_max = 1.0
        plot = False
        if 'flux_min' in kwargs:
            flux_min = kwargs['flux_min']
        if 'flux_max' in kwargs:
            flux_max = kwargs['flux_max']
        if 'plot' in kwargs:
            if (kwargs['plot'] == True):
                plot = True
        lc_flux = (self.flux - flux_min)/(flux_max-flux_min)
        if (np.all(mask) == None):
            preliminar_occ = self.occ_detect()
            tmax = preliminar_occ['emersion_time']+1.00*preliminar_occ['occultation_duration']
            tmin = preliminar_occ['immersion_time']-1.00*preliminar_occ['occultation_duration']
            chord = preliminar_occ['occultation_duration']
            mask = np.invert((self.time > tmin-(chord/2)) & (self.time < tmax+(chord/2)))
        norm_time = (self.time - self.time.min())/(self.time.max()-self.time.min())
        if (poly_deg != None):
            n = poly_deg
            p, err = fit_pol(norm_time[mask],lc_flux[mask],n)
            flux_poly_model = np.zeros(len(norm_time))
            for ii in np.arange(n+1):
                flux_poly_model = flux_poly_model + p[ii]*(norm_time**(n-ii))
            if (plot == True):
                pl.plot(norm_time[mask],lc_flux[mask],'k.-')
                pl.plot(norm_time[mask],flux_poly_model[mask],'r-')
                pl.title('Polynomial degree = {}'.format(n),fontsize=15)
                pl.show()       
        if (poly_deg == None):
            n = 0
            p, err = fit_pol(norm_time[mask],lc_flux[mask],n)
            flux_poly_model = np.zeros(len(norm_time))
            for ii in np.arange(n+1):
                flux_poly_model += p[ii]*(norm_time**(n-ii))
            if (plot == True):
                pl.plot(norm_time[mask],lc_flux[mask],'k.-')
                pl.plot(norm_time[mask],flux_poly_model[mask],'r-')
                pl.title('Polynomial degree = {}'.format(n),fontsize=15)
                pl.show()
            flux_poly_model_old = flux_poly_model.copy()
            for nn in np.arange(1,10):
                p, err = fit_pol(norm_time[mask],lc_flux[mask],nn)
                flux_poly_model_new = np.zeros(len(norm_time))
                for ii in np.arange(nn+1):
                    flux_poly_model_new += p[ii]*(norm_time**(nn-ii))
                F = np.var(flux_poly_model[mask]-lc_flux[mask])/np.var(flux_poly_model_new[mask]-lc_flux[mask])
                if (F > 1.05):
                    flux_poly_model = flux_poly_model_new.copy()
                    n = nn
                    if (plot == True):
                        pl.plot(norm_time[mask],lc_flux[mask],'k.-')
                        pl.plot(norm_time[mask],flux_poly_model[mask],'r-')
                        pl.title('Polynomial degree = {}'.format(nn),fontsize=15)
                        pl.show()
                else:
                    print('Normalization using a {} degree polynom'.format(n))
                    print('There is no improvement with a {} degree polynom'.format(n+1))
                    break
        self.flux = lc_flux/flux_poly_model
        self.normalizer_flux = flux_poly_model
        self.normalizer_mask = mask
        return

    def reset_flux(self):
        '''
        Reset flux for original values
        '''
        try:
            self.flux = self.flux_obs
        except:
            raise ValueError('Reset is only possible when a LightCurve is instatiated with time and flux.')
        return

    
    def occ_model(self,immersion_time, emersion_time, opacity, mask, npt_star=12, time_resolution_factor=10,flux_min=0,flux_max=1):
        """ Returns the modelled light curve considering fresnel difraction, star diameter and intrumental response.
        ----------
        Parameters
        ----------
        immersion_time (int, float): Immersion time, in seconds.                (input)
        emersion_time  (int, float): Emersion time, in seconds.                 (input)
        opacity (int, float): Opacity, opaque = 1.0, transparent = 0.0.         (input)
        mask (array with Booleans): Mask with True values to be computed        (input)
        npt_star  (int): Number of subdivisions for computing the star size's effects, default equal to 12. (auto)
        time_resolution_factor (int,float): Steps for fresnel scale used for modelling the light curve,     (auto)
            default equals to 10 steps for fresnel scale.
        flux_min (int,float): Bottom flux (only object), default equal to 0.0              (auto)
        flux_max (int,float): Base flux (object plus star), default equal to 1.0           (auto)
        ----------
        Returns
        ----------
        self.flux_inst       (array): Modelled Instrumental light flux.
        self.time_model      (array): Modelled timing.
        self.flux_star       (array): Modelled light flux considering fresnel difraction and star's diameter.
        self.flux_fresnel    (array): Modelled light flux considering fresnel difraction.
        self.model_geometric (array): Modelled light flux considering a box model.
        """
        # Computing the fresnel scale
        lamb  = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist  = self.dist*u.au.to('km')
        vel = np.absolute(self.vel)
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist,lamb+dlamb/2.0)
        fresnel_scale   = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/vel,self.exptime]))/time_resolution_factor
        #
        #Creating a high resolution curve to compute fresnel difraction, stellar diameter and instrumental integration 
        time_model = np.arange(time_obs.min()-5*self.exptime,time_obs.max()+5*self.exptime,time_resolution)
        #
        #Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x   = time_model*vel    
        x01 = immersion_time*vel
        x02 = emersion_time*vel
        #
        #Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = self.__bar_fresnel(x,x01,x02,fresnel_scale_1,opacity)
        flux_fresnel_2 = self.__bar_fresnel(x,x01,x02,fresnel_scale_2,opacity)
        flux_fresnel   = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star      = flux_fresnel.copy()
        if (self.d_star > 0):
            #Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao   = self.d_star/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            #Computing stellar diameter only near the immersion or emersion times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star,npt_star)*resolucao
            coeff = np.sqrt(np.absolute(self.d_star**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_1,opacity)
                flux2 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_2,opacity)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii]   = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs)) 
        for i in range(len(time_obs)):
            event_model  = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        self.model[mask] = flux_inst*(flux_max - flux_min) + flux_min
        self.time_model = time_model
        self.model_star = flux_star*(flux_max - flux_min) + flux_min
        self.model_fresnel = flux_fresnel*(flux_max - flux_min) + flux_min
        ev_model = (time_model > immersion_time) & (time_model < emersion_time)
        flux_box = np.ones(len(time_model))
        flux_box[ev_model] = (1-opacity)**2
        flux_box = flux_box*(flux_max - flux_min) + flux_min
        self.model_geometric = flux_box
        self.baseflux = flux_max
        self.bottomflux = flux_min
        return


  
    def occ_lcfit(self, **kwargs):
        """ Monte Carlo chi square fit for occultations lightcurve.
        ----------
        Parameters
        ----------
        tmin (int,float): Minimum time to consider in the fit prcedure                  
        tmax (int,float): Maximum time to consider in the fit prcedure                  
        flux_min (int,float): Bottom flux (only object), default equal to 0.0          
        flux_max (int,float): Base flux (object plus star), default equal to 1.0       
        immersion_time  (int, float): Initial guess for immersion time, in seconds.                       
        emersion_time   (int, float): Initial guess for emersion time, in seconds.                        
        opacity  (int, float): Initial guess for opacity, opaque = 1.0, transparent = 0.0, default equal to 0.0.   
        delta_t (int, float): Interval to fit immersion or emersion time               
        dopacity   (int, float): Interval to fit opacity, default equal to 0, no fit.  
        loop       (int): Number of tests to be done, default equal to 10000.          
        ----------
        Returns
        ----------
        chi2 (ChiSquare): ChiSquare object
        """
        #
        if not hasattr(self, 'flux'):
            raise ValueError('Fit curve is only possible when a LightCurve is instatiated with time and flux.')
        delta_t = 2*self.cycle
        loop = 10000
        t_i = np.zeros(loop)
        t_e = np.zeros(loop)
        tmax = self.time.max()
        tmin = self.time.min()
        immersion_time = tmin - self.exptime
        do_immersion = False
        emersion_time = tmax + self.exptime
        do_emersion = False
        opacity = 1.0
        delta_opacity = 0.0
        do_opacity = False
        if 'loop' in kwargs:
            loop = kwargs['loop']
        if ('immersion_time' not in kwargs) and ('emersion_time' not in kwargs):
            preliminar_occ = self.occ_detect()
            immersion_time = preliminar_occ['immersion_time']
            do_immersion = True
            emersion_time = preliminar_occ['emersion_time']
            do_emersion = True
            delta_t = 5*preliminar_occ['time_err']
            tmax = emersion_time+2*preliminar_occ['occultation_duration']
            tmin = immersion_time-2*preliminar_occ['occultation_duration']
            if 2*preliminar_occ['occultation_duration'] < 10*self.cycle:
                tmax = emersion_time + 10*self.cycle
                tmin = immersion_time - 10*self.cycle
        if 'tmax' in kwargs:
            tmax = kwargs['tmax']
        if 'tmin' in kwargs:
            tmin = kwargs['tmin']
        if 'delta_t' in kwargs:
            delta_t = kwargs['delta_t']
        if 'immersion_time' in kwargs:
            immersion_time = kwargs['immersion_time']
            do_immersion = True
        t_i = immersion_time + delta_t*(2*np.random.random(loop) - 1)
        if 'emersion_time' in kwargs:
            emersion_time = kwargs['emersion_time']
            do_emersion = True
        t_e = emersion_time  + delta_t*(2*np.random.random(loop) - 1)
        mask = (self.time >= tmin) & (self.time <= tmax)
        mask_sigma = (((self.time >= tmin) & (self.time < immersion_time - self.exptime)) +
                      ((self.time > emersion_time + self.exptime) & (self.time <= tmax)))
        sigma = self.flux[mask_sigma].std(ddof=1)
        if 'sigma' in kwargs:
            sigma = kwargs['sigma']
        if 'opacity' in kwargs:
            opacity = kwargs['opacity']
        if 'dopacity' in kwargs:
            delta_opacity = kwargs['dopacity']
            do_opacity = True
        opas = opacity + delta_opacity*(2*np.random.random(loop) - 1)
        opas[opas>1.], opas[opas<0.] = 1.0, 0.0
        flux_min = 0
        flux_max = 1
        if 'flux_min' in kwargs:
            flux_min = kwargs['flux_min']
        if 'flux_max' in kwargs:
            flux_max = kwargs['flux_max']
        #
        tflag = np.zeros(loop)
        tflag[t_i > t_e] = t_i[t_i > t_e]
        t_i[t_i > t_e]    = t_e[t_i > t_e]
        t_e[t_i > t_e]    = tflag[t_i > t_e]
        chi2 = 999999*np.ones(loop)
        #tcontrol_f0 = datetime.now()
        for i in range(loop):
            model_test = self.__occ_model(t_i[i],t_e[i],opas[i],mask,flux_min=flux_min,flux_max=flux_max)
            chi2[i] = np.sum((self.flux[mask] -  model_test)**2)/(sigma**2)
        #tcontrol_f3 = datetime.now()
        #print('Elapsed time: {:.3f} seconds.'.format((tcontrol_f3 - tcontrol_f0).total_seconds()))
        kkargs = {}
        if do_immersion:
            kkargs['immersion'] = t_i
        if do_emersion:
            kkargs['emersion'] = t_e
        if do_opacity:
            kkargs['opacity'] = opas
        chisquare = ChiSquare(chi2, len(self.flux[mask]), **kkargs)
        onesigma = chisquare.get_nsigma(1)
        if 'immersion' in onesigma:
            self._immersion = self.tref + onesigma['immersion'][0]*u.s
            self.immersion_err = onesigma['immersion'][1]
            immersion_time = onesigma['immersion'][0]
        if 'emersion' in onesigma:
            self._emersion = self.tref + onesigma['emersion'][0]*u.s
            self.emersion_err = onesigma['emersion'][1]
            emersion_time = onesigma['emersion'][0]
        if 'opacity' in onesigma:
            opacity = onesigma['opacity'][0]
        # Run occ_model() to save best parameters in the Object.
        self.occ_model(immersion_time,emersion_time,opacity,mask,flux_min=flux_min,flux_max=flux_max)
        self.lc_sigma = sigma
        self.chisquare = chisquare
        return chisquare

    def plot_lc(self):
        '''
        Plot the light curve
        '''
        if (np.any(self.flux) != None):
            pl.close()
            pl.plot(self.time,self.flux,'k.-',label='Obs.',zorder=0)
            if (np.any(self.model) != None):
                pl.plot(self.time,self.model,'r-',label='Model',zorder=2)
                pl.scatter(self.time,self.model, s=50, facecolors='none', edgecolors='r',zorder=3)
            pl.tight_layout()
            pl.xlabel('Time [seconds]',fontsize=20)
            pl.ylabel('Relative Flux',fontsize=20)
            pl.legend()
        else:
            raise ValueError('Plotting the light curve is only possible when the Object LightCurve is instatiated with time and flux')
        return

    def plot_model(self):
        '''
        Plot the modelled light curve
        '''
        if np.all(self.time_model) != None:
            pl.plot(self.time_model,self.model_geometric,'c-',label='Geometric',zorder=1)
            pl.plot(self.time_model,self.model_fresnel,'b-',label='Fresnel',zorder=1)
            pl.plot(self.time_model,self.model_star,'g-',label='Star diam.',zorder=1)
            pl.tight_layout()
            pl.xlabel('Time [seconds]',fontsize=20)
            pl.ylabel('Relative Flux',fontsize=20)
            pl.legend()
        else:
            raise ValueError('Plotting the model light curve is only possible after the model [LightCurve.occ_model()] or the fit [LightCurve.occ_lcfit()]')
        return

    def to_log(self,namefile=None):
        """ Save the light curve log to a file

        Parameters:
            namefile (str): Filename to save the log
        """
        if (namefile == None):
            namefile = self.name.replace(' ','_')+'.log'
        f = open(namefile, 'w')
        f.write(self.__str__())
        f.close()
        

    def to_file(self, namefile=None):
        """ Save the light curve to a file

        Parameters:
            namefile (str): Filename to save the data
        """
        ###### Observational data
        if (namefile == None):
            folder=''
            file = self.name.replace(' ','_')+'.dat'
        else:
            folder = os.path.dirname(namefile)
            file = os.path.basename(namefile)
        data = np.array([(self.time*u.s + self.tref).jd,self.time,self.flux,self.model,self.flux-self.model])
        colunm_names = ['Time JD','Time relative to {} UTC in seconds'.format(self.tref.iso),'Observational Flux','Modelled Flux','Residual O-C']
        np.savetxt(os.path.join(folder,file), data.T, fmt='%11.8f')
        f = open(os.path.join(folder,file) + '.label', 'w')
        for i, name in enumerate(colunm_names):
            f.write('Column {}: {}\n'.format(i+1,name))
        f.close()
        ###### Complete Model
        if (np.all(self.time_model) != None):
            data_model = np.array([(self.time_model*u.s + self.tref).jd,self.time_model,self.model_geometric,self.model_fresnel,self.model_star])
            colunm_names_model = ['Model time JD','Model time relative to {} UTC in seconds'.format(self.tref.iso),'Geometric Model','Model with Fresnel diffraction','Model with star diameter']
            np.savetxt(os.path.join(folder,'model_'+file), data_model.T, fmt='%11.8f')
            f = open(os.path.join(folder,'model_'+file)+'.label', 'w')
            for i, name in enumerate(colunm_names_model):
                f.write('Column {}: {}\n'.format(i+1,name))
            f.close()



    
    def occ_detect(self, maximum_duration=None, dur_step=None, snr_limit=None, \
                  n_detections=None, plot=False):
        """
        Detect automatically the occultation event in the light curve
        
        Parameters
        ----------
        maximum_duration (float): Maximum duration of the occultation event (default is the light curve's time span). (optional)
        dur_step (float): Step size to sweep occultation duration event (default value is 1/2 of sampling). (optional)
        snr_limit (float): Minimum occultation SNR. (optional)
        n_detections (int): N best detections regardless from SNR. n_detections is superseded by snr_limit. (optional)
        plot (boolean): True if output plots are desired. (optional)
            
        Returns
        -------
        OrderedDict
            An ordered dictionary of :attr:`name`::attr:`value` pairs for each Parameter.
            
        Examples
        --------
        >>> lc = sora.LightCurve(time=time, flux=flux, exptime=0.0, name='lc_example')
        >>> params = lc.occ_detect()
        >>> params
        {'rank': 1, 
        'occultation_duration': 40.1384063065052, 
        'central_time': 7916.773870512843, 
        'immersion_time': 7896.7046673595905, 
        'emersion_time': 7936.843073666096, 
        'time_err': 0.05011036992073059, 
        'depth': 0.8663887801707082, 
        'depth_err': 0.10986223384336465, 
        'baseline': 0.9110181732552853, 
        'baseline_err': 0.19045768512595365, 
        'snr': 7.886138392251848,
        'occ_mask': array([False, False, False, ..., False, False, False])}
        """

        if not hasattr(self, 'flux'):
            raise ValueError('time and flux must be instantiated to use '\
                             'occ_detect function.')

        # duration of the light curve
        time_span = self.time[-1]-self.time[0]
        if maximum_duration and (maximum_duration > time_span):
            warnings.warn('Occultation duration (maximum_duration={0}) ' \
                             'exceeds the time series lenght ({1:0.5f}).' \
                          ' maximum_duration reset to the time series lenght.'\
                          .format(maximum_duration, time_span))
            maximum_duration = time_span
        if not maximum_duration:
            maximum_duration = time_span

        if not dur_step:
            dur_step = self.cycle/2
        
        if dur_step < self.cycle/2:
            warnings.warn('The given dur_step is oversampled by a factor ' \
                          'of {0:0.1f} and has been reset to half a cycle ' \
                          '({1:0.4f}).' \
                          .format((self.cycle/2.)/dur_step, self.cycle/2.))
            dur_step = self.cycle/2
       
        duration_grid = np.arange(dur_step, maximum_duration, dur_step)
        # initial occultation mask (all data points)
        mask = np.ones(len(self.time), dtype=bool)
        # inital detection rank
        rank = 1

        if snr_limit:
            # minimum SNR accepted in a detection for multiple search
            snr_value = snr_limit+1
            occ0 = self.__run_bls(time_span, duration_grid)
            mask *= ~occ0['occ_mask']
            while (snr_value > snr_limit):
                rank += 1
                occ1 = self.__run_bls(time_span, duration_grid, mask=mask, \
                                      rank=rank)
                if occ1['snr'] > snr_limit:
                    snr_value = occ1['snr']
                    mask *= ~occ1['occ_mask']                    
                    occ0 = self.__summarize_bls(occ0,occ1)
                else:
                    snr_value = snr_limit
            
            if plot:
                self.__plot_occ_detect(occ0)                    
            return occ0
        elif n_detections:
            # search the n best fits
            occ0 = self.__run_bls(time_span, duration_grid)
            mask *= ~occ0['occ_mask']    
            for i in range(n_detections-1):
                rank += 1
                occ1 = self.__run_bls(time_span, duration_grid, mask=mask, \
                                      rank=rank)
                snr_value = occ1['snr']
                mask *= ~occ1['occ_mask']                    
                occ0 = self.__summarize_bls(occ0,occ1)
            
            if plot:
                self.__plot_occ_detect(occ0)
            return occ0
        else:
            # search only the first best fit
            occ0 = self.__run_bls(time_span, duration_grid)
           
            if plot:
                self.__plot_occ_detect(occ0)
                
            return occ0
        
    def __plot_occ_detect(self, occ):
        n = np.size(occ['rank'])
        if n > 1:
            # case for multiple detections
            pl.plot(self.time, self.flux, 'k.-')
            mask = np.zeros(len(self.time), dtype=bool)
            for i in range(n):
                trues = np.sum(occ['occ_mask'][i])
                pl.plot(self.time[occ['occ_mask'][i]], np.repeat(np.mean(self.flux[occ['occ_mask'][i]]),trues), '.', label='Rank: '+str(i+1))
                mask += occ['occ_mask'][i]
 
            falses = list(mask).count(False)
            pl.plot(self.time[~mask], np.repeat(np.mean(self.flux[~mask]),falses), 'r.', label='Baseline')
            pl.xlabel('Time [seconds]')
            pl.ylabel('Relative Flux')
            pl.legend()    
            
        else:
            # case for single occultation
            trues = list(occ['occ_mask']).count(True)
            falses = list(occ['occ_mask']).count(False)
            pl.plot(self.time, self.flux, 'k.-')
            pl.plot(self.time[occ['occ_mask']], np.repeat(np.mean(self.flux[occ['occ_mask']]),trues), '.', label='Occultation')
            pl.plot(self.time[~occ['occ_mask']], np.repeat(np.mean(self.flux[~occ['occ_mask']]),falses), 'r.', label='Baseline')
            pl.xlabel('Time [seconds]')
            pl.ylabel('Relative Flux')
            pl.legend()
            
    
    def __run_bls(self, per_grid, dur_grid, mask=None, rank=None):
        """
        Private function to find the best box fit suitable to the data
        """

        # object with no occultation mask
        mskmodel = BoxLeastSquares(self.time, self.flux, dy=self.dflux)
        # if there is no dflux array, reset it to None in case of 
        # using a mask
        if self.dflux is None:
            dfluxmask = None
        else:
            dfluxmask = self.dflux[mask]
            
        # object with occultation mask
        if np.sum(mask):
            model = BoxLeastSquares(self.time[mask], self.flux[mask], \
                                    dy=dfluxmask)
        else:
            model = mskmodel   
        
        r = model.power(per_grid, dur_grid, objective='snr', method='fast')
        # statistics of the BLS fit
        stats = model.compute_stats(r.period, r.duration, r.transit_time)
        # occultation mask of the event with respect to all data
        occ_mask = mskmodel.transit_mask(self.time, r.period, r.duration, \
                                         r.transit_time)
        # parameters computation for clarity purposes
        occultation_duration = r.duration[0]
        central_time = stats['transit_times'][0]
        immersion_time = stats['transit_times'][0] - r.duration[0]/2
        emersion_time = stats['transit_times'][0] + r.duration[0]/2
        time_err = np.median(self.time[1:-1]-self.time[0:-2])/2
        depth = np.mean(self.flux[~occ_mask])-np.mean(self.flux[occ_mask])
        depth_err = np.std(self.flux[occ_mask],ddof=1)
        baseline = np.mean(self.flux[~occ_mask])
        baseline_err = np.std(self.flux[~occ_mask],ddof=1)
        # If there is only one measurement during the occultation it will
        # use the baseline_err to compute SNR, otherwise it will use depth_err
        if np.sum(occ_mask) < 2:
            snr = depth/baseline_err
        else:
            snr = depth/depth_err


        # define rank
        if rank:
            rank = rank
        else:
            rank = 1
        
        return {'rank' : rank,
                'occultation_duration' : occultation_duration, 
                'central_time' : central_time, 
                'immersion_time' : immersion_time, 
                'emersion_time' : emersion_time,
                'time_err' : time_err,
                'depth' : depth, 'depth_err' : depth_err,
                'baseline' : baseline, 'baseline_err' : baseline_err,
                'snr' : snr,'occ_mask' : occ_mask}
    
    def __summarize_bls(self, dict1, dict2):
        ''' Private function to merge dictionaries returned by BLS and 
            to keep values of common keys in list.
        '''
        dict3 = {}
        for key, value in dict1.items():
            if key == 'occ_mask':
                sz = int(np.size(dict1[key])/np.size(dict2[key]))
                if sz > 1:
                    l = [None]*(sz+1)
                    for i in range(sz):
                        l[i] = dict1[key][i]
                    l[i+1] = dict2[key]
                    dict3[key] = l
                else:
                    dict3[key] = [dict1[key],dict2[key]]
            else:
                dict3[key] = np.append(dict1[key],dict2[key])
        return dict3

    def __bar_fresnel(self,X,X01,X02,fresnel_scale,opacity):
        """ Returns the modelled light curve considering fresnel difraction.
        ----------
        Parameters
        ----------
        X   (array): Array with time values converted in km using the event velocity.
        X01 (int, float): Immersion time converted in km using the event velocity.
        X02 (int, float): Emersion time converted in km using the event velocity.
        fresnel_scale (int, float): Fresnel scale.
        opacity (int, float): Opacity, opaque = 1.0, transparent = 0.0
        ----------
        Returns
        ----------
        flux_fresnel (array):
        """
        # Converting from km to units of fresnel scale
        x   = X/fresnel_scale
        x01 = X01/fresnel_scale
        x02 = X02/fresnel_scale
        # Fresnel difraction parameters 
        x1 = x - x01
        x2 = x - x02
        s1,c1 = scsp.fresnel(x1)
        s2,c2 = scsp.fresnel(x2)
        cc = c1 - c2
        ss = s1 - s2
        r_ampli = - (cc+ss)*(opacity/2.)
        i_ampli =   (cc-ss)*(opacity/2.)
        # Determining the flux considering fresnel difraction
        flux_fresnel = (1.0 + r_ampli)**2 + (i_ampli)**2
        return flux_fresnel


    def __occ_model(self,immersion_time, emersion_time, opacity, mask, npt_star=12, time_resolution_factor=10,flux_min=0.0,flux_max=1.0):
        """ Private function returns the modelled light curve considering fresnel difraction, star diameter and intrumental response, intended for fitting inside the self.occ_lcfit().
        ----------
        Parameters
        ----------
        immersion_time (int, float): Immersion time, in seconds.                (input)
        emersion_time  (int, float): Emersion time, in seconds.                 (input)
        opacity (int, float): Opacity, opaque = 1.0, transparent = 0.0.       (input)
        mask (array with Booleans): Mask with True values to be computed        (input)
        npt_star  (int): Number of subdivisions for computing the star size's effects, default equal to 12. (auto)
        time_resolution_factor (int,float): Steps for fresnel scale used for modelling the light curve,     (auto)
            default equals to 10 steps for fresnel scale.
        flux_min (int,float): Bottom flux (only object), default equal to 0.0              (auto)
        flux_max (int,float): Base flux (object plus star), default equal to 1.0           (auto)
        ----------
        Returns
        ----------
        flux_inst       (array): Modelled Instrumental light flux.
        """
        # Computing the fresnel scale
        lamb  = self.lambda_0*u.micrometer.to('km')
        dlamb = self.delta_lambda*u.micrometer.to('km')
        dist  = self.dist*u.au.to('km')
        vel = np.absolute(self.vel)
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist,lamb+dlamb/2.0)
        fresnel_scale   = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/vel,self.exptime]))/time_resolution_factor
        self.model_resolution = time_resolution
        #
        #Creating a high resolution curve to compute fresnel difraction, stellar diameter and instrumental integration 
        time_model = np.arange(time_obs.min()-5*self.exptime,time_obs.max()+5*self.exptime,time_resolution)
        #
        #Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x   = time_model*vel    
        x01 = immersion_time*vel
        x02 = emersion_time*vel
        #
        #Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = self.__bar_fresnel(x,x01,x02,fresnel_scale_1,opacity)
        flux_fresnel_2 = self.__bar_fresnel(x,x01,x02,fresnel_scale_2,opacity)
        flux_fresnel   = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star      = flux_fresnel.copy()
        if (self.d_star > 0):
            #Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao   = self.d_star/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            #Computing stellar diameter only near the immersion or emersion times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star,npt_star)*resolucao
            coeff = np.sqrt(np.absolute(self.d_star**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_1,opacity)
                flux2 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_2,opacity)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii]   = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs)) 
        for i in range(len(time_obs)):
            event_model  = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        return flux_inst*(flux_max - flux_min) + flux_min

    def __str__(self):
        """ String representation of the LightCurve Object
        """
        output  = 'Light curve name: {}\n'.format(self.name)
        output += 'Initial time: {} UTC\n'.format(self.initial_time.iso)
        output += 'End time:     {} UTC\n'.format(self.end_time.iso)
        output += 'Duration:     {:.3f} minutes\n'.format((self.end_time - self.initial_time).value*u.d.to('min'))
        output += 'Time offset:  {:.3f} seconds\n\n'.format(self.dt)
        try:
            output += 'Exposure time:    {:.4f} seconds\n'.format(self.exptime)
            output += 'Cycle time:       {:.4f} seconds\n'.format(self.cycle)
            output += 'Num. data points: {}\n\n'.format(len(self.time))
        except:
            output += 'Object LightCurve was not instantiated with time and flux.\n\n'
        try:
            output += 'Bandwidth:            {:.3f} +/- {:.3f} microns\n'.format(self.lambda_0,self.delta_lambda)
            output += 'Used shadow velocity: {:.3f} km/s\n'.format(self.vel)
            output += 'Fresnel scale:        {:.3f} seconds or {:.2f} km\n'.format(self.fresnel_scale/self.vel,self.fresnel_scale)
            output += 'Stellar size effect:  {:.3f} seconds or {:.2f} km\n'.format(self.d_star/self.vel,self.d_star)
        except:
            output += '\nThere is no occultation associated with this light curve.\n'
        try:
            output += 'Inst. response:       {:.3f} seconds or {:.2f} km\n'.format(self.exptime, self.exptime*self.vel)
            output += 'Dead time effect:     {:.3f} seconds or {:.2f} km\n'.format(self.cycle-self.exptime,(self.cycle-self.exptime)*self.vel)
            output += 'Model resolution:     {:.3f} seconds or {:.2f} km\n'.format(self.model_resolution, self.model_resolution*self.vel)
            output += 'Modelled baseflux:    {:.3f}\n'.format(self.baseflux)
            output += 'Modelled bottomflux:  {:.3f}\n'.format(self.bottomflux)
            output += 'Light curve sigma:    {:.3f}\n\n'.format(self.lc_sigma)
        except:
            output += '\nObject LightCurve model was not fitted.\n\n'
        try:
            output += 'Immersion time: {} UTC +/- {:.3f} seconds\n'.format(self.immersion.iso,self.immersion_err)
            output += 'Emersion time:  {} UTC +/- {:.3f} seconds\n\n'.format(self.emersion.iso,self.emersion_err)
        except:
            output += 'Immersion and emersion times were not fitted or instantiated.\n\n'

        try:
            output += 'Monte Carlo chi square fit.\n\n' + self.chisquare.__str__() + '\n'
        except:
            pass
        return output


    def __del__(self):
        """ When this object is deleted, it removes the name from the Class name list.
        """
        try:
            self.__names.remove(self.__name)
        except:
            pass
