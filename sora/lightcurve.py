import numpy as np
import pylab as pl
import astropy.units as u
from astropy.timeseries import BoxLeastSquares
import scipy.special as scsp
from datetime import datetime
from .extra import ChiSquare

def calc_fresnel(distance,lambida):
    """ Returns the fresnel scale.
    ----------
    Parameters
    ----------
    distance  (int, float, array): distances, in km.
    lambida   (int, float, array): wavelength, in km.
    ----------
    """
    return np.sqrt(lambida*distance/2)


class LightCurve():
    '''
    Docstring
    Define a Light Curve
    '''
    def __init__(self, time, flux, dflux=None, **kwargs):
        if len(flux) != len(time):
            raise ValueError('The time and the flux should have the same length')
        self.flux = flux
        self.dflux = dflux
        self.time = time
        self.exptime = None
        self.model = np.ones(len(self.time))
        self.lambda_0 = 0.70 #microns #### *u.micrometer.to('km')
        self.delta_lambda = 0.30 #microns #### *u.micrometer.to('km')
        self.sigma = None
        self.tref = None
        self.dist = None
        self.d_star = None
        self.vel = None
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
        
    def set_exposure(self,exp):
        '''
        Set the exposure time utilized
        Inputs:
        exp = float, in seconds
        '''
        self.exptime = np.float(exp)
        return
    
    def set_vel(self,vel):
        '''
        Set the occultation velocity
        Inputs:
        vel = float, in km/s
        '''
        if type(vel) == u.quantity.Quantity:
            vel = vel.to(u.km/u.s).value
        elif type(vel) == [float,int]:
            pass
        else:
            raise TypeError('vel must be a float or an Astropy Unit object')
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
            raise TypeError('dist must be a float or an Astropy Unit object')
        self.dist = dist

    def set_diam(self,diam):
        '''
        Set the star diameter
        Inputs:
        diam = float, in km
        '''
        if type(diam) == u.quantity.Quantity:
            diam = diam.to(u.km).value
        elif type(diam) in [float,int]:
            pass
        else:
            raise TypeError('diam must be a float or an Astropy Unit object')
        self.d_star = diam
    
    def sigma_noise(self,sigma=None):
        '''
        Set the standard deviation of the light-curve (sigma)
        Inputs:
        sigma = float, in flux
        '''
        if sigma is None:
            cut = np.absolute(self.flux-1) < np.std(self.flux)
            for loop in range(10):
                cut = np.absolute(self.flux-1) < 3*np.std(self.flux[cut])
            self.sigma = np.std(self.flux[cut])
        else:
            self.sigma = sigma
        return    
    
    def magnitude_drop(self,star,obj):
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
    
    def occ_model(self,t_ingress,t_egress,opa_ampli,mask,npt_star=12,time_resolution_factor=10):
        """ Returns the modelled light curve considering fresnel difraction, star diameter and intrumental response.
        ----------
        Parameters
        ----------
        t_ingress (int, float): Ingrees time, in seconds.                       (input)
        t_egress  (int, float): Egress time, in seconds.                        (input)
        opa_ampli (int, float): Opacity, opaque = 1.0, transparent = 0.0.       (input)
        mask (array with Booleans): Mask with True values to be computed        (input)
        npt_star  (int): Number of subdivisions for computing the star size's effects, default equal to 12. (auto)
        time_resolution_factor (int,float): Steps for fresnel scale used for modelling the light curve,     (auto)
            default equals to 10 steps for fresnel scale.
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
        time_obs = self.time[mask]
        fresnel_scale_1 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale_2 = calc_fresnel(dist,lamb-dlamb/2.0)
        fresnel_scale   = (fresnel_scale_1 + fresnel_scale_2)/2.0
        time_resolution = (np.min([fresnel_scale/self.vel,self.exptime]))/time_resolution_factor
        #
        #Creating a high resolution curve to compute fresnel difraction, stellar diameter and instrumental integration 
        time_model = np.arange(time_obs.min()-5*self.exptime,time_obs.max()+5*self.exptime,time_resolution)
        #
        #Changing X: time (s) to distances in the sky plane (km), considering the tangential velocity (vel in km/s)
        x   = time_model*self.vel    
        x01 = t_ingress*self.vel
        x02 = t_egress*self.vel
        #
        #Computing fresnel diffraction for the case where the star size is negligenciable
        flux_fresnel_1 = self.__bar_fresnel(x,x01,x02,fresnel_scale_1,opa_ampli)
        flux_fresnel_2 = self.__bar_fresnel(x,x01,x02,fresnel_scale_2,opa_ampli)
        flux_fresnel   = (flux_fresnel_1 + flux_fresnel_2)/2.
        flux_star      = flux_fresnel.copy()
        if (self.d_star > 0):
            #Computing fresnel diffraction for the case where the star size is not negligenciable
            resolucao   = self.d_star/npt_star
            flux_star_1 = np.zeros(len(time_model))
            flux_star_2 = np.zeros(len(time_model))
            #Computing stellar diameter only near the ingrees or egrees times
            star_diam = (np.absolute(x - x01) < 3*self.d_star) + (np.absolute(x - x02) < 3*self.d_star)
            p = np.arange(-npt_star,npt_star)*resolucao
            coeff = np.sqrt(np.absolute(self.d_star**2 - p**2))
            for ii in np.where(star_diam == True)[0]:
                xx = x[ii] + p
                flux1 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_1,opa_ampli)
                flux2 = self.__bar_fresnel(xx,x01,x02,fresnel_scale_2,opa_ampli)
                flux_star_1[ii] = np.sum(coeff*flux1)/coeff.sum()
                flux_star_2[ii] = np.sum(coeff*flux2)/coeff.sum()
                flux_star[ii]   = (flux_star_1[ii] + flux_star_2[ii])/2.
        flux_inst = np.zeros(len(time_obs)) 
        for i in range(len(time_obs)):
            event_model  = (time_model > time_obs[i]-self.exptime/2.) & (time_model < time_obs[i]+self.exptime/2.)
            flux_inst[i] = (flux_star[event_model]).mean()
        self.model[mask] = flux_inst
        self.time_model = time_model
        self.model_star = flux_star
        self.model_fresnel = flux_fresnel
        ev_model = (time_model > t_ingress) & (time_model < t_egress)
        flux_box = np.ones(len(time_model))
        flux_box[ev_model] = (1-opa_ampli)**2
        self.model_geometric = flux_box
        return self.model
    
    def occ_lcfit(self,mask, t_ingress, t_egress, opa_ampli=1, dt_ingress=0, dt_egress=0, dopacity=0, loop=10000):
        """ Brute force chi square fit for occultations lightcurve.
        ----------
        Parameters
        ----------
        mask (array with Booleans): Mask with True values to be computed                (input)
        t_ingress  (int, float): Ingrees time, in seconds.                              (input)
        t_egress   (int, float): Egress time, in seconds.                               (input)
        opa_ampli  (int, float): Opacity, opaque = 1.0, transparent = 0.0, default = 1. (input)
        dt_ingress (int, float): Interval to fit ingress time, default equal to 0, no fit. (auto)
        dt_egress  (int, float): Interval to fit egress time, default equal to 0, no fit.  (auto)
        dopacity   (int, float): Interval to fit opacity, default equal to 0, no fit.      (auto)
        loop       (int): Number of tests to be done, default equal to 10000.              (auto)
        ----------
        Returns
        ----------
        chi2 (array): Tested chi squared values.
        t_i  (array): Tested ingress values
        t_e  (array): Tested egress values
        opa  (array): Tested opacity values
        """
        #
        t_i = t_ingress + dt_ingress*(2*np.random.random(loop) - 1)
        t_e = t_egress  + dt_egress*(2*np.random.random(loop) - 1)
        opa = opa_ampli + dopacity*(2*np.random.random(loop) - 1)
        opa[opa>1.], opa[opa<0.] = 1.0, 0.0
        #
        chi2 = 999999*np.ones(loop)
        tcontrol_f0 = datetime.now()
        for i in range(loop):
           model_test = self.occ_model(t_i[i],t_e[i],opa[i],mask)
           chi2[i] = np.sum((self.flux[mask] -  model_test[mask])**2)/(self.sigma**2)
        tcontrol_f3 = datetime.now()
        print('Elapsed time: {:.3f} seconds.'.format((tcontrol_f3 - tcontrol_f0).total_seconds()))
        chisquare = ChiSquare(chi2, len(self.flux[mask]), immersion=t_i, emersion=t_e, opacity=opa)
        return chisquare

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
    
    def DetectOcc(self, maximum_duration=None, dur_step=None, snr_limit=None, \
                  n_detections=None):
        """
        Detect automatically the occultation event in the light curve
        
        Parameters
        ----------
        maximum_duration: float, optional
            Maximum duration of the occultation event (default is 1/4th of the
            light curve's time span).
        dur_step: float, optionl
            Step size to sweep occultation duration event (default value is 1/2
            of sampling).
        snr_limit: float, optional
            Minimum occultation SNR.
        n_detections: int, optional
            N best detections regardless from SNR. n_detections is superseded
            by snr_limit.
            
        Returns
        -------
        OrderedDict
            An ordered dictionary of :attr:`name`::attr:`value` pairs for each
            Parameter.
            
        Examples
        --------
        >>> params = lc.DetectOcc()
        >>> params
        {'occultation_duration': 0.0004645648878067732,
         'central_time': 2457852.5916293273,
         'imersion_time': 2457852.5913970447,
         'emersion_time': 2457852.59186161,
         'time_err': 5.799811333417892e-07,
         'depth': 0.8663887801707082,
         'depth_err': 0.10972550419008305,
         'baseline': 0.9110181732552853,
         'baseline_err': 0.1904360360568157,
         'snr': 91.21719495827487,
         'occ_mask': array([False, False, False, ..., False, False, False])}
        """

        # duration of the light curve
        time_span = self.time[-1]-self.time[0]
        if maximum_duration and (maximum_duration > time_span):
            raise ValueError('Occultation duration (maximum_duration={0}) ' \
                             'exceeds the time series lenght ({1:0.5f}).' \
                             .format(maximum_duration, time_span))
        if not maximum_duration:
            maximum_duration = time_span*0.25
        if not dur_step:
            dur_step = np.median(self.time[1:-1]-self.time[0:-2])/2      
        duration_grid = np.arange(dur_step, maximum_duration, dur_step)
        # initial occultation mask (all data points)
        mask = np.ones(len(self.time), dtype=bool)

        if snr_limit:
            # minimum SNR accepted in a detection for multiple search
            snr_value = snr_limit+1
            occ0 = self._runBLS(time_span, duration_grid)
            mask *= ~occ0['occ_mask']
            while (snr_value > snr_limit):
                occ1 = self._runBLS(time_span, duration_grid, mask=mask)
                if occ1['snr'] > snr_limit:
                    snr_value = occ1['snr']
                    mask *= ~occ1['occ_mask']                    
                    occ0 = self._summarizeBLS(occ0,occ1)
                else:
                    snr_value = snr_limit
            return occ0
        elif n_detections:
            # search the n best fits
            occ0 = self._runBLS(time_span, duration_grid)
            mask *= ~occ0['occ_mask']    
            for i in range(n_detections-1):
                occ1 = self._runBLS(time_span, duration_grid, mask=mask)
                snr_value = occ1['snr']
                mask *= ~occ1['occ_mask']                    
                occ0 = self._summarizeBLS(occ0,occ1)
            return occ0
        else:
            # search only the first best fit
            return self._runBLS(time_span, duration_grid)
        
        
    def _runBLS(self, per_grid, dur_grid, mask=None):
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
        imersion_time = stats['transit_times'][0] - r.duration[0]/2
        emersion_time = stats['transit_times'][0] + r.duration[0]/2
        time_err = dur_grid[1] - dur_grid[0]
        depth = np.mean(self.flux[~occ_mask])-np.mean(self.flux[occ_mask])
        depth_err = np.std(self.flux[occ_mask])
        baseline = np.mean(self.flux[~occ_mask])
        baseline_err = np.std(self.flux[~occ_mask])
        snr = (depth/baseline_err)*np.sqrt(np.sum(occ_mask))

        return {'occultation_duration' : occultation_duration, 
                     'central_time' : central_time, 
                     'imersion_time' : imersion_time, 
                     'emersion_time' : emersion_time,
                     'time_err' : time_err,
                     'depth' : depth, 'depth_err' : depth_err,
                     'baseline' : baseline, 'baseline_err' : baseline_err,
                     'snr' : snr,'occ_mask' : occ_mask}
    
    def _summarizeBLS(self, dict1, dict2):
        ''' Private function to merge dictionaries returned by BLS and 
            keep values of common keys in list.
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

    def __bar_fresnel(self,X,X01,X02,fresnel_scale,opa_ampli):
         """ Returns the modelled light curve considering fresnel difraction.
         ----------
         Parameters
         ----------
         X   (array): Array with time values converted in km using the event velocity.
         X01 (int, float): Ingrees time converted in km using the event velocity.
         X02 (int, float): Egress time converted in km using the event velocity.
         fresnel_scale (int, float): Fresnel scale.
         opa_ampli     (int, float): Opacity, opaque = 1.0, transparent = 0.0
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
         r_ampli = - (cc+ss)*(opa_ampli/2.)
         i_ampli =   (cc-ss)*(opa_ampli/2.)
         # Determining the flux considering fresnel difraction
         flux_fresnel = (1.0 + r_ampli)**2 + (i_ampli)**2
         return flux_fresnel


