import numpy as np
import pylab as pl
import astropy.units as u


class LightCurve():
    '''
    Docstring
    Define a Light Curve
    '''
    def __init__(self, time, flux, dflux=None, **kwargs):
        if len(flux) != len(time):
            print('The time and the flux should have the same length')
            return
        self.flux = flux
        self.dflux = dflux
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

