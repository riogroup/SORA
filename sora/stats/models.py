from numpy import pi, arctan2, sin, cos, sqrt



def ellipse(parameters, x_values, y_values):
    '''
    Returns a

    [extended_summary]

    Parameters
    ----------
    parameters : Parameters object
        Parameters that describe the ellipse
    x_values : `np.ndarray`, `list`
        X axis array of values
    y_values : `np.ndarray`, `list`
        Y axis array of values

    Returns
    -------
    [x_model, y_model] : `list`
        [description]
    '''
    p = parameters.valuesdict()
    
    b = p['equatorial_radius'] - p['equatorial_radius']*p['oblateness']
    phi = p['position_angle'] * np.pi/180.0
    angle = theta + phi
    radial_model = ( p['equatorial_radius'] * b )/np.sqrt( ( p['equatorial_radius'] * np.sin( angle ) )**2 + ( b * np.cos( angle ) )**2 )
    x_model = p['center_f'] + radial_model*np.cos( theta )
    y_model = p['center_g'] + radial_model*np.sin( theta )
    
    return [x_model, y_model]


def ellipseError(parameters, f, g, uncertainty, ellipse_error=0):
    '''
    ellipseError [summary]S

    [extended_summary]

    Parameters
    ----------
    parameters : [type]
        [description]
    f : [type]
        [description]
    g : [type]
        [description]
    uncertainty : [type]
        [description]
    ellipse_error : int, optional
        [description], by default 0

    Returns
    -------
    [type]
        [description]
    '''
    f_model, g_model = ellipse(parameters, f, g)
    return (( (f - f_model)**2 + (g - g_model)**2 )/(uncertainty**2 + ellipse_error**2) )