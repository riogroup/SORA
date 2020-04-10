### Useful classes
class colors():
    '''
    Docstring
    '''
    __cores = {
    'positive': 'blue',
    'negative': 'green',
    'error': 'red'
    }
    def __init__(self):
        pass

    @property
    def positive_color(self):
        return self.__cores['positive']

    @positive_color.setter
    def positive_color(self, value):
        self.__cores['positive']=value

    @property
    def negative_color(self):
        return self.__cores['negative']

    @negative_color.setter
    def negative_color(self, value):
        self.__cores['negative']=value

    @property
    def error_bar(self):
        return self.__cores['error']

    @error_bar.setter
    def error_bar(self, value):
        self.__cores['error']=value


### useful functions
def test_attr(attr, typ, name):
    """
    This function tests if the attribute "attr" belongs to the type "typ",
    if not it gives an error message informing the name of the variable
    """
    try:
        return typ(attr)
    except:
        raise TypeError('"{}" keyword must be a {} object'.format(name, typ))


praia_occ_head = (" Notes:\n"
                " C/A: geocentric closest approach, in arcsec\n"
                " P/A: Planet position angle wrt to star at C/A, in deg\n"
                " vel: velocity in plane of sky, in km/sec, positive= prograde, negative= retrograde\n"
                " R: instrumental magnitude in UCAC2 system\n"
                " J, H, K: 2MASS magnitudes (50.0 = not in 2MASS)\n"
                " R*, J*, H*, K* are normalized magnitudes to a common\n"
                " shadown velocity of 20 km/sec by the relationship:\n"
                " Mag* = Mag_actual + 2.5*log10[velocity/20 (km/sec)]\n"
                " Delta: Planet range to Earth, AU\n"
                " long: East longitude of sub-planet point, deg, positive towards East\n"
                " loc. t.= UT + long: local solar time at sub-planet point, hh:mm\n"
                " --------------------------------------------------------------------------------------------------------------------\n"
                " Selection criteria:\n"
                " Maximum geocentrique closest approach considered:   {max_ca:.3f}\n"
                " Day light exclusion (range loc. t.): 12.0hs to 12.0 hs (t1 = t2 -> no exclusions)\n"
                " List below has:           {size}  entries\n"
                " Reference ephemeris: {ephem}\n"
                " Offset applied to ephemeris off_ra(mas) = A * (t-t0) + B \n"
                " Offset applied to ephemeris off_de(mas) = C * (t-t0) + D \n"
                " t0 =    2005.0000000000000       yrs\n"
                " A =    0.0000000000000000       (mas/yr)\n"
                " B =    0.0000000000000000       (mas)\n"
                " C =    0.0000000000000000       (mas/yr)\n"
                " D =    0.0000000000000000       (mas)\n"
                " pm = proper motion applied? (ok, no)\n"
                " ct = uc (UCAC2); 2m (2MASS); fs (field star); g2 (Gaia-DR2)\n"
                " f = multiplicity flag\n"
                " 0 - no multiple entries per star in astrometry\n"
                " 1 - single position from 2 or more uc/2m entries\n"
                " 2 - single position from 1 uc/2m entry only\n"
                " 3 - fs position from entry with more N contributions\n"
                " 4 - fs position from entry with best (x,y) error\n"
                " 5 - fs position from entry with brightest R mag.\n"
                " 6 - fs position from average over all entries\n"
                " (see details in Assafin et al. 2009)\n"
                " E_ra, E_de: error of star position (mas); (9999 = no estimation)\n"
                " pmra, pmde: star proper motion (mas); (9999 = no proper motion)\n"
                " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------\n"
                "  d  m year  h:m:s UT     ra___dec___J2000_candidate     ra_dec_J2000_target_geocen     C/A    P/A     vel  Delta  G*   J*   H*   K*   long  loc. t.  off_ra   off_de pm ct f E_ra E_de pmra pmde\n"
                " ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------"
                 )


ow_occ_head = (" Planete: {name}: Star GAIA-DR2, {ephem}\n"
            " Notes:\n"
            " C/A: geocentric closest approach, in arcsec\n"
            " P/A: Planet's position angle wrt to star at C/A, in deg\n"
            " vel: velocity in plane of sky, in km/sec, positive= prograde, negative= retrograde\n"
            " G*, K*: quite rough R and K magnitudes from 1.2m ESO Swiss telescope,\n"
            "         **corrected to a standard velocity of 20 km/sec**\n"
            " Delta: Planet's range to Earth, AU\n"
            " long: East longitude of sub-planet point, deg, positive towards East\n"
            " loc. t.= UT + long: **rough** local solar time at sub-planet point, hh:mm\n"
            " ---------------------------------------------------------------------------------------------------------------------------------------\n"
            " Selection criteria:\n"
            " Maximum corrected magnitudes R* and K* considered: 100.0\n"
            " Maximum geocentrique closest approach considered: {max_ca:.3f} \n"
            " List below has:          {size}  entries\n"
            " Offset applied to ephemeris (mas):     0.0     0.0\n"
            " {radius}                                                Radius (km)\n"
            " {ow_des}          {ow_des}          {name}             Provisional Designation, Number, Name\n"
            " ----------------------------------------------------------------------------------------------------------------------------------------------------\n"
            " d  m year  h:m:s UT    ra___dec___J2000_candidate    ra___dec___J2000_ephemeris     C/A    P/A     vel   Delta   G*   K*   long  loc.t.  e_ra  e_dec\n"
            " ----------------------------------------------------------------------------------------------------------------------------------------------------\n"
              )
