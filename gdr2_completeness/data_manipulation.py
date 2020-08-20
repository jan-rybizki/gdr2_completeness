import numpy as np
import urllib.request
import os
from healpy import ang2pix, pix2ang
from astropy.coordinates import SkyCoord

def hpx2lb(hpx):
    '''
    Function to convert between HEALpix level 6 to Galactic l,b values
    '''
    ra,dec = pix2ang(64,hpx,nest = True, lonlat=True)
    coods = SkyCoord(ra,dec,unit = 'deg', frame = 'icrs')
    l,b = coods.galactic.l.value, coods.galactic.b.value
    return(l,b)

def g_grp4specific_grvs_color(grvs,col):
    """
    returns the g and grp (input for selection function) for specific grvs and g-grp values
    Only takes scalar magnitudes!
    """
    g = np.linspace(-2,18,2001)
    grp = g - col
    grvs_array = grvs_from_g_grp_array(g,grp)
    cut = np.where(np.abs(grvs_array - grvs) == np.min(np.abs(grvs_array - grvs)))
    return(g[cut][0],grp[cut][0])

def grvs_from_g_grp(g,rp):
    """
    returns grvs from g and grp
    """
    if g-rp < 1.4:
        grvs = grvs_for_blue_sources(g,rp)
    else:
        grvs = grvs_for_red_sources(g,rp)
    return(grvs)    

def grvs_from_g_grp_array(g,rp):
    """
    returns grvs from g and grp
    """
    grvs = np.zeros_like(g)    
    blue = (g-rp < 1.4)
    red = (g-rp >= 1.4)
    grvs[blue] = grvs_for_blue_sources(g[blue],rp[blue])
    grvs[red] = grvs_for_red_sources(g[red],rp[red])
    return(grvs)    


def grvs_for_blue_sources(g,rp):
    """
    This comes from the GDR2 release paper if g-grp < 1.4
    """
    grp = g-rp
    a0 = 0.042319
    a1 = -0.65124
    a2 = 1.0215
    a3 = -1.3947
    a4 = 0.53768
    grvs = rp + a0 + a1 * grp + a2 * grp**2 + a3 * grp**3 + a4 * grp**4  
    return(grvs)

def grvs_for_red_sources(g,rp):
    '''
    for g-grp >= 1.4
    '''
    grp = g-rp
    a0 = 132.32
    a1 = -377.28
    a2 = 402.32
    a3 = -190.97
    a4 = 34.026
    grvs = rp + a0 + a1 * grp + a2 * grp**2 + a3 * grp**3 + a4 * grp**4  
    return(grvs)


class rvs_selection_function(object):
    '''
    This objects helps to query the Gaia DR2 radial velocity selection function
    '''
    def __init__(self, from_rp_to_g_completeness=True):
        '''
        Initialisation. Loads the data from the compressed npy files.
        if from_rp_to_g_completeness is True, then the array containing
        the internal GDR2 RP completeness with respect to G sources is loaded
        and will be used for completeness calculation.
        '''
        x = np.load("../data/3d_cube_level_6_completeness_function.npz")
        self.from_rp_to_g_completeness = from_rp_to_g_completeness
        self.rvs = x['arr_0']
        self.dr2 =x['arr_1']
        self.xrange =x['arr_2']
        self.yrange =x['arr_3']
        self.lower =x['arr_4']
        self.upper =x['arr_5']
        self.compl =x['arr_6']
        self.hplvl =x['arr_7']
        self.nall =x['arr_8']
        self.nrvs =x['arr_9']

        if from_rp_to_g_completeness:
            y = np.load("../data/g_rp_internal_completeness.npz")
            self.data_rp = y["arr_0"]
            self.data_all = y["arr_1"]
            self.bins_all = y["arr_2"]
            
    def query(self, g, rp, l, b, extra = 'none'):
        '''
        This queries the completeness function for single values
        INPUT:
           g = G magnitude
           rp = RP magnitude
           l = Galactic Longitude
           b = Galactic Latitude
        OUTPUT:
           completeness = function between 0 and 1
        '''
        # converting g and grp magnitude to grvs
        col = g-rp
        if col < 1.4:
            grvs = grvs_for_blue_sources(g,rp)
        else:
            grvs = grvs_for_red_sources(g,rp)
        
        # calculating the respective healpix from l and b
        coods = SkyCoord(l, b, unit='deg', frame='galactic')
        ra, dec = coods.icrs.ra.value, coods.icrs.dec.value 
        hpx = ang2pix(64,ra,dec,nest = True,lonlat = True)

        # finding the right mag and color bin
        icol = np.where(np.abs(self.xrange-col)==np.min(np.abs(self.xrange-col)))[0][0]
        imag = np.where(np.abs(self.yrange-grvs)==np.min(np.abs(self.yrange-grvs)))[0][0]
        
        # Warnings for col and mag out of range
        if col < 0.05:
            print('Warning: G-G_RP is too blue, limit is 0.05 mag, your source has %.2f mag' %(col))
        elif col > 1.75:
            print('Warning: G-G_RP is too red, limit is 1.75 mag, your source has %.2f mag' %(col))
        if grvs < 2.9:
            print('Warning: G_RVS is too bright, limit is 2.9 mag, your source has %.2f mag' %(grvs))
        elif grvs > 14.1:
            print('Warning: G_RVS is too faint, limit is 14.1 mag, your source has %.2f mag' %(grvs))
        
        # querying the array
        completeness = self.compl[hpx,imag,icol]
        
        # adding rp completeness if from_rp_to_g_completeness is True
        if self.from_rp_to_g_completeness:
            imagg = np.where(np.abs(self.bins_all-g)==np.min(np.abs(self.bins_all-g)))[0][0]
            rp_completeness = np.divide(self.data_rp[hpx,imagg],self.data_all[hpx,imagg])
            if rp_completeness == np.nan:
                rp_completeness = 1.0            
            completeness *= rp_completeness
            
            # Warning for G out of range
            if g < 1.7:
                print('Warning: G is too bright, limit is 1.7 mag, your source has %.2f mag' %(g))
            elif g > 15.1:
                print('Warning: G is too faint, limit is 15.1 mag, your source has %.2f mag' %(g))
        if extra == 'none':
            return(completeness)
        elif extra == 'rvs_count':
            return(completeness,self.rvs[hpx,imag,icol])
        elif extra == 'all_count':
            return(completeness,self.dr2[hpx,imag,icol])
        elif extra == 'col_bins':
            return(completeness,self.xrange)
        elif extra == 'rvsmag_bins':
            return(completeness,self.yrange)
        elif extra == 'gmag_bins':
            return(completeness,self.bins_all)
        elif extra == 'hpx_level':
            return(completeness,self.hplvl[hpx,imag,icol])
        elif extra == 'rvs_count_averaged':
            return(completeness,self.nrvs[hpx,imag,icol])
        elif extra == 'all_count_averaged':
            return(completeness,self.nall[hpx,imag,icol])
        elif extra == 'rp_count_grp_internal_completeness':
            return(completeness,self.data_rp[hpx,imagg])
        elif extra == 'all_count_grp_internal_completeness':
            return(completeness,self.data_all[hpx,imagg])
        
        elif extra == 'completeness_upper':
            if self.from_rp_to_g_completeness:
                upper = self.upper[hpx,imag,icol]*rp_completeness
                return(completeness,upper)
            else:
                return(completeness,self.upper[hpx,imag,icol])
        elif extra == 'completeness_lower':
            if self.from_rp_to_g_completeness:
                lower = self.lower[hpx,imag,icol] * rp_completeness
                return(completeness,lower)
            else:
                return(completeness,self.lower[hpx,imag,icol])

        
    def query_array(self, g, rp, l, b, extra = 'none'):
        '''
        This queries the completeness function in array form
        INPUT:
           g = G magnitude
           rp = RP magnitude
           l = Galactic Longitude
           b = Galactic Latitude
        OUTPUT:
           completeness = function between 0 and 1
        '''
        # converting g and grp magnitude to grvs
        grvs = np.zeros_like(g)
        col = g-rp
        blue_cut = (col<1.4) 
        grvs[blue_cut] = grvs_for_blue_sources(g[blue_cut],rp[blue_cut])
        red_cut = (col>=1.4) 
        grvs[red_cut] = grvs_for_red_sources(g[red_cut],rp[red_cut])        
        
        
        # calculating the respective healpix from l and b
        coods = SkyCoord(l, b, unit='deg', frame='galactic')
        ra, dec = coods.icrs.ra.value, coods.icrs.dec.value 
        hpx = ang2pix(64,ra,dec,nest = True,lonlat = True)

        # finding the right mag and color bin
        icol = np.zeros(len(col),dtype=np.int32)
        for i in range(len(col)):
            icol[i] = np.where(np.abs(self.xrange-col[i])==np.min(np.abs(self.xrange-col[i])))[0][0]
        imag = np.zeros_like(icol)
        for i in range(len(col)):
            imag[i] = np.where(np.abs(self.yrange-grvs[i])==np.min(np.abs(self.yrange-grvs[i])))[0][0]
        
        # Warnings for col and mag out of range
        if np.any(col < 0.05):
            print('Warning: G-G_RP is too blue, limit is 0.05 mag, your source has %.2f mag' %(np.min(col)))
        elif np.any(col > 1.75):
            print('Warning: G-G_RP is too red, limit is 1.75 mag, your source has %.2f mag' %(np.max(col)))
        if np.any(grvs < 2.9):
            print('Warning: G_RVS is too bright, limit is 2.9 mag, your source has %.2f mag' %(np.min(grvs)))
        elif np.any(grvs > 14.1):
            print('Warning: G_RVS is too faint, limit is 14.1 mag, your source has %.2f mag' %(np.max(grvs)))
        
        # querying the array
        completeness = np.zeros(len(g))
        for i in range(len(g)):
            completeness[i] = self.compl[hpx[i],imag[i],icol[i]]
        
        # adding rp completeness if from_rp_to_g_completeness is True
        if self.from_rp_to_g_completeness:
            imagg = np.zeros_like(imag)
            for i in range(len(g)):
                imagg[i] = np.where(np.abs(self.bins_all-g[i])==np.min(np.abs(self.bins_all-g[i])))[0][0]
            rp_completeness = np.zeros(len(g))
            for i in range(len(g)):
                rp_completeness[i] = np.divide(self.data_rp[hpx[i],imagg[i]],self.data_all[hpx[i],imagg[i]])
            nan = np.isnan(rp_completeness)
            rp_completeness[nan] = 1.0                        
            completeness *= rp_completeness
            
            # Warning for G out of range
            if np.any(g < 1.7):
                print('Warning: G is too bright, limit is 1.7 mag, your source has %.2f mag' %(np.min(g)))
            elif np.any(g > 15.1):
                print('Warning: G is too faint, limit is 15.1 mag, your source has %.2f mag' %(np.max(g)))
        if extra == 'none':
            return(completeness)

        elif extra == 'rvs_count':
            rvs_count = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                rvs_count[i] = self.rvs[hpx[i],imag[i],icol[i]]
            return(completeness,rvs_count)
 
        elif extra == 'all_count':
            all_count = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                all_count[i] = self.dr2[hpx[i],imag[i],icol[i]]
            return(completeness,all_count)
        
        elif extra == 'col_bins':
            return(completeness,self.xrange)
        
        elif extra == 'rvsmag_bins':
            return(completeness,self.yrange)
        
        elif extra == 'gmag_bins':
            return(completeness,self.bins_all)
        
        elif extra == 'hpx_level':
            hplvl = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                hplvl[i] = self.hplvl[hpx[i],imag[i],icol[i]]
            return(completeness,hplvl)
        
        elif extra == 'rvs_count_averaged':
            rvs_count_a = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                rvs_count_a[i] = self.nrvs[hpx[i],imag[i],icol[i]]
            return(completeness,rvs_count_a)
        
        elif extra == 'all_count_averaged':
            all_count_a = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                all_count_a[i] = self.nall[hpx[i],imag[i],icol[i]]
            return(completeness,all_count_a)
        
        elif extra == 'rp_count_grp_internal_completeness':
            rvs_count_grp = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                rvs_count_grp[i] = self.data_rp[hpx[i],imagg[i]]
            return(completeness,rvs_count_grp)
        
        elif extra == 'all_count_grp_internal_completeness':
            all_count_grp = np.zeros(len(g),dtype = np.int32)
            for i in range(len(g)):
                all_count_grp[i] = self.data_all[hpx[i],imagg[i]]
            return(completeness,all_count_grp)
        
        elif extra == 'completeness_upper':
            upper = np.zeros(len(g))
            for i in range(len(g)):
                upper[i] = self.upper[hpx[i],imag[i],icol[i]]
            if self.from_rp_to_g_completeness:
                upper *= rp_completeness
                return(completeness,upper)
            else:
                return(completeness,upper)

        elif extra == 'completeness_lower':
            lower = np.zeros(len(g))
            for i in range(len(g)):
                lower[i] = self.lower[hpx[i],imag[i],icol[i]]
            if self.from_rp_to_g_completeness:
                lower *= rp_completeness
                return(completeness,lower)
            else:
                return(completeness,lower)



def download_completeness_maps(folder, also_dr2_all_starcounts = False):
    if not os.path.exists(folder):
        print('folder did not exist before, is created now')
        os.makedirs(folder)
    urllib.request.urlretrieve("https://keeper.mpdl.mpg.de/f/4cace26d760646faa6f0/?dl=1", folder + 'dr2_bins.npy')
    urllib.request.urlretrieve("https://keeper.mpdl.mpg.de/f/1c4451ef5faa4aad8c78/?dl=1", folder + '15_18.npy')
    urllib.request.urlretrieve("https://keeper.mpdl.mpg.de/f/c724ef1785e644efb80b/?dl=1", folder + '12_15.npy')
    urllib.request.urlretrieve("https://keeper.mpdl.mpg.de/f/96938a21f932487fb0df/?dl=1", folder + '8_12.npy')
    if also_dr2_all_starcounts:
        urllib.request.urlretrieve("https://keeper.mpdl.mpg.de/f/76165998e7e7470193a2/?dl=1", folder + 'dr2_starcounts.npy')

def create_completeness_matrix(folder):
    '''
    uses the files from Ronald Drimmel (DR2 completeness from xmatch with 2MASS) 
    to create a file that has the same structure as the completeness counts
    '''
    bright = np.load(folder + '8_12.npy')
    mid = np.load(folder + '12_15.npy')
    faint = np.load(folder + '15_18.npy')
    bins = np.load(folder + 'dr2_bins.npy')
    cdr2 = np.zeros(shape=(len(bright),len(bins)))
    for item in bins:
        if item >= 8 and item < 12:
            cdr2[:,np.where(bins==item)[0][0]] = np.array(bright)
        elif item >= 12 and item <= 15:
            cdr2[:,np.where(bins==item)[0][0]] = np.array(mid)
        elif item > 15 and item <= 18:
            cdr2[:,np.where(bins==item)[0][0]] = np.array(faint)
    return(cdr2,bins)

def decrease_hpx(data,ntimes):
    for i in range(ntimes):
        data = data[0::4,:] + data[1::4,:] + data[2::4,:] + data[3::4,:]
    return(data)
def total_mag_function(data):
    return(np.sum(data,axis=0))

def density_sky(data):
    return(np.sum(data,axis=1))
def last_complete_bin(data,bins):
    max_values = np.max(data[:,:], axis = 1)
    mag_limit = np.zeros(len(data))
    for i,item in enumerate(data[:,:]):
        #print(i,item)
        cut = np.where(item == max_values[i])
        mag_limit[i] = bins[cut[0][0]]
    return(mag_limit)
def get_mag_bin(data,bins,mag):
    assert(mag in bins)
    cut = np.where(bins==mag)[0]
    #print(cut)
    return(data[:,cut][:,0])

def percentile_maglim(data,bins,p):
    '''
    this gives the magnitude at which p from (0,1) percent of the sources has been seen
    '''
    nHpx = len(data)
    maglim = np.zeros(nHpx)
    for i in range(nHpx):
        #if i%10000 == 0:
        #    print(i,nHpx)
        cut = np.where(np.abs(np.divide(np.cumsum(data[i]),sum(data[i])) - p) == np.min(np.abs(np.divide(np.cumsum(data[i]),sum(data[i])) - p)))[0][0]
        maglim[i] = bins[cut]
    return(maglim)

def plot_mollweide_log(data):
    import healpy as hp
    from matplotlib.colors import LogNorm
    import matplotlib.pylab as plt
    norm = LogNorm()
    map_mollweide = data        
    total = np.nansum(map_mollweide)
    hp.mollview(map_mollweide, cbar = True, min=None, max=None, nest = True,norm = norm, coord= "CG", unit = 'starcount per hpx',notext =True)
    plt.title("total starcount=%d" %(total))        
    plt.show()
    
def plot_mollweide_linear(data):
    import healpy as hp
    import matplotlib.pylab as plt

    map_mollweide = data
    hp.mollview(map_mollweide, cbar = True, nest = True, coord= "CG", unit = 'completeness',notext =True, max = 1)
    plt.show()    
def plot_mollweide_linear_mag(data,mag):
    import healpy as hp
    import matplotlib.pylab as plt

    map_mollweide = data
    hp.mollview(map_mollweide, cbar = True, nest = True, coord= "CG", unit = 'completeness',notext =True,min = 0, max = 1)
    plt.title(r"Completeness RVS sample with $\varpi>0$ and $\varpi/\sigma_\varpi>20$  at G=%.1f" %(mag))
    plt.savefig('plots/per_mag_bin/rvs_tutorial/rvs_completeness_%.1f.png' %(mag))
    plt.clf()
    plt.close()

