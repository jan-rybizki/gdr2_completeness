import numpy as np
import urllib.request
import os

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

def plot_mollweide_log(data):
    import healpy as hp
    from matplotlib.colors import LogNorm
    import matplotlib.pylab as plt
    norm = LogNorm()
    map_mollweide = data        
    total = sum(map_mollweide)
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

