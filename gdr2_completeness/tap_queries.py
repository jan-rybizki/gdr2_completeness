import time
import pyvo
import numpy as np
import os
import glob
from numpy.lib.recfunctions import stack_arrays

def stack_healpix_files(folder):
    '''
    This is for random downloaded format: stack files to one file and delete them afterwards
    '''
    files = glob.glob1(folder,'*')
    array_list = []
    for item in files: 
        array_list.append(np.load(folder + item))
        os.remove(folder + item)
    result = stack_arrays(array_list,usemask = False, autoconvert=True)
    np.save(folder + 'result',result)
    return(result)

def convert_to_starcount_cube(folder):
    '''
    function that creates the matrix that contains the starcounts per hpx and mag bin
    '''

    def combine_healpix_files(folder):
        file_names = glob.glob1(folder,'*')
        array_list = []
        for file in file_names: 
            array_list.append(np.load(folder + file))
        result = stack_arrays(array_list, usemask=False, autoconvert=True)
        print('Length of the stacked download:',len(result),' And the dtype:', result.dtype)
        return(result) 

    def make_result(data,all_bins,all_hpx):
        healpix_level = np.arange(13)
        number_of_healpixes = np.array([4**(item) * 12 for item in healpix_level])
        if len(all_hpx) not in number_of_healpixes:
            hpx_number = np.min(number_of_healpixes[number_of_healpixes>len(all_hpx)])
            result = np.zeros(shape = (hpx_number,len(all_bins)), dtype = int)
        else:
            result = np.zeros(shape = (len(all_hpx),len(all_bins)), dtype = int)
        for i,item in enumerate(all_bins):
            temp = data[np.where(data['bin']==item)]
            bins = temp['hpx'].astype(np.int32)
            result[bins,i] = temp['ct']
        return(result)

    def save_and_clean(result,all_bins,folder):        
        number_files = len(glob.glob1(folder,'*'))
        np.save(folder+'result',result)
        np.save(folder+'all_bins', all_bins)
        for i in np.arange(number_files):
            os.remove(folder + '%d.npy' %(i))
    
    x = combine_healpix_files(folder)
    all_bins = np.unique(x['bin'])
    all_hpx = np.unique(x['hpx'])
    print('healpixels with entries:',len(all_hpx))
    print('the magnitude bins:',all_bins)
    result = make_result(x,all_bins,all_hpx)
    save_and_clean(result,all_bins,folder)


def get_results(QUERIES, folder, verbose = True):
    """
    This function has been taken from a tutorial by Markus Demleitner
       Input:
          QUERIES: an object that contains, names, urls and the query_text of a tap query
          verbose: weather to print out additional information
    
    """
    job_results = []  
    jobs = set()
    for short_name, access_url, query in QUERIES:
        if not os.path.isfile(folder + short_name + '.npy'):
            # in async, you first create a job:
            job = pyvo.dal.TAPService(access_url).submit_job(query.format(**locals()), maxrec=16000000)
            # then start it. This immediately returns.
            job.run()
            # we keep note of the jobs we started -- we'll watch them later.
            jobs.add((short_name, job))
            if verbose:
                print(short_name + ' added to the job list')
        else:
            if verbose:
                print(short_name + ' is already queried')

    try:
        total_jobs = len(list(jobs))
        jobs_old = len(list(jobs))
        if verbose:
            print('total ADQL queries to download: %d' %(total_jobs))
        while jobs:
            jobs_now = len(list(jobs))
            if jobs_now != jobs_old and verbose:
                print(total_jobs - jobs_now,'/',total_jobs,'done')
            jobs_old = len(list(jobs))
            # we do the list(.) so we can remove jobs with impunity
            for short_name, job in list(jobs):
                # async jobs are in phases; they're done (or failed) when
                # they're neither queued nor executing.
                #print(short_name, job.phase)
                if job.phase not in ('QUEUED', 'EXECUTING'):
                    t = np.array(job.fetch_result().table)
                    np.save(folder + short_name,t)
                    jobs.remove((short_name, job))
                    job.delete()
            # wait a bit before doing the next round of polling
            time.sleep(10)
    finally:
        if jobs:
            print('failed:', jobs)
            for short_name, job in list(jobs):
                job.delete()
        time.sleep(1)

def gaia_hpx_factor(healpix_number = 1):
    """
    returns the number by which to divide the source_id in order to get a hpx number of a specific hpx level
    INPUT:
       healpix_number: the healpix level, ranging from 0 to 12, an integer
    OUTPUT:
       the gaia source id factor to get a specific hpx dicretization
    """
    if healpix_number == -1:
        return(6917528997577384321)
    else:
        return(np.power(2,35)*np.power(4,12-healpix_number))

def number_of_healpixels(healpix_number = 1):
    """
    returns the number of pixels for a specific level
    """
    if healpix_number == -1:
        return(1)
    else:
        return(np.power(4,healpix_number)*12)

## need flag for single hpx query
def tap_query_gdr2_hpx_sliced(service = "CDS", hpx_level = 1, folder = 'data/',
                              Select_what = 'COUNT(*)', under_condition = 'AND phot_g_mean_mag < 12', join_text = '', verbose = True, test_1st_hpx_only = True):
    """
    A function that splits GDR2 ADQL queries into hpx chunks using the source_id and saves them in individual files as they arrive from the server into 'folder'.
    
    Input:
       service: which TAP service should be used. For choice there is one of the following:
           CDS - this is the vizier service
           ARI - this is the ARI Gaia service
           ESA - this is the ESA service (did not seem to work with pyvo so far)
           GDR2light - this is the GAVO service for Gaia DR2 light (only the main columns are included, bad if you need quality flags)
           GDR2mock - this is also hosted by GAVO and contains the GDR2 mock data which is a model Galaxy with a gmag limit of 20.7
       hpx_level: into how many chunks the catalogue should be chopped.-1 = 1, 0 = 12, 1 = 48, 2 = 192 and so on
       folder: Where to store the downloaded data
       Select_what: the ADQL syntax for which rows one is interested, eg. 'parallax, ROUND(phot_g_mean_mag,1) as gmag'
       under_condition = the selection criteria, starting with 'AND' because it follows the healpix condition, for example 'AND phot_g_mean_mag<17 AND astrometric_excess_noise < 1'
       join_text: An additional entry after the table selection before the WHERE condition which can be used to specify crossmatches to other tables
       verbose: Since the queries can get timely you might want to get some feedback how far the download has proceeded
       test_1st_hpx_only: In order to not overload the services we should first test our query on a single hpx and in case it works and it is what we want we can run the other hpx by setting this to 'False'.
    """

    hpx_gaia_factor = gaia_hpx_factor(hpx_level)
    hpx_number = number_of_healpixels(hpx_level)

    if verbose:
        print("Healpix level, gaia_factor, #healpixels")
        print(hpx_level,hpx_gaia_factor, hpx_number)

    if test_1st_hpx_only:
        print('Only the first healpix will be queried. If you want to download all set "test_1st_hpx_only" to "False".')
        hpx_number = 1

    SERVICES = {
      'CDS': ('http://tapvizier.u-strasbg.fr/TAPVizieR/tap','"I/345/gaia2"'),
      'ARI': ('http://gaia.ari.uni-heidelberg.de/tap','gaiadr2.gaia_source'),
      'ESA': ('http://gea.esac.esa.int/tap-server/tap','gaiadr2.gaia_source'),
      'GDR2light': ('http://dc.zah.uni-heidelberg.de/tap','gaia.dr2light'),
      'GDR2mock': ('http://dc.zah.uni-heidelberg.de/tap','gdr2mock.main'),
    'GDR3mock': ('http://dc.zah.uni-heidelberg.de/tap','gdr3mock.main'),}
    access_url, table_name = SERVICES[service]

    if not os.path.exists(folder):
        print('folder did not exist before, is created now')
        os.makedirs(folder)
    else:
        print('folder existed and existing files will not be queried again')
        
    Query_text = """ SELECT %s
    FROM
      %s
      %s
    WHERE source_id BETWEEN %d AND %d
    %s"""
        
    list_of_border_source_ids = []
    for i in np.arange(hpx_number+1):
        list_of_border_source_ids.append(i*hpx_gaia_factor)

    QUERIES = []
    for i,item in enumerate(list_of_border_source_ids):
        if i == len(list_of_border_source_ids)-1:
            continue
        QUERIES.append(("%d" %(i),access_url,Query_text %(Select_what,table_name,join_text,item,list_of_border_source_ids[i+1]-1,under_condition)))
      
    counter = 1
    while len([n for n in os.listdir(folder) if os.path.isfile(folder + n)]) < hpx_number and counter < 4:
        if verbose:
            print('number of attempt:', counter)
            print('number of files in folder:',len([name for name in os.listdir(folder) if os.path.isfile(name)]))
            print('number of hpx chunks:', hpx_number)
        get_results(QUERIES,folder,verbose)
        counter += 1
    if len([n for n in os.listdir(folder) if os.path.isfile(folder + n)]) == hpx_number and verbose:
        print('all files have been successfully downloaded to', folder)
    if len([n for n in os.listdir(folder) if os.path.isfile(folder + n)]) < hpx_number and verbose:
        print('after %d attempts only' %(counter),len([n for n in os.listdir(folder) if os.path.isfile(folder + n)]), 'of', hpx_number, 'have been downloaded')
        print('the server might be overloaded at the moment, try at a later time or maybe go to another service')
