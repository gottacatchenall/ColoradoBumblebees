from chelsa_cmip6.GetClim import chelsa_cmip6, CmipClimat, ChelsaClimat, BioClim, DeltaChangeClim

import os
import pathlib
import xarray as xr
import rioxarray as rio 

print('Warning: this script relies on having an internet connection and tends to take a while...')

class MyChelsaClimat:
    def __init__(self, dicts):
        for var in ['pr', 'tas', 'tasmax', 'tasmin']:
            print("getting variable: " + var)
            setattr(self, var, dicts[var])

def get_chelsa():
    """ 
    Calculate chelsa cmip 6 climatological normals and bioclimatic variables
    
    :param activity_id: the activity_id according to CMIP6
    :param table_id: the table id according to CMIP6
    :param experiment_id: the experiment_id according to CMIP6
    :param instituion_id: the instituion_id according to CMIP6
    :param source_id: the source_id according to CMIP6
    :param member_id: the member_id according to CMIP6
    :param output: output directory, string
    """
    activity_id='ScenarioMIP'
    table_id='Amon'
    institution_id='NOAA-GFDL'
    source_id='GFDL-ESM4'
    member_id='r1i1p1f1'
    SSPS = ["ssp126", "ssp245", "ssp370"]

    REFERENCE_START = '2000-01-15'
    REFERENCE_END = '2014-12-15'
    #START_YEARS = ["20%i0-01-15" % i for i in range(2,10)]
    #END_YEARS = ["20%i9-12-15" % i for i in range(2,10)]
    START_YEARS = ["20%i0-01-15" % i for i in range(7,10)]
    END_YEARS = ["20%i9-12-15" % i for i in range(7,10)]


    EXTENT = {
        'bottom': 34.0,
        'top': 44.0,
        'left': -110.5,
        'right': -103.5
    }
    xmin, xmax, ymin, ymax = EXTENT["left"], EXTENT["right"], EXTENT["bottom"], EXTENT["top"]

    # Needs to be uncommented to run the first time 
    #ch_climat = ChelsaClimat(xmin, xmax, ymin, ymax)
    #for var in ['pr', 'tas', 'tasmax', 'tasmin']:
    #    getattr(ch_climat, var).to_netcdf("CHELSA_%s_baseline.nc" % var)
    
    dicts = {var: xr.open_dataset("CHELSA_%s_baseline.nc" % var) for var in ['pr', 'tas', 'tasmax', 'tasmin']}
    my_ch_climat = MyChelsaClimat(dicts)

    baseline_dir = "CHELSA/baseline/"
    pathlib.Path(baseline_dir).mkdir(parents=True, exist_ok=True)


    pathlib.Path("CHELSA").mkdir(parents=True, exist_ok=True)
    for ssp in SSPS:
        pathlib.Path("CHELSA/%s" % ssp).mkdir(parents=True, exist_ok=True)
        for x in range(len(START_YEARS)):
            fefps, fefpe = START_YEARS[x], END_YEARS[x]

            startyear, endyear = fefps.split("-")[0], fefpe.split("-")[0]

            output_dir = "CHELSA/%s/%s-%s/" % (ssp, startyear, endyear)
            pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

            print('start downloading CMIP data:')
            cm_climat = CmipClimat(activity_id, table_id,
                                    ssp,
                                    institution_id, source_id,
                                    member_id, REFERENCE_START,
                                    REFERENCE_END, fefps,
                                    fefpe)

            print('applying delta change:')
            dc = DeltaChangeClim(my_ch_climat, cm_climat, REFERENCE_START,
                                    REFERENCE_END, fefps,
                                    fefpe)

            print('start building climatologies data:')
            biohist = BioClim(dc.hist_pr, dc.hist_tas, dc.hist_tasmax, dc.hist_tasmin)
            biofutr = BioClim(dc.futr_pr, dc.futr_tas, dc.futr_tasmax, dc.futr_tasmin)

            print('saving bioclims:')
                
            for n in range(1, 20):
                name = output_dir + str('bio' + str(n)) + '.csv'
                a = getattr(biofutr, 'bio' + str(n))()
                a['bio'+str(n)].to_pandas().to_csv(name)

get_chelsa()


