import astropy.units as u
from pyia import GaiaData
import matplotlib.pyplot as plt
import pandas as pd
from astropy.coordinates import SkyCoord
import matplotlib.cm as cm
import numpy as np
import os
import random
import time

WFS = {'nxSubaps': 12,
       'pxlsPerSubap': 10,
       'subapFOV': 10,
       }

light = {
    'wavelength': 550e-9, # Wavelength of light to test
    'int_frequency_m': 440e-9, # Bandpass for Gaia G band in m
    'int_frequency_nm': 440 # Bandpass for Gaia G band in nm
}

PLANCK = 6.626e-34
C = 3e8

class cameras:
    def __init__(self,name,n_pixels,qe,bits,fw,dc,rn,sens):
        self.name = name
        self.n_pixels = n_pixels
        self.qe = qe
        self.bits = bits
        self.fw = fw
        self.dc = dc
        self.rn = rn
        self.sens = sens
    
    def saturation(self):
        return (self.n_pixels)*((2**self.bits)*self.sens)/self.qe
    
class telescopes:
    def __init__(self, name, diameter, fov):
        self.diameter = diameter
        self.fov = fov # In degrees

    def area(self):
        return (np.pi*(self.diameter/2)**2)*10000 # Converting to cm2 for CGS units
    
class shwf_sesnor:
    def __init__(self,n_subaps, subap_fov):
        self.n_subaps = n_subaps
        self.subap_fov = subap_fov

    def pxl_per_subap(self,camera_pixels):
        return camera_pixels/self.n_subaps
        
    
ASI6200 = cameras('ASI6200',4096*2160, 0.9,16,1504500,0.02,5, 0.8)
Marana42B6 = cameras('Marana42B6',128*128, 0.95,16,55000,0.15,1.6, 0.8)
telescope_1_2 = telescopes("1.2 M, FOV: 0.4",1.2,0.4)
default_shwf = shwf_sesnor(64, 10)

def fetch_data(min_mag = 7, max_mag = 12, save_data = True, filepath = ""):

    # Define query to fetch the desired star data
    query = 'SELECT l, b, phot_g_mean_mag FROM gaiadr3.gaia_source WHERE phot_g_mean_mag BETWEEN {} AND {}'.format(min_mag,max_mag)

    # Specify the filename to save the data
    filename = filepath + 'star_data_{}_{}.csv'.format(min_mag,max_mag)

    if save_data == True:
        # Check if the data file already exists
        try:
            # Load the data from the file using Pandas
            df = pd.read_csv(filename)
            print('Data fetched from local file.')
        except FileNotFoundError:
            # Fetch the data using from_query() and save it to the file
            print('No local data found, fetching from Gaia catologue')
            g = GaiaData.from_query(query)

            # Create a DataFrame from the extracted columns
            data = {'l': g.l, 'b': g.b,'mag': g.phot_g_mean_mag}
            df = pd.DataFrame(data)
            
            # Save the DataFrame to a CSV file
            df.to_csv(filename, index=False)
            print('Data fetched and saved to file.')
    else:
            print('No local data found, fetching from Gaia catologue')
            g = GaiaData.from_query(query)

            # Create a DataFrame from the extracted columns
            data = {'l': g.l, 'b': g.b,'mag': g.phot_g_mean_mag}
            df = pd.DataFrame(data)
                   

    return df


def create_folder(folder_path):
    if not os.path.exists(folder_path):
        os.makedirs(folder_path)
        print(f"Folder '{folder_path}' created successfully.")
    else:
        print(f"Folder '{folder_path}' already exists.")

def convert_to_aitoff(star_data):
    
    # Convert the coordinates to the Aitoff projection
    coords = SkyCoord(star_data['l'], star_data['b'], unit = (u.deg, u.deg) ,frame='galactic')
    l_rad = coords.l.wrap_at('180d').radian
    b_rad = coords.b.radian

    return l_rad, b_rad


def plot_aitoff(star_data):
    l_rad, b_rad = convert_to_aitoff(star_data)
    mag = star_data['mag']
    # Plot the star map in Aitoff projection
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='aitoff')
    ax.grid(True)

    # Adjust the marker size and color as desired
    ax.scatter(l_rad, b_rad, s=0.1, c=mag, cmap = 'inferno')

    # Create a ScalarMappable object for colorbar
    sm = cm.ScalarMappable(cmap='inferno')
    sm.set_array(mag)

    # Add colorbar
    ax.cbar = plt.colorbar(sm, shrink = 0.7)
    ax.cbar.set_label('Magnitude')

    plt.title('Star Map of Stars with Magnitudes between {} and {}'.format(np.floor(star_data['mag'].min()),np.floor(star_data['mag'].max())))
    plt.show()

def plot_inset_graph(star_data, telescope, poi_l = 0, poi_b = 0):
    fig, ax = plt.subplots()
    coords = SkyCoord(star_data['l'], star_data['b'], unit = (u.deg, u.deg) ,frame='galactic')
    l_deg_wrapped = coords.l.wrap_at('180d')

    ax.scatter(l_deg_wrapped, star_data['b'], s=0.1, c=star_data['mag'], cmap = 'viridis')

    # Create a ScalarMappable object for colorbar
    sm = cm.ScalarMappable(cmap='viridis')
    sm.set_array(star_data['mag'])

    # Add colorbar
    cbar = plt.colorbar(sm, ax=ax)
    cbar.set_label('Magnitude')

    plt.xlabel('l (degrees)')
    plt.ylabel('b (degrees)')
    plt.title('Postional Map of Stars with Magnitudes between {} and {}'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max())))
    plt.grid(True)

    # Add inset, define position and size
    axin = ax.inset_axes([1.25, 0.7, 0.3, 0.25])

    # Plot inset points
    axin.scatter(l_deg_wrapped, star_data['b'], s = 6, c = star_data['mag'], cmap = 'viridis')

    # Inset limits
    x_lims = (poi_l,poi_l + telescope.fov)
    y_lims = (poi_b, poi_b + telescope.fov)
    axin.set_xlim(x_lims[0], x_lims[1])
    axin.set_ylim(y_lims[0], y_lims[1])
    axin.set_xlabel('l (degrees)')
    axin.set_ylabel('b (degrees)')
    ax.indicate_inset_zoom(axin)
    plt.tight_layout()

    plt.show()

def plot_mag_hist(star_data):
    rounded_mag = []

    for _mag in star_data['mag']:
        rounded_mag.append(np.floor(_mag))

    fig = plt.figure()
    labels, counts = np.unique(rounded_mag, return_counts=True)
    plt.bar(labels, counts, align='center')
    plt.xlabel('Star Magnitude')
    plt.ylabel('Counts')
    plt.title('Number of stars by magnitude in data set')
    plt.gca().set_xticks(labels)
    plt.show()


def e_ph(wavelength = light['wavelength']):
    return PLANCK * (C/wavelength) * 1e7 # Converting to erg for CGS units


def stellar_flux(star_data):
    f0 = 3229e-23 # Vega zero-point flux in the G band with units of erg⋅s−1⋅cm−2⋅Hz−1
    star_flux_mag = {
        'mag':[],
        'flux':[]
    }
    
    for mags in star_data['mag']:
        star_flux_mag['mag'].append(mags)
        star_flux_mag['flux'].append((10**(mags/-2.5) * f0*3e8)/(light['int_frequency_nm']*light['int_frequency_m'])) 
 
    
    return star_flux_mag

def calc_mag(exp_t, telescope, photon_flux):
    f0 =  ((3229e-23) * u.erg / u.cm**2 / u.s / u.Hz).to(u.erg / u.cm**2 / u.s / u.nm,
                equivalencies=u.spectral_density(5500 * u.AA)).value # Vega zero-point flux in the G band with units of erg⋅s−1⋅cm−2⋅nm-1
    telescope_area = telescope.area()
    energy_flux = (photon_flux*e_ph())/(exp_t*telescope_area*(light['int_frequency_nm']))
    return -2.5*np.log10(energy_flux/f0)

def n_ph(star_data, exp_t, telescope):
    star_params = stellar_flux(star_data)
    star_params['n_ph'] = []
    photon_energy = e_ph()
    telescope_area = telescope.area()

    for fluxes in star_params['flux']:
        star_params['n_ph'].append((fluxes*exp_t*telescope_area*(light['int_frequency_nm']))/(photon_energy))
    
    return star_params

def calc_sky_background(exp_t, telescope, bg_mag_equivalent = 12):
    f0 = 3229e-23
    flux = (10**(bg_mag_equivalent/-2.5) * f0*3e8)/(light['int_frequency_nm']*light['int_frequency_m'])
    return ((flux)*exp_t*telescope.area()*(light['int_frequency_nm']))/(e_ph())

def calc_snr_standard(star_data, exp_t, telescope, camera):

    star_params = n_ph(star_data, exp_t, telescope)
    star_params['ideal_snr'] = []
    star_params['actual_snr'] = []
    star_params['n_ph_per_pix'] = []
    
    for photons in star_params['n_ph']:
        dark_noise = camera.rn**2 + (camera.dc*exp_t)
        star_params['ideal_snr'].append(np.sqrt((photons)))
        star_params['actual_snr'].append((photons*camera.qe)/np.sqrt(camera.n_pixels*(dark_noise) + (calc_sky_background(exp_t,telescope)*camera.qe) + (photons*camera.qe)))
        star_params['n_ph_per_pix'].append(photons/camera.n_pixels)

    return star_params

def calc_standard_noise_contribution(star_data, exp_t, telescope, camera):
        star_params = n_ph(star_data, exp_t, telescope)

        noise_contribution_standard = {
             'Dark Current': [],
             'Read Out Noise': [],
             'Shot Noise':[],
             'Background Noise':[]
        }
        
        for photons in star_params['n_ph']:
             sky_background = (calc_sky_background(exp_t,telescope))
             total_noise = (photons*camera.qe) + (camera.n_pixels*camera.rn**2) + (camera.n_pixels*(camera.dc*exp_t)) + (sky_background* camera.qe)
             noise_contribution_standard['Shot Noise'].append((photons*camera.qe)/total_noise)
             noise_contribution_standard['Read Out Noise'].append((camera.n_pixels*camera.rn**2)/total_noise)
             noise_contribution_standard['Dark Current'].append((camera.n_pixels*camera.dc*exp_t)/total_noise)
             noise_contribution_standard['Background Noise'].append((sky_background* camera.qe)/total_noise)

        return noise_contribution_standard


def min_photon_standard(exp_t,telescope, camera, snr = 5):
    return (camera.qe*snr**2)+((camera.qe*snr)*(np.sqrt((4*calc_sky_background(exp_t,telescope)*camera.qe)+(4*camera.n_pixels*camera.dc*exp_t)+(4*camera.n_pixels*camera.rn**2)+(snr**2))))/(2*camera.qe**2)

def calc_shwf_snr_standard(star_data, exp_t, telescope, camera, shwf):

    star_params = n_ph(star_data, exp_t, telescope)
    star_params['ideal_snr'] = []
    star_params['actual_snr'] = []
    star_params['n_ph_per_pix'] = []
    pix_per_subap = shwf.pxl_per_subap(camera.n_pixels)
    
    for photons in star_params['n_ph']:
        photons = photons/shwf.n_subaps
        dark_noise = camera.rn**2 + (camera.dc*exp_t)
        star_params['ideal_snr'].append(np.sqrt(photons))
        star_params['actual_snr'].append((photons*camera.qe)/np.sqrt(pix_per_subap*(dark_noise) + ((calc_sky_background(exp_t,telescope)/shwf.n_subaps)*camera.qe) + (photons*camera.qe)))
        star_params['n_ph_per_pix'].append(photons/pix_per_subap)

    return star_params

def calc_shwf_noise_contribution(star_data, exp_t, telescope, camera, shwf):
        star_params = n_ph(star_data, exp_t, telescope)
        pix_per_subap = shwf.pxl_per_subap(camera.n_pixels)

        noise_contribution_standard = {
             'Dark Current': [],
             'Read Out Noise': [],
             'Shot Noise':[],
             'Background Noise':[]
        }
        
        for photons in star_params['n_ph']:
             photons = photons/shwf.n_subaps
             sky_background = calc_sky_background(exp_t,telescope)/shwf.n_subaps
             total_noise = (photons*camera.qe) + (pix_per_subap*camera.rn**2) + (pix_per_subap*(camera.dc*exp_t)) + (sky_background* camera.qe)
             noise_contribution_standard['Shot Noise'].append((photons*camera.qe)/total_noise)
             noise_contribution_standard['Read Out Noise'].append((pix_per_subap*camera.rn**2)/total_noise)
             noise_contribution_standard['Dark Current'].append((pix_per_subap*camera.dc*exp_t)/total_noise)
             noise_contribution_standard['Background Noise'].append((sky_background* camera.qe)/total_noise)

        return noise_contribution_standard



def plot_snr_standard(star_data, exp_t, telescope, camera, log_xscale = False, log_yscale = False, save_plots = False):
    start_time = time.time()
    fig, ax = plt.subplots()

    snr_data = calc_snr_standard(star_data, exp_t, telescope,camera)

    ax.scatter(snr_data['mag'], snr_data['actual_snr'], label = 'Theoretical SNR')
    plt.axvline(x=calc_mag(exp_t, telescope,camera.saturation()), color='r', linestyle='--', label='Saturation Point')
    plt.axvline(x=calc_mag(exp_t, telescope,min_photon_standard(exp_t,telescope,camera)), color='r', linestyle='solid', label='5-SNR Mag Threshold')
    ax.set_title('Theoretical SNR for Star Magnitudes of {} - {} using a {} Camera for a {} m telescope with a {} second exposure'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()),camera.name, telescope.diameter,exp_t))
    ax.set_xlabel('Star Mag')
    ax.set_ylabel('SNR')
    ax.legend()

    if log_xscale == True:
        ax.set_xscale('log')

    if log_yscale == True:
        ax.set_yscale('log')

    if save_plots == True:
        plt.close()

        folder_path = 'C:/Users/Muhamad/Documents/Python Scripts/AO/plots/{}_{}_{}_{}'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()), telescope.diameter,exp_t)
        create_folder(folder_path)
        fig.savefig(folder_path + '/standard')

    else:
        plt.show()

    end_time = time.time()  # Stop the timer
    runtime = end_time - start_time


def plot_standard_noise_contributions(star_data, exp_t, telescope, camera, save_plots = False):
    start_time = time.time()
    noise = calc_standard_noise_contribution(star_data, exp_t, telescope, camera)

    fig, ax = plt.subplots()
    
    for key in noise:
        ax.scatter(star_data['mag'],noise[key], label = key)

    ax.set_title('Noise Contribution for Standard SNR model for Star Magnitudes of {} - {} using a {} Camera for a {} m telescope with a {} second exposure'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()), camera.name, telescope.diameter,exp_t))
    ax.set_xlabel('Star Mag')
    ax.set_ylabel('Contribution')
    ax.legend()

    if save_plots == True:
        plt.close()

        folder_path = 'C:/Users/Muhamad/Documents/Python Scripts/AO/plots/{}_{}_{}_{}'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()), telescope.diameter,exp_t)
        create_folder(folder_path)
        fig.savefig(folder_path + '/noise_standard')

    else:
        plt.show()

    end_time = time.time()  # Stop the timer
    runtime = end_time - start_time



def plot_shwf_snr_standard(star_data, exp_t, telescope, camera, shwf, log_xscale = False, log_yscale = False, save_plots = False):
    start_time = time.time()
    fig, ax = plt.subplots()

    snr_data = calc_shwf_snr_standard(star_data, exp_t, telescope,camera,shwf)

    ax.scatter(snr_data['mag'], snr_data['actual_snr'], label = 'Theoretical SNR')
    plt.axvline(x=calc_mag(exp_t, telescope,camera.saturation()), color='r', linestyle='--', label='Saturation Point')
    plt.axvline(x=calc_mag(exp_t, telescope,min_photon_standard(exp_t,telescope,camera)), color='r', linestyle='solid', label='5-SNR Mag Threshold')
    ax.set_title('Theoretical SNR for Star Magnitudes of {} - {} using a {} Camera for a {} m telescope with a {} second exposure and {} subaperatures'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()),camera.name, telescope.diameter,exp_t, shwf.n_subaps))
    ax.set_xlabel('Star Mag')
    ax.set_ylabel('SNR')
    ax.legend()

    if log_xscale == True:
        ax.set_xscale('log')

    if log_yscale == True:
        ax.set_yscale('log')

    if save_plots == True:
        plt.close()

        folder_path = 'C:/Users/Muhamad/Documents/Python Scripts/AO/plots/{}_{}_{}_{}'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()), telescope.diameter,exp_t)
        create_folder(folder_path)
        fig.savefig(folder_path + '/standard')

    else:
        plt.show()

    end_time = time.time()  # Stop the timer
    runtime = end_time - start_time



def plot_shwf_noise_contributions(star_data, exp_t, telescope, camera, shwf, save_plots = False):
    start_time = time.time()
    noise = calc_shwf_noise_contribution(star_data, exp_t, telescope, camera, shwf)

    fig, ax = plt.subplots()
    
    for key in noise:
        ax.scatter(star_data['mag'],noise[key], label = key)

    ax.set_title('Noise Contribution for Standard SNR model for Star Magnitudes of {} - {} using a {} Camera for a {} m telescope with a {} second exposure'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()), camera.name, telescope.diameter,exp_t))
    ax.set_xlabel('Star Mag')
    ax.set_ylabel('Contribution')
    ax.legend()

    if save_plots == True:
        plt.close()

        folder_path = 'C:/Users/Muhamad/Documents/Python Scripts/AO/plots/{}_{}_{}_{}'.format(np.floor(star_data['mag'].min()),np.ceil(star_data['mag'].max()), telescope.diameter,exp_t)
        create_folder(folder_path)
        fig.savefig(folder_path + '/noise_standard')

    else:
        plt.show()

    end_time = time.time()  # Stop the timer
    runtime = end_time - start_time


def search_stars_old(star_data, telescope, num_searches, mag_thresh):
    # Define Zones
    zones = [
        {'name': 'Galactic Centre', 'b_min': -25, 'b_max': 25},
        {'name': 'Galactic North', 'b_min': 25, 'b_max': 90},
        {'name': 'Galactic South', 'b_min': -90, 'b_max': -25}
    ]

    rows_list = {
        'Galactic Centre': {'l':[], 'b':[], 'mag':[]},
        'Galactic North': {'l':[], 'b':[], 'mag':[]},
        'Galactic South': {'l':[], 'b':[], 'mag':[]}
    }

    # Loop through each zone
    for zone in zones:


        # Perform the specified number of searches
        for _ in range(num_searches):
            # Generate random coordinates within the zone
            l = random.uniform(-180, 180)
            b = random.uniform(zone['b_min'], zone['b_max'])


            # Filter stars within the field of view
            for row in star_data.iterrows():
                if (row[1]['l'] >= l - telescope.fov and row[1]['l'] <= l + telescope.fov) and (row[1]['b'] >= b - telescope.fov and row[1]['b'] <= b + telescope.fov) and (row[1]['mag'] < mag_thresh):
                    rows_list[zone['name']]['l'].append(row[1]['l'])
                    rows_list[zone['name']]['b'].append(row[1]['b'])
                    rows_list[zone['name']]['mag'].append(row[1]['mag'])

            
    gc_df = pd.DataFrame(rows_list['Galactic Centre'])
    gn_df = pd.DataFrame(rows_list['Galactic North'])
    gs_df = pd.DataFrame(rows_list['Galactic South'])


    rounded_mag_gc = []
    rounded_mag_gn = []
    rounded_mag_gs = []

    for _mag in gc_df['mag']:
        rounded_mag_gc.append(np.floor(_mag))

    for _mag in gn_df['mag']:
        rounded_mag_gn.append(np.floor(_mag))

    for _mag in gs_df['mag']:
        rounded_mag_gs.append(np.floor(_mag))

    fig = plt.figure()
    labels_gc, counts_gc = np.unique(rounded_mag_gc, return_counts=True)
    labels_gn, counts_gn = np.unique(rounded_mag_gn, return_counts=True)
    labels_gs, counts_gs = np.unique(rounded_mag_gs, return_counts=True)
    plt.bar(labels_gc, counts_gc, align='center', color = 'r', label = 'Galactic Centre')
    plt.bar(labels_gn, counts_gn, align='center', color = 'g', label = 'Galactic North')
    plt.bar(labels_gs, counts_gs, align='center', color = 'b', label = 'Galactic South')
    plt.xlabel('Star Magnitude')
    plt.ylabel('Counts')
    plt.title('Number of stars found in a random search for each region')
    plt.legend()
    #plt.gca().set_xticks(labels)
    plt.show()


    print('Number of Stars found with a threshold magnitude greater than {:0.2f} after {} random searches with a FOV of {} degrees squared'.format(mag_thresh, num_searches, telescope.fov))
    print('Galactic Centre:', len(gc_df['mag']))
    print('Galactic North', len(gn_df['mag']))
    print('Galactic South:', len(gs_df['mag']))

def search_stars(star_data, telescope, mag_thresh, num_searches = 10):
    start_time = time.time()
    # Define Zones
    zones = [
        {'name': 'Galactic Centre', 'b_min': -25, 'b_max': 25},
        {'name': 'Galactic North', 'b_min': 25, 'b_max': 90},
        {'name': 'Galactic South', 'b_min': -90, 'b_max': -25}
    ]

    rows_list = {
        'Galactic Centre': {'l':[], 'b':[], 'mag':[]},
        'Galactic North': {'l':[], 'b':[], 'mag':[]},
        'Galactic South': {'l':[], 'b':[], 'mag':[]}
    }

    # Convert star_data DataFrame to arrays
    star_l = star_data['l'].values
    star_b = star_data['b'].values
    star_mag = star_data['mag'].values

    # Loop through each zone
    for zone in zones:
        # Perform the specified number of searches
        for _ in range(num_searches):
            # Generate random coordinates within the zone
            l = random.uniform(-180, 180)
            b = random.uniform(zone['b_min'], zone['b_max'])

            # Filter stars within the field of view
            fov_mask = ((star_l >= l - telescope.fov) & (star_l <= l + telescope.fov) & 
                        (star_b >= b - telescope.fov) & (star_b <= b + telescope.fov) &
                        (star_mag < mag_thresh))

            filtered_l = star_l[fov_mask]
            filtered_b = star_b[fov_mask]
            filtered_mag = star_mag[fov_mask]

            rows_list[zone['name']]['l'].extend(filtered_l)
            rows_list[zone['name']]['b'].extend(filtered_b)
            rows_list[zone['name']]['mag'].extend(filtered_mag)

    gc_df = pd.DataFrame(rows_list['Galactic Centre'])
    gn_df = pd.DataFrame(rows_list['Galactic North'])
    gs_df = pd.DataFrame(rows_list['Galactic South'])

    rounded_mag_gc = np.floor(gc_df['mag'])
    rounded_mag_gn = np.floor(gn_df['mag'])
    rounded_mag_gs = np.floor(gs_df['mag'])

    fig = plt.figure()
    labels_gc, counts_gc = np.unique(rounded_mag_gc, return_counts=True)
    labels_gn, counts_gn = np.unique(rounded_mag_gn, return_counts=True)
    labels_gs, counts_gs = np.unique(rounded_mag_gs, return_counts=True)
    plt.bar(labels_gc, counts_gc, align='center', color='r', label='Galactic Centre')
    plt.bar(labels_gn, counts_gn, align='center', color='g', label='Galactic North', width = 0.50)
    plt.bar(labels_gs, counts_gs, align='center', color='b', label='Galactic South', width = 0.25)
    plt.xlabel('Star Magnitude')
    plt.ylabel('Counts')
    plt.title('Number of Stars found with a threshold magnitude greater than {:0.2f} after {} random searches with an FOV of {} degrees squared'.format(mag_thresh, num_searches, telescope.fov))
    plt.legend()
    plt.show()


    end_time = time.time()  # Stop the timer
    runtime = end_time - start_time

    print('Galactic Centre:', len(gc_df['mag']))
    print('Galactic North', len(gn_df['mag']))
    print('Galactic South:', len(gs_df['mag']))

