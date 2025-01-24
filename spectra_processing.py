#! /usr/bin/env python3

'''
Written by Oleksandra Rebrysh.
'''

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob
import os

from scipy.interpolate import CubicSpline
from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter, SimplexLSQFitter
from astropy.modeling import models
from astropy import units as u
import lineid_plot

def LoadData(filename, skiprows = 14):

    df = pd.read_csv(filename, skiprows = skiprows, sep = '\t', header = None)
    #print(df)

    df.columns = ['wavelength', 'intensity']

    # Convert all data to strings first, then replace ',' with '.'
    df['wavelength'] = df['wavelength'].astype(str).str.replace(',', '.')

    # Remove commas and convert to numeric values for each relevant column
    df['wavelength'] = pd.to_numeric(df['wavelength'],downcast='float')

    # Convert all data to strings first, then replace ',' with '.'
    df['intensity'] = df['intensity'].astype(str).str.replace(',', '.')

    # Convert the column back to numeric, coercing any invalid values to NaN
    df['intensity'] = pd.to_numeric(df['intensity'],  downcast='float')


    return df


def LoadAllSpectra(path, common_part):
    # Use glob to find all files matching the pattern
    pattern = path + common_part + '*.txt'
    files = glob.glob(pattern)

    # Load each file into a DataFrame and store in a list
    spectra = []
    for filename in files:
        print(filename)
        df = LoadData(filename)

        spectra.append(df)

    return spectra

def LoadAttenuation(path, filename):

    file_path = os.path.join(path,filename)
    _, file_extension = os.path.splitext(file_path)

    if file_extension in ['.xlsx', '.xls', '.xlsm']:
        df = pd.read_excel(file_path, sheet_name='Attenuation Data', usecols=[2,3], skiprows = 1)
    elif file_extension in ['.txt', '.csv']:
        df = pd.read_csv(file_path, skiprows = 1, sep = '\t', header = None)
    #print(df)
    df = df.dropna()
    df.columns = ['wavelength', 'intensity']
    # Remove commas and convert to numeric values for each relevant column
    df['wavelength'] = pd.to_numeric(df['wavelength'])
    # print(df['wavelength'])
    df['intensity'] = pd.to_numeric(df['intensity'], downcast='float')
    # print(df['intensity'])

    # plt.plot(df['wavelength'], df['intensity'])
    # plt.show()

    return df

def CorrectAttenuation(spectra, attenuation):


    attenuation_wvl = attenuation.loc[:, 'wavelength'].values
    attenuation_coef = attenuation.loc[:, 'intensity'].values
    # Convert attenuation to transmission efficiency
    transmission = np.power(10, -(attenuation_coef*1e-3) / 10)

    corrected_spectra = []
    for data in spectra:

        corrected_df = data.copy()
        wavelength = data.loc[:, 'wavelength'].values
        intensity = data.loc[:, 'intensity'].values
        minwavelength = np.min(attenuation_wvl)
        maxwavelength = np.max(wavelength)
        # print(minwavelength)

        mask = (wavelength >= minwavelength) & (wavelength <= maxwavelength)
        wavelength_truncated = wavelength[mask]
        intensity_truncated = intensity[mask]
        # print(intensity_truncated.shape)

        transmission_interpolate = CubicSpline(attenuation_wvl, transmission)
        transmission_corsize = transmission_interpolate(wavelength_truncated)

        # print(intensity)
        corrected_df.loc[mask, 'intensity'] = intensity_truncated / transmission_corsize
        corrected_spectra.append(corrected_df)


        # Plot original and corrected spectra
        # fig, ax = plt.subplots(figsize=(10, 8))
        # ax.plot(corrected_df['wavelength'], corrected_df['intensity'], label='Corrected Spectra', color='green')
        # ax.plot(wavelength, intensity, label='Original Spectra', color='red')
        # ax.set_xlabel('Wavelength')
        # ax.set_ylabel('Intensity')
        # ax.legend()
        # ax.grid(True)
        # plt.title("Original vs Corrected Spectra")
        # plt.show()


    return corrected_spectra

def CalcSNRData(sum_data, sigma_dark, sigma_background,subtract_background = False, subtract_dark = False,plot_snr=False):
    #required_columns = ['wavelength', 'intensity']

    ## Check if input data frames are valid
    #for df in [background, dark, *rawdata]:
        #if not isinstance(df, pd.DataFrame) or set(df.columns) != set(required_columns):
            #raise ValueError("Invalid input data frame format or columns.")

    ## Check if input data frames have the same length
    #background_length = len(background)
    #dark_length = len(dark)
    #if any(len(data) != background_length or len(data) != dark_length for data in rawdata):
        #raise ValueError("Input data frames must have the same length.")




    wavelength = sum_data.loc[:, 'wavelength'].values
    intensity = sum_data.loc[:, 'intensity'].values

    # Add all background and dark contributions for each spectrum
    noise_aux=0

    processed_df = sum_data.copy()
    if subtract_background:
        noise_aux = sigma_background['intensity'].values**2

    # print(noise_aux)
    if subtract_dark:
        noise_aux += sigma_dark['intensity'].values**2
    # print(noise_aux)
    noise_aux += processed_df['intensity'].values
    # print(noise_aux)
    snr = processed_df['intensity'].values/np.sqrt(noise_aux)
    #processed_df['intensity'] = processed_df['intensity'].clip(lower=0)  # Ensure non-negative intensities

    if plot_snr:
        fig, ax = plt.subplots(2, 1, figsize=(10, 8))

        ax[0].plot(wavelength, snr,label= 'SNR')

        ax[0].minorticks_on()
        ax[0].xaxis.set_minor_locator(plt.MultipleLocator(10))
        ax[0].grid(True, 'both')
        ax[0].set_xlabel(f'Wavelength, nm')
        ax[0].set_ylabel(f'SNR')
        ax[0].set_title('Signal to Noise Ratio')
        ax[0].legend()


        print(sigma_dark)
        ax[1].plot(wavelength,sigma_background['intensity'].values,label='Background Noise')
        ax[1].plot(wavelength,sigma_dark['intensity'].values,label='Dark Noise')
        ax[1].plot(wavelength, np.sqrt(processed_df['intensity'].values),label= 'Spectra Noise')
        ax[1].plot(wavelength, np.sqrt(noise_aux),label= 'Total Noise')

        ax[1].minorticks_on()
        ax[1].xaxis.set_minor_locator(plt.MultipleLocator(10))
        ax[1].grid(True, 'both')
        ax[1].set_xlabel(f'Wavelength, nm')
        ax[1].set_ylabel(f'Intensity, counts')
        ax[1].set_title('Noise Contributors')
        ax[1].legend()


        plt.show()



    return snr

def AverageSpectra(spectra):
    '''
    Please have your spectra in the following form:
    spectra = [spectrum1,
                spectrum2]
    '''

    combined_spectra = pd.concat(spectra)
    average_spectra = combined_spectra.groupby('wavelength', as_index=False)['intensity'].mean()

    return average_spectra

def ProcessDataEach(rawdata, background, dark, flat, subtract_background = False, subtract_dark = False, divide_flat = False):
    #required_columns = ['wavelength', 'intensity']

    ## Check if input data frames are valid
    #for df in [background, dark, *rawdata]:
        #if not isinstance(df, pd.DataFrame) or set(df.columns) != set(required_columns):
            #raise ValueError("Invalid input data frame format or columns.")

    ## Check if input data frames have the same length
    #background_length = len(background)
    #dark_length = len(dark)
    #if any(len(data) != background_length or len(data) != dark_length for data in rawdata):
        #raise ValueError("Input data frames must have the same length.")

    processed_data = []

    # Subtract background and dark from each raw data spectrum
    for data in rawdata:
        processed_df = data.copy()
        if subtract_background:
            processed_df['intensity'] -= background['intensity'].values
        if subtract_dark:
            processed_df['intensity'] -= dark['intensity'].values
        if divide_flat:
            processed_df['intensity'] /= flat['intensity'].values/np.max(flat['intensity'])
        processed_data.append(processed_df)

    return processed_data

def Build_Instrument_Response(data, reference, plot=False):
    #required_columns = ['wavelength', 'intensity']

    #get wavelengths of our data and of our reference
    wavelength = data['wavelength'].values
    intensity  =data['intensity'].values
    wavelength_stand = reference['wavelength'].values
    intensity_stand = reference['intensity'].values

    #Mininum and maximum wavelength of Ocean View
    minwavelength = 400.0
    maxwavelength = np.max(wavelength)

    #Mask reference data to have values only inside our ocean range
    mask = (wavelength_stand >= minwavelength) & (wavelength_stand <= maxwavelength)

    #mask our data to have values only above 350 mm
    mask_2 =  (wavelength >= minwavelength) 
    # print(mask, mask.shape)
    # print(mask_2, mask_2.shape)

    wavelength_truncated = wavelength_stand[mask]
    wavelength = wavelength[mask_2]
    intensity_truncated = intensity_stand[mask]
    intensity = intensity[mask_2]
    # print(intensity_truncated.shape)



    #Interpolate the reference data
    interpolate = CubicSpline(wavelength_truncated, intensity_truncated)

    #Get the values of the reference data at our ocean wavelenghts
    intensity_standard = interpolate(wavelength)

    if plot==True:
        fig1, axs1 = plt.subplots(2, 1, figsize=(12, 8))
        PlotData(axs1[0], reference, 'Reference Data Raw', 'Raw data')
        axs1[1].plot(wavelength_truncated,intensity_truncated)
        plt.show()

    #Save this data on a pandas dataframe
    d = {'wavelength': wavelength, 'intensity': intensity}
    norma_data_df = pd.DataFrame(data=d)
    # print(norma_data_df['wavelength'])

    d1 = {'wavelength': wavelength, 'intensity': intensity_standard}
    norma_ref_df = pd.DataFrame(data=d1)

    # Calculate scale factor and normalize observed counts
    ourdata_sum = np.sum(intensity)
    standard_sum = np.sum(intensity_standard)
    scale_factor = standard_sum / ourdata_sum


    norma_data_df['intensity'] *= scale_factor


    ins_response = norma_data_df.copy()
    
    ins_response['intensity'] = norma_data_df['intensity'].values/norma_ref_df['intensity'].values
    

    if plot==True:
        fig, axs = plt.subplots(2, 1, figsize=(12, 8))
        PlotData(axs[0], norma_ref_df, 'Reference Data Normalized', 'Reference Spectrophotometric Standard')
        PlotData(axs[0], norma_data_df, 'Data Normalized', 'Ocean Data Spectrophotometric Standard')
        PlotData(axs[1], ins_response, 'Data/Ref', 'Instrument Response Function')
        plt.show()


    return ins_response

def ProcessData(rawdata, background, dark):
    # Ensure the dataframes have the expected columns
    required_columns = ['wavelength', 'intensity']

    for df in [rawdata, background, dark]:
        if df.shape[1] != 2:
            print(f'ERROR EXCEPTION: Data should have exactly two columns')
            return None
        df.columns = required_columns

    # Ensure the lengths of the dataframes match
    if not (rawdata.shape[0] == background.shape[0] == dark.shape[0]):
        print(f'ERROR EXCEPTION: Data should have the same length')
        return None

    # Process the data
    processed_data = rawdata.copy()
    processed_data['intensity'] = rawdata['intensity'] - background['intensity'] - dark['intensity']

    return processed_data

def AddSpectra(spectra):

    '''
    Please have your spectra in the following form:
    spectra = [spectrum1,
                spectrum2]
    '''

    combined_spectra = pd.concat(spectra)
    summed_spectra = combined_spectra.groupby('wavelength', as_index=False)['intensity'].sum()
    #print(summed_spectra)
    return summed_spectra

def PlotData(ax, data, name, labelname):
    #print(data)
    wavelength = data.loc[:, 'wavelength']
    intensity = data.loc[:, 'intensity']
    ax.plot(wavelength, intensity, label = labelname)
    ax.minorticks_on()
    ax.xaxis.set_minor_locator(plt.MultipleLocator(10))
    ax.grid(True, 'both')
    ax.set_xlabel(f'Wavelength, nm')
    ax.set_ylabel(f'Intensity, counts')
    ax.set_title(name)
    ax.legend()

    return

def DefSpectrumForLines(spectr, minwavelength, maxwavelength):

    wavelength = spectr.loc[:, 'wavelength'].values
    intensity = spectr.loc[:, 'intensity'].values

    mask = (wavelength >= minwavelength) & (wavelength <= maxwavelength)
    wavelength_truncated = wavelength[mask]
    intensity_truncated = intensity[mask]
    # spectrum = Spectrum1D(flux=intensity_truncated*u.adu, spectral_axis=wavelength_truncated*u.nm)


    return wavelength_truncated, intensity_truncated

def FitContinuum(degree, wavelength, intensity, image = False):

    chebyshev_model = models.Chebyshev1D(degree=degree)
    fitter = LinearLSQFitter()
    g1_fit = fitter(chebyshev_model, wavelength, intensity)
    spec_norm = (intensity - g1_fit(wavelength))/np.linalg.norm(intensity - g1_fit(wavelength))

    if image == True:
        fig, ax = plt.subplots(2, 1, figsize=(10, 8))
        ax[0].plot(wavelength, intensity, label='Original Spectrum')
        ax[0].plot(wavelength, g1_fit(wavelength), label='Fitted Continuum')
        ax[0].set_title('Continuum Fitting')
        ax[0].legend()
        ax[0].grid(True)

        # ax[1].plot(wavelength_truncated, spec_norm.flux)
        ax[1].plot(wavelength, spec_norm)
        ax[1].set_title('Continuum Normalized Spectrum')
        ax[1].grid(True)

    return spec_norm

def BlackBodyFit(temperature, wavelength, intensity, image = False):

    bb = models.BlackBody(temperature=temperature*u.K, scale = 1)
    fitter = SimplexLSQFitter()
    intensity_flux = intensity / np.max(intensity)
    fit = fitter(bb, wavelength*u.nm, intensity_flux)
    # plt.plot(wavelength*u.nm, bb(wavelength*u.nm))
    spec_norm = (intensity_flux - fit(wavelength*u.nm).value)/np.linalg.norm(intensity_flux - fit(wavelength*u.nm).value)

    if image == True:
        fig, ax = plt.subplots(2, 1, figsize=(10, 8))
        ax[0].plot(wavelength, intensity_flux, label='Original Spectrum')
        ax[0].plot(wavelength, fit(wavelength*u.nm), label='Fitted Continuum')
        ax[0].set_title('Continuum Fitting')
        ax[0].legend()
        ax[0].grid(True)

        # ax[1].plot(wavelength_truncated, spec_norm.flux)
        ax[1].plot(wavelength, spec_norm)
        ax[1].set_title('Continuum Normalized Spectrum')
        ax[1].grid(True)


def DetectAndPlotLines(spectra, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = None, save = False, path = None, name = None):

    detected_lines = {}
    telluric_lines_d = {}

    # Step 4: Check for each line if its depth/strength exceeds the threshold
    for line_name, line_wavelength in spectral_lines.items():
        index = np.argmin(np.abs(wavelength_truncated - line_wavelength))
        flux = spectra[index]
        depth = - flux
        if depth >= threshold:
            if flux < spectra[index - 10] and flux < spectra[index + 10]:
                detected_lines[line_name] = line_wavelength

    for line_wavelength in telluric_lines["TL"]:
        index = np.argmin(np.abs(wavelength_truncated - line_wavelength))
        flux = spectra[index]
        depth = - flux
        if depth >= 0.02:
            if flux < spectra[index - 1] and flux < spectra[index + 1]:
                telluric_lines_d[line_wavelength] = line_wavelength

    fig, ax = lineid_plot.plot_line_ids(
        wavelength_truncated,
        spectra,
        list(detected_lines.values()) + list(telluric_lines_d.values()),  # Flux values of detected lines
        list(detected_lines.keys()) + ["TL"]*len(telluric_lines_d),    # Names of detected lines
        max_iter = 100)
    fig.set_size_inches(12, 6)
    num_detected_lines = len(detected_lines)  # Count of detected lines
    num_telluric_lines = len(telluric_lines_d)  # Count of telluric lines

    for index in range(num_detected_lines+1, num_detected_lines + num_telluric_lines+1):
        line = ax.lines[index]
        line.set_color("red")
        line.set_linestyle("--")

    for text in ax.texts:
        if text.get_text() == "TL":  # Target labels with "TL"
            text.set_visible(False)
        if text.get_text() != "TL":  # Ignore telluric lines
            text.set_y(text.get_position()[1] + 0.005)

    ax.set_xlabel('Wavelength, nm')
    ax.set_ylabel('Normalised flux')
    tl = mpl.lines.Line2D([], [], color='r', linestyle = '--')
    ax.legend([tl], ['Telluric line'])
    ax.set_title(title, y = 1.2, fontweight="bold")
    if save:
        plt.savefig(path+name, dpi = 600)


    return
