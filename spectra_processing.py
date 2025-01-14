#! /usr/bin/env python3

'''
Written by Oleksandra Rebrysh.
'''

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import glob

from astropy.modeling.fitting import LevMarLSQFitter, LinearLSQFitter
from astropy.modeling import models
from astropy import units as u
import lineid_plot

def LoadData(filename):

    df = pd.read_csv(filename, skiprows = 14, sep = '\t', header = None)
    #print(df)

    df.columns = ['wavelength', 'intensity']
    # Remove commas and convert to numeric values for each relevant column
    df['wavelength'] = pd.to_numeric(df['wavelength'].str.replace(',', '.'))
    df['intensity'] = pd.to_numeric(df['intensity'].str.replace(',', '.'), downcast='float')

    return df

def LoadAllSpectra(path, common_part):
    # Use glob to find all files matching the pattern
    pattern = path + common_part + '*.txt'
    files = glob.glob(pattern)

    # Load each file into a DataFrame and store in a list
    spectra = []
    for filename in files:
        df = LoadData(filename)
        spectra.append(df)

    return spectra

def AverageSpectra(spectra):
    '''
    Please have your spectra in the following form:
    spectra = [spectrum1,
                spectrum2]
    '''

    combined_spectra = pd.concat(spectra)
    average_spectra = combined_spectra.groupby('wavelength', as_index=False)['intensity'].mean()

    return average_spectra

def ProcessDataEach(rawdata, background, dark, subtract_background = False, subtract_dark = False):
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
        #processed_df['intensity'] = processed_df['intensity'].clip(lower=0)  # Ensure non-negative intensities
        processed_data.append(processed_df)

    return processed_data

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

    wavelength = sum_spectra.loc[:, 'wavelength'].values
    intensity = sum_spectra.loc[:, 'intensity'].values

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

def DetectAndPlotLines(spectra, wavelength, spectral_lines, telluric_lines, threshold, title = None, save = False, path = None, name = None):

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

#if __name__ == '__main__':

    #path = ''
    #savefold = ''

    #spectral_lines = {
    #"Hα": 656.3,
    #"Hβ": 486.1,
    #"Hγ": 434.0,
    #"Hδ": 410.2,
    #"He II": 420.0,
    #"He II": 454.1,
    #"He I": 447.1,
    #"He I": 402.6,
    #"He I": 667.8,
    #"Fe I": 495.8,
    #"Fe I": 466.8,
    #"Fe I": 438.4,
    #"Ca I": 420.8,
    #"Fe I": 527.0,
    #"Fe II": 516.90,
    #"Mg I": 518.0,
    #"Na I D1": 589.00,
    #"Na I D2": 589.60,
    #"Ca II H": 396.85,
    #"Ca II K": 393.37,
    #"Ca II IR 1": 849.80,
    #"Ca II IR 2": 854.20,
    #"Ca II IR 3": 866.20,
    #"[O I] 1": 630.0,
    #"[O I] 2": 636.4,
    #"C II": 426.7,
    #"Si II": 412.8,
    #"Si II": 634.7,
    #"Si II": 637.1,
    #"Mg II": 448.1,
    #"O I": 898.8,
    #"O I": 822.7,
    #"O I": 759.4,
    #"O I": 686.7,
    #"O I": 627.7,
    #"O I": 777.1,
    #"O I": 777.4,
    #"O I": 777.5,
    #"He I": 587.6,
    #"Ti II": 336.1,
    #"Ni I": 299.4,
    #"TiO": 476.1,
    #"TiO": 495.4,
    #"TiO": 516.7
    #}

    #telluric_lines = {"TL":
    #[687.8,
     #718.5,
     #719.4,
     #725.0,
    #759.3,  # O2 (760.0 nm region)
    #760.5,
    #761.0,  # O2 (760.5 nm region)
    #762.0,  # O2 (762.0 nm region)
    #764.0,  # H2O (763.0-764.0 nm region)
    #820.5,  # H2O (820.0-821.0 nm region)
    #822.2,  # H2O (822.0 nm region)
    #935.0,  # H2O (934.0-935.0 nm region)
    #940.0,  # H2O (940.0 nm region)
    #942.0,  # H2O (942.0 nm region)
    #940.5,  # H2O (940.5 nm region)
    #943.0,  # H2O (943.0 nm region)
    #946.0,  # H2O (946.0 nm region)
    #953.0,  # H2O (953.0 nm region)
    #960.0,  # H2O (960.0 nm region)
#]}
