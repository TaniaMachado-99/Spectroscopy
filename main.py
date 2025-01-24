from spectra_processing import *

if __name__ == '__main__':
 
    savefold = '/Users/taniamachado/Documents/PhD/Oceanview/Asiago_Jan16-18/plots/'
    path = '/Users/taniamachado/Documents/PhD/Oceanview/Asiago_Jan16-18/'
    path_at = '/Users/taniamachado/Documents/PhD/Oceanview/'
    file_at = 'FG_UGA_LGA_UCA_LCA_UEA_LEA_Attenuation_Data_2.xlsx'

    spectral_lines = {
    "Hα": 656.3,
    "Hβ": 486.1,
    "Hγ": 434.0,
    "Hδ": 410.2,
    "He II": 420.0,
    "He II": 454.1,
    "He I": 447.1,
    "He I": 402.6,
    "He I": 667.8,
    "Fe I": 495.8,
    "Fe I": 466.8,
    "Fe I": 438.4,
    "Ca I": 420.8,
    "Fe I": 527.0,
    "Fe II": 516.90,
    "Mg I": 518.0,
    "Na I D1": 589.00,
    "Na I D2": 589.60,
    "Ca II H": 396.85,
    "Ca II K": 393.37,
    "Ca II IR 1": 849.80,
    "Ca II IR 2": 854.20,
    "Ca II IR 3": 866.20,
    "[O I] 1": 630.0,
    "[O I] 2": 636.4,
    "C II": 426.7,
    "Si II": 412.8,
    "Si II": 634.7,
    "Si II": 637.1,
    "Mg II": 448.1,
    "O I": 898.8,
    "O I": 822.7,
    "O I": 759.4,
    "O I": 686.7,
    "O I": 627.7,
    "O I": 777.1,
    "O I": 777.4,
    "O I": 777.5,
    "He I": 587.6,
    "Ti II": 336.1,
    "Ni I": 299.4,
    "TiO": 476.1,
    "TiO": 495.4,
    "TiO": 516.7
    }

    telluric_lines = {"TL":
    [687.8,
     718.5,
     719.4,
     725.0,
    759.3,  # O2 (760.0 nm region)
    760.5,
    761.0,  # O2 (760.5 nm region)
    762.0,  # O2 (762.0 nm region)
    764.0,  # H2O (763.0-764.0 nm region)
    820.5,  # H2O (820.0-821.0 nm region)
    822.2,  # H2O (822.0 nm region)
    935.0,  # H2O (934.0-935.0 nm region)
    940.0,  # H2O (940.0 nm region)
    942.0,  # H2O (942.0 nm region)
    940.5,  # H2O (940.5 nm region)
    943.0,  # H2O (943.0 nm region)
    946.0,  # H2O (946.0 nm region)
    953.0,  # H2O (953.0 nm region)
    960.0,  # H2O (960.0 nm region)
]}

    threshold = 0.02

    #=======================================================================================#
    #                                  Zeta Cas Spectrophometric Standard                   #
    #=======================================================================================#
    #Reference Data Read
    ref_spectra = LoadData(path_at + 'zetaCas_standard.txt',skiprows=1)
    ref_spectra['wavelength']/=10

    #######
    
    path_aux = path+'zetaCas/'
    spectra = LoadAllSpectra(path_aux, 'zetaCas_data_HDX017711')
    # print(spectra)

    background = LoadAllSpectra(path_aux, 'zetaCas_background_HDX017711')

    dark = LoadAllSpectra(path+'dark_16012025/', 'dark_16s_HDX017711')
    ###print(spectra[0])
    ###print(background[0])
    #flat = LoadAllSpectra('/Users/taniamachado/Documents/PhD/Oceanview/Asiago_Jan16-18/Flats_16012025/', 'flat_16s_HDX017711')


    fiber_at = LoadAttenuation(path_at, file_at)

    data_atcor = CorrectAttenuation(spectra, fiber_at)

    sum_spectra_raw = AddSpectra(spectra)
    background_av = AverageSpectra(background)
    ###print(background_av)
    dark_av = AverageSpectra(dark)
    ###print(dark_av)
    #flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    #flat_av = AverageSpectra(flat_subtr)



    background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    subtr = ProcessDataEach(spectra, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    ####print(f'subtr: \n {subtr}')

    sum_spectra = AddSpectra(subtr)
    inst_response = Build_Instrument_Response(sum_spectra, ref_spectra, plot=True)
    wavelength = sum_spectra.loc[:, 'wavelength'].values
    intensity = sum_spectra.loc[:, 'intensity'].values
    
    wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)
    
    spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)

    DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'zetaCas, x [m(V) = x]', save = False, path = savefold, name = "zetaCas_spectrum_lines.png")

    fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    PlotData(axs[0], spectra[0], 'zetaCas', 'zetaCas raw data')
    PlotData(axs[0], background_av, 'zetaCas', 'background')
    PlotData(axs[0], dark_av, 'zetaCas', 'dark')
    #PlotData(axs[0], flat_av, 'Alkaid, x [m(V) = x]', 'flat')
    PlotData(axs[1], sum_spectra, None, 'zetaCas processed data')
    plt.tight_layout()
    plt.savefig(savefold+'/zetaCas.png', dpi = 600)

    fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    PlotData(axs, sum_spectra, 'zetaCas, x [m(V) = x]', 'zetaCas processed data')
    plt.savefig(savefold+'/zetaCas_spectrum.png', dpi = 600)
    plt.show()
    '''
    #=======================================================================================#
    #                                  Castor                                               #
    #=======================================================================================#
    '''

    spectra = LoadAllSpectra(path, 'Castor_data_HDX017711')
    fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    PlotData(axs, spectra[0], 'Castor', 'Castor')
    #plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Castor_spectrum.png', dpi = 600)
    plt.show()
    '''


    #=======================================================================================#
    #                                  Alkaid                                               #
    #=======================================================================================#
    '''
    path_aux = path+'Alkaid/'
    spectra = LoadAllSpectra(path_aux, 'Alkaid_spectrum_HDX017711')
    # print(spectra)

    background = LoadAllSpectra(path_aux, 'Alkaid_background_HDX017711')

    dark = LoadAllSpectra(path+'dark/', 'dark_HDX017711')
    ###print(spectra[0])
    ###print(background[0])
    flat = LoadAllSpectra('/Users/taniamachado/Documents/PhD/Oceanview/Asiago_Jan16-18/Flats_16012025/', 'flat_16s_HDX017711')


    fiber_at = LoadAttenuation(path_at, file_at)

    data_atcor = CorrectAttenuation(spectra, fiber_at)

    sum_spectra_raw = AddSpectra(spectra)
    background_av = AverageSpectra(background)
    ###print(background_av)
    dark_av = AverageSpectra(dark)
    ###print(dark_av)
    flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    flat_av = AverageSpectra(flat_subtr)

    fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    PlotData(axs, flat_av, 'Flat', 'Alkaid Flat')

    background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    subtr = ProcessDataEach(spectra, None, dark_av, flat_av, subtract_background = False, subtract_dark = True, divide_flat = True)
    ####print(f'subtr: \n {subtr}')

    sum_spectra = AddSpectra(subtr)
    wavelength = sum_spectra.loc[:, 'wavelength'].values
    intensity = sum_spectra.loc[:, 'intensity'].values

    wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)

    spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)

    DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Alkaid, x [m(V) = x]', save = False, path = savefold, name = "Alkaid_spectrum_lines.png")

    fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    PlotData(axs[0], spectra[0], 'Alkaid', 'Alkaid raw data')
    PlotData(axs[0], background_av, 'Alkaid', 'background')
    PlotData(axs[0], dark_av, 'Alkaid', 'dark')
    PlotData(axs[0], flat_av, 'Alkaid, x [m(V) = x]', 'flat')
    PlotData(axs[1], sum_spectra, None, 'Alkaid processed data')
    plt.tight_layout()
    plt.savefig(savefold+'/Alkaid.png', dpi = 600)

    # fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    # PlotData(axs, sum_spectra, 'Alkaid, x [m(V) = x]', 'Alkaid processed data')
    # plt.savefig(savefold+'/Alkaid_spectrum.png', dpi = 600)
    # plt.show()
    #=======================================================================================#
    #                                  Capella                                               #
    #=======================================================================================#
    
    path_aux = path+'Capella/'
    spectra = LoadAllSpectra(path_aux, 'Capella_data_HDX017711')
    # print(spectra)

    background = LoadAllSpectra(path_aux, 'Capella_background_HDX017711')

    dark = LoadAllSpectra(path+'dark/', 'dark_HDX017711')
    ###print(spectra[0])
    ###print(background[0])
    flat = LoadAllSpectra(path+'flat/', 'flat_HDX017711')

    fiber_at = LoadAttenuation(path_at, file_at)

    data_atcor = CorrectAttenuation(spectra, fiber_at)

    sum_spectra_raw = AddSpectra(spectra)
    background_av = AverageSpectra(background)
    ###print(background_av)
    dark_av = AverageSpectra(dark)
    ###print(dark_av)
    flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    flat_av = AverageSpectra(flat_subtr)

    background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    subtr = ProcessDataEach(spectra, None, dark_av, flat_av, subtract_background = False, subtract_dark = True, divide_flat = True)
    ####print(f'subtr: \n {subtr}')

    sum_spectra = AddSpectra(subtr)
    wavelength = sum_spectra.loc[:, 'wavelength'].values
    intensity = sum_spectra.loc[:, 'intensity'].values

    wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)

    spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)

    DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Capella, α Aurigae [m(V) = 0.08]', save = False, path = savefold, name = "Capella_spectrum_lines.png")

    fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    PlotData(axs[0], spectra[0], 'Capella', 'Capella raw data')
    PlotData(axs[0], background_av, 'Capella', 'background')
    PlotData(axs[0], dark_av, 'Capella', 'dark')
    PlotData(axs[0], flat_av, 'Capella, α Aurigae [m(V) = 0.08]', 'flat')
    PlotData(axs[1], sum_spectra, None, 'Capella processed data')
    plt.tight_layout()
    plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Capella.png', dpi = 600)

    fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    PlotData(axs, sum_spectra, 'Capella, α Aurigae [m(V) = 0.08]', 'Capella processed data')
    plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Capella_spectrum.png', dpi = 600)
    plt.show()
    '''
    #=======================================================================================#
    #                                  HR6212                                               #
    #=======================================================================================#

    #spectra = LoadAllSpectra(path, 'HR6212_spectrum_HDX017711')
    ###print(spectra)

    #background = LoadAllSpectra(path, 'HR6212_background_HDX017711')

    #dark = LoadAllSpectra(path, 'dark_HDX017711')
    ##print(spectra[0])
    ##print(background[0])

    #sum_spectra_raw = AddSpectra(spectra)
    #background_av = AverageSpectra(background)
    ##print(background_av)
    #dark_av = AverageSpectra(dark)
    ##print(dark_av)
    #background_subtr = ProcessDataEach(background, None, dark_av, subtract_background = False, subtract_dark = True)
    #background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, subtract_background = True, subtract_dark = True)
    ###print(f'subtr: \n {subtr}')

    #sum_spectra = AddSpectra(subtr)
    #wavelength = sum_spectra.loc[:, 'wavelength'].values
    #intensity = sum_spectra.loc[:, 'intensity'].values

    #wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)

    #spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)

    #DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'HR6212, ζ Herculis [m(V) = 2.81]', save = True, path = savefold, name = "HR6212_spectrum_lines.png")

    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], 'HR6212', 'HR6212 raw data')
    # PlotData(axs[0], background_av, 'HR6212', 'background')
    # PlotData(axs[0], dark_av, 'HR6212, ζ Herculis [m(V) = 2.81]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'HR6212 processed data')
    # plt.tight_layout()
    # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/HR6212.png', dpi = 600)

    #fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    #PlotData(axs, sum_spectra, 'HR6212, ζ Herculis [m(V) = 2.81]', 'HR6212 processed data')
    #plt.savefig(f'D:\PhD\oceanview\Spectra_stars\\HR6212_spectrum.png', dpi = 600)
    #plt.show()

    #=======================================================================================#
    #                                  Arcturus                                               #
    #=======================================================================================#
    #

    #spectra = LoadAllSpectra(path, 'Arcturus_spectrum_HDX017711')
    #print(spectra)

    #background = LoadAllSpectra(path, 'Arcturus_background_HDX017711')

    #dark = LoadAllSpectra(path, 'dark_3s_HDX017711')

    #flat = LoadAllSpectra(path, 'flat_HDX017711')
    ###print(spectra[0])
    ####print(background[0])
    #fiber_at = LoadAttenuation(path_at, file_at)

    #data_atcor = CorrectAttenuation(spectra, fiber_at)

    ## sum_spectra_raw = AddSpectra(spectra)
    #background_av = AverageSpectra(background)
    ###print(background_av)
    #dark_av = AverageSpectra(dark)

    #flat_subtr = ProcessDataEach(flat, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    #flat_av = AverageSpectra(flat_subtr)
    ###print(dark_av)
    #background_subtr = ProcessDataEach(background, None, dark_av, None, subtract_background = False, subtract_dark = True, divide_flat = False)
    #background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(data_atcor, background_subtr_av, dark_av, flat_av, subtract_background = True, subtract_dark = True, divide_flat = True)
    ###print(f'subtr: \n {subtr}')

    #sum_spectra = AddSpectra(subtr)

    #wavelength = sum_spectra.loc[:, 'wavelength'].values
    #intensity = sum_spectra.loc[:, 'intensity'].values


    #wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)
    ## spec_norm = BlackBodyFit(4200, wavelength_truncated, intensity_truncated, True)
    #spec_norm = FitContinuum(16, wavelength_truncated, intensity_truncated, image = True)

    #DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Arcturus, α Boötis [m(V) = -0.05]', save = False, path = savefold, name = "Arcturus_spectrum_lines.png")
    # #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # # PlotData(axs[0], spectra[0], 'Arcturus', 'Arcturus raw data')
    # PlotData(axs[0], background_av, 'Arcturus', 'background')
    # PlotData(axs[0], flat_av, 'Arcturus', 'flat')
    # PlotData(axs[0], dark_av, 'Arcturus, α Boötis [m(V) = -0.05]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'Arcturus processed data')
    # plt.tight_layout()
    # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Arcturus.png', dpi = 600)

    #fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    #PlotData(axs, sum_spectra, 'Arcturus, α Boötis [m(V) = -0.05]', 'Arcturus processed data')
    #plt.savefig(f'D:\PhD\oceanview\Spectra_stars\\Arcturus_spectrum.png', dpi = 600)
    #plt.show()
    #
    # #=======================================================================================#
    # #                                  Deneb                                              #
    # #=======================================================================================#
    #
    #spectra = LoadAllSpectra(path, 'Deneb_spectrum_HDX017711')
    ##print(spectra)

    #background = LoadAllSpectra(path, 'Deneb_background_HDX017711')

    #dark = LoadAllSpectra(path, 'dark_HDX017711')
    ###print(spectra[0])
    ####print(background[0])

    #sum_spectra_raw = AddSpectra(spectra)
    #background_av = AverageSpectra(background)
    ##print(background_av)
    #dark_av = AverageSpectra(dark)
    ##print(dark_av)
    #background_subtr = ProcessDataEach(background, None, dark_av, subtract_background = False, subtract_dark = True)
    #background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, subtract_background = True, subtract_dark = True)
    ####print(f'subtr: \n {subtr}')

    #sum_spectra = AddSpectra(subtr)

    #wavelength = sum_spectra.loc[:, 'wavelength'].values
    #intensity = sum_spectra.loc[:, 'intensity'].values

    #wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)

    #spec_norm = FitContinuum(12, wavelength_truncated, intensity_truncated, image = True)

    #DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Deneb, α Cygni [m(V) = 1.25]', save = True, path = savefold, name = "Deneb_spectrum_lines.png")
    ##
    ## fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    ## PlotData(axs[0], spectra[0], None, 'Deneb raw data')
    ## PlotData(axs[0], background_av, None, 'background')
    ## PlotData(axs[0], dark_av, 'Deneb, α Cygni [m(V) = 1.25]', 'dark')
    ## PlotData(axs[1], sum_spectra, None, 'Deneb processed data')
    ## plt.tight_layout()
    ## plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Deneb.png', dpi = 600)

    #fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    #PlotData(axs, sum_spectra, 'Deneb, α Cygni [m(V) = 1.25]', 'Deneb processed data')
    #plt.savefig(f'D:\PhD\oceanview\Spectra_stars\\Deneb_spectrum.png', dpi = 600)
    #plt.show()
    #
    # #=======================================================================================#
    # #                                  alpha Oph                                              #
    # #=======================================================================================#
    #
    #spectra = LoadAllSpectra(path, 'alphaOph_spectrum_HDX017711')
    ##print(spectra)

    #background = LoadAllSpectra(path, 'alphaOph_background_HDX017711')

    #dark = LoadAllSpectra(path, 'dark_HDX017711')
    ###print(spectra[0])
    ##print(background[0])

    #sum_spectra_raw = AddSpectra(spectra)
    #background_av = AverageSpectra(background)
    ##print(background_av)
    #dark_av = AverageSpectra(dark)
    ##print(dark_av)
    #background_subtr = ProcessDataEach(background, None, dark_av, subtract_background = False, subtract_dark = True)
    #background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, subtract_background = True, subtract_dark = True)
    ##print(f'subtr: \n {subtr}')

    #sum_spectra = AddSpectra(subtr)

    #wavelength = sum_spectra.loc[:, 'wavelength'].values
    #intensity = sum_spectra.loc[:, 'intensity'].values

    #wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)

    #spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)

    #DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Rasalhague, α Ophiuchi [m(V) = 2.07]', save = True, path = savefold, name = "Rasalhague_spectrum_lines.png")

    ## fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    ## PlotData(axs[0], spectra[0], None, 'Rasalhague raw data')
    ## PlotData(axs[0], background_av, None, 'background')
    ## PlotData(axs[0], dark_av, 'Rasalhague, α Ophiuchi [m(V) = 2.07]', 'dark')
    ## PlotData(axs[1], sum_spectra, None, 'Rasalhague processed data')
    ## plt.tight_layout()
    ## plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/alphaOph.png', dpi = 600)

    #fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    #PlotData(axs, sum_spectra, 'Rasalhague, α Ophiuchi [m(V) = 2.07]', 'Rasalhague processed data')
    #plt.savefig(f'D:\PhD\oceanview\Spectra_stars\\alphaOph_spectrum.png', dpi = 600)
    #plt.show()
    #
    # #=======================================================================================#
    # #                                  Alkaid                                             #
    # #=======================================================================================#
    #
    #spectra = LoadAllSpectra(path, 'Alkaid_spectrum_HDX017711')
    ####print(spectra)

    #background = LoadAllSpectra(path, 'Alkaid_background_HDX017711')

    #dark = LoadAllSpectra(path, 'dark_HDX017711')
    ###print(spectra[0])
    ####print(background[0])

    #sum_spectra_raw = AddSpectra(spectra)
    #background_av = AverageSpectra(background)
    ####print(background_av)
    #dark_av = AverageSpectra(dark)
    ####print(dark_av)
    #background_subtr = ProcessDataEach(background, None, dark_av, subtract_background = False, subtract_dark = True)
    #background_subtr_av = AverageSpectra(background_subtr)
    #subtr = ProcessDataEach(spectra, background_subtr_av, dark_av, subtract_background = True, subtract_dark = True)
    ####print(f'subtr: \n {subtr}')

    #sum_spectra = AddSpectra(subtr)

    #wavelength = sum_spectra.loc[:, 'wavelength'].values
    #intensity = sum_spectra.loc[:, 'intensity'].values

    #wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)

    #spec_norm = FitContinuum(15, wavelength_truncated, intensity_truncated, image = True)

    #DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'Alkaid, η Ursae Majoris [m(V) = 1.86]', save = True, path = savefold, name = "Alkaid_spectrum_lines.png")


    ### Identify absorption lines

    # absorption_lines = DetectLines(spec_norm, 0.1, False)
    #
    # nist_lines = NISTquery(wavelength_truncated)
    #
    # tolerance = 0.07 # nm
    # MatchElements(tolerance, absorption_lines, nist_lines, spec_norm)
    #
    # fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    # PlotData(axs[0], spectra[0], None, 'Alkaid raw data')
    # PlotData(axs[0], background_av, None, 'background')
    # PlotData(axs[0], dark_av, 'Alkaid, η Ursae Majoris [m(V) = 1.86]', 'dark')
    # PlotData(axs[1], sum_spectra, None, 'Alkaid processed data')
    # plt.tight_layout()
    # plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/Alkaid.png', dpi = 600)

    #fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    #PlotData(axs, sum_spectra, 'Alkaid, η Ursae Majoris [m(V) = 1.86]', 'Alkaid processed data')
    #plt.savefig(f'D:\PhD\oceanview\Spectra_stars\\Alkaid_spectrum.png', dpi = 600)
    #plt.show()
    #
    # #=======================================================================================#
    # #                                  delta Sge  (no background data)                      #
    # #=======================================================================================#
    #
    # spectra = LoadAllSpectra(path, 'deltaSag_spectrum_HDX017711')
    ##print(spectra)

    #background = LoadAllSpectra(path, 'deltaSag_background_HDX017711')

    # dark = LoadAllSpectra(path, 'dark_HDX017711')
    #print(spectra[0])
    ##print(background[0])

    # sum_spectra_raw = AddSpectra(spectra)
    #background_av = AverageSpectra(background)
    ##print(background_av)
    # dark_av = AverageSpectra(dark)
    ##print(dark_av)
    #background_subtr = ProcessDataEach(background, None, dark_av, subtract_background = False, subtract_dark = True)
    #background_subtr_av = AverageSpectra(background_subtr)
    # subtr = ProcessDataEach(spectra, None, dark_av, subtract_background = False, subtract_dark = True)
    ##print(f'subtr: \n {subtr}')

    # sum_spectra = AddSpectra(subtr)
    #
    # wavelength = sum_spectra.loc[:, 'wavelength'].values
    # intensity = sum_spectra.loc[:, 'intensity'].values
    #
    # wavelength_truncated, intensity_truncated = DefSpectrumForLines(sum_spectra, 360, 1100)
    #
    # spec_norm = FitContinuum(20, wavelength_truncated, intensity_truncated, image = True)
    #
    # DetectAndPlotLines(spec_norm, wavelength_truncated, spectral_lines, telluric_lines, threshold, title = 'δ Sagittae [m(V) = 4.21 (3.82)]', save = True, path = savefold, name = "DeltaSge_spectrum_lines.png")
    ##
    ## fig, axs = plt.subplots(2, 1, figsize=(8, 8))
    ## PlotData(axs[0], spectra[0], None, 'delta Sag raw data')
    ## PlotData(axs[0], background_av, None, 'background')
    ## PlotData(axs[0], dark_av, 'Kaus Media, δ Sagittarii [m(V) = 2.7]', 'dark')
    ## PlotData(axs[1], sum_spectra, None, 'delta Sag processed data (only dark subtracted)')
    ## plt.tight_layout()
    ## plt.savefig(f'/media/astronomer/Transcend/PhD/oceanview/Spectra_stars/deltaSag.png', dpi = 600)

    #fig, axs = plt.subplots(1, 1, figsize=(12, 8))
    #PlotData(axs, sum_spectra, 'δ Sagittae [m(V) = 4.21 (3.82)]', 'delta Sge processed data (only dark subtracted)')
    #plt.savefig(f'D:\PhD\oceanview\Spectra_stars\\delSge_spectrum.png', dpi = 600)
    plt.show()

