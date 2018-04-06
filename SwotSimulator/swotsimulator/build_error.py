# import params as p
import numpy
import swotsimulator.rw_data as rw_data
import swotsimulator.const as const
import swotsimulator.mod_tools as mod_tools
from math import pi, sqrt
import netCDF4
import sys
from scipy.ndimage.filters import gaussian_filter
import logging

# Define logging level for debug purposes
logger = logging.getLogger(__name__)


def reconstruct_2D_error(x_ac, err_out, dict_noise):
    nac = numpy.shape(x_ac)[0]
    if 'phase' in dict_noise.keys():
        ac_l = numpy.mat(x_ac[:int(nac / 2)])
        ac_r = numpy.mat(x_ac[int(nac / 2):])
        phase1d = dict_noise['phase']
        err_out.phase[:, :int(nac / 2)] = numpy.mat(phase1d[:, 0]).T * ac_l
        err_out.phase[:, int(nac / 2):] = numpy.mat(phase1d[:, 1]).T * ac_r
    if 'roll' in dict_noise.keys():
        ac = numpy.mat(x_ac)
        roll1d = dict_noise['roll']
        err_out.roll[:, :] = numpy.mat(roll1d).T * ac
    if 'baseline_dilation' in dict_noise.keys():
        ac2 = numpy.mat((x_ac)**2)
        baseline_dilation1d = numpy.mat(dict_noise['baseline_dilation'])
        err_out.baseline_dilation[:, :] = baseline_dilation1d.T * ac2
    if 'timing' in dict_noise.keys():
        timing1d = dict_noise['timing']
        ones_ac = numpy.mat(numpy.ones((int(nac/2))))
        err_out.timing[:, :int(nac / 2)] = numpy.mat(timing1d[:, 0]).T*ones_ac
        err_out.timing[:, int(nac / 2):] = numpy.mat(timing1d[:, 1]).T*ones_ac
    return None



class error():
    '''Class error define all the possible errors that can be computed using
    SWOT simulator.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in
    file file_coeff, the random realisations are read directly using
    load_coeff.  The corresponding errors on a swath can be computed using
    make_error. '''
    def __init__(self, p, roll=None, ssb=None, wet_tropo=None, phase=None,
                 baseline_dilation=None, karin=None, timing=None, SSH=None,
                 wt=None):
        self.roll = roll
        self.ssb = ssb
        self.wet_tropo = wet_tropo
        self.phase = phase
        self.baseline_dilation = baseline_dilation
        self.karin = karin
        self.timing = timing
        self.SSH_error = SSH
        self.wt = wt
        self.ncomp1d = getattr(p, 'ncomp1d', 2000)
        p.ncomp1d = self.ncomp1d
        self.ncomp2d = getattr(p, 'ncomp2d', 2000)
        p.ncomp2d = self.ncomp2d
        self.nrand = getattr(p, 'nrandkarin', 1000)
        p.nrandkarin = self.nrand

    def init_error(self, p, nac):
        '''Initialization of errors: Random realisation of errors are
        computed using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of
        each random realisation.
        By default, there are ncomp1d=2000 random realisations for
        the instrumental errors (1d spectrum) and ncomp2d=2000 random
        realisations for the geophysical errors (2d spectrum)
        and nrandkarin*x_ac km of random number for KaRIN noise.'''
        # # Run reprodictible: generate or load nrand random numbers:
        if p.karin is True:
            self.A_karin_l = numpy.random.normal(0.0, numpy.float64(1),
                                                 (self.nrand, nac))
            self.A_karin_r = numpy.random.normal(0.0, numpy.float64(1),
                                                 (self.nrand, nac))
        # - If one instrumental error needs to be computed, load frequencies
        #   and spectrum from the instrumental netcdf file:
        if p.phase or p.roll or p.baseline_dilation or p.timing:
            if p.lambda_cut is None:
                p.lambda_cut = 20000
            file_instr = rw_data.file_instr(file=p.file_inst_error)
            spatial_frequency = []
            file_instr.read_var(spatial_frequency=spatial_frequency)
            # - Cut frequencies larger than Nyquist frequency and cut long
            #   wavelengths (larger than p.lambda_max)
            cut_min = 1. / float(2 * p.delta_al)
            cut_max = 1. / p.lambda_max
            ind = numpy.where((file_instr.spatial_frequency < cut_min)
                              & (file_instr.spatial_frequency > 0)
                              & (file_instr.spatial_frequency > cut_max))[0]
            freq = file_instr.spatial_frequency[ind]
            f_cut = 1. / p.lambda_cut
            ind_cut = numpy.where(freq < f_cut)
            min_ind_cut = numpy.min(numpy.where(freq >= f_cut))
            if p.roll is True:
                # - Read and compute roll power spectrum, wavelength longer
                #   than p.lambda_cut are supposed to be corrected
                PSroll = []
                PSgyro = []
                file_instr.read_var(rollPSD=PSroll, gyroPSD=PSgyro)
                PSTroll = file_instr.rollPSD[ind] + file_instr.gyroPSD[ind]
                PSTroll[ind_cut] = PSTroll[min_ind_cut]
                # - Compute random coefficients using the power spectrum
                if p.savesignal is True:
                    gencoef = mod_tools.gen_rcoeff_signal1d(freq, PSTroll,
                                                            2 * p.delta_al,
                                                            p.lambda_max,
                                                            p.npseudoper,
                                                            p.len_repeat)
                    self.A_roll, self.phi_roll = gencoef
                else:
                    gencoef = mod_tools.gen_coeff_signal1d(freq, PSTroll,
                                                           self.ncomp1d)
                    self.A_roll, self.phi_roll, self.fr_roll = gencoef
            if p.phase is True:
                # - Read Phase power spectrum, wavelength longer than
                #   p.lambda_cut are supposed to be corrected
                phasePSD = []
                file_instr.read_var(phasePSD=phasePSD)
                PSphase = file_instr.phasePSD[ind]
                PSphase[ind_cut] = PSphase[min_ind_cut]
                # - Compute left and right random coefficients using
                #   the power spectrum
                if p.savesignal is True:
                    gencoef = mod_tools.gen_rcoeff_signal1d(freq, PSphase,
                                                            2 * p.delta_al,
                                                            p.lambda_max,
                                                            p.npseudoper,
                                                            p.len_repeat)
                    self.A_phase_l, self.phi_phase_l = gencoef
                    gencoef = mod_tools.gen_rcoeff_signal1d(freq, PSphase,
                                                            2 * p.delta_al,
                                                            p.lambda_max,
                                                            p.npseudoper,
                                                            p.len_repeat)
                    self.A_phase_r, self.phi_phase_r = gencoef
                else:
                    gencoef = mod_tools.gen_coeff_signal1d(freq, PSphase,
                                                           self.ncomp1d)
                    self.A_phase_l, self.phi_phase_l, self.fr_phase_l = gencoef
                    gencoef = mod_tools.gen_coeff_signal1d(freq, PSphase,
                                                           self.ncomp1d)
                    self.A_phase_r, self.phi_phase_r, self.fr_phase_r = gencoef
            if p.baseline_dilation is True:
                # - Read baseline dilation power spectrum, wavelength longer
                #   than p.lambda_cut are supposed to be corrected
                dilationPSD = []
                file_instr.read_var(dilationPSD=dilationPSD)
                PSbd = file_instr.dilationPSD[ind]
                PSbd[ind_cut] = PSbd[min_ind_cut]
                # - Compute random coefficients using the power spectrum
                if p.savesignal is True:
                    gencoef = mod_tools.gen_rcoeff_signal1d(freq, PSbd,
                                                            2 * p.delta_al,
                                                            p.lambda_max,
                                                            p.npseudoper,
                                                            p.len_repeat)
                    self.A_bd, self.phi_bd = gencoef
                else:
                    gencoef = mod_tools.gen_coeff_signal1d(freq, PSbd,
                                                           self.ncomp1d)
                    self.A_bd, self.phi_bd, self.fr_bd = gencoef
            if p.timing is True:
                # - Read timing power spectrum, wavelength longer than
                #   p.lambda_cut are supposed to be corrected
                timingPSD = []
                file_instr.read_var(timingPSD=timingPSD)
                PStim = file_instr.timingPSD[ind]
                PStim[ind_cut] = PStim[min_ind_cut]
                # - Compute random coefficients using the power spectrum
                if p.savesignal is True:
                    gencoef = mod_tools.gen_rcoeff_signal1d(freq, PStim,
                                                            2 * p.delta_al,
                                                            p.lambda_max,
                                                            p.npseudoper,
                                                            p.len_repeat)
                    self.A_tim_l, self.phi_tim_l = gencoef
                    gencoef = mod_tools.gen_rcoeff_signal1d(freq, PStim,
                                                            2 * p.delta_al,
                                                            p.lambda_max,
                                                            p.npseudoper,
                                                            p.len_repeat)
                    self.A_tim_r, self.phi_tim_r = gencoef
                else:
                    gencoef = mod_tools.gen_coeff_signal1d(freq, PStim,
                                                           self.ncomp1d)
                    self.A_tim_l, self.phi_tim_l, self.fr_tim_l = gencoef
                    gencoef = mod_tools.gen_coeff_signal1d(freq, PStim,
                                                           self.ncomp1d)
                    self.A_tim_r, self.phi_tim_r, self.fr_tim_r = gencoef
            if p.wet_tropo is True:
                # - Define power spectrum of error in path delay
                #   due to wet tropo
                f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
                # - Global mean wet tropo power spectrum in
                #   cm**2/(cycle/km) for L >= 100 km
                PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
                # - Wet tropo power spectrum in cm**2/(cycle/km)
                #   for L < 100 km
                indf = numpy.where(f > 10**-2)
                PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
                # - Compute random coefficients in 2D using the previously
                #   defined power spectrum
                gencoef = mod_tools.gen_coeff_signal2d(f, PSwt, self.ncomp2d)
                self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gencoef
                # - Define radiometer error power spectrum for a beam
                #   High frequencies are cut to filter the associated error:
                #   during the reconstruction of the wet trop signal
                # f=numpy.arange(1./3000.,1./float(20.),1./3000.)
                PSradio = 9.5 * 10**-5 * f**-1.79
                PSradio[numpy.where((f < 1./1000.))] = \
                    9.5 * 10**-5 * (10**-3)**-1.79
                indf = numpy.where((f > 0.0023) & (f <= 0.0683))
                PSradio[indf] = 0.036 * f[indf]**-0.814
                PSradio[numpy.where(f > 0.0683)] = 0.32
                # - Compute random coefficients (1D) for the radiometer error
                #   power spectrum for right and left beams
                gencoef = mod_tools.gen_coeff_signal1d(f, PSradio,
                                                       self.ncomp2d)
                self.A_radio_r, self.phi_radio_r, self.fr_radio_r = gencoef
                self.A_radio_l, self.phi_radio_l, self.fr_radio_l = gencoef

    def load_coeff(self, p):
        '''Load existing random realisations that has been stored in file_coeff.
        The outputs are the amplitude, the phase and the frequency of
        each random realisation.
        There are ncomp1d random realisations for instrumental and ncomp2d
        for geophysical errors and nrandkarin *x_ac random number for
        KaRin noise. '''
        try:
            fid = netCDF4.Dataset(p.file_coeff, 'r')
        except IOError:
            logger.error('There was an error opening file '
                         '{}'.format(p.file_coeff))
            sys.exit(1)
        if p.karin is True:
            _A_karin = fid.variables['A_karin_l'][:, :]
            self.A_karin_l = numpy.array(_A_karin).squeeze()
            _A_karin = fid.variables['A_karin_r'][:, :]
            self.A_karin_r = numpy.array(_A_karin).squeeze()
            if numpy.shape(self.A_karin_l)[0] != self.nrand:
                logger.error('{1} dimensions are different from nrandkarin={2}'
                             '\n remove {1} or adjust nrandkarin number in '
                             'parameter file'.format(p.file_coeff, self.nrand))
                sys.exit(1)
        if p.roll is True:
            self.A_roll = numpy.array(fid.variables['A_roll'][:]).squeeze()
            if numpy.shape(self.A_roll)[0] != self.ncomp1d:
                logger.error('{1} dimensions are different from ncomp1d={2}\n'
                             'remove {1} or adjust ncomp1d number in parameter'
                             'file'.format(p.file_coeff, self.ncomp1d))
                sys.exit(1)
            self.phi_roll = numpy.array(fid.variables['phi_roll'][:]).squeeze()
            self.fr_roll = numpy.array(fid.variables['fr_roll'][:]).squeeze()
        if p.phase is True:
            _tmpphase = fid.variables['A_phase_l'][:]
            self.A_phase_l = numpy.array(_tmpphase).squeeze()
            _tmpphase = fid.variables['phi_phase_l'][:]
            self.phi_phase_l = numpy.array(_tmpphase).squeeze()
            _tmpphase = fid.variables['fr_phase_l'][:]
            self.fr_phase_l = numpy.array(_tmpphase).squeeze()
            _tmpphase = fid.variables['A_phase_r'][:]
            self.A_phase_r = numpy.array(_tmpphase).squeeze()
            _tmpphase = fid.variables['phi_phase_r'][:]
            self.phi_phase_r = numpy.array(_tmpphase).squeeze()
            _tmpphase = fid.variables['fr_phase_r'][:]
            self.fr_phase_r = numpy.array(_tmpphase).squeeze()
        if p.baseline_dilation is True:
            self.A_bd = numpy.array(fid.variables['A_bd'][:]).squeeze()
            self.phi_bd = numpy.array(fid.variables['phi_bd'][:]).squeeze()
            self.fr_bd = numpy.array(fid.variables['fr_bd'][:]).squeeze()
        if p.timing is True:
            self.A_tim_l = numpy.array(fid.variables['A_tim_l'][:]).squeeze()
            _tmptim = fid.variables['phi_tim_l'][:]
            self.phi_tim_l = numpy.array(_tmptim).squeeze()
            self.fr_tim_l = numpy.array(fid.variables['fr_tim_l'][:]).squeeze()
            self.A_tim_r = numpy.array(fid.variables['A_tim_r'][:]).squeeze()
            _tmptim = fid.variables['phi_tim_r'][:]
            self.phi_tim_r = numpy.array(_tmptim).squeeze()
            self.fr_tim_r = numpy.array(fid.variables['fr_tim_r'][:]).squeeze()
        if p.wet_tropo is True:
            self.A_wt = numpy.array(fid.variables['A_wt'][:]).squeeze()
            if numpy.shape(self.A_wt)[0] != self.ncomp2d:
                logger.error('{1} dimensions are different from ncomp1d={2}\n'
                             'remove {1} or adjust ncomp1d number in parameter'
                             'file'.format(p.file_coeff, self.ncomp2d))
                sys.exit(1)
            self.phi_wt = numpy.array(fid.variables['phi_wt'][:]).squeeze()
            self.frx_wt = numpy.array(fid.variables['frx_wt'][:]).squeeze()
            self.fry_wt = numpy.array(fid.variables['fry_wt'][:]).squeeze()
            _tmpradio = fid.variables['A_radio_r'][:]
            self.A_radio_r = numpy.array(_tmpradio).squeeze()
            _tmpradio = fid.variables['phi_radio_r'][:]
            self.phi_radio_r = numpy.array(_tmpradio).squeeze()
            _tmpradio = fid.variables['fr_radio_r'][:]
            self.fr_radio_r = numpy.array(_tmpradio).squeeze()
            _tmpradio = fid.variables['A_radio_l'][:]
            self.A_radio_l = numpy.array(_tmpradio).squeeze()
            _tmpradio = fid.variables['phi_radio_l'][:]
            self.phi_radio_l = numpy.array(_tmpradio).squeeze()
            _tmpradio = fid.variables['fr_radio_l'][:]
            self.fr_radio_l = numpy.array(_tmpradio).squeeze()

        fid.close()
        return None

    def make_error(self, sgrid, cycle, SSH_true, p):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, the phase between the two signals,
        the timing error, the roll of the satellite, the sea surface bias,
        the distorsion of the mast,
        the karin noise due to the sensor itself. '''
        # - Errors array are the same size as the swath size
        nal, nac = numpy.shape(SSH_true)
        # ind_al=numpy.arange(0,nal)
        if p.karin is True:
            # - Formula of karin noise as a function of x_ac (smile shape)
            # - Load Karin noise from file:
            karin_data = rw_data.file_karin(file=p.karin_file)
            karin_data.read_karin(p.swh)
            sigma_karin = numpy.interp(numpy.abs(sgrid.x_ac),
                                       karin_data.x_ac, karin_data.hsdt)
            size_grid = sqrt(numpy.float64(p.delta_al * p.delta_ac))
            sigma_karin = sigma_karin / size_grid
            # - Compute random karin error
            Ai = (((numpy.float64(sgrid.x_al) + float(cycle * sgrid.al_cycle))
                  / p.delta_al) % self.nrand).astype('int')
            for j in range(0, nac):
                self.karin[:, j] = (sigma_karin[j]) * self.A_karin_r[Ai, j]

        if p.phase is True:
            # - Compute left and right phase angles using random
            #   coefficients or signals previously initialized
            if p.savesignal is True:
                xx = (numpy.float64(sgrid.x_al[:]) + float(cycle
                      * sgrid.al_cycle)) % (p.len_repeat)
                theta_l = mod_tools.gen_signal1d(xx, self.A_phase_l,
                                                 self.phi_phase_l,
                                                 2 * p.delta_al,
                                                 p.lambda_max, p.npseudoper)
                theta_r = mod_tools.gen_signal1d(xx, self.A_phase_r,
                                                 self.phi_phase_r,
                                                 2 * p.delta_al,
                                                 p.lambda_max, p.npseudoper)
            else:
                theta_l = numpy.zeros((nal))
                theta_r = numpy.zeros((nal))
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2.*pi*float(self.fr_phase_l[comp])
                                  * (numpy.float64(sgrid.x_al[:]) +
                                  + float(cycle * sgrid.al_cycle))) % (2.*pi)
                    theta_l[:] = (theta_l[:] + 2*self.A_phase_l[comp]
                                  * numpy.cos(phase_x_al[:]
                                  + self.phi_phase_l[comp]))
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2. * pi * float(self.fr_phase_r[comp])
                                  * (numpy.float64(sgrid.x_al[:])
                                  + float(cycle*sgrid.al_cycle))) % (2.*pi)
                    theta_r[:] = (theta_r[:] + 2 * self.A_phase_r[comp]
                                  * numpy.cos(phase_x_al[:]
                                  + self.phi_phase_r[comp]))
            # - Compute the associated phase error on the swath in m
            _to_km = (1 / (const.Fka * 2*pi / const.C * const.B)
                       * (1 + const.sat_elev / const.Rearth) *pi/180. * 10**3)
            phase1d = numpy.concatenate(([_to_km * theta_l[:].T],
                                        [_to_km * theta_r[:].T]),
                                             axis=0)
            self.phase1d = phase1d.T
            del theta_l, theta_r
        if p.roll is True:
            # - Compute roll angle using random coefficients or signals
            # previously initialized with the power spectrum
            if p.savesignal is True:
                xx = (numpy.float64(sgrid.x_al[:]) + float(cycle
                      * sgrid.al_cycle)) % (p.len_repeat)
                theta = mod_tools.gen_signal1d(xx, self.A_roll, self.phi_roll,
                                               2 * p.delta_al, p.lambda_max,
                                               p.npseudoper)
            else:
                theta = numpy.zeros((nal))
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2. * pi * float(self.fr_roll[comp])
                                  * (numpy.float64(sgrid.x_al[:])
                                  + float(cycle * sgrid.al_cycle))) % (2.*pi)
                    theta[:] = (theta[:] + 2 * self.A_roll[comp]
                                * numpy.cos(phase_x_al[:]
                                + self.phi_roll[comp]))
            # - Compute the associated roll  error on the swath in m
            self.roll1d = ((1 + const.sat_elev/const.Rearth)
                          * theta[:]*pi/180./3600. * 10**3)
            del theta
        if p.baseline_dilation is True:
            # - Compute baseline dilation using random coefficients or signals
            # previously initialized with the power spectrum
            if p.savesignal is True:
                xx = (numpy.float64(sgrid.x_al[:]) + float(cycle
                      * sgrid.al_cycle)) % (p.len_repeat)
                dil = mod_tools.gen_signal1d(xx, self.A_bd, self.phi_bd,
                                             2 * p.delta_al, p.lambda_max,
                                             p.npseudoper)
            else:
                dil = numpy.zeros((nal))
                # self.baseline_dilation=numpy.empty((nal,nac))
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2. * pi * float(self.fr_bd[comp])
                                  * (numpy.float64(sgrid.x_al[:])
                                  + float(cycle*sgrid.al_cycle))) % (2.*pi)
                    dil[:] = (dil[:] + 2 * self.A_bd[comp]
                              * numpy.cos(phase_x_al[:]
                              + self.phi_bd[comp]))
            # - Compute the associated baseline dilation  error
            #   on the swath in m
            _to_km = ((1 + const.sat_elev / const.Rearth) # * 10**-6 * 10**6
                      /(const.sat_elev * const.B))
            self.baseline_dilation1d = _to_km * dil[:]
            del dil
        if p.timing is True:
            # - Compute timing delay using random coefficients previously
            #   initialized with the power spectrum
            if p.savesignal is True:
                xx = (numpy.float64(sgrid.x_al[:]) + float(cycle
                      * sgrid.al_cycle)) % (p.len_repeat)
                tim_l = mod_tools.gen_signal1d(xx, self.A_tim_l,
                                               self.phi_tim_l,
                                               2 * p.delta_al, p.lambda_max,
                                               p.npseudoper)
                tim_r = mod_tools.gen_signal1d(xx, self.A_tim_r,
                                               self.phi_tim_r,
                                               2 * p.delta_al, p.lambda_max,
                                               p.npseudoper)
            else:
                tim_l = numpy.zeros((nal))
                tim_r = numpy.zeros((nal))
            # - Compute the associated phase error on the swath in m
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2.  * pi * float(self.fr_tim_r[comp])
                                  * (numpy.float64(sgrid.x_al[:])
                                  + float(cycle*sgrid.al_cycle))) % (2.*pi)
                    tim_r[:] = (tim_r[:] + 2 * self.A_tim_r[comp]
                                * numpy.cos(phase_x_al[:]
                                + self.phi_tim_r[comp]))
                for comp in range(0, self.ncomp1d):
                    phase_x_al = (2. * pi * float(self.fr_tim_l[comp])
                                  * (numpy.float64(sgrid.x_al[:])
                                  + float(cycle*sgrid.al_cycle))) % (2.*pi)
                    tim_l[:] = (tim_l[:] + 2 * self.A_tim_l[comp]
                                * numpy.cos(phase_x_al[:]
                                + self.phi_tim_l[comp]))
            # - Compute the corresponding timing error on the swath in m
            _to_m = const.C/2 * 10**(-12)
            timing1d = numpy.concatenate(([_to_m * tim_l[:].T],
                                         [_to_m * tim_r[:].T]), axis=0)
            self.timing1d = timing1d.T
        if p.wet_tropo is True:
            # - Initialization of radiometer error in right and left beam
            err_radio_r = numpy.zeros((nal))
            err_radio_l = numpy.zeros((nal))
            # - Initialization of swath matrices and large swath matrices
            #   (which include wet tropo data around the nadir and outside
            #   the swath)
            #   x_ac_large and wt_large are necessary to compute the gaussian
            #   footprint of a beam on the nadir or near the edge of the swath
            x, y = numpy.meshgrid(sgrid.x_ac, numpy.float64(sgrid.x_al[:])
                                  + float(cycle * sgrid.al_cycle))
            start_x = -2. * p.sigma / float(p.delta_ac) + sgrid.x_ac[0]
            stop_x = (2. * p.sigma / float(p.delta_ac) + sgrid.x_ac[-1]
                      + p.delta_ac)
            x_ac_large = numpy.arange(start_x, stop_x, p.delta_ac)
            naclarge = numpy.shape(x_ac_large)[0]
            wt_large = numpy.zeros((nal, naclarge))
            y_al = numpy.float64(sgrid.x_al[:]) + float(cycle*sgrid.al_cycle)
            x_large, y_large = numpy.meshgrid(x_ac_large, y_al)
            # - Compute path delay error due to wet tropo and radiometer error
            #   using random coefficient initialized with power spectrums
            for comp in range(0, self.ncomp2d):
                phase_x_al = (2. * pi * (float(self.frx_wt[comp])
                              * (numpy.float64(x)) + float(self.fry_wt[comp])
                              * numpy.float64(y))) % (2.*pi)
                self.wt = (self.wt + self.A_wt[comp]
                           * numpy.cos(phase_x_al + self.phi_wt[comp])*10**-2)
                phase_x_al_large = (2. * pi * (float(self.frx_wt[comp])
                                    * (numpy.float64(x_large))
                                    + float(self.fry_wt[comp])
                                    * numpy.float64(y_large))) % (2.*pi)
                wt_large = (wt_large + self.A_wt[comp]
                            * numpy.cos(phase_x_al_large
                            + self.phi_wt[comp])*10**-2)
                phase_x_al = (2. * pi * float(self.fr_radio_r[comp])
                              * (numpy.float64(sgrid.x_al[:])
                              + float(cycle*sgrid.al_cycle))) % (2.*pi)
                err_radio_r = (err_radio_r + 2*self.A_radio_r[comp]
                               * numpy.cos(phase_x_al
                               + self.phi_radio_r[comp])*10**-2)
                phase_x_al = (2. * pi * float(self.fr_radio_l[comp])
                              * (numpy.float64(sgrid.x_al[:])
                              + float(cycle*sgrid.al_cycle))) % (2.*pi)
                err_radio_l = (err_radio_l + 2*self.A_radio_l[comp]
                               * numpy.cos(phase_x_al
                               + self.phi_radio_l[comp])*10**-2)
            # - Compute Residual path delay error after a 1-beam radiometer
            #   correction
            if p.nbeam == 1 or p.nbeam == 'both':
                beam = numpy.zeros((nal))
                diff_h1 = numpy.zeros((nal, nac))
                # indac = numpy.where((sgrid.x_ac<2.*p.sigma)
                #                   & (sgrid.x_ac>-2.*p.sigma))[0]
                # - Find across track indices in the gaussian footprint of
                #   2.*p.sigma
                indac = numpy.where((x_ac_large < 2.*p.sigma)
                                    & (x_ac_large > -2.*p.sigma))[0]
                for i in range(0, nal):
                    # - Find along track indices in the gaussian footprint of
                    #   2.*p.sigma
                    delta_x_al = sgrid.x_al[:] - sgrid.x_al[i]
                    indal = numpy.where((delta_x_al <= (2 * p.sigma))
                                        & (delta_x_al > (-2*p.sigma)))[0]
                    # indal=numpy.where(((sgrid.x_al[:]-sgrid.x_al[i])<=(2*
                    # p.sigma))&((sgrid.x_al[:]-sgrid.x_al[i])>(-2*p.sigma)))[0]
                    slice_al = slice(min(indal), max(indal) + 1)
                    slice_ac = slice(min(indac), max(indac) + 1)
                    x, y = numpy.meshgrid(x_ac_large[slice_al],
                                          sgrid.x_al[slice_ac] - sgrid.x_al[i])
                    # x,y = numpy.meshgrid(sgrid.x_al[indal],sgrid.x_ac[indac])
                    # - Compute path delay on gaussian footprint
                    G = 1./(2.*pi*p.sigma**2) * numpy.exp(-(x**2. + y**2.)
                                                          / (2.*p.sigma**2))
                    beam[i] = (sum(sum(G * wt_large[slice_al, slice_ac]))
                               / sum(sum(G)) + err_radio_l[i])
                # - Filtering beam signal to cut frequencies higher than 125 km
                beam = gaussian_filter(beam, 30. / p.delta_al)
                # - Compute residual path delay
                for i in range(0, nal):
                    diff_h1[i, :] = self.wt[i, :] - beam[i]
                self.wet_tropo1[:, :] = diff_h1[:, :]  # en 2d
                if p.nadir is True:
                    self.wet_tropo1nadir = (wt_large[:, int(naclarge/2.)]
                                            - beam[:])
                # self.wet_tropo1[:, (nac-1)/2]=0 #numpy.nan
            # - Compute Residual path delay error after a 2-beams radiometer
            #   correction
            if p.nbeam == 2 or p.nbeam == 'both':
                beam_r = numpy.zeros((nal))
                beam_l = numpy.zeros((nal))
                diff_h2 = numpy.zeros((nal, nac))
                diff_h2nadir = numpy.zeros((nal))
                # - Find righ and leftacross track indices in the gaussian
                #   footprint of 2.*p.sigma
                ind_r = x_ac_large + p.beam_pos_r
                indac_r = numpy.where((ind_r < 2.*p.sigma)
                                      & (ind_r > -2.*p.sigma))[0]
                ind_l = x_ac_large + p.beam_pos_l
                indac_l = numpy.where((ind_l < 2.*p.sigma)
                                      & (ind_l > -2.*p.sigma))[0]
                for i in range(0, nal):
                    # - Find along track indices in the gaussian footprint
                    #   of 2.*p.sigma
                    delta_al = sgrid.x_al[:] - sgrid.x_al[i]
                    indal = numpy.where((delta_al <= (2*p.sigma))
                                        & (delta_al > (-2*p.sigma)))[0]
                    # indal=numpy.where(((sgrid.x_al[:]-sgrid.x_al[i])<=(2*
                    slice_al = slice(min(indal), max(indal) + 1)
                    slice_acr = slice(min(indac_r), max(indac_r) + 1)
                    slice_acl = slice(min(indac_l), max(indac_l) + 1)
                    x, y = numpy.meshgrid(x_ac_large[slice_acr],
                                          sgrid.x_al[slice_al] - sgrid.x_al[i])
                    # - Compute path delay on left and right gaussian footprint
                    G = 1. / (2.*pi*p.sigma**2) * numpy.exp(-(x**2. + y**2.)
                                                            / (2.*p.sigma**2))
                    beam_r[i] = (sum(sum(G*wt_large[slice_al, slice_acr]))
                                 / sum(sum(G))+err_radio_r[i])
                    beam_l[i] = (sum(sum(G*wt_large[slice_al, slice_acl]))
                                 / sum(sum(G)) + err_radio_l[i])
                # - Filtering beam signal to cut frequencies higher than 125 km
                beam_r = gaussian_filter(beam_r, 30. / p.delta_al)
                beam_l = gaussian_filter(beam_l, 30. / p.delta_al)
                for i in range(0, nal):
                    # - Compute residual path delay (linear combination of left
                    #   and right path delay)
                    pol = numpy.polyfit([p.beam_pos_l, p.beam_pos_r],
                                        [beam_l[i], beam_r[i]], 1)
                    pbeam = numpy.polyval(pol, -sgrid.x_ac)
                    diff_h2[i, :] = self.wt[i, :] - pbeam[:]
                    diff_h2nadir[i] = (wt_large[i, int(naclarge/2.)]
                                       - pbeam[int(nac/2.)])
                self.wet_tropo2[:, :] = diff_h2[:, :]  # en 2d
                if p.nadir:
                    self.wet_tropo2nadir = diff_h2nadir  # en 1d
                # self.wet_tropo2[:, (nac-1)/2]=0 #numpy.nan
            if p.nadir is True:
                self.wtnadir = wt_large[:, int(naclarge/2.)]
            if (not p.nbeam == 1) and (not p.nbeam == 2) \
               and (not p.nbeam == 'both'):
                logger.error("\n nbeam = {} \n".format(p.nbeam))
                logger.error("wrong number of beam, nbeam should be either 1"
                             "or 2 or 'both'")
                sys.exit(1)
        if p.ssb is True:
            logger.error("No SSB error implemented yet")
        return None

    def reconstruct_2D(self, p, x_ac):
        dict_noise = {}
        if p.phase is True:
            dict_noise['phase'] = self.phase1d
        if p.roll is True:
            dict_noise['roll'] = self.roll1d
        if p.baseline_dilation is True:
            dict_noise['baseline_dilation'] = self.baseline_dilation1d
        if p.timing is True:
            dict_noise['timing'] = self.timing1d
        reconstruct_2D_error(x_ac, self, dict_noise)
        return None

    def make_SSH_error(self, SSH_true, p):
        '''Compute observed SSH adding all the computed error to the model SSH.
        If residual path delay errors after 2-beams and 1-beam radiometer
        correction are both computed (nbeam='both'), only the path delay error
        after 1-beam radiometer correction is considered in the
        observed SSH. '''
        numpy.seterr(invalid='ignore')
        self.SSH = SSH_true
        if p.karin is True:
            self.SSH = self.SSH + self.karin
        if p.timing is True:
            self.SSH = self.SSH + self.timing
        if p.roll is True:
            self.SSH = self.SSH + self.roll
        if p.baseline_dilation is True:
            self.SSH = self.SSH + self.baseline_dilation
        if p.phase is True:
            self.SSH = self.SSH + self.phase
        if p.wet_tropo is True:
            if p.nbeam == 1 or p.nbeam == 'both':
                self.SSH = self.SSH + self.wet_tropo1
            else:
                self.SSH = self.SSH + self.wet_tropo2
        if p.file_input is not None:
            self.SSH[numpy.where(SSH_true == p.model_nan)] = p.model_nan

    def save_coeff(self, p, nac):
        '''Save random realisations to enable runs to be reproducible.
        The ncomp1d random phase phi, amplitude A and frequency fr in 1D and
        ncomp2d random phase phi, amplitude A and frequencies frx and fry in 2D
        are saved in file_coeff for each error (except KaRIN) and can be
        loaded using load_coeff.
        Random numbers on a grid nrandkarin km long and x_ac large are stored
        for KaRIN noise.
        '''
        # - Open Netcdf file in write mode
        fid = netCDF4.Dataset(p.file_coeff, 'w')
        fid.description = "Random coefficients computed for SWOT simulator"

        # - Create dimensions
        fid.createDimension('nrand1d', self.ncomp1d)
        fid.createDimension('nrand2d', self.ncomp2d)
        fid.createDimension('nkarin', self.nrand)
        fid.createDimension('x_ac', nac)
        # - Create and write Variables
        if p.karin is True:
            var = fid.createVariable('A_karin_l', 'f4', ('nkarin', 'x_ac'))
            var[:, :] = self.A_karin_l
            var = fid.createVariable('A_karin_r', 'f4', ('nkarin', 'x_ac'))
            var[:, :] = self.A_karin_r
        if p.roll is True:
            var = fid.createVariable('A_roll', 'f4', ('nrand1d', ))
            var[:] = self.A_roll
            var = fid.createVariable('phi_roll', 'f4', ('nrand1d', ))
            var[:] = self.phi_roll
            var = fid.createVariable('fr_roll', 'f4', ('nrand1d', ))
            var[:] = self.fr_roll
        if p.phase is True:
            var = fid.createVariable('A_phase_r', 'f4', ('nrand1d', ))
            var[:] = self.A_phase_r
            var = fid.createVariable('phi_phase_r', 'f4', ('nrand1d', ))
            var[:] = self.phi_phase_r
            var = fid.createVariable('fr_phase_r', 'f4', ('nrand1d', ))
            var[:] = self.fr_phase_r
            var = fid.createVariable('A_phase_l', 'f4', ('nrand1d', ))
            var[:] = self.A_phase_l
            var = fid.createVariable('phi_phase_l', 'f4', ('nrand1d', ))
            var[:] = self.phi_phase_l
            var = fid.createVariable('fr_phase_l', 'f4', ('nrand1d', ))
            var[:] = self.fr_phase_l
        if p.baseline_dilation is True:
            var = fid.createVariable('A_bd', 'f4', ('nrand1d', ))
            var[:] = self.A_bd
            var = fid.createVariable('phi_bd', 'f4', ('nrand1d', ))
            var[:] = self.phi_bd
            var = fid.createVariable('fr_bd', 'f4', ('nrand1d', ))
            var[:] = self.fr_bd
        if p.timing is True:
            var = fid.createVariable('A_tim_l', 'f4', ('nrand1d', ))
            var[:] = self.A_tim_l
            var = fid.createVariable('phi_tim_l', 'f4', ('nrand1d', ))
            var[:] = self.phi_tim_l
            var = fid.createVariable('fr_tim_l', 'f4', ('nrand1d', ))
            var[:] = self.fr_tim_l
            var = fid.createVariable('A_tim_r', 'f4', ('nrand1d', ))
            var[:] = self.A_tim_r
            var = fid.createVariable('phi_tim_r', 'f4', ('nrand1d', ))
            var[:] = self.phi_tim_r
            var = fid.createVariable('fr_tim_r', 'f4', ('nrand1d', ))
            var[:] = self.fr_tim_r
        if p.wet_tropo is True:
            var = fid.createVariable('A_wt', 'f4', ('nrand2d', ))
            var[:] = self.A_wt
            var = fid.createVariable('phi_wt', 'f4', ('nrand2d', ))
            var[:] = self.phi_wt
            var = fid.createVariable('frx_wt', 'f4', ('nrand2d', ))
            var[:] = self.frx_wt
            var = fid.createVariable('fry_wt', 'f4', ('nrand2d', ))
            var[:] = self.fry_wt
            var = fid.createVariable('A_radio_r', 'f4', ('nrand2d', ))
            var[:] = self.A_radio_r
            var = fid.createVariable('phi_radio_r', 'f4', ('nrand2d', ))
            var[:] = self.phi_radio_r
            var = fid.createVariable('fr_radio_r', 'f4', ('nrand2d', ))
            var[:] = self.fr_radio_r
            var = fid.createVariable('A_radio_l', 'f4', ('nrand2d', ))
            var[:] = self.A_radio_l
            var = fid.createVariable('phi_radio_l', 'f4', ('nrand2d', ))
            var[:] = self.phi_radio_l
            var = fid.createVariable('fr_radio_l', 'f4', ('nrand2d', ))
            var[:] = self.fr_radio_l
        fid.close()
        return None


class errornadir():
    '''Class errornadir defines the error on the nadir.
    Random realisation of errors can be initialized using init_error.
    If the random realisations have already been computed and stored in file
    file_coeff, the random realisations are read directly using load_coeff.
    The correspondg errors on a swath can be computed using make_error. '''
    def __init__(self, p, nadir=None, wet_tropo1=None, wt=None):
        self.nadir = nadir
        self.wet_tropo1 = wet_tropo1
        self.wt = wt
        self.nrand = getattr(p, 'nrandkarin', 1000)
        p.nrandkarin = self.nrand
        self.ncomp2d = getattr(p, 'ncomp2d', 2000)
        p.ncomp2d = self.ncomp2d
        self.ncomp1d = getattr(p, 'ncomp1d', 2000)
        p.ncomp1d = self.ncomp1d

    def init_error(self, p):
        '''Initialization of errors: Random realisation of errors are computed
        using a known power spectrum.
        The outputs are the amplitude, the phase and the frequency of each
        random realisation.
        By default, there are ncomp2d=2000 random realisations for the
        wet tropo and ncomp1d=2000 random realisations for the nadir 1d
        spectrum error.'''
        # Run reprodictible: generate or load nrand random numbers:
        # - Compute random coefficients in 1D for the nadir error
        wnoise = getattr(p, 'wnoise', 100)
        p.wnoise = wnoise
        # - Define the sepctrum of the nadir instrument error
        # self.A=numpy.random.normal(0.0,sqrt(p.wnoise)
        # /numpy.float64(sqrt(2*p.delta_al)), (self.nrand))*0.01
        f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
        PSD = 8 + 1.05 * 10**(-4) * f**(-2.2)
        indf = numpy.where(f < 0.00023627939582672978)
        PSD[indf] = 10**4
        # Convert spectrum in m2/cy
        PSD = PSD * 10**(-4)
        if p.savesignal is True:
            gencoef = mod_tools.gen_rcoeff_signal1d(f, PSD, 2 * p.delta_al,
                                                    p.lambda_max, p.npseudoper,
                                                    p.len_repeat)
            self.A, self.phi = gencoef
        else:
            gencoef = mod_tools.gen_coeff_signal1d(f, PSD, self.ncomp1d)
            self.A, self.phi, self.f = gencoef
        if p.wet_tropo is True:
            # - Define power spectrum of error in path delay due to wet tropo
            f = numpy.arange(1./3000., 1./float(2.*p.delta_al), 1./3000.)
            # - Global mean wet tropo power spectrum in cm**2/(cycle/km)
            #   for L >= 100 km
            PSwt = 3.156 * 10**-5 * f**(-8./3.)  # *10**-4
            # - Wet tropo power spectrum in cm**2/(cycle/km) for L < 100 km
            indf = numpy.where(f > 10**-2)
            PSwt[indf] = 1.4875 * 10**-4 * f[indf]**(-2.33)  # *10**-4
            # - Compute random coefficients in 2D using the previously defined
            #   power spectrum
            gencoef = mod_tools.gen_coeff_signal2d(f, PSwt, self.ncomp2d)
            self.A_wt, self.phi_wt, self.frx_wt, self.fry_wt = gencoef
            # - Define radiometer error power spectrum for a beam
            #   High frequencies are cut to filter the associated error during
            #   the reconstruction of the wet trop signal
            # f=numpy.arange(1./3000.,1./float(20.),1./3000.)
            PSradio = 9.5 * 10**-5 * f**-1.79
            PSradio[numpy.where((f < 1./1000.))] = 9.5 * 10**-5*(10**-3)**-1.79
            indf = numpy.where((f > 0.0023) & (f <= 0.0683))
            PSradio[indf] = 0.036 * f[indf]**-0.814
            PSradio[numpy.where(f > 0.0683)] = 0.32
            # - Compute random coefficients (1D) for the radiometer error power
            #   spectrum for right and left beams
            gencoef = mod_tools.gen_coeff_signal1d(f, PSradio, self.ncomp2d)
            self.A_radio, self.phi_radio, self.fr_radio = gencoef
        return None

    def load_coeff(self, p):
        '''Load existing random realisations that has been stored in
        nadir+file_coeff. The outputs are the amplitude,
        the phase and the frequency of each random realisation.
        There are ncomp random realisations.'''
        try:
            fid = netCDF4.Dataset('{}_nadir.nc'.format(p.file_coeff[:-3]), 'r')
        except IOError:
            logger.error('There was an error opening the file nadir '
                         '{}_nadir.nc'.format(p.file_coeff[:-3]))
            sys.exit(1)
        self.A = numpy.array(fid.variables['A'][:]).squeeze()
        self.f = numpy.array(fid.variables['f'][:]).squeeze()
        self.phi = numpy.array(fid.variables['phi'][:]).squeeze()
        if p.wet_tropo is True:
            self.A_wt = numpy.array(fid.variables['A_wt'][:]).squeeze()
            if numpy.shape(self.A_wt)[0] != self.ncomp2d:
                logger.error('{1} dimensions are different from ncomp2d={2}\n'
                             'remove {1} or adjust ncomp2d number in parameter'
                             'file'.format(p.file_coeff, self.ncomp2d))
                sys.exit(1)
            self.phi_wt = numpy.array(fid.variables['phi_wt'][:]).squeeze()
            self.frx_wt = numpy.array(fid.variables['frx_wt'][:]).squeeze()
            self.fry_wt = numpy.array(fid.variables['fry_wt'][:]).squeeze()
            self.A_radio = numpy.array(fid.variables['A_radio'][:]).squeeze()
            _tmpradio = fid.variables['phi_radio'][:]
            self.phi_radio = numpy.array(_tmpradio).squeeze()
            _tmpradio = fid.variables['fr_radio'][:]
            self.fr_radio = numpy.array(_tmpradio).squeeze()
        fid.close()
        return None

    def make_error(self, orb, cycle, SSH_true, p):
        ''' Build errors corresponding to each selected noise
        among the effect of the wet_tropo, and the instrumental error
        '''
        nal = numpy.shape(SSH_true)[0]
        # - Compute random noise of 10**2 cm**2/(km/cycle)
        # - Compute the correspond error on the nadir in m
        if p.savesignal is True:
            xx = (numpy.float64(orb.x_al[:]) + float(cycle
                  * orb.al_cycle)) % p.len_repeat
            errnadir = mod_tools.gen_signal1d(xx, self.A, self.phi,
                                              2 * p.delta_al, p.lambda_max,
                                              p.npseudoper)
        else:
            errnadir = numpy.zeros((nal))
            for comp in range(0, self.ncomp1d):
                phase_x_al = (2. * pi * float(self.f[comp])
                              * (numpy.float64(orb.x_al[:])
                              + float(cycle*orb.al_cycle))) % (2.*pi)
                errnadir[:] = (errnadir[:] + 2*self.A[comp]
                               * numpy.cos(phase_x_al[:]+self.phi[comp]))
        # - Compute the correspond timing error on the swath in m
        self.nadir = errnadir[:]
        if p.wet_tropo:
            # - Initialization of radiometer error in right and left beam
            err_radio = numpy.zeros((nal))
            # - Initialization of swath matrices and large swath matrices
            #   (which include wet tropo data around the nadir and outside the
            #   swath)
            #   x_ac_large and wt_large are necessary to compute the gaussian
            #   footprint of a beam on the nadir or near the edge of the swath
            x_ac_large = numpy.arange(-2. * p.sigma/float(p.delta_ac),
                                      2.*p.sigma/float(p.delta_ac)+p.delta_ac,
                                      p.delta_ac)
            wt_large = numpy.zeros((numpy.shape(orb.x_al[:])[0],
                                   numpy.shape(x_ac_large)[0]))
            y_al = numpy.float64(orb.x_al[:]) + float(cycle * orb.al_cycle)
            x_large, y_large = numpy.meshgrid(x_ac_large, y_al)
            # - Compute path delay error due to wet tropo and radiometer error
            #   using random coefficient initialized with power spectrums
            for comp in range(0, self.ncomp2d):
                phase_x_al_large = (2. * pi * (float(self.frx_wt[comp])
                                    * (numpy.float64(x_large))
                                    + float(self.fry_wt[comp])
                                    * numpy.float64(y_large))) % (2.*pi)
                wt_large = (wt_large + self.A_wt[comp]
                            * numpy.cos(phase_x_al_large
                            + self.phi_wt[comp])*10**-2)
                phase_x_al = (2. * pi * float(self.fr_radio[comp])
                              * (numpy.float64(orb.x_al[:])
                              + float(cycle*orb.al_cycle))) % (2.*pi)
                err_radio = (err_radio + 2*self.A_radio[comp]
                             * numpy.cos(phase_x_al + self.phi_radio[comp])
                             * 10**-2)
            # - Compute Residual path delay error after a 1-beam radiometer
            #   correction
            beam = numpy.zeros((nal))
            diff_h1 = numpy.zeros((nal))
            # indac=numpy.where((sgrid.x_ac<2.*p.sigma)
            # & (sgrid.x_ac>-2.*p.sigma))[0]
            # - Find across track indices in the gaussian footprint of
            #   2.*p.sigma
            indac = numpy.where((x_ac_large < 2.*p.sigma)
                                & (x_ac_large > -2.*p.sigma))[0]
            for i in range(0, nal):
                # - Find along track indices in the gaussian footprint of
                #   2.*p.sigma
                delta_al = orb.x_al[:] - orb.x_al[i]
                indal = numpy.where((delta_al <= (2*p.sigma))
                                    & (delta_al > -2*p.sigma))[0]
                # indal=numpy.where(((sgrid.x_al[:]-sgrid.x_al[i])<=
                slice_ac = slice(min(indac), max(indac) + 1)
                slice_al = slice(min(indal), max(indal) + 1)
                x, y = numpy.meshgrid(x_ac_large[slice_ac],
                                      orb.x_al[slice_al] - orb.x_al[i])
                # - Compute path delay on gaussian footprint
                G = 1. / (2.*pi*p.sigma**2) * numpy.exp(-(x**2.+y**2.)
                                                        / (2.*p.sigma**2))
                beam[i] = (sum(sum(G*wt_large[slice_al, slice_ac]))
                           / sum(sum(G)) + err_radio[i])
            # - Filtering beam signal to cut frequencies higher than 125 km
            beam = gaussian_filter(beam, 30./p.delta_al)
            # - Compute residual path delay
            diff_h1 = wt_large[:, int(numpy.shape(wt_large)[1]/2.)] - beam
            self.wet_tropo1 = diff_h1
            self.wt = wt_large[:, int(numpy.shape(wt_large)[1]/2.)]
        return None

    def save_coeff(self, p):
        '''Save random realisations to enable runs to be reproducible.
        The ncomp1d random phase phi, amplitude A and frequency fr for
        1D spectrum and ncomp2d random phase phi, amplitude A and frequencies
        frx and fry for 2D spectrum are saved in nadirfile_coeff for each error
        and can be loaded using load_coeff.
        '''
        # - Open Netcdf file in write mode
        fid = netCDF4.Dataset('{}_nadir.nc'.format(p.file_coeff[:-3]), 'w')
        fid.description = "Random coefficients from orbit simulator"

        # - Create dimensions
        fid.createDimension('nrand1d', self.ncomp1d)
        fid.createDimension('nrand2d', self.ncomp2d)
        # - Create and write Variables
        var = fid.createVariable('A', 'f4', ('nrand1d', ))
        var[:] = self.A
        var = fid.createVariable('f', 'f4', ('nrand1d', ))
        var[:] = self.f
        var = fid.createVariable('phi', 'f4', ('nrand1d', ))
        var[:] = self.phi

        # var = fid.createVariable('phi', 'f4', ('ninstr',))
        # var[:] = self.phi
        # var = fid.createVariable('fr', 'f4', ('ninstr',))
        # var[:] = self.fr
        if p.wet_tropo is True:
            var = fid.createVariable('A_wt', 'f4', ('nrand2d', ))
            var[:] = self.A_wt
            var = fid.createVariable('phi_wt', 'f4', ('nrand2d', ))
            var[:] = self.phi_wt
            var = fid.createVariable('frx_wt', 'f4', ('nrand2d', ))
            var[:] = self.frx_wt
            var = fid.createVariable('fry_wt', 'f4', ('nrand2d', ))
            var[:] = self.fry_wt
            var = fid.createVariable('A_radio', 'f4', ('nrand2d', ))
            var[:] = self.A_radio
            var = fid.createVariable('phi_radio', 'f4', ('nrand2d', ))
            var[:] = self.phi_radio
            var = fid.createVariable('fr_radio', 'f4', ('nrand2d', ))
            var[:] = self.fr_radio
        fid.close()
        return None
