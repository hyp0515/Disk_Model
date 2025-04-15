import os
import numpy as np
from radmc3dPy.image import *
from radmc3dPy.analyze import *
from radmc3dPy.data import *

class generate_simulation:
    
    def __init__(self, parms,
                 channel=False,
                 pv=False,
                 conti=False,
                 sed=False,
                 line_spectrum=False
                 ):
        
        self.save_out = getattr(parms, 'save_out', True)
        self.save_npz = getattr(parms, 'save_npz', True)
        if (self.save_out is False) and (self.save_npz is False):
            self.save_out = True
        
        if channel is True:
            try:
                self.generate_cube(parms, parms.channel_cube_parms)
                print('Image cube is generated')
            except:
                print('Wrong parameters to generate image cube')
        
        if pv is True:
            try:
                self.generate_cube(parms, parms.pv_cube_parms)
                print('Image cube is generated')
            except:
                print('Wrong parameters to generate image cube')
        
        if sed is True:
            try:
                self.generate_sed(parms.sed_parms)
                print('SED is generated')
                
            except:
                print('Wrong parameters to generate SED')
                
        if line_spectrum is True:
            try:
                self.generate_line_spectrum(parms, parms.spectrum_parms)
                print('Line spectrum is generated')
            except:
                print('Wrong parameters to generate line spectra')
        
        if conti is True:
            try:
                self.generate_continuum(parms.conti_parms)
                print('Dust continuum is generated')
            except:
                print('Wrong parameters to generate dust continuum')

    def generate_cube(self, parms, cube_parms, **kwargs):
                    
        """
        incl               : Inclination angle of the disk
        line               : Transistion level (see 'molecule_ch3oh.inp')
        v_width            : Range of velocity to simulate
        nlam               : Number of velocities
        npix               : Number of map's pixels
        nodust             : If False, dust effect is included
        scat               : If True and nodust=False, scattering is included. (Time-consuming)
        extracted_gas      : If True, spectral line is extracted (I_{dust+gas}-I_{dust})
        sizeau             : Map's span
        """
        
        condition = {}
        for c in ['nodust', 'scat', 'extract_gas']:
            if c in kwargs:
                condition[c] = kwargs[c]
            else:
                condition[c] = getattr(parms.condition_parms, c, False)
        nodust = condition['nodust']
        scat = condition['scat']
        extract_gas = condition['extract_gas']
        
        npix        = getattr(cube_parms,        'npix',   100)
        sizeau      = getattr(cube_parms,      'sizeau',   100)
        incl        = getattr(cube_parms,        'incl',    73)
        line        = getattr(cube_parms,        'line',   240)
        v_width     = getattr(cube_parms,     'v_width',    10)
        nlam        = getattr(cube_parms,        'nlam',    10)
        vkms        = getattr(cube_parms,        'vkms',     0)
        posang      = getattr(cube_parms,      'posang',    45)

        prompt = f'npix {npix} sizeau {sizeau} incl {incl} posang {-posang} iline {line} vkms {vkms} widthkms {v_width} linenlam {nlam}'
        
        if nlam > 15:
            type_note = 'pv'
        else:
            type_note = 'channel'
        
        
        if extract_gas is False:
            if nodust is True:
                prompt = prompt + ' noscat nodust'
                f      = '_nodust'
            elif nodust is False:
                if scat is True:
                    f      = '_scat'
                elif scat is False:
                    prompt = prompt +' noscat'
                    f      = '_noscat'
                else:
                    pass
            else:
                pass
                    
            os.system(f"radmc3d image "+prompt)
            

            self.cube = readImage('image.out')
            if self.save_npz is True:
                if 'fname' in kwargs.keys():
                    self.save_npzfile(self.cube, cube_parms, fname=kwargs['fname'], f=f, note=type_note)
                else:
                    self.save_npzfile(self.cube, cube_parms, f=f, note=type_note)

            
            if self.save_out is True:
                if 'fname' in kwargs.keys():
                    self.save_outfile(cube_parms, fname=kwargs['fname'], f=f, note=type_note)
                else:
                    self.save_outfile(cube_parms, f=f, note=type_note)
            else:
                pass
                    
        elif extract_gas is True:
            os.system(f"radmc3d image "+prompt)

            self.cube = readImage('image.out')
            if self.save_npz is True:
                if 'fname' in kwargs.keys():
                    self.save_npzfile(self.cube, cube_parms, fname=kwargs['fname'], f='_scat', note=type_note)
                else:
                    self.save_npzfile(self.cube, cube_parms, f='_scat', note=type_note)
            else:
                pass
            
            if self.save_out is True:
                if 'fname' in kwargs.keys():
                    self.save_outfile(cube_parms, fname=kwargs['fname'], f='_scat', note=type_note)
                else:
                    self.save_outfile(cube_parms, f='_scat', note=type_note)
            else:
                pass

            if self.cube.nwav%2 == 0:
                mid_wav = 0.5 * (self.cube.wav[self.cube.nwav//2] + self.cube.wav[(self.cube.nwav//2)+1])
            else:
                mid_wav = self.cube.wav[(self.cube.nwav//2)+1]\
            
            os.system(f"radmc3d image npix {npix} sizeau {sizeau} incl {incl} posang {-posang} lambda {mid_wav} noline")

            self.conti = readImage('image.out')
            if self.save_npz is True:
                if 'fname' in kwargs.keys():
                    self.save_npzfile(self.conti, cube_parms, fname=kwargs['fname'], f='_conti', note=type_note)
                else:
                    self.save_npzfile(self.conti, cube_parms, f='_conti', note=type_note)
            else:
                pass
            
            if self.save_out is True:
                if 'fname' in kwargs.keys():
                    self.save_outfile(cube_parms, fname=kwargs['fname'], f='_conti', note=type_note)
                else:
                    self.save_outfile(cube_parms, f='_conti', note=type_note)
            else:
                pass


            self.cube_list = [self.cube, self.conti]
            if self.save_npz is True:
                if 'fname' in kwargs.keys():
                    self.save_npzfile(self.cube_list, cube_parms, fname=kwargs['fname'], f='_extracted', note=type_note)
                else:
                    self.save_npzfile(self.cube_list, cube_parms, f='_extracted', note=type_note)
            else:
                pass
            
        else:
            pass
            
    def generate_continuum(self, conti_parms, **kwargs):
        
        type_note = 'conti'
        
        incl   = getattr(conti_parms,   'incl',   70)
        wav    = getattr(conti_parms,    'wav', 1300)
        npix   = getattr(conti_parms,   'npix',  200)
        sizeau = getattr(conti_parms, 'sizeau',  100)
        posang = getattr(conti_parms, 'posang',  45)

        prompt = f'radmc3d image npix {npix} sizeau {sizeau} incl {incl} lambda {wav} posang {-posang} noline'
        
        f = '_scat'
        
        if conti_parms.scat is False:
            prompt = prompt + ' noscat'
            f = '_noscat'
        os.system(prompt)
        
        
        self.conti = readImage('image.out')
        if self.save_npz is True:
            if 'fname' in kwargs.keys():
                self.save_npzfile(self.conti, conti_parms, fname=kwargs['fname'], f=f, note=type_note)
            else:
                self.save_npzfile(self.conti, conti_parms, f=f, note=type_note)
        
        if self.save_out is True:
            if 'fname' in kwargs.keys():
                self.save_outfile(conti_parms, fname=kwargs['fname'], f=f, note=type_note)
            else:
                self.save_outfile(conti_parms, f=f, note=type_note)
        else:
            pass
 
    def generate_sed(self, sed_parms, **kwargs):
        
        type_note = 'sed'
        
        incl     = getattr(sed_parms,     'incl',     70)
        freq_min = getattr(sed_parms, 'freq_min',    5e1)
        freq_max = getattr(sed_parms, 'freq_max',    1e3)
        nlam     = getattr(sed_parms,     'nlam',    100)
        wav_max  = ((cc*1e-2)/(freq_min*1e9))*1e+6
        wav_min  = ((cc*1e-2)/(freq_max*1e9))*1e+6
        
        
        prompt = f"radmc3d spectrum incl {incl} lambdarange {wav_min} {wav_max} nlam {nlam} noline"
        f = '_scat'
        if sed_parms.scat is False:
            prompt = prompt + ' noscat'
            f = '_noscat'
        os.system(prompt)
        if (sed_parms.read_sed is True) or (kwargs['read_sed'] is True):
            self.spectrum = readSpectrum('spectrum.out')
            if self.save_npz is True:
                if 'fname' in kwargs.keys():
                    self.save_npzfile(self.spectrum, sed_parms, fname=kwargs['fname'], f=f, note=type_note)
                else:
                    self.save_npzfile(self.spectrum, sed_parms, f=f, note=type_note)
        else:
            pass
        
        if self.save_out is True:
            if 'fname' in kwargs.keys():
                self.save_outfile(sed_parms, fname=kwargs['fname'], f=f, note=type_note)
            else:
                self.save_outfile(sed_parms, f=f, note=type_note)
        else:
            pass

    def generate_line_spectrum(self, parms, spectrum_parms, **kwargs):
        """
        incl               : Inclination angle of the disk
        line               : Transistion level (see 'molecule_ch3oh.inp')
        v_width            : Range of velocity to simulate
        nlam               : Number of velocities
        npix               : Number of map's pixels
        nodust             : If False, dust effect is included
        scat               : If True and nodust=False, scattering is included. (Time-consuming)
        extracted_gas      : If True, spectral line is extracted (I_{dust+gas}-I_{dust})
        sizeau             : Map's span
        """
        
        type_note = 'spectrum'
        
        
        condition = {}
        for c in ['nodust', 'scat', 'extract_gas']:
            if c in kwargs:
                condition[c] = kwargs[c]
            else:
                condition[c] = getattr(parms.condition_parms, c, False)
        nodust = condition['nodust']
        scat = condition['scat']
        extract_gas = condition['extract_gas']
        
        incl        = getattr(spectrum_parms,        'incl',    70)
        line        = getattr(spectrum_parms,        'line',   240)
        v_width     = getattr(spectrum_parms,     'v_width',    10)
        nlam        = getattr(spectrum_parms,        'nlam',    10)
        vkms        = getattr(spectrum_parms,        'vkms',     0)
        
        prompt = f"radmc3d spectrum incl {incl} iline {line} vkms {vkms} widthkms {v_width} linenlam {nlam}"
        
        
        if extract_gas is False:
            if nodust is True:
                prompt = prompt + ' noscat nodust'
                f      = '_nodust'
            elif nodust is False:
                if scat is True:
                    f      = '_scat'
                elif scat is False:
                    prompt = prompt +' noscat'
                    f      = '_noscat'
                else:
                    pass
            else:
                pass
                    
            os.system(prompt)
            
            if (spectrum_parms.read_spectrum is True) or (kwargs['read_spectrum'] is True):
                self.spectrum = readSpectrum('spectrum.out')
                if self.save_npz is True:
                    if 'fname' in kwargs.keys():
                        self.save_npzfile(self.spectrum, spectrum_parms, fname=kwargs['fname'], f=f, note=type_note)
                    else:
                        self.save_npzfile(self.spectrum, spectrum_parms, f=f, note=type_note)
            else:
                pass
            
            if self.save_out is True:
                if 'fname' in kwargs.keys():
                    self.save_outfile(spectrum_parms, fname=kwargs['fname'], f=f, note=type_note)
                else:
                    self.save_outfile(spectrum_parms, f=f, note=type_note)
            else:
                pass
                    
        elif extract_gas is True:
            os.system(prompt)
            if spectrum_parms.read_spectrum is True:
                self.spectrum = readSpectrum('spectrum.out')
                if self.save_npz is True:
                    if 'fname' in kwargs.keys():
                        self.save_npzfile(self.spectrum, spectrum_parms, fname=kwargs['fname'], f='_scat', note=type_note)
                    else:
                        self.save_npzfile(self.spectrum, spectrum_parms, f='_scat', note=type_note)
                else:
                    pass
                
                if self.save_out is True:
                    if 'fname' in kwargs.keys():
                        self.save_outfile(spectrum_parms, fname=kwargs['fname'], f='_scat', note=type_note)
                    else:
                        self.save_outfile(spectrum_parms, f='_scat', note=type_note)
                else:
                    pass
            else:
                pass
            
            os.system(f"radmc3d spectrum incl {incl} lambdarange {self.spectrum[0,0]} {self.spectrum[-1,0]} nlam {nlam} noline")
            
            if spectrum_parms.read_spectrum is True:
                self.conti = readSpectrum('spectrum.out')
                if self.save_npz is True:
                    if 'fname' in kwargs.keys():
                        self.save_npzfile(self.conti, spectrum_parms, fname=kwargs['fname'], f='_conti', note=type_note)
                    else:
                        self.save_npzfile(self.conti, spectrum_parms, f='_conti', note=type_note)
                else:
                    pass
                
                if self.save_out is True:
                    if 'fname' in kwargs.keys():
                        self.save_outfile(spectrum_parms, fname=kwargs['fname'], f='_conti', note=type_note)
                    else:
                        self.save_outfile(spectrum_parms, f='_conti', note=type_note)
                else:
                    pass
            else:
                pass
            
            if spectrum_parms.read_spectrum is True:
                
                self.spectrum_list = [self.spectrum, self.conti]
                if self.save_npz is True:
                    if 'fname' in kwargs.keys():
                        self.save_npzfile(self.spectrum_list, spectrum_parms, fname=kwargs['fname'], f='_extracted', note=type_note)
                    else:
                        self.save_npzfile(self.spectrum_list, spectrum_parms, f='_extracted', note=type_note)
                else:
                    pass
            else:
                pass
        else:
            pass

    def save_outfile(self, parms, **kwargs):
        """
        This will save whole information of simulation.
        The file type will be '*.out', which the storage may be MB-scale
        """

        dir = getattr(parms, 'dir', './test/')
        
        if 'fname' in kwargs.keys():
            fname = kwargs['fname']
        else:
            fname = getattr(parms, 'fname', 'test')
        
        if 'f' in kwargs.keys():
            f = kwargs['f']
        else:
            f = ''
            
        if 'note' in kwargs.keys():
            head = kwargs['note'] + '_'
        else:
            head = ''
        
        os.makedirs(dir+'outfile/', exist_ok=True)
        if (kwargs['note'] == 'sed') or (kwargs['note'] == 'spectrum'):
            os.system('mv spectrum.out '+dir+'outfile/'+head+fname+f+'.out')
        else:
            os.system('mv image.out '+dir+'outfile/'+head+fname+f+'.out')

    def save_npzfile(self, data, parms, **kwargs):
        """
        This will only save image data regardless of other information
        The file tyep will be '*.npz', which the storage may be kB-scale
        """
        
        dir   = getattr(parms, 'dir', './test/')
        
        if 'fname' in kwargs.keys():
            fname = kwargs['fname']
        else:
            fname = getattr(parms, 'fname', 'test')
        
        
        if 'f' in kwargs.keys():
            f = kwargs['f']
        else:
            f = ''
        
        if 'note' in kwargs.keys():
            head = kwargs['note'] + '_'
        else:
            head = ''
        
        os.makedirs(dir+'npzfile/', exist_ok=True)
        
        if (kwargs['note'] == 'sed') or (kwargs['note'] == 'spectrum'):
            np.savez(dir+'npzfile/'+head+fname+f+'.npz',
                    lam = data[:, 0],
                    fnu = data[:, 1])
        else:
            if isinstance(data, list):
                line_cube = data[0]
                dust_conti_value = data[1].imageJyppix
                dust_conti_value = np.tile(dust_conti_value[:, :, np.newaxis], (1, 1, line_cube.nwav))

                np.savez(dir+'npzfile/'+head+fname+f+'.npz',
                        imageJyppix = line_cube.imageJyppix - dust_conti_value,
                        x           = line_cube.x,
                        nx          = line_cube.nx,
                        y           = line_cube.y,
                        ny          = line_cube.ny,
                        wav         = line_cube.wav,
                        freq        = line_cube.freq,
                        nwav        = len(line_cube.wav),
                        nfreq       = len(line_cube.freq),)
            else:
                np.savez(dir+'npzfile/'+head+fname+f+'.npz',
                        imageJyppix = data.imageJyppix,
                        x           = data.x,
                        nx          = data.nx,
                        y           = data.y,
                        ny          = data.ny,
                        wav         = data.wav,
                        freq        = data.freq,
                        nwav        = len(data.wav),
                        nfreq       = len(data.freq),)
        
