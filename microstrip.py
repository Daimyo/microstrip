# -*- coding: utf-8 -*-

#This class is based on an article of Frank Schnieder and Wolfgang Heinrich
#"Model of thin-Film Microstrip Line for Circuit Design"
# IEEE Transactions on Microwave Theory And Techniques, vol 49, n° 1, January 2001

#Copyright (C) 2013 Dumur Étienne

#This program is free software; you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation; either version 2 of the License, or
#(at your option) any later version.

#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#GNU General Public License for more details.

#You should have received a copy of the GNU General Public License along
#with this program; if not, write to the Free Software Foundation, Inc.,
#51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

import numpy as np
import scipy.constants as cst

class microstrip:
    
    def __init__(self, epsilon_r = 11.68, tan_delta = 7e-4, kappa = 3.53e50, w = 5e-6, h = 500e-6, t = 25e-9, w_g = 5e-3):
        '''Class allowing the calculation of the parameters of a microstrip line.
            
            Input:
                - epsilon_r (float) : Relative permitivity of the substrat in farad per meter.
                - tan_delta (float) : Loss tangent.
                - kappa     (float) : Conductivity of the metal layer in Siemens per meter.
                
                - w         (float) : Width of the central line in meter.
                - h         (float) : Height of the substrat separation in meter.
                - t         (float) : Thickness of the metal layers in meter.
                - w_g       (float) : Width of the ground plane in meter.
        '''
        
        self._epsilon_r = epsilon_r
        self._tan_delta = tan_delta
        self._kappa     = kappa
        self._w   = w
        self._h   = h
        self._t   = t
        self._w_g = w_g
        
        self._lambda = 16e-9



    #################################################################################
    #
    #
    #                                    set methods
    #
    #
    #################################################################################



    def set_width(self, w):
        self._w = float(w)



    def set_height(self, h):
        self._h = float(h)



    def set_thickness(self, t):
        self._t = float(t)



    def set_relative_permitivity(self, epsilon_r):
        self._epsilon_r = float(epsilon_r)



    def set_electrical_conductivity(self, kappa):
        self._kappa = float(kappa)



    def set_tan_delta(self, tanDelta):
        self._tan_delta = float(tanDelta)


    #################################################################################
    #
    #
    #                                    Usefull methods
    #
    #
    #################################################################################

    def _omega(self, f):
        return 2.*cst.pi*f

    def _cotanh(self, a):
        return np.cosh(a)/np.sinh(a)

    #################################################################################
    #
    #
    #                                    Methods
    #
    #
    #################################################################################

    def _F1(self):
        return 6. + (2.*np.pi - 6.)*np.exp(-(30.666*self._h/self._w)**0.7528)

    def _eta0(self):
        return np.sqrt(cst.mu_0/cst.epsilon_0)

    def _ZL0(self, w):
        return self._eta0()*np.log(self._F1()*self._h/w + np.sqrt(1. + 4.*self._h**2./w**2.))/2./np.pi

    def _a(self):
        return 1. +1./49.*np.log((self._w**4./self._h**4. + self._w**2./self._h**2./2704.)/(self._w**4./self._h**4. + 0.432)) + np.log(1. + self._w**3./self._h**3./5929.741)/18.7

    def _b(self):
        return 0.564*((self._epsilon_r - 0.9)/(self._epsilon_r + 3.))**0.053

    def _epsilon_r_eff_0(self, w):
        return (self._epsilon_r + 1.)/2. +  (self._epsilon_r - 1.)*(1. + 10.*self._h/w)**(-self._a()*self._b())/2.

    def _w_eq_0(self):
        temp = np.sqrt(6.517*self._w/self._h)
        temp = np.cosh(temp)/np.sinh(temp)
        
        return self._w + self._t*np.log(1. + 4.*np.exp(1.)/(self._t*temp**2./self._h ))/np.pi

    def _w_eq_Z(self):
        return self._w + (self._w_eq_0() - self._w)*(1. + 1./np.cosh(np.sqrt(self._epsilon_r - 1.)))/2.

    def _ZL(self):
        return self._ZL0(self._w_eq_Z())/np.sqrt(self._epsilon_r_eff_0(self._w_eq_Z()))

    def _epsilon_r_eff(self):
        return self._epsilon_r_eff_0(self._w_eq_Z())*self._ZL0(self._w_eq_0())**2./self._ZL0(self._w_eq_Z())**2.

    def _Ca_prim(self):
        return 1./cst.c/self._ZL0(self._w_eq_0())

    def _C_epsilon_prim(self):
        return ((self._epsilon_r_eff() - 1.)/(self._epsilon_r - 1.))*self._epsilon_r*self._Ca_prim()

    def _A(self):
        return 1. + self._h*(1. + 1.25*np.log(2.*self._h/self._t)/np.pi)/self._w_eq_0()

    def _R_S(self, f):
        return np.sqrt(np.pi*f*cst.mu_0/self._kappa)

    def _alpha_c(self, f):
        if self._w/self._h <= 1. :
            return 0.1589*self._A()*self._R_S(f)*(32. - self._w_eq_0()**2./self._h**2.)/self._h/self._ZL()/(32. + self._w_eq_0()**2./self._h**2.)
        else:
            return 7.0229e-6*self._A()*self._R_S(f)*self._ZL()*self._epsilon_r_eff()*(self._w_eq_0()/self._h + 0.667*self._w_eq_0()/self._h/(self._w_eq_0()/self._h + 1.444))/self._h

    def _R_se_prim(self, f):
        return 2.*self._ZL()*self._alpha_c(f)

    def _R_prim_se(self, f):
        return 2.*self.Z(self)*self._alpha_c(f)

    def _L_a_prim(self):
        return 1./cst.c**2./self._Ca_prim()

    def _L_i_prim(self, f):
        return self._R_se_prim(f)/self._omega(f)

    def _L_se_prim(self, f):
        return self._L_a_prim() + self._L_i_prim()

    def _R_DC_prim_w(self):
        return 1./self._kappa/self._w/self._t

    def _R_DC_prim_g(self):
        return 1./self._kappa/self._w_g*self._t

    def _R_DC_prim(self):
        return self._R_DC_prim_w() + self._R_DC_prim_g()

    def _K4(self, z):
        return z**4.*(np.log(z) - 25./12.)/24.

    def _K_s(self, a, b):
        return (4.*(self._K4(a) +self._K4(1j*b)) - 2.*(self._K4(a + 1j*b) + self._K4(a - 1j*b))).real + np.pi*a*b**3./3.

    def _K_m(self, a, b, c, d, h):
        temp1 = [-a/2.   , a/2.]
        temp2 = [-1j*b/2., 1j*b/2.]
        temp3 = [-c/2.   , c/2.]
        temp4 = [-1j*d/2., 1j*d/2.]
        
        end = 0.
        for z1 in temp1:
            for z2 in temp2:
                for z3 in temp3:
                    for z4 in temp4:
                        end += self._K4(z1 + z2 - z3 - z4 - 1j*h)
        
        return -end.real

    def _L_DC_prim(self):
        return -cst.mu_0*(self._K_s(self._w, self._t)/self._w**2. - 2./self._w/self._w_g*self._K_m(self._w, self._t, self._w_g, self._t, self._h + self._t) + self._K_s(self._w_g, self._t)/self._w_g**2.)/2./np.pi/self._t**2.

    def _f_0(self):
        return 2.*self._R_DC_prim_w()*self._R_DC_prim_g()/cst.mu_0/(self._R_DC_prim_w() + self._R_DC_prim_g())

    def _f_se(self):
        return (1.6 + (10.*self._t/self._w)/(1. + self._w/self._h))/cst.pi/cst.mu_0/self._kappa/self._t**2.

    #################################################################################
    #
    #
    #                                    Dynamic correction
    #
    #
    #################################################################################

    def _P4(self):
        return 1. + 2.751 * (1. - np.exp(-(self._epsilon_r/15.916)**8.))

    def _P3(self, f):
        return 0.0363*np.exp(-4.6*self._w/self._h)*(1. - np.exp(-(f*self._h/3.87*1e9*1e-2)**4.97))

    def _P2(self):
        return 0.33622*(1. - np.exp(-0.03442*self._epsilon_r))

    def _P1(self, f):
        return 0.27488 + self._w/self._h*(0.6315 + 0.525/(1. + 0.157 * f * self._h/1e9/1e-2)**20.) - 0.065683*np.exp(-8.7513*self._w/self._h)

    def _P(self, f):
        return self._P1(f)*self._P2()*((0.1844 + self._P3(f)*self._P4())*10.*f*self._h/1e9/1e-2)**1.5763

    def _F_epsilon(self, f):
        return self._epsilon_r/self._epsilon_r_eff_0(self._w) - ((self._epsilon_r/self._epsilon_r_eff_0(self._w) - 1.))/(1. + self._P(f))

    def _F_Z(self, f):
        return (self._epsilon_r_eff_0(self._w)*self._F_epsilon(f) - 1.)/(self._epsilon_r_eff_0(self._w) - 1.)/np.sqrt(self._F_epsilon(f))

    def _F_L(self, f):
        return self._F_Z(f)*np.sqrt(self._F_epsilon(f))

    def _F_C(self, f):
        return np.sqrt(self._F_epsilon(f))/self._F_Z(f)

    #################################################################################
    #
    #
    #                                    Superconducting correction (only for the inductance)
    #
    #
    #################################################################################

#    def _beta(self):
#        return 1. + self._t/self._h

#    def _p(self):
#        return 2.*self._beta()**2. - 1. + np.sqrt((2.*self._beta()**2. - 1.)**2. - 1.)

#    def _eta(self):
#        return np.pi*self._w/2./self._h*np.sqrt(self._p()) + (self._p() + 1.)/2.*(1. + np.log(4./(np.sqrt(self._p()) + 1.))) -np.sqrt(self._p())*np.log(np.sqrt(self._p()) + 1.) - (np.sqrt(self._p()) - 1.)**2./2.*np.log(np.sqrt(self._p()) - 1.)
##        return np.sqrt(self._p())*(np.pi*self._w/2./self._h + (self._p() + 1.)*(1. + np.log( 4./(self._p() - 1.)))/2./np.sqrt(self._p()) - 2.*np.arctanh(np.sqrt(self._p())))

#    def _r_bo(self):
#        if self._eta() >= self._p() :
#            delta = self._eta()
#        else:
#            delta = self._p()
#        return self._eta() + (self._p() + 1.)*np.log(delta)/2.

#    def _r_b(self):
#        if self._w/self._h >= 5.:
#            return self._r_bo()
#        else:
#            return self._r_bo() + np.sqrt((self._r_bo() - 1)*(self._r_bo() - self._p())) + (self._p() + 1.)*np.arctanh(np.sqrt((self._r_bo() - self._p())/(self._r_bo() - 1.))) - 2.*np.sqrt(self._p())*np.arctanh(np.sqrt((self._r_bo() - self._p())/(self._p() * (self._r_bo() - 1.)))) + np.pi*self._w*np.sqrt(self._p())/2./self._h 

#    def _r_a(self):
#        return np.exp(-1. - np.pi*self._w/2./self._h + np.log(4.*self._p()) - (np.sqrt(self._p()) + 1.)**2./2./np.sqrt(self._p())*np.log(np.sqrt(self._p()) + 1) + (np.sqrt(self._p()) - 1.)**2./2./np.sqrt(self._p())*np.log(np.sqrt(self._p()) - 1.))
##        return np.exp(-1. - np.pi*self._w/2./self._h - (self._p() + 1.)*np.arctanh(np.sqrt(self._p()))/np.sqrt(self._p()) - np.log((self._p() - 1.)/(4.*self._p())))

#    def _K(self, w, h, t):
#        return h*2.*np.log(2.*self._r_b()/self._r_a())/w/np.pi

#    def L(self):
#        
#        return cst.mu_0*(self._lambda*(self._cotanh(self._t/self._lambda) + 2.*np.sqrt(self._p())/np.sinh(self._t/self._lambda)/self._r_b()) +self._lambda*self._cotanh(self._t/self._lambda) )/self._w/self._K(self._w, self._h, self._t)

    #################################################################################
    #
    #
    #                                    Final result
    #
    #
    #################################################################################

    def get_conductance_per_unit_length(self, f):
        '''Return the length conductance of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Conductance per unit length in Siemens per meter
        '''
        
        return self._omega(f) * self._C_epsilon_prim()*self._tan_delta

    def get_capacitance_per_unit_length(self, f):
        '''Return the length capacitance of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Capacitance per unit length in Farrad per meter
        '''
        
        return (self._epsilon_r_eff() * self._Ca_prim())*self._F_C(f)

    def get_inductance_per_unit_length(self, f):
        '''Return the inductance per unit length of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Inductance per unit length in Henrys per meter
        '''
        
        return (self._L_a_prim() + self._L_i_prim(self._f_se())/(1. + np.sqrt(f/self._f_se())) + (self._L_DC_prim() - self._L_a_prim() + self._L_i_prim(self._f_se()))/np.sqrt(1. + f**2./self._f_0()**2.))*self._F_L(f)

#    def get_superconducting_inductance_per_unit_length(self):
#        '''Return the superconducting inductance per unit length of the microstrip line
#                - Input :
#                    - Frequency (float | list | numpy.ndarray) in Hertz
#                
#                - Output :
#                    - Inductance per unit length in Henrys per meter
#        '''
#        
##        return cst.mu_0*2.*self._lambda**2./self._t/self._w
#        pass
#        return cst.mu_0/

    def get_resistance_per_unit_length(self, f):
        '''Return the length resistance of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Resistance per unit length in Ohms per meter
        '''
        
        A = self._R_se_prim(self._f_se())*(np.sqrt(f/self._f_se()) + np.sqrt(1. + f**2./self._f_se()**2.))/(1. + np.sqrt(f/self._f_se()))
        B = (self._R_se_prim(self._f_se()) - self._R_DC_prim())/np.sqrt(1. + f**2./self._f_0()**2.)
        C = 1.+ 0.2*np.log(1. + self._f_se()/f)/(1. + self._w/self._h)
        return self._R_DC_prim() + (A - B - self._R_DC_prim())/C

    def get_characteristic_impedance(self, f):
        '''Return the absolute value of the characteristic impedance of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Characteristic impedance in Ohms
        '''
        
        return abs(np.sqrt((self.get_resistance_per_unit_length(f) + 1j*self.get_inductance_per_unit_length(f)*self._omega(f))/(self.get_conductance_per_unit_length(f) + 1j*self.get_capacitance_per_unit_length(f)*self._omega(f))))

    def get_complex_wave_vector(self, f):
        '''Return the absolute value of the complex wave vector coefficient of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Absolute value of the complex wave vector coefficient
        '''
        
        return abs(np.sqrt((self.get_resistance_per_unit_length(f) + 1j*self.get_inductance_per_unit_length(f)*self._omega(f))*(self.get_conductance_per_unit_length(f) + 1j*self.get_capacitance_per_unit_length(f)*self._omega(f))))

    def get_attenuation(self, f):
        '''Return the attenuation coefficient of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Attenuation coefficient
        '''
        return np.sqrt((self.get_resistance_per_unit_length(f) + 1j*self._omega(f)*self.get_inductance_per_unit_length(f))*(self.get_conductance_per_unit_length(f) + 1j*self._omega(f)*self.get_capacitance_per_unit_length(f))).real

    def get_wave_vector(self, f):
        '''Return the wave vector coefficient of the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Beta coefficient
        '''
        return np.sqrt((self.get_resistance_per_unit_length(f) + 1j*self._omega(f)*self.get_inductance_per_unit_length(f))*(self.get_conductance_per_unit_length(f) + 1j*self._omega(f)*self.get_capacitance_per_unit_length(f))).imag

    def get_velocity(self, f):
        '''Return the velocity of the wave in the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Velocity in unit of c (speed of light)
        '''
        
        return 1./np.sqrt(self.get_capacitance_per_unit_length(f) * self.get_inductance_per_unit_length(f))/cst.c

    def get_effective_permitivity(self, f):
        '''Return the effective permitivity felt from the microstrip line
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Effective permitivity in farrad per meter
        '''
        
        return (self.get_wave_vector(f)*cst.c/self._omega(f))**2.

    #################################################################################
    #
    #
    #                                    Print methods
    #
    #
    #################################################################################

    def print_parameters(self):
        '''
            Summarize all parameters of the microstrip object.
        '''
        
        print '------------------------------------------------------------'
        print '            Parameters'
        print ''
        print '    Central line width:        '+str(self._w*1e6)+'        µm'
        print '    Substrat height :        '+str(self._h*1e6)+'        µm'
        print '    Thickness:            '+str(self._t*1e6)+'        µm'
        print '    Ground plane width:        '+str(self._w_g*1e6)+'        µm'
        print ''
        print '    Relative permitivity:        '+str(self._epsilon_r)+'        F/m'
        print '    Loss tangente:            '+str(self._tan_delta)
        print '    Electrical conductivity:    '+str(self._kappa)+'    S/m'

    def print_results(self, frequency):
        '''
            Summarize all results of the microstrip object.
        '''
        
        print '------------------------------------------------------------'
        print '            Results'
        print ''
        print '    Inductance per unit length:    '+str(self.get_inductance_per_unit_length(frequency))+'    H/m'
        print '    Capacitance per unit length:    '+str(self.get_capacitance_per_unit_length(frequency))+'    F/m'
        print '    Resistance per unit length:    '+str(self.get_resistance_per_unit_length(frequency))+'    ohm/m'
        print '    Conductance per unit length:    '+str(self.get_conductance_per_unit_length(frequency))+'    S/m'
        print ''
        print '    Attenuation:            '+str(self.get_attenuation(frequency))+'    /m'
        print '    Wave vector:            '+str(self.get_wave_vector(frequency))+'    /m'
        print '    Velocity:            '+str(self.get_velocity(frequency))+'    c'
        print '    Characteristic impedance:    '+str(self.get_characteristic_impedance(frequency))+'    ohm'

    #################################################################################
    #
    #
    #                                    Find methods
    #
    #
    #################################################################################

    def find_length_lambda_over_two_resonance(self, f0):
        '''Return the length of the cavity needed to obtain a lambda over two resonance at the frequency f0
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Length in meter
        '''
        
        return 1./np.sqrt(self.get_capacitance_per_unit_length(f0)*self.get_inductance_per_unit_length(f0))/2./f0

    def find_length_lambda_over_four_resonance(self, f0):
        '''Return the length of the cavity needed to obtain a lambda over four resonance at the frequency f0
                - Input :
                    - Frequency (float | list | numpy.ndarray) in Hertz
                
                - Output :
                    - Length in meter
        '''
        
        return 1./np.sqrt(self.get_capacitance_per_unit_length(f0)*self.get_inductance_per_unit_length(f0))/4./f0



    #################################################################################
    #
    #
    #                                    Get equivalent electrical component
    #
    #
    #################################################################################

    def get_equivalent_resistance_lambda_over_two_resonance(self, f):
        l = self.find_length_lambda_over_two_resonance(f)
        return self.get_characteristic_impedance(f)/self.get_attenuation(f)/l



    def get_equivalent_capacitance_lambda_over_two_resonance(self, f):
        return 1./4./f/self.get_characteristic_impedance(f)



    def get_equivalent_conductance_lambda_over_two_resonance(self, f):
        return self.get_characteristic_impedance(f)/np.pi**2./f
