import textwrap
import sys
from engineering_notation import EngNumber
import argparse
from argparse import RawTextHelpFormatter
import time
import os
import math
import subprocess
from subprocess import Popen, PIPE
import cmath


class psikic:

    def __init__(self):
        # See figure1: The layout of a square planar spiral coil and its design parameters.
        # Inductance Formula for Square Planar Spiral Inductors with Rectangular
        # Conductor Cross Section
        # H. A. Aebischer
        # ADVANCED ELECTROMAGNETICS, VOL. 8, NO. 4, SEPTEMBER 2019
        #
        # Note: units of length are millmetres (mm)
        #
        # N, number of turns or windings, N>=2
        # Ai, innermost mid-conductor side length
        # w, winding distance or pitch (w = s + g)
        # g, gap or spacing between windings
        # s, conductor width
        # h, conductor height or thickness, hidden in Figure 1 of the paper.
        self.parseCommandLineArgs()

    def kicadFootprintStartAndEndCoordinates(self):
        Ai = self.Ai
        s = self.s
        g = self.g
        start_x = [eval('-Ai/2 + s + g'), eval('-Ai/2 + s + g')]
        start_y = [0]
        end_x = []
        end_y = []
        sign = -1
        count = -1
        for i in range(2, self.n):
            if i % 4 == 2:
                count += 1
            if i % 2 == 0:
                sign *= -1
            start_x.append(eval('({0})*Ai/2 + ({0})*({1})*(s + g)'.format(sign, count)))

        sign = -1
        count = -1
        for i in range(1, self.n):
            if i % 4 == 1:
                count += 1
            if i % 2 == 1:
                sign *= -1
            start_y.append(eval('{0}*Ai/2 + ({0})*{1}*(s+g)'.format(sign, count)))

        sign = -1
        count = -1
        end_y = []
        for i in range(self.n - 1):
            if i % 4 == 0:
                count += 1
            if i % 2 == 0:
                sign *= -1
            end_y.append(eval('{0}*Ai/2 + ({0})*{1}*(s+g)'.format(sign, count)))
        end_y.append(0)

        end_x = [eval('-Ai/2 + s + g')]
        for i in range(1, self.n):
            end_x.append(end_y[i - 1])

        size_x = []
        size_y = []

        for i in range(len(start_x)):
            size_x.append(end_x[i] - start_x[i])
            size_y.append(end_y[i] - start_y[i])

        return start_x, start_y, end_x, end_y, size_x, size_y

    def formCourtYard(self, x, y):
        courtYard = ''
        courtYard += '  (fp_line (start -{0} -{1}) (end {0} -{1}) (layer F.CrtYd) (width 0.12))\n'.format(x + self.w,
                                                                                                          y + self.w)
        courtYard += '  (fp_line (start -{0} {1}) (end {0} {1}) (layer F.CrtYd) (width 0.12))\n'.format(x + self.w,
                                                                                                        y + self.w)
        courtYard += '  (fp_line (start -{0} -{1}) (end -{0} {1}) (layer F.CrtYd) (width 0.12))\n'.format(x + self.w,
                                                                                                          y + self.w)
        courtYard += '  (fp_line (start {0} {1}) (end {0} -{1}) (layer F.CrtYd) (width 0.12))\n'.format(x + self.w,
                                                                                                        y + self.w)
        return courtYard

    def formForwardSilk(self, x, y):
        silk = ''
        silk += '  (fp_line (start -{0} -{1}) (end {0} -{1}) (layer F.SilkS) (width 0.12))\n'.format(x + 1 * self.w,
                                                                                                     y + 1 * self.w)
        silk += '  (fp_line (start -{0} {1}) (end {0} {1}) (layer F.SilkS) (width 0.12))\n'.format(x + 1 * self.w,
                                                                                                   y + 1 * self.w)
        silk += '  (fp_line (start -{0} -{1}) (end -{0} {1}) (layer F.SilkS) (width 0.12))\n'.format(x + 1 * self.w,
                                                                                                     y + 1 * self.w)
        silk += '  (fp_line (start {0} {1}) (end {0} -{1}) (layer F.SilkS) (width 0.12))\n'.format(x + 1 * self.w,
                                                                                                   y + 1 * self.w)
        return silk

    def formKicadPCBFootprint(self, start_x, start_y, end_x, end_y):
        L = float(self.L)
        header_text = \
            """
            (module planar_inductance (layer F.Cu) (tedit 5E123001)
                (fp_text reference REF** (at 0 0.5) (layer F.SilkS)
                    (effects (font (size 1 1) (thickness 0.15)))
                )
                (fp_text value planar_inductance_{0} (at 0 {8}) (layer F.Fab)
                    (effects (font (size 1 1) (thickness 0.15)))
                )
                (fp_text user N={1}_Ai={2}_s={3}_g={4}_h={5}__{0} (at {6} {7}) (layer Cmts.User)
                    (effects (font (size 1 1) (thickness 0.15)))
                )
            """.format(str(EngNumber(L, precision=4)) + 'H', self.N, self.Ai, self.s, self.g, self.h, 0,
                       end_y[-2] - 2 * self.w, -1 * (end_y[-2] - 2 * self.w))
        footprint = ''
        footprint += header_text
        footprint += self.formCourtYard(abs(start_x[-2]), abs(start_y[-2]))
        footprint += self.formForwardSilk(abs(start_x[-2]), abs(start_y[-2]))
        line_format = '(fp_line (start {} {}) (end {} {}) (layer F.Cu) (width {}))'
        for i in range(len(start_x)):
            footprint += line_format.format(start_x[i], start_y[i], end_x[i], end_y[i], self.s)

        footprint += \
            """
           (pad 2 thru_hole circle (at {:.4f} {:.4f}) (size 0.6 0.6) (drill 0.3) (layers *.Cu)
             (zone_connect 2))
           (pad 1 thru_hole circle (at {:.4f} {:.4f}) (size 0.6 0.6) (drill 0.3) (layers *.Cu)
             (zone_connect 2)) 
            """.format(start_x[0], start_y[0], end_x[-1], end_y[-1])

        footprint += ')'

        filename = "planar_inductor_{}_{}_{}_{}.kicad_mod".format(self.Ai, self.s, self.g, self.N)
        f = open(filename, "w+")
        f.write(footprint)
        f.close()

        return footprint

    def generateFasthenryInputFile(self, start_x, start_y, end_x, end_y):
        fasthenry = ''
        title = '** Planar Inductor **\n'
        fasthenry += title
        units = '* The following line names millimeters as the length units for the rest\n'
        units += '* of the file.\n'
        units += '.Units MM\n\n'
        fasthenry += units
        defaults = ''
        defaults += '* Make z=0 the default z coordinate and copper the default conductivity.\n'
        defaults += '* Note that the conductivity is in units of 1/(mm*Ohms), and not 1/(m*Ohms)\n'
        defaults += '* since the default units are millimetres\n'
        defaults += '.Default z=0 sigma=5.9595e4\n\n'
        fasthenry += defaults

        nodes = ''
        nodes += '* The nodes of the planar inductor (z=0 is the default)\n'
        nodeCoordinates_x = []
        nodeCoordinates_y = []
        for i in range(len(start_x)):
            nodeCoordinates_x.append('{}'.format(start_x[i]))
            nodeCoordinates_y.append('{}'.format(-1 * start_y[i]))
        nodeCoordinates_x.append('{}'.format(start_x[-1]))
        nodeCoordinates_y.append('{}'.format(0))

        for node in range(len(nodeCoordinates_x)):
            nodes += '{:10} x={:10} y={:10}\n'.format('N' + str(node + 1), nodeCoordinates_x[node],
                                                      nodeCoordinates_y[node])

        nodes += '\n\n* The segments connecting the nodes\n'
        for node in range(self.n):
            nodes += 'E{0} N{0} N{1} w={2} h={3} nhinc={4} nwinc={4}\n'.format(node + 1, node + 2, self.s, self.h, 5)
        nodes += '\n*  Define one input \'port\' of the network\n'
        nodes += '.external N{} N{}\n'.format(1, self.n + 1)
        nodes += '\n* Frequency range of interest\n'
        nodes += '.freq fmin=200000 fmax=200000 ndec=1\n\n'
        nodes += '.end\n'

        fasthenry += nodes

        filename = "planar_inductor.inp"
        f = open(filename, "w+")
        f.write(fasthenry)
        f.close()

        return fasthenry

    def calculateImpedanceOfSquarePlanarInductorUsingFormula(self):
        """
        L_SquarePlanarSpiral to calculate the DC inductance of square planar
        spiral coils with the help of equation (25) is given in appendix 2.
        """
        # from MATLAB code for L_SquarePlanarSpiral (N, A, w, s, h)
        A = self.A / 1000
        Ai = self.Ai / 1000
        N = self.N
        s = self.s / 1000
        g = self.g / 1000
        w = s + g
        h = self.h / 1000
        print('Ai={}, A={}, N={}, s={}, g={}, w={}, h={}'.format(Ai, A, N, s, g, w, h))
        print('Paste the following into the MATLAB script:')
        print('L_SquarePlanarSpiral({}, {}, {}, {},{})\n'.format(N, A, w, s, h))
        if self.N > 20:
            print('Formula is only valid for N in the range 2 to 20.')
            return float('nan'), float('nan')

        mu2 = 2e-7
        rho = ((N - 1) * w + s) / (A - (N - 1) * w)

        l = 4 * N * Ai + (4 * N ** 2 - 5 * N) * w
        dp4 = 4 * N * (N + 1) / 3 * w
        c1 = '{}'.format(float("nan"))
        c2 = '{}'.format(float("nan"))

        if N == 2:
            if rho > 0.36001:
                print('Invalid parameter combination.')
                return float('nan'), float('nan')
            c1 = '0*{0}'
            c2 = '0*{0}'
        if N in range(3, 8):
            if rho > 0.52001:
                print('Invalid parameter combination.')
                return float('nan'), float('nan')
            c1 = '-(7.2*{0} + 0.35) / (2.78*{0} + 1)'
            c2 = '-12.8*{0}**2 + 11*{0} + 0.80'
        if N in range(8, 13):
            if rho > 0.78001:
                print('Invalid parameter combination.')
                return float('nan'), float('nan')
            c1 = '-(2.1*{0} + 0.17) / (0.75*{0} + 1)'
            c2 = '-(0.59*{0} + 1.2) / (-0.90*{0} + 1)'
        if N in range(13, 21):
            if rho > 86001:
                print('Invalid parameter combination.')
                return float('nan'), float('nan')
            c1 = '-(1.5*{0} + 0.11) / (0.98*{0} + 1)'
            c2 = '-(2.85*{0} + 2.48) /(-0.75*{0} + 1)`'

        L_partial = mu2 * l * (math.log(l / N / (s + h)) - 0.19315)
        Mm = mu2 * l * 0.46716
        Mp = mu2 * l * (math.log(math.sqrt(1 + (l / dp4) ** 2) + l / dp4) - math.sqrt(1 + (dp4 / l) ** 2) + dp4 / l)
        L0 = L_partial - N * Mm + (N - 1) * Mp
        L = L0 * (1 - (eval(c1.format(rho)) * N + eval(c2.format(rho))) / 100)

        return L, rho

    def determineImpedanceOfSquarePlanarInductorUsingFasthenry(self):
        inp = "planar_inductor.inp"
        if os.path.exists(inp):
            logfile = 'logfile.txt'
            output_to_logfile = open(logfile, 'w+')
            if os.path.exists('Zc.mat'):
                os.remove('Zc.mat')
            p = Popen(["fasthenry", inp], stdout=output_to_logfile, stderr=subprocess.PIPE, universal_newlines=True)
            p.wait()

            if os.path.exists('Zc.mat'):
                print('Zc.mat successfully created.')
                lines = []
                f = open('Zc.mat', 'r')
                for line in f:
                    lines.append(line)
                f.close()
                z = lines[-1].split()
                z = complex(float(z[0]), float(z[1][:-1]))
                print('real part of impedance = {}'.format(z.real))
                print('imag part of impedance = {}'.format(z.imag))
                inductance = z.imag / (2 * math.pi * 200000)
                self.L = inductance
                print('Fasthenry determined Inductance of the planar inductor to be: {}H'.format(
                    EngNumber(self.L, precision=3)))
            else:
                print('Zc.mat NOT generated!')
                sys.exit()
        else:
            print("The file 'planar_inductor.inp' the input file to Fasthenry doesn't exist!")
            sys.exit()

    def generateOctaveScript(self):
        script = \
            """
            % Inductance Formula for Square Planar Spiral Inductors with Rectangular 
            % Conductor Cross Section
            % H. A. Aebischer
            % ADVANCED ELECTROMAGNETICS, VOL. 8, NO. 4, SEPTEMBER 2019
            %
            % source code of a MATLAB® function to implement the formula given in the
            % paper.
            % Besides the inductance L, the function also returns the filling factor ρ.
            %
            % NOTE NOTE NOTE A is the outside length see the diagram in the paper!!!!
    
            function [L, rho] = L_SquarePlanarSpiral (N, A, w, s, h)
                L = NaN;
                mu2 = 2e-7;
                rho = ((N - 1)*w + s) / (A - (N - 1)*w);
                Ai = A - 2*(N - 1)*w;
                l = 4*N*Ai + (4*N^2 - 5*N)*w;
                dp4 = 4*N * (N + 1)/3*w;
    
                switch N,
                    case 2;
                     if rho > 0.36001,
                        disp ('Invalid parameter combination.');
                        return;
                     end
                     c1 = @(rho) 0;
                     c2 = @(rho) 0;
    
                    case {3, 4, 5, 6, 7};
                     if rho > 0.52001,
                        disp ('Invalid parameter combination.');
                        return;
                     end
                     c1 = @(rho) -(7.2*rho + 0.35) / (2.78*rho + 1);
                     c2 = @(rho) -12.8*rho^2 + 11*rho + 0.80;
    
                    case {8, 9, 10, 11, 12};
                     if rho > 0.78001,
                        disp ('Invalid parameter combination.');
                        return;
                     end
                     c1 = @(rho) -(2.1*rho + 0.17) / (0.75*rho + 1);
                     c2 = @(rho) -(0.59*rho + 1.2) / (-0.90*rho + 1);
    
                    case {13, 14, 15, 16, 17, 18, 19, 20};
                     if rho > 0.86001,
                        disp ('Invalid parameter combination.');
                        return;
                     end
                     c1 = @(rho) -(1.5*rho + 0.11) / (0.98*rho + 1);
                     c2 = @(rho) -(2.85*rho + 2.48) /(-0.75*rho + 1);
                end 
    
                L_partial = mu2*l*(log(l/N/(s + h)) - 0.19315);
                Mm = mu2*l*0.46716;
                Mp = mu2*l*(log(sqrt(1 + (l/dp4)^2) + l/dp4)...
                 - sqrt(1 + (dp4/l)^2) + dp4/l);
                L0 = L_partial - N*Mm + (N - 1)*Mp;
                L = L0 * (1 - (c1(rho)*N + c2(rho))/100); 
            """
        filename = 'L_SquarePlanarSpiral.m'
        f = open(filename, 'w+')
        f.write(script)
        f.close()

    def runOctaveScript(self):
        if self.N >= 2 and self.N <= 20:
            self.generateOctaveScript()
        else:
            print("Octave script is only valid for N>=2 and N<=20!")
            return

        if os.path.exists('./L_SquarePlanarSpiral.m'):
            output_to_logfile = open('octave.log', 'w+')
            script = ''
            script += 'N={};\n'.format(self.N)
            script += 'A={};\n'.format(self.Ai + 2 * (self.N - 1) * (self.w))
            script += 'w={};\n'.format(self.w)
            script += 's={};\n'.format(self.s)
            script += 'h={};\n'.format(self.h)
            script += '[L, rho] = L_SquarePlanarSpiral (N, A, w, s, h);\n'
            script += 'fprintf("%e,%f",L, rho);'
            f = open('doit.m', 'w+')
            f.write(script)
            f.close()
            if os.path.exists('./doit.m'):
                p = Popen(["octave", 'doit.m'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
                p.wait()
                output, errors = p.communicate()
                if output != None:
                    L, rho = output.split(',')
                    if 'NaN' in L:
                        print('{}'.format('** Invalid parameter combination reported from Octave script. **'))
                        print('** Calculated value of rho = {} **\n'.format(rho))
                    else:
                        print('\nCalculated value of L = {}H'.format(EngNumber(float(L) / 1000, precision=3)))
                        print('Calculated value of rho = {:.4f}\n'.format(float(rho)))
                elif errors != None:
                    print(errors)
                    sys.exit()
                else:
                    print('No calculated output available!')
                    sys.exit()
            else:
                print('doit.m does not exist!')
                sys.exit()
        else:
            print('L_SquarePlanarSpiral.m does not exist!')
            sys.exit()

        return

    def parseCommandLineArgs(self):

        tag = textwrap.dedent( \
            """
 _______  _______ _________ _       _________ _______ 
(  ____ )(  ____ {0}__   __/| \    /{0}__   __/(  ____ {1}
| (    )|| (    \/   ) (   |  \  / /   ) (   | (    \/
| (____)|| (_____    | |   |  (_/ /    | |   | |      
|  _____)(_____  )   | |   |   _ (     | |   | |      
| (            ) |   | |   |  ( \ \    | |   | |      
| )      /\____) |___) (___|  /  \ \___) (___| (____/{1}
|/       \_______)\_______/|_/    \/\_______/(_______/
            """).format("\\\\", "\\")
        print(tag)

        examples = textwrap.dedent( \
            """
              examples:
                %(prog)s -N {0:>2} -Ai {4} -s {8} -g {9} -t {10}
                %(prog)s -N {1:>2} -Ai {5} -s {8} -g {9} -t {10}
                %(prog)s -N {2:>2} -Ai {6} -s {8} -g {9} -t {10}
                %(prog)s -N {3:>2} -Ai {7} -s {8} -g {9} -t {10}
            """.format(2, 5, 10, 15, 47.2, 38.8, 24.8, 10.8, 0.7, 0.7, 0.035))

        parser = argparse.ArgumentParser(description='KiCad PCB Footprint Utility.', \
                                         epilog=examples, \
                                         formatter_class=RawTextHelpFormatter)
        parser.add_argument('-N', required=True, help=u'number of turns or windings, N \u2265 2'),
        parser.add_argument('-Ai', required=True, help='innermost mid-conductor side length, mm')
        parser.add_argument('-s', required=True, help='conductor width, mm')
        parser.add_argument('-g', required=True, help='gap or spacing between windings, mm')
        parser.add_argument('-t', required=True, help='trace thickness or height, mm')

        args = vars(parser.parse_args())
        self.N = eval(args['N'])
        self.Ai = eval(args['Ai'])
        self.g = eval(args['g'])
        self.s = eval(args['s'])
        self.h = eval(args['t'])

        if self.s < 0.127:
            print('The track width has to be greater than 0.127mm or 5mil.')

        if self.g < 0.127:
            print('The track spacing has to be greater than 0.127mm or 5mil.')

        constraint = self.s / 2 + self.g
        if constraint < 0.35:
            print('s/2 + g > 0.35mm')
            print('{} + {} = {}'.format(self.s / 2, self.g, constraint))
            sys.exit()

        self.w = self.s + self.g
        self.n = 4 * self.N + 1  # number of sides
        self.A = self.Ai + 2 * (self.N - 1) * self.w

        print('\nN={}'.format(args['N']))
        print('Ai={}'.format(args['Ai']))
        print('s={}'.format(args['s']))
        print('g={}'.format(args['g']))
        print('t={}'.format(args['t']))
        print('A={}'.format(self.A))
        print('w={}\n'.format(self.w))


#  instantiate the kicad class

k = psikic()
k.runOctaveScript()
start_x, start_y, end_x, end_y, size_x, size_y = k.kicadFootprintStartAndEndCoordinates()
fasthenry = k.generateFasthenryInputFile(start_x, start_y, end_x, end_y)
k.determineImpedanceOfSquarePlanarInductorUsingFasthenry()
footprint = k.formKicadPCBFootprint(start_x, start_y, end_x, end_y)



