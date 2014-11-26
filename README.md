microstrip
===

microstrip is a class allowing the calculation of microstrip Waveguide characteristic values following geometric parameters, conductivity, relative permitivity and, loss tangent.

All method are well documented into the class.


Example:

`````python

import microstrip

line = microstrip.microstrip()

line.set_thickness(0.5e-6)

line.get_inductance_per_unit_length(2e9)
//4.26e-7

`````


This class is based on an article of Frank Schnieder and Wolfgang Heinrich
"Model of thin-Film Microstrip Line for Circuit Design"
 IEEE Transactions on Microwave Theory And Techniques, vol 49, n° 1, January 2001

