import numpy as np
from pylab import *
z, xe, tau = np.loadtxt('tau.txt', unpack = True)

'''
plot(xe, z); axvline(0.5); axvline(0.25); axvline(0.75)
plot(tau, z);
ylabel(r'Redshift $z$', fontsize = 12)
xlabel(r'$\tau$ or $x(z)$', fontsize = 12)
show()
'''
plot(z, tau); axvline(8); show()

plot(z, xe); axvline(8); axhline(0.5);
plot(z, tau); axvline(8); 
xlabel(r'Redshift $z$', fontsize = 12)
ylabel(r'$\tau$ or $x(z)$', fontsize = 12)
show()


pars.Reion.set_tau(param_dict['tau'])
camb.get_zre_from_tau(pars,param_dict['tau'])

pars.Reion.set_tau(param_dict['tau'], delta_redshift = 10.)
camb.get_zre_from_tau(pars,param_dict['tau'])