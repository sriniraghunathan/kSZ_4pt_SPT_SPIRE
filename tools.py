import numpy as np, os, sys
from pylab import *

def convert_amber_deltaz90_deltaz50(delta_z_90_or_50, a = 0.387, b = 0.009, invert = False):
    """
    if not invert: #convert delta_z_50 to delta_z_90
        delta_z_50 = a * delta_z_90 + b
    else: #convert delta_z_90 to delta_z_50
        delta_z_90 = (delta_z_50 - b)/a 

    Parameters
    ----------
    delta_z_90_or_50: Either delta_z_90 or  delta_z_50 value. Either float or np.array    
    a: slope. Default = 0.387; check Fig. 1 of https://arxiv.org/pdf/2203.04337.pdf
    b: intercept. Default = 0.009; check Fig. 1 of https://arxiv.org/pdf/2203.04337.pdf

    Returns
    ----------    
    delta_z_50: float or array if invert = False
    delta_z_90: float or array if invert = True
    """
    if not invert:
        return a * delta_z_90_or_50 + b
    else:
        return (delta_z_90_or_50 - b)/a

def get_ksz_4pt_parameterisation(delta_z_90, z_re, pivot_deltaz = 4., pivot_zre = 8.):
    """

    Parameters
    ----------
    delta_z_90: float
    z_re: float
    

    Returns
    ----------
    ksz_4pt_amber_fit: float
    """
    fit_params = [2.72567876e-05, 1.73720621, 2.51875476]
    return fit_params[0] * (delta_z_90 / pivot_deltaz)**fit_params[1] * (z_re / pivot_zre)**fit_params[2]

def get_latex_param_str(param, use_H = False):
    """
    function that returns latex str for parameters.

    Parameters
    ----------
    param: parameter symbol.
    use_H: Changed h to H0. Default is False

    Returns
    -------
    equivalent latex str of the parameter.
    """
    params_str_dic= {\
    'norm_YszM': r'${\rm log}(Y_{\ast})$', 'logYstar': r'${\rm log}(Y_{\ast})$', 'alpha_YszM': r'$\alpha_{_{Y}}$',\
    'beta_YszM': r'$\beta_{_{Y}}$', 'gamma_YszM': r'$\gamma_{_{Y}}$', \
    'sigma_logYszM': r'$\sigma_{_{\rm logY}}$', 'alpha_sigma_logYszM': r'$\alpha_{\sigma}$', 'gamma_sigma_logYszM': r'$\gamma_{\sigma}$', \
    'alpha': r'$\eta_{\rm v}$', 'sigma_8': r'$\sigma_{\rm 8}$', \
    'one_minus_hse_bias': r'$1-b_{\rm HSE}$', 'omega_m': r'$\Omega_{\rm m}$', 'thetastar': r'$\theta_{\ast}$',\
    'h':r'$h$', 'h0':r'$h$', 'm_nu':r'$\sum m_{\nu}$ [eV]', 'ombh2': r'$\Omega_{b}h^{2}$', 'omch2': r'$\Omega_{c}h^{2}$', 'w0': r'$w_{0}$', 'wa': r'$w_{a}$', \
    'tau': r'$\tau_{\rm re}$', 'As': r'$A_{\rm s}$', 'ns': r'$n_{\rm s}$', 'neff': r'N$_{\rm eff}$', \
    'rho_snr_mlens': r'$\rho_{\rm SNR-lens}$', \
    'slope_vir': r'$A_{\rm v}$', 'intercept_vir': r'$B_{\rm v}$', \
    'r': r'$r$', \
    'A_phi_sys': r'$A_{\phi}^{\rm sys}$', 'alpha_phi_sys': r'$\alpha_{\phi}^{\rm sys}$', \
    'alpha_radio': r'$\alpha_{\rm rad}$', 'alpha_radio_sigma': r'$\sigma(\alpha_{\rm rad})$', \
    'Aksz': r'$A_{\rm kSZ}$', \
    'Aksz_h': r'$A_{\rm kSZ}^{h}$', 'alphaksz_h': r'$\alpha_{\rm kSZ}^{h}$', \
    'zmid': r'$z_{\rm re}^{\rm mid}$', 'zdur': r'$\Delta z_{\rm re}$', \
    'Acibtsz': r'$A_{\rm CIB+tSZ}$', \
    #reionisation stuff
    'delta_z': r'$\Delta z_{\rm re, 50}$', 'zre': 'z_{\rm re}^{\rm mid}', 'ampcmbfg': r'$A_{\rm CMB-FG}$', \
    }

    if use_H:
        params_str_dic['h0'] = r'$H_{0}$'

    return params_str_dic[param]

def get_ini_param_dict(fpath = 'data/params_planck_r_0.0_2018_cosmo.txt'):
    """
    read params file and initialise cosmology

    Parameters
    ----------
    fpath: parameter file
    

    Returns
    ----------
    param_dict: dict containing parameter names and values.
    """
    try:
        params = np.recfromtxt(fpath, delimiter = '=', encoding = 'utf-8')
    except:
        params = np.recfromtxt(fpath, delimiter = '=')
    param_dict = {}
    for rec in params:
        val = rec[1].strip()##.decode("utf-8")
        try:
            val = val.decode("utf-8")
        except:
            pass
        try:
            if val.find('.')>-1:
                val = float(val)
            else:
                val = int(val)
        except:
            val = str(val)

        if val == 'None':
            val = None
        paramname = rec[0].strip()
        try:
            paramname = paramname.decode("utf-8")
        except:
            pass
        param_dict[paramname] = val

    return param_dict

def parent_get_amber_ksz_4pt_using_scaling(reqd_delta_z_50_arr, reqd_z_re_arr, basline_amber_ksz_clkk, pivot_delta_z_90 = 4., pivot_zre = 8.):
    """
    return ksz 4pt as fn(delta_z_90, z_re)

    Parameters
    ----------
    reqd_delta_z_50_arr: np.array
        Duration of reionisation.
    reqd_z_re_arr: np.array
        Redshift midpoint of reionisation.
    basline_amber_ksz_clkk: np.array
        baseline AMBER kSZ 4-pt CL_KK.
    pivot_delta_z_90: float
        Pivotal duration of reionisation in delta_z_90.
        Default is 4.
    pivot_zre: float
        Pivotal redshift midpoint of reionisation.
        Default is 8.


    Returns
    ----------    
    amber_ksz_clkk_arr: array
        CL_KK calcualted by scaling the baseline basline_amber_ksz_clkk.
    """
    reqd_delta_z_90_arr = convert_amber_deltaz90_deltaz50(reqd_delta_z_50_arr, invert = True)
    
    if len(reqd_z_re_arr) == 1:
        reqd_z_re_arr = np.tile(reqd_z_re_arr[0], len(reqd_delta_z_90_arr))

    baseline_ksz_4pt_amber_fit = get_ksz_4pt_parameterisation(pivot_delta_z_90, pivot_zre)
    amber_ksz_clkk_arr = []
    for tmpcntr,(reqd_delta_z_90_val, reqd_z_re_val) in enumerate(zip(reqd_delta_z_90_arr, reqd_z_re_arr)):
        curr_ksz_4pt_amber_fit = get_ksz_4pt_parameterisation(reqd_delta_z_90_val, reqd_z_re_val)
        scaling_value = curr_ksz_4pt_amber_fit / baseline_ksz_4pt_amber_fit
        curr_amber_ksz_clkk = basline_amber_ksz_clkk * scaling_value

        amber_ksz_clkk_arr.append( curr_amber_ksz_clkk )

    return amber_ksz_clkk_arr

def format_axis(ax, fx, fy, maxxloc=None, maxyloc = None):
    """
    function to format axis fontsize.


    Parameters
    ----------
    ax: subplot axis.
    fx: fontsize for xaxis.
    fy: fontsize for yaxis.
    maxxloc: total x ticks.
    maxyloc: total y ticks.

    Returns
    -------
    formatted axis "ax".
    """
    for label in ax.get_xticklabels(): label.set_fontsize(fx)
    for label in ax.get_yticklabels(): label.set_fontsize(fy)
    if maxyloc is not None:
        ax.yaxis.set_major_locator(MaxNLocator(nbins=maxxloc))
    if maxxloc is not None:
        ax.xaxis.set_major_locator(MaxNLocator(nbins=maxxloc))

    return ax

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height])#,axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    #subax.xaxis.set_tick_params(labelsize=x_labelsize)
    #subax.yaxis.set_tick_params(labelsize=y_labelsize)
    return subax    

def truncate_colormap(cmap, minval=0., maxval=0.6, n=255):
    '''
    https://stackoverflow.com/a/18926541
    '''
    if isinstance(cmap, str):
        cmap = plt.get_cmap(cmap)
    new_cmap = mpl.colors.LinearSegmentedColormap.from_list(
        'trunc({n},{a:.2f},{b:.2f})'.format(n=cmap.name, a=minval, b=maxval),
        cmap(np.linspace(minval, maxval, n)))
    return new_cmap   


def get_planck_zmid_from_tau(tau=0.0544):
    """
    return zmin given tau. Uses tanh reionisation model.


    Parameters
    ----------
    tau: float
        optical depth.
        Default is 0.0544

    Returns
    -------
    z_re: float
        redshift midpoint of reionisation.
    """
    param_dict = get_ini_param_dict()
    param_dict['tau'] = tau
    import camb
    pars = camb.CAMBparams(max_l_tensor = param_dict['max_l_tensor'], max_eta_k_tensor = param_dict['max_eta_k_tensor'])
    pars.set_accuracy(AccuracyBoost = param_dict['AccuracyBoost'], lAccuracyBoost = param_dict['lAccuracyBoost'], lSampleBoost = param_dict['lSampleBoost'],\
        DoLateRadTruncation = param_dict['do_late_rad_truncation'])
    pars.set_cosmology(thetastar=param_dict['thetastar'], ombh2=param_dict['ombh2'], omch2=param_dict['omch2'], nnu = param_dict['neff'], mnu=param_dict['mnu'], \
        omk=param_dict['omk'], tau=param_dict['tau'], YHe = param_dict['YHe'], Alens = param_dict['Alens'], \
        num_massive_neutrinos = param_dict['num_nu_massive'])
    pars.set_tau(param_dict['tau'])

    z_re = camb.get_zre_from_tau(pars,param_dict['tau'])

    return z_re


