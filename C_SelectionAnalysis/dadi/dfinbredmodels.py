"""
Population demographic models edited for inbreeding inclusion
"""
import numpy

from dadi import Numerics, PhiManip, Integration
from dadi.Spectrum_mod import Spectrum
from dadi.PortikModels.portik_models_2d import *

def growth_inbreeding(params, ns, pts):
    """
    Exponential growth beginning some time ago.

    params = (nu,T)
    ns = (n1,)

    nu: Ratio of contemporary to ancient population size
    T: Time in the past at which growth began (in units of 2*Na 
       generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T, F = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: numpy.exp(numpy.log(nu) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

def two_epoch_inbreeding(params, ns, pts):
    """
    Instantaneous size change some time ago.

    params = (nu,T)
    ns = (n1,)

    nu: Ratio of contemporary to ancient population size
    T: Time in the past at which size change happened (in units of 2*Na 
       generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu,T, F = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    
    phi = Integration.one_pop(phi, xx, T, nu)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

def bottlegrowth_1d_inbreeding(params, ns, pts):
    """
    Instantanous size change followed by exponential growth.

    params = (nuB,nuF,T)
    ns = (n1,)

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contemporary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    n1: Number of samples in resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T, F = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
    phi = Integration.one_pop(phi, xx, T, nu_func)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F,), (2,))
    return fs

def snm_2d_inbred(params, ns, pts):
    """
    ns = (n1,n2)

    Standard neutral model, populations never diverge.
    """
    F1 = params[0]
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,), (F1,), (2,))
    return fs
snm_2d_inbred.__param_names__ = []
snm = snm_2d_inbred

def bottlegrowth_2d_inbred(params, ns, pts):
    """
    params = (nuB,nuF,T)
    ns = (n1,n2)

    Instantanous size change followed by exponential growth with no population
    split.

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T,F1,F2 = params
    return bottlegrowth_split_mig_inbred((nuB,nuF,0,T,0,F1,F2), ns, pts)
bottlegrowth_2d_inbred.__param_names__ = ['nuB', 'nuF', 'T']
bottlegrowth = bottlegrowth_2d_inbred

def bottlegrowth_split_inbred(params, ns, pts):
    """
    params = (nuB,nuF,T,Ts)
    ns = (n1,n2)

    Instantanous size change followed by exponential growth then split.

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    Ts: Time in the past at which the two populations split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,T,Ts,F1, F2 = params
    return bottlegrowth_split_mig_inbred((nuB,nuF,0,T,Ts,F1,F2), ns, pts)
bottlegrowth_split_inbred.__param_names__ = ['nuB', 'nuF', 'T', 'Ts']

def bottlegrowth_split_mig_inbred(params, ns, pts):
    """
    params = (nuB,nuF,m,T,Ts)
    ns = (n1,n2)

    Instantanous size change followed by exponential growth then split with
    migration.

    nuB: Ratio of population size after instantanous change to ancient
         population size
    nuF: Ratio of contempoary to ancient population size
    m: Migration rate between the two populations (2*Na*m).
    T: Time in the past at which instantaneous change happened and growth began
       (in units of 2*Na generations) 
    Ts: Time in the past at which the two populations split.
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuB,nuF,m,T,Ts,F1,F2 = params

    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)

    if T >= Ts:
        nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
        phi = Integration.one_pop(phi, xx, T-Ts, nu_func)

        phi = PhiManip.phi_1D_to_2D(xx, phi)
        nu0 = nu_func(T-Ts)
        nu_func = lambda t: nu0*numpy.exp(numpy.log(nuF/nu0) * t/Ts)
        phi = Integration.two_pops(phi, xx, Ts, nu_func, nu_func, m12=m, m21=m)
    else:
        phi = PhiManip.phi_1D_to_2D(xx, phi)
        phi = Integration.two_pops(phi, xx, Ts-T, 1, 1, m12=m, m21=m)
        nu_func = lambda t: nuB*numpy.exp(numpy.log(nuF/nuB) * t/T)
        phi = Integration.two_pops(phi, xx, T, nu_func, nu_func, m12=m, m21=m)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    return fs
bottlegrowth_split_mig_inbred.__param_names__ = ['nuB', 'nuF', 'm', 'T', 'Ts']

def split_mig_inbred(params, ns, pts):
    """
    params = (nu1,nu2,T,m)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m: Migration rate between populations (2*Na*m)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T,m,F1,F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m, m21=m)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    return fs
split_mig_inbred.__param_names__ = ['nu1', 'nu2', 'T', 'm']

def split_mig_mscore(params):
    """
    ms core command for split_mig.
    """
    nu1,nu2,T,m = params

    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-ej %(T)f 2 1 -en %(T)f 1 1"

    sub_dict = {'nu1':nu1, 'nu2':nu2, 'm12':2*m, 'm21':2*m, 'T': T/2}

    return command % sub_dict
split_mig_mscore.__param_names__ = ['nu1', 'nu2', 'T', 'm']

def split_asym_mig_inbred(params, ns, pts):
    """
    params = (nu1,nu2,T,m12,m21)
    ns = (n1,n2)

    Split into two populations of specifed size, with asymetric migration .

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2 (2*Na*m21)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,T,m12,m21,F1,F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    phi = Integration.two_pops(phi, xx, T, nu1, nu2, m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    return fs
split_asym_mig_inbred.__param_names__ = ['nu1', 'nu2', 'T', 'm12', 'm21']

def split_delay_mig_inbred(params, ns, pts):
    """
    params = (nu1,nu2,Tpre,Tmig,m12,m21)
    ns = (n1,n2)

    Split into two populations of specifed size, with migration after some time has passed post split.

    nu1: Size of population 1 after split.
    nu2: Size of population 2 after split.
    Tpre: Time in the past after split but before migration (in units of 2*Na generations) 
    Tmig: Time in the past after migration starts (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2 (2*Na*m21)
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nu1,nu2,Tpre,Tmig,m12,m21,F1,F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = Integration.two_pops(phi, xx, Tpre, nu1, nu2, m12=0, m21=0)
    phi = Integration.two_pops(phi, xx, Tmig, nu1, nu2, m12=m12, m21=m21)
    
    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    return fs
split_delay_mig_inbred.__param_names__ = ['nu1', 'nu2', 'Tpre', 'Tmig', 'm12', 'm21']

def IM_inbred(params, ns, pts):
    """
    ns = (n1,n2)
    params = (s,nu1,nu2,T,m12,m21)

    Isolation-with-migration model with exponential pop growth.

    s: Size of pop 1 after split. (Pop 2 has size 1-s.)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    s,nu1,nu2,T,m12,m21,F1,F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_func = lambda t: s * (nu1/s)**(t/T)
    nu2_func = lambda t: (1-s) * (nu2/(1-s))**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)
    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    return fs
IM_inbred.__param_names__ = ['s', 'nu1', 'nu2', 'T', 'm12', 'm21']

def IM_mscore(params):
    """
    ms core command for IM.
    """
    s,nu1,nu2,T,m12,m21 = params

    alpha1 = numpy.log(nu1/s)/T
    alpha2 = numpy.log(nu2/(1-s))/T
    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-eg 0 1 %(alpha1)f -eg 0 2 %(alpha2)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-ej %(T)f 2 1 -en %(T)f 1 1"

    sub_dict = {'nu1':nu1, 'nu2':nu2, 'alpha1':2*alpha1, 'alpha2':2*alpha2,
                'm12':2*m12, 'm21':2*m21, 'T': T/2}

    return command % sub_dict
IM_mscore.__param_names__ = ['s', 'nu1', 'nu2', 'T', 'm12', 'm21']

def IM_pre_inbred(params, ns, pts):
    """
    params = (nuPre,TPre,s,nu1,nu2,T,m12,m21)
    ns = (n1,n2)

    Isolation-with-migration model with exponential pop growth and a size change
    prior to split.

    nuPre: Size after first size change
    TPre: Time before split of first size change.
    s: Fraction of nuPre that goes to pop1. (Pop 2 has size nuPre*(1-s).)
    nu1: Final size of pop 1.
    nu2: Final size of pop 2.
    T: Time in the past of split (in units of 2*Na generations) 
    m12: Migration from pop 2 to pop 1 (2*Na*m12)
    m21: Migration from pop 1 to pop 2
    n1,n2: Sample sizes of resulting Spectrum
    pts: Number of grid points to use in integration.
    """
    nuPre,TPre,s,nu1,nu2,T,m12,m21,F1,F2 = params

    xx = Numerics.default_grid(pts)

    phi = PhiManip.phi_1D(xx)
    phi = Integration.one_pop(phi, xx, TPre, nu=nuPre)
    phi = PhiManip.phi_1D_to_2D(xx, phi)

    nu1_0 = nuPre*s
    nu2_0 = nuPre*(1-s)
    nu1_func = lambda t: nu1_0 * (nu1/nu1_0)**(t/T)
    nu2_func = lambda t: nu2_0 * (nu2/nu2_0)**(t/T)
    phi = Integration.two_pops(phi, xx, T, nu1_func, nu2_func,
                               m12=m12, m21=m21)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx,xx), (F1,F2), (2,2))
    return fs
IM_pre_inbred.__param_names__ = ['nuPre', 'TPre', 's', 'nu1', 'nu2', 'T', 'm12', 'm21']

def IM_pre_mscore(params):
    """
    ms core command for IM_pre.
    """
    nuPre,TPre,s,nu1,nu2,T,m12,m21 = params
    
    nu1_0 = nuPre*s
    nu2_0 = nuPre*(1-s)
    alpha1 = numpy.log(nu1/nu1_0)/T
    alpha2 = numpy.log(nu2/nu2_0)/T
    command = "-n 1 %(nu1)f -n 2 %(nu2)f "\
            "-eg 0 1 %(alpha1)f -eg 0 2 %(alpha2)f "\
            "-ma x %(m12)f %(m21)f x "\
            "-ej %(T)f 2 1 -en %(T)f 1 %(nuP)f "\
            "-en %(TP)f 1 1"

    sub_dict = {'nu1':nu1, 'nu2':nu2, 'alpha1':2*alpha1, 'alpha2':2*alpha2,
                'm12':2*m12, 'm21':2*m21, 'T': T/2, 'nuP':nuPre, 'TP':(T+TPre)/2}

    return command % sub_dict
IM_pre_mscore.__param_names__ = ['nuPre', 'TPre', 's', 'nu1', 'nu2', 'T', 'm12', 'm21']

def iso_inbred(params, ns, pts):
    T, nu1, nu2, F1, F2 = params
    xx = Numerics.default_grid(pts)
    phi = PhiManip.phi_1D(xx)
    phi = PhiManip.phi_1D_to_2D(xx, phi)
    phi = phi.Integration.two_pops(phi, xx, T, nu1, nu2)

    fs = Spectrum.from_phi_inbreeding(phi, ns, (xx, xx), (F1, F2), (2, 2))
    return fs