"""
hizea.prospector
================

Code to do SED-fitting with Prospector.

"""
import os, time
import numpy as np

def logmass2mass(logmass=11.0, **extras):
    return 10**logmass

def load_model(template_library='delayed-tau', redshift=0.0, verbose=False):
    """
    http://dfm.io/python-fsps/current/stellarpop_api/#api-reference
    https://github.com/moustakas/siena-astrophysics/blob/master/research/redmapper/redmapper-stellar-mass.py#L125-L197    
    
    """
    from prospect.models import priors
    from prospect.models.sedmodel import SedModel
    from prospect.models.templates import TemplateLibrary
    from prospect.models.transforms import dustratio_to_dust1

    def base_delayed_tau():
        model_params = TemplateLibrary['parametric_sfh']

        # Initialize with sensible numbers.
        model_params['tau']['init'] = 10.0
        model_params['tage']['init'] = 1.0
        model_params['logzsol']['init'] = 0.2

        model_params['logmass'] = {'N': 1, 'isfree': True, 'init': 11.0,
                                   'prior': priors.TopHat(mini=10.0, maxi=12.0),
                                   'units': '$M_{\odot}$'}

        model_params['mass']['isfree'] = False
        model_params['mass']['init'] = 10**model_params['logmass']['init']
        model_params['mass']['prior'] = None
        model_params['mass']['depends_on'] = logmass2mass
        
        # Adjust the prior ranges.
        model_params['tau']['prior'] = priors.LogUniform(mini=0.1, maxi=30.0)
        model_params['tage']['prior'] = priors.LogUniform(mini=0.01, maxi=10.0)
        model_params['logzsol']['prior'] = priors.TopHat(mini=-0.5, maxi=0.2)

        #print('HACK!!!!!!!!!!!!!')
        model_params['tau']['isfree'] = False
        model_params['tage']['isfree'] = False
        model_params['logzsol']['isfree'] = False
        model_params['dust2']['isfree'] = False

        return model_params

    if template_library == 'delayed-tau':
        # Underlying delayed tau model.
        model_params = base_delayed_tau()

    if template_library == 'bursty':
        # Underlying delayed tau model.
        model_params = base_delayed_tau()
        
        # Add bursts
        model_params.update(TemplateLibrary['burst_sfh'])
        
        #model_params['tburst']['isfree'] = True
        #model_params['tburst']['init'] = 8.0
        #model_params['tburst']['prior'] = priors.TopHat(mini=0.0, maxi=10.0)

        model_params['fburst']['isfree'] = True
        model_params['fburst']['init'] = 0.1
        model_params['fburst']['prior'] = priors.TopHat(mini=0.0, maxi=1.0)

        model_params['fage_burst']['isfree'] = True
        model_params['fage_burst']['init'] = 0.9
        model_params['fage_burst']['prior'] = priors.TopHat(mini=0.5, maxi=1.0)

    # Add dust emission (with fixed dust SED parameters).
    #model_params.update(TemplateLibrary['dust_emission'])

    model_params['dust2']['init'] = 1.0 # diffuse dust
    model_params['dust2']['prior'] = priors.TopHat(mini=0.0, maxi=4.0)

    ## Add more dust flexibility.
    #model_params['dust_type'] = {'N': 1, 'isfree': False, 'init': 0, 'units': 'dust model'}
    #model_params['dust_index'] = {'N': 1, 'isfree': False, 'init': -0.7,
    #                              'units': 'power-law index', 'prior': None}
    #
    #model_params['dust1'] = {'N': 1, 'isfree': False, 'init': 0.0, 'prior': None,
    #                         'units': 'optical depth towards young stars',
    #                         'depends_on': dustratio_to_dust1}
    #model_params['dust_ratio'] = {'N': 1, 'isfree': True, 'init': 1.0,
    #                              'prior': priors.TopHat(mini=1.0, maxi=10.0),
    #                              'units': 'dust1/dust2 ratio (optical depth to young stars vs diffuse)'}

    ## Add nebular emission.
    #model_params.update(TemplateLibrary['nebular'])
    ##model_params['add_neb_continuum']['init'] = False
    #model_params['gas_logu']['init'] = -1.0 # harder radiation field [default is -2.0]

    # Fixed redshift.
    model_params['zred']['init'] = redshift
    model_params['zred']['isfree'] = False 

    # Change the IMF from Kroupa to Salpeter.
    #model_params['imf_type']['init'] = 0
        
    # Now instantiate the model using this new dictionary of parameter specifications
    model = SedModel(model_params)
    if verbose:
        print(model)

    return model

def load_sps(zcontinuous=1, verbose=False):
    """zcontinuous - interpolate between metallicity values.

    """
    from prospect.sources import CSPSpecBasis
    t0 = time.time()
    sps = CSPSpecBasis(zcontinuous=zcontinuous)
    if verbose:
        print('Loading SPS models took {:.2f} sec'.format(time.time()-t0))
    return sps

def lnprobfn(theta, model, obs, sps, nested=False, verbose=False):
    """
    Given a parameter vector, a dictionary of observational data 
    a model object, and an sps object, return the ln of the posterior. 
    This requires that an sps object (and if using spectra 
    and gaussian processes, a GP object) be instantiated.

    """
    from prospect.likelihood import lnlike_spec, lnlike_phot, write_log

    # Calculate prior probability and exit if not within prior
    # Also if doing nested sampling, do not include the basic priors, 
    # since the drawing method includes the prior probability
    lnp_prior = model.prior_product(theta, nested=nested)
    if not np.isfinite(lnp_prior):
        return -np.infty
        
    # Generate "mean" model
    t1 = time.time()
    spec, phot, mfrac = model.mean_model(theta, obs, sps=sps)
    d1 = time.time() - t1
 
    # Calculate likelihoods
    t2 = time.time()
    lnp_spec = lnlike_spec(spec, obs=obs)
    lnp_phot = lnlike_phot(phot, obs=obs)
    d2 = time.time() - t2
    if verbose:
        write_log(theta, lnp_prior, lnp_spec, lnp_phot, d1, d2)

    return lnp_prior + lnp_phot + lnp_spec

def chivecfn(theta, model, obs, sps, verbose):
    """A version of lnprobfn that returns the simple uncertainty 
    normalized residual instead of the log-posterior, for use with 
    least-squares optimization methods like Levenburg-Marquardt.
    
    It's important to note that the returned chi vector does not 
    include the prior probability.
    
    """
    from prospect.likelihood import chi_spec, chi_phot

    lnp_prior = model.prior_product(theta)
    if not np.isfinite(lnp_prior):
        return -np.infty

    # Generate mean model
    t1 = time.time()
    try:
        spec, phot, x = model.mean_model(theta, obs, sps=sps)
    except(ValueError):
        return -np.infty
    d1 = time.time() - t1

    chispec = chi_spec(spec, obs)
    chiphot = chi_phot(phot, obs)
    
    return np.concatenate([chispec, chiphot])

def run_minimize(model, obs, sps, seed=1, nmin=10, xtol=5e-16, ftol=5e-16,
                 maxnfev=500, min_method='trf', nproc=1, verbose=True, pool=None):
    """Simple maximum likelihood fitting.

    """
    import multiprocessing
    from scipy.optimize import least_squares
    from prospect import fitting

    if verbose:
        print('Minimization method: {}, number of initial tries: {}, number of cores: {}.'.format(
            min_method, nmin, nproc))

    # Start minimization from "nmin" separate places, the first based on
    # the "init" values of each parameter and the rest drawn from the prior.
    qinit = fitting.minimizer_ball(model.theta.copy(), nmin,
                                   model, seed=seed)

    pool = multiprocessing.Pool(nproc)

    min_opts = {'xtol': xtol, 'ftol': ftol, 'max_nfev': maxnfev}

    t0 = time.time()
    with np.errstate(invalid='ignore'):
        res = list(pool.map([least_squares(chivecfn, np.array(qq), method=min_method,
                                           x_scale='jac', **min_opts,
                                           args=(model, obs, sps, verbose)) for qq in qinit]))
    dt = time.time() - t0

    pool.close()
    
    chisq = [np.sum(rr.fun**2) for rr in res]
    chi2min = np.argmin(chisq)
    theta_best = fitting.reinitialize(res[best].x, model)

    if verbose:
        print('{} minimization took {:.3f} seconds.'.format(min_method, dt))
        print('Chi2min = {:.3f} '.format(chi2min))
        print('Best-fitting paramaters:')
        for k, t in zip(model.theta_labels(), theta_best):
            print('  {} = {}'.format(k, t))

    return theta_best, chi2min

def run_dynesty(run_params, model, nested=False):

    from prospect import fitting
    from dynesty.dynamicsampler import stopping_function, weight_function

    def prior_transform(u):        
        return model.prior_transform(u)

    t0 = time.time()  # time it


    out = fitting.run_dynesty(obs, model, sps, lnprobfn=lnprobfn, 
                              stop_function=stopping_function,
                              wt_function=weight_function,
                              logl_args=(model, obs, sps, nested, verbose),
                              **run_params)

    print('done sampling in {:.3f} min'.format((time.time() - t0)/60))

    return out
