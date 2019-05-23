"""
hizea.prospector
================

Code to do SED-fitting with Prospector.

"""
import os, time
import numpy as np

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
