"""
hizea.prospector
================

Code to do SED-fitting with Prospector.

"""
import os, time
import numpy as np

def load_obs(snr=10):
    """
    Load the photometry for J2118.
    
    From Christy:
      W1      14.8739    0.0154054
      W2      14.1342    0.0283280
      W3      10.4872    0.0699660
      W4      8.02651     0.179578
    
    """
    import sedpy
    from prospect.utils.obsutils import fix_obs    

    phot = dict(
        FUV=(1.82491e-09, 1.15115e+19),
        NUV=(8.06441e-09, 2.25963e+19),
        u=(1.27666e-08, 1.53319e+18),
        g=(1.59991e-08, 7.47418e+18),
        r=(3.15573e-08, 2.18797e+18),
        i=(3.63049e-08, 1.53877e+18),
        z=(4.14564e-08, 2.71207e+17),
        ch1=(9.25971e-08, 3.91873e+17),
        ch2=(9.67009e-08, 3.54276e+17),
        W1=(9.40651e-08, 5.61366e+17),
        W2=(1.02882e-07, 1.38784e+17),
        W3=(5.44324e-07, 8.12757e+14),
        W4=(1.38524e-06, 1.90498e+13))    

    galex = ['galex_FUV', 'galex_NUV']
    sdss = ['sdss_{}0'.format(b) for b in ['u','g','r','i','z']]
    spitzer = ['spitzer_irac_ch{}'.format(n) for n in ['1','2']]
    wise = ['wise_w{}'.format(n) for n in ['1','2', '3', '4']]
    filternames = galex + sdss + spitzer + wise

    obs = {}
    obs['redshift'] = 0.535
    obs["filters"] = sedpy.observate.load_filters(filternames)

    obs["maggies"] = np.array([phot[filt][0] for filt in phot.keys()])
    obs["maggies_unc"] = np.array([1/np.sqrt(phot[filt][1]) for filt in phot.keys()])

    # mask out W4
    #obs["phot_mask"] = np.array(['w4' in f.name for f in obs["filters"]])    
    
    # Create a handy vector of effective wavelengths (optional) 
    obs["phot_wave"] = [f.wave_effective for f in obs["filters"]]
    obs["wavelength"] = None # spectral wavelength
    obs["spectrum"] = None
    obs['unc'] = None  # spectral uncertainties are given here
    obs['mask'] = None
    obs = fix_obs(obs)

    run_params = {}
    run_params['redshift'] = obs['redshift']
    run_params["verbose"] = True

    run_params["nmin"] = 5
    run_params['ftol'] = 3e-16 
    run_params['maxfev'] = 500
    run_params['xtol'] = 3e-16
    run_params["min_method"] = 'levenberg_marquardt'
    
    return obs, run_params

def load_sps(zcontinuous=1):
    """zcontinuous - interpolate between metallicity values.

    """
    from prospect.sources import CSPSpecBasis
    sps = CSPSpecBasis(zcontinuous=zcontinuous)
    return sps

def load_model(template_library='delayed-tau', redshift=0.0):
    
    from prospect.models import priors
    from prospect.models.sedmodel import SedModel
    from prospect.models.templates import TemplateLibrary

    if template_library == 'delayed-tau':
        model_params = TemplateLibrary["parametric_sfh"]

        # Initial values
        model_params["zred"]['init'] = redshift
        model_params["mass"]["init"] = 3e11
        model_params["tau"]["init"] = 10.0
        model_params["tage"]["init"] = 1.0
        model_params["dust2"]["init"] = 1.0
        model_params["logzsol"]['init'] = 0.1
        
        model_params["zred"]['isfree'] = False # fixed redshift

        # Prior ranges
        model_params['mass']['prior'] = priors.LogUniform(mini=1e10, maxi=1e12)
        model_params['tau']['prior'] = priors.LogUniform(mini=1, maxi=1e2)
        model_params['tage']['prior'] = priors.LogUniform(mini=0.001, maxi=10)
        model_params['dust2']['prior'] = priors.TopHat(mini=0.0, maxi=2.0)
        model_params['logzsol']['prior'] = priors.TopHat(mini=-1.0, maxi=0.3)

        # If we are going to be using emcee, it is useful to provide a 
        # minimum scale for the cloud of walkers (the default is 0.1)
        model_params["mass"]["disp_floor"] = 1e9
        model_params["tau"]["disp_floor"] = 1.0
        model_params["tage"]["disp_floor"] = 0.05
        model_params['dust2']['disp_floor'] = 0.1
        model_params['logzsol']['disp_floor'] = 0.05

        # Add dust emission (with fixed dust SED parameters)
        model_params.update(TemplateLibrary["dust_emission"])
        
    # Now instantiate the model using this new dictionary of parameter specifications
    model = SedModel(model_params)

    return model

def lnprobfn(theta, nested=False, verbose=False):
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

def max_likelihood(run_params, model, obs, sps, seed=1, verbose=True):
    """Simple maximum likelihood fitting.

    """
    from prospect import fitting
    from scipy.optimize import least_squares

    # We'll start minimization from "nmin" separate places, 
    # the first based on the "init" values of each parameter and the 
    # rest drawn from the prior.  This can guard against local minima.
    ts = time.time()  # time it
    pinitial = fitting.minimizer_ball(model.initial_theta.copy(),
                                      run_params['nmin'], model,
                                      seed=seed)
    guesses = []
    for i, pinit in enumerate(pinitial): #loop over initial guesses
        res = least_squares(chivecfn, np.array(pinit), 
                            method='dogbox', x_scale='jac',
                            xtol=run_params["xtol"], ftol=run_params["ftol"], 
                            max_nfev=run_params["maxfev"],
                            args=(model, obs, sps, verbose))
        guesses.append(res)

    # Calculate chi-square of the results, and choose the best one
    # fitting.reinitialize moves the parameter vector away from edges of the prior.
    chisq = [np.sum(r.fun**2) for r in guesses]
    best = np.argmin(chisq)
    theta_best = fitting.reinitialize(guesses[best].x, model,
                                      edge_trunc=run_params.get('edge_trunc', 0.1))
    initial_prob = None
    pdur = time.time() - ts

    # output results
    if verbose:
        print('done {0} in {1}s'.format(run_params['min_method'], pdur))
        print('best {0} chi-sq: {1}'.format(run_params['min_method'], chisq[best]))
        print('best guess paramaters:')
        for k, t in zip(model.theta_labels(), theta_best):
            print('  {} = {}'.format(k, t))

    return theta_best

