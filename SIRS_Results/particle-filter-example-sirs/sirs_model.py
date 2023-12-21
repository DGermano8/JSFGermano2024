import scipy.stats
import numpy as np
import pypfilt                  # type: ignore
from pypfilt.model import Model
from pypfilt.obs import Univariate, Obs
import pdb
import JSF_Solver_BasePython as JSF
from pypfilt.io import time_field


# --------------------------------------------------------------------
# Define the process models
#
# - SIS_ODE :: ODE
# - SIS_CTMC :: CTMC
# - SIS_Hybrid :: Hybrid
# --------------------------------------------------------------------
class SIRS_XXX(Model):

    def field_types(self, ctx):
        """
        """
        return [('S', np.float_),
                ('I', np.float_),
                ('R', np.float_),
                ('N', np.float_),
                ('betaCoef', np.float_),
                ('gammaCoef', np.float_),
                ('omegaCoef', np.float_),
                ('muCoef', np.float_),
                ('kappaCoef', np.float_),
                ('intervention_time', np.float_),
                ('intervention_magnitude', np.float_)]

    def init(self, ctx, vec):
        """
        """
        prior = ctx.data['prior']
        self.num_particles = prior['I'].shape[0]

        vec['S'] = prior['S']
        vec['I'] = prior['I']
        vec['R'] = prior['R']
        vec['N'] = prior['S'] + prior['I'] + prior['R']
        vec['betaCoef'] = prior['betaCoef']
        vec['gammaCoef'] = prior['gammaCoef']
        vec['omegaCoef'] = prior['omegaCoef']
        vec['muCoef'] = prior['muCoef']
        vec['kappaCoef'] = prior['kappaCoef']
        vec['intervention_time'] = prior['intervention_time']
        vec['intervention_magnitude'] = prior['intervention_magnitude']
    
    def can_smooth(self):
        """
        """
        return {'betaCoef', 'gammaCoef', 'omegaCoef', 'muCoef', 'kappaCoef'}
    
       

class SIRS_ODE(SIRS_XXX):
    """
    """
    
    def update(self, ctx, time_step, is_forecast, prev, curr):
        """
        (Destructively) update the state vector `curr`.
        """
        deriv = np.zeros(curr['betaCoef'].shape)

        # only print time modulo 1.0
        if( time_step.start % 1.0 == 0.0 ):
            print('time = ', time_step.start, end='\r')
        # print('time = ', time_step.start)
    
        for p_ix in range(self.num_particles):
            curr['intervention_time'][p_ix] = prev['intervention_time'][p_ix]
            curr['intervention_magnitude'][p_ix] = prev['intervention_magnitude'][p_ix]
            
            beta_scaling_factor = 1.0
            if(time_step.start >= curr['intervention_time'][p_ix]):
                beta_scaling_factor = curr['intervention_magnitude'][p_ix]
            
            # print('time = ', time_step.start, 'beta_fact = ', beta_scaling_factor, 'fact_start = ', curr['intervention_time'][p_ix] ,end='\r')

            curr['betaCoef'][p_ix] = prev['betaCoef'][p_ix]
            curr['gammaCoef'][p_ix] = prev['gammaCoef'][p_ix]
            curr['omegaCoef'][p_ix] = prev['omegaCoef'][p_ix]
            curr['muCoef'][p_ix] = prev['muCoef'][p_ix]
            curr['kappaCoef'][p_ix] = prev['kappaCoef'][p_ix]

            # print('partile ', p_ix, ' state (', prev['S'][p_ix], prev['I'][p_ix], prev['R'][p_ix], ')')
            # print('partile ', p_ix, ' params (', prev['betaCoef'][p_ix], prev['gammaCoef'][p_ix], prev['omegaCoef'][p_ix], prev['muCoef'][p_ix], prev['kappaCoef'][p_ix], ')')

            # Derivative for the SIRS model
            dSdt = -(beta_scaling_factor*prev['betaCoef'][p_ix]/7.0) * prev['S'][p_ix] * prev['I'][p_ix] / prev['N'][p_ix] +(1.0/(365.0*prev['omegaCoef'][p_ix])) * prev['R'][p_ix] - (1.0/(365.0*prev['muCoef'][p_ix])) * prev['S'][p_ix] + (1.0/(365.0*prev['kappaCoef'][p_ix]))* prev['N'][p_ix]
            dIdt = (beta_scaling_factor*prev['betaCoef'][p_ix]/7.0) * prev['S'][p_ix] * prev['I'][p_ix] / prev['N'][p_ix] - (prev['gammaCoef'][p_ix]/7.0) * prev['I'][p_ix] - (1.0/(365.0*prev['muCoef'][p_ix])) * prev['I'][p_ix]
            dRdt = (prev['gammaCoef'][p_ix]/7.0) * prev['I'][p_ix] - (1.0/(365.0*prev['omegaCoef'][p_ix])) * prev['R'][p_ix] - (1.0/(365.0*prev['muCoef'][p_ix])) * prev['R'][p_ix]
            
            curr['S'][p_ix] = prev['S'][p_ix] + time_step.dt * dSdt
            curr['I'][p_ix] = prev['I'][p_ix] + time_step.dt * dIdt
            curr['R'][p_ix] = prev['R'][p_ix] + time_step.dt * dRdt
            curr['N'][p_ix] = curr['S'][p_ix] + curr['I'][p_ix] + curr['R'][p_ix]

            # if(curr['R'][p_ix] < 0.0):
            #     pdb.set_trace()
            assert curr['S'][p_ix] >= 0.0
            assert curr['I'][p_ix] >= 0.0
            assert curr['R'][p_ix] >= 0.0

    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 10**-10
        return np.ceil(hist['I'] < extinction_thresh)

            
    

class SIRS_JSF(SIRS_XXX):
    """
    """

    # num_particles = -1
    threshold = 10**3
    _nu_reactants = [[1, 1, 0],
                     [0, 1, 0],
                     [1, 0, 0],
                     [0, 1, 0],
                     [0, 0, 1],
                     [1, 0, 0],
                     [0, 1, 0],
                     [0, 0, 1],
                     [0, 0, 1]]
    _nu_products = [[0,2,0],
                    [0,0,1],
                    [2,0,0],
                    [1,1,0],
                    [1,0,1],
                    [0,0,0],
                    [0,0,0],
                    [0,0,0],
                    [1,0,0]]
    _nu = [[a - b for a, b in zip(r1, r2)]
           for r1, r2 in zip(_nu_products, _nu_reactants)]
    _stoich = {'nu': _nu,
               'DoDisc': [0, 0, 0],
               'nuReactant': _nu_reactants,
               'nuProduct': _nu_products}

    def _rates(self, x, theta, time):
        """
        """
        s = x[0]
        i = x[1]
        r = x[2]
        m_beta = theta[0]
        m_gamma = theta[1]
        m_wane = theta[2]
        m_mu = theta[3]
        m_kappa = theta[4]
        return [(m_beta/7.0)*(s*i)/(s+i+r),
                (m_gamma/7.0)*i,
                (1.0/(365.0*m_kappa))*(s),
                (1.0/(365.0*m_kappa))*(i),
                (1.0/(365.0*m_kappa))*(r),
                (1.0/(365.0*m_mu))*(s),
                (1.0/(365.0*m_mu))*(i),
                (1.0/(365.0*m_mu))*(r),
                (1.0/(365.0*m_wane))*(r)]
    

    def update(self, ctx, time_step, is_forecast, prev, curr):

        _my_opts = {'EnforceDo': [0, 0, 0],
                    'dt': 0.01,
                    'SwitchingThreshold': [self.threshold,
                                           self.threshold,
                                           self.threshold]}
        
        print('time = ', time_step.start, end='\r')

        for p_ix in range(self.num_particles):
            
            ptcl = prev[p_ix].copy()

            beta_scaling_factor = 1.0
            if(time_step.start >= ptcl['intervention_time']):
                beta_scaling_factor = ptcl['intervention_magnitude']
            
            # print('time = ', time_step.start, 'beta_fact = ', beta_scaling_factor, 'fact_start = ', curr['intervention_time'][p_ix] ,end='\r')

            x0 = [ptcl['S'], ptcl['I'], ptcl['R']]
            theta = [beta_scaling_factor*ptcl['betaCoef'], ptcl['gammaCoef'], ptcl['omegaCoef'], ptcl['muCoef'], ptcl['kappaCoef']]

            xs, ts = JSF.JumpSwitchFlowSimulator(
                x0,
                lambda x, time: self._rates(x, theta, time),
                self._stoich,
                time_step.dt,
                _my_opts
            )

            
            # print('time = ', time_step.start, 'partile = ', p_ix, ' state = (', xs[0][-1], xs[1][-1], xs[2][-1], ')')

            ptcl['S'] = xs[0][-1]
            ptcl['I'] = xs[1][-1]
            ptcl['R'] = xs[2][-1]
            curr[p_ix] = ptcl

    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 1
        return np.ceil(hist['I'] < extinction_thresh)
      

class BackcastStateCIs(pypfilt.summary.BackcastPredictiveCIs):
    """
    Summary statistic for the smoothing problem.
    """
    def n_rows(self, ctx, forecasting):
        n_obs_models = len(ctx.component['obs'])
        n_backcast_times = ctx.summary_count()
        return len(self._BackcastPredictiveCIs__probs) * n_backcast_times * n_obs_models
    
class PrExtinction(pypfilt.summary.Table):
    """
    predict the probability of extinction at time ext_time
    """
    def field_types(self, ctx, obs_list, name):
        self.__model = ctx.component['model']
        self.__time = ctx.component['time']
        time = time_field('time')
        return [time, ('pr', np.float64)]
    def n_rows(self, ctx, forecasting):
        return ctx.summary_count()
    def add_rows(self, ctx, fs_time, window, insert_fn):
        for snapshot in window:
            mask = self.__model.is_extinct(snapshot.state_vec)
            seeded_weights = snapshot.weights * mask
            insert_fn((snapshot.time, np.sum(seeded_weights)))


# --------------------------------------------------------------------
# Define the observation models
# --------------------------------------------------------------------

class UniformObservation(Univariate):
    """
    Observation without error.
    """
    def distribution(self, ctx, snapshot):
        expected_value = snapshot.state_vec['I']
        return scipy.stats.randint(low=np.round(expected_value),
                                   high=np.round(expected_value+10))

class GaussianStateObservation(Univariate):
    """
    Observation with Gaussian noise.
    """
    def distribution(self, ctx, snapshot):
        expected_value = snapshot.state_vec['I']
        obs = scipy.stats.norm(loc=expected_value, scale=5)
        return obs

class NoisyStateObservation(Univariate):
    """
    Observation with Poisson noise.
    """
    def distribution(self, ctx, snapshot):
        expected_value = snapshot.state_vec['I']
        return scipy.stats.poisson(mu=expected_value)

class NegativeBinomailObservation(Univariate):
    """
    Observation with nagative binomial noise.
    """
    def distribution(self, ctx, snapshot):
        expected_value = snapshot.state_vec['I']
        dispersion_param = 100.0
        variance_value = expected_value+expected_value*expected_value/dispersion_param 

        # if np.any(expected_value == 0):
        #     pdb.set_trace()
        p_val = expected_value/variance_value
        n_val = expected_value*expected_value/(variance_value-expected_value)

        zero_mask = expected_value == 0
        n_val[zero_mask] = 1.0
        p_val[zero_mask] = 1.0
        return scipy.stats.nbinom(n=n_val, p=p_val)

