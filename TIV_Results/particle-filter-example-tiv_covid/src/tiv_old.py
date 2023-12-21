import scipy.stats              # type: ignore
import numpy as np
import pypfilt                  # type: ignore
from pypfilt.model import Model # type: ignore
from pypfilt.obs import Univariate, Obs # type: ignore
import pdb
import src.JSF_Solver_BasePython as JSF
import statsmodels.stats.weightstats as smws

from pypfilt.io import time_field

class RefractoryCellFullModel_JSF(Model):
    """
    """
    num_particles = -1
    threshold = 50
    #    T, E, I, R, V, F
    _nu_reactants = [
        [1, 0, 0, 0, 1, 0],
        [1, 0, 0, 0, 0, 1],
        [0, 0, 0, 1, 0, 0],
        [0, 1, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 1, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 1]
    ]
    _nu_products = [
        [0, 1, 0, 0, 1, 0],
        [0, 0, 0, 1, 0, 1],
        [1, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 1, 0],
        [0, 0, 0, 0, 0, 0],
        [0, 0, 1, 0, 0, 1],
        [0, 0, 0, 0, 0, 0]
    ]

    _nu = [[a - b for a, b in zip(r1, r2)]
           for r1, r2 in zip(_nu_products, _nu_reactants)]
    _stoich = {'nu': _nu,
               'DoDisc': [0, 0, 0, 0, 0, 0],
               'nuReactant': _nu_reactants,
               'nuProduct': _nu_products}
    def _rates(self, x, theta, time):
        """
        """
        t = x[0]
        e = x[1]
        i = x[2]
        r = x[3]
        v = x[4]
        f = x[5]
        m_beta = theta[0]
        m_phi = theta[1]
        m_rho = theta[2]
        m_k = theta[3]
        m_delta = theta[4]
        m_pi = theta[5]
        m_c = theta[6]
        m_s = theta[7]
        m_mu = theta[8]
        return [m_beta*(t*v),
                m_phi*f*t,
                m_rho*r,
                m_k*e,
                m_delta*i,
                m_pi*i,
                m_c*v,
                m_s*i,
                m_mu*f]

    def field_types(self, ctx):
        """
        """
        return [
            ('lnV0', np.float_),
            ('beta', np.float_),
            ('phi', np.float_),
            ('rho', np.float_),
            ('k', np.float_),
            ('delta', np.float_),
            ('pi', np.float_),
            ('c', np.float_),
            ('s', np.float_),
            ('mu', np.float_),
            ('T', np.float_),
            ('E', np.float_),
            ('I', np.float_),
            ('R', np.float_),
            ('V', np.float_),
            ('F', np.float_),
            ('next_event', np.int_),
            ('next_time', np.float_),
        ]
    def init(self, ctx, vec):
        """
        """
        prior = ctx.data['prior']
        self.num_particles = prior['T0'].shape[0]
        
        vec['lnV0'] = prior['lnV0']
        vec['beta'] = prior['beta']
        vec['phi'] = prior['phi']
        vec['rho'] = prior['rho']
        vec['k'] = prior['k']
        vec['delta'] = prior['delta']
        vec['pi'] = prior['pi']
        vec['c'] = prior['c']
        vec['s'] = prior['s']
        vec['mu'] = prior['mu']
        vec['T'] = prior['T0']
        vec['E'] = prior['E0']
        vec['I'] = prior['I0']
        vec['R'] = prior['R0']
        vec['V'] = np.round(np.exp(vec['lnV0']))
        vec['F'] = prior['F0']

    def update(self, ctx, time_step, is_forecast, prev, curr):
        """
        """
        assert time_step.dt==1
        _my_opts = {'EnforceDo': [0, 0, 0, 0, 0, 0],
                    'dt': 0.0001,
                    'SwitchingThreshold': [self.threshold,
                                           self.threshold,
                                           self.threshold,
                                           self.threshold,
                                           self.threshold,
                                           self.threshold]}
        for p_ix in range(self.num_particles):

            ptcl = prev[p_ix].copy()
            x0 = [ptcl['T'], ptcl['E'], ptcl['I'], ptcl['R'],  ptcl['V'], ptcl['F']]
            theta = [
                ptcl['beta'],
                ptcl['phi'],
                ptcl['rho'],
                ptcl['k'],
                ptcl['delta'],
                ptcl['pi'],
                ptcl['c'],
                ptcl['s'],
                ptcl['mu'],
            ]

            # pdb.set_trace()
            xs, ts = JSF.JumpSwitchFlowSimulator(
                x0,
                lambda x, time: self._rates(x, theta, time),
                self._stoich,
                time_step.dt,
                _my_opts
            )
            
            print('time = ', time_step.start, ' partile = ', p_ix, end='\r')

            # if (xs[0][-1] != np.round(xs[0][-1]) and xs[0][-1] <= (self.threshold-1)) or (xs[1][-1] != np.round(xs[1][-1]) and xs[1][-1] <= (self.threshold-1)) or (xs[2][-1] != np.round(xs[2][-1]) and xs[2][-1] <= (self.threshold-1)) or (xs[3][-1] != np.round(xs[3][-1]) and xs[3][-1] <= (self.threshold-1)):
            #     pdb.set_trace()
            # # if any of the states are negative, go into debug mode
            # if (xs[0][-1] < 0 or xs[1][-1] < 0 or xs[2][-1] < 0 or xs[3][-1] < 0):
            #     pdb.set_trace()

            ptcl['T'] = xs[0][-1]
            ptcl['E'] = xs[1][-1]
            ptcl['I'] = xs[2][-1]
            ptcl['R'] = xs[3][-1]
            ptcl['V'] = xs[4][-1]
            ptcl['F'] = xs[5][-1]

            curr[p_ix] = ptcl
    def can_smooth(self):
        """
        """
        return {'lnV0', 'beta', 'phi', 'rho', 'k', 'delta', 'pi', 'c', 's', 'mu'}
    
    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 1
        return np.ceil( (hist['V'] < extinction_thresh)*(hist['I'] < extinction_thresh)*(hist['E'] < extinction_thresh) )




class TargetCellLimited_JSF(Model):
    """
    """

    num_particles = -1
    threshold = 50
    #    T, E, I, V
    _nu_reactants = [
        [1, 0, 0, 1],
        [0, 1, 0, 0],
        [0, 0, 1, 0],
        [0, 0, 1, 0],
        [0, 0, 0, 1]
    ]
    _nu_products = [
        [0, 1, 0, 1],
        [0, 0, 1, 0],
        [0, 0, 0, 0],
        [0, 0, 0, 1],
        [0, 0, 0, 0]
    ]

    _nu = [[a - b for a, b in zip(r1, r2)]
           for r1, r2 in zip(_nu_products, _nu_reactants)]
    _stoich = {'nu': _nu,
               'DoDisc': [0, 0, 0, 0],
               'nuReactant': _nu_reactants,
               'nuProduct': _nu_products}
    
    def _rates(self, x, theta, time):
        """
        """
        t = x[0]
        e = x[1]
        i = x[2]
        v = x[3]
        m_beta = theta[0]
        m_k = theta[1]
        m_delta = theta[2]
        m_pi = theta[3]
        m_c = theta[4]
        return [m_beta*(t*v),
                m_k*e,
                m_delta*i,
                m_pi*i,
                m_c*v]

    def field_types(self, ctx):
        """
        """
        return [
            ('lnV0', np.float_),
            ('beta', np.float_),
            ('k', np.float_),
            ('delta', np.float_),
            ('pi', np.float_),
            ('c', np.float_),
            ('T', np.float_),
            ('E', np.float_),
            ('I', np.float_),
            ('V', np.float_),
            ('next_event', np.int_),
            ('next_time', np.float_),
        ]

    def init(self, ctx, vec):
        """
        """
        prior = ctx.data['prior']
        self.num_particles = prior['T0'].shape[0]

        vec['beta'] = prior['beta']
        vec['k'] = prior['k']
        vec['lnV0'] = prior['lnV0']
        vec['c'] = prior['c']
        vec['delta'] = prior['delta']
        vec['pi'] = prior['pi']
        vec['T'] = prior['T0']
        vec['I'] = prior['I0']
        vec['V'] = np.round(np.exp(vec['lnV0']))

    def update(self, ctx, time_step, is_forecast, prev, curr):
        """
        """
        assert time_step.dt==1
        _my_opts = {'EnforceDo': [0, 0, 0, 0],
                    'dt': 0.001,
                    'SwitchingThreshold': [self.threshold,
                                           self.threshold,
                                           self.threshold,
                                           self.threshold]}

        for p_ix in range(self.num_particles):

            # assert prev[p_ix]['V'] >= 0
            # if prev[p_ix]['V'] < 0:
            #     pdb.set_trace()

            # if prev[p_ix]['V'] != np.round(prev[p_ix]['V']) and prev[p_ix]['V'] < (self.threshold-1):
            #     pdb.set_trace()
            ptcl = prev[p_ix].copy()
            x0 = [ptcl['T'], ptcl['E'], ptcl['I'], ptcl['V']]
            theta = [
                ptcl['beta'],
                ptcl['k'],
                ptcl['delta'],
                ptcl['pi'],
                ptcl['c'],
            ]

            # pdb.set_trace()
            xs, ts = JSF.JumpSwitchFlowSimulator(
                x0,
                lambda x, time: self._rates(x, theta, time),
                self._stoich,
                time_step.dt,
                _my_opts
            )
            
            # print('partile ', p_ix, ' state (', xs[0][-1], xs[1][-1], xs[2][-1], ')')
            print('time = ', time_step.start, ' partile = ', p_ix, end='\r')

            if (xs[0][-1] != np.round(xs[0][-1]) and xs[0][-1] <= (self.threshold-1)) or (xs[1][-1] != np.round(xs[1][-1]) and xs[1][-1] <= (self.threshold-1)) or (xs[2][-1] != np.round(xs[2][-1]) and xs[2][-1] <= (self.threshold-1)) or (xs[3][-1] != np.round(xs[3][-1]) and xs[3][-1] <= (self.threshold-1)):
                pdb.set_trace()

            # if any of the states are negative, go into debug mode
            if (xs[0][-1] < 0 or xs[1][-1] < 0 or xs[2][-1] < 0 or xs[3][-1] < 0):
                pdb.set_trace()

            ptcl['T'] = xs[0][-1]
            ptcl['E'] = xs[1][-1]
            ptcl['I'] = xs[2][-1]
            ptcl['V'] = xs[3][-1]
            curr[p_ix] = ptcl
    
    def can_smooth(self):
        """
        """
        return {'lnV0', 'beta', 'k', 'delta', 'pi', 'c'}
    
    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 1
        return np.ceil( (hist['V'] < extinction_thresh)*(hist['I'] < extinction_thresh)*(hist['E'] < extinction_thresh) )


class TIV_JSF(Model):
    """
    """

    num_particles = -1
    threshold = 50
    _nu_reactants = [
        [1, 0, 1],
        [0, 1, 0],
        [0, 1, 0],
        [0, 0, 1]
    ]
    _nu_products = [
        [0, 1, 1],
        [0, 0, 0],
        [0, 1, 1],
        [0, 0, 0],
    ]

    _nu = [[a - b for a, b in zip(r1, r2)]
           for r1, r2 in zip(_nu_products, _nu_reactants)]
    _stoich = {'nu': _nu,
               'DoDisc': [0, 0, 0],
               'nuReactant': _nu_reactants,
               'nuProduct': _nu_products}
    
    def _rates(self, x, theta, time):
        """
        """
        t = x[0]
        i = x[1]
        v = x[2]
        m_beta = theta[0]
        m_p = theta[1]
        m_c = theta[2]
        m_gamma = theta[3]
        return [m_beta*(t*v),
                m_gamma*i,
                m_p*i,
                m_c*v]

    def field_types(self, ctx):
        """
        """
        return [
            ('lnV0', np.float_),
            ('beta', np.float_),
            ('p', np.float_),
            ('c', np.float_),
            ('gamma', np.float_),
            ('T', np.float_),
            ('I', np.float_),
            ('V', np.float_),
            ('next_event', np.int_),
            ('next_time', np.float_),
        ]

    def init(self, ctx, vec):
        """
        """
        prior = ctx.data['prior']
        self.num_particles = prior['T0'].shape[0]

        vec['beta'] = prior['beta']
        vec['p'] = prior['p']
        vec['lnV0'] = prior['lnV0']
        vec['c'] = prior['c']
        vec['gamma'] = prior['gamma']
        vec['T'] = prior['T0']
        vec['I'] = prior['I0']
        vec['V'] = np.round(np.exp(vec['lnV0']))

    def update(self, ctx, time_step, is_forecast, prev, curr):
        """
        """
        assert time_step.dt==1
        _my_opts = {'EnforceDo': [0, 0, 0],
                    'dt': 0.001,
                    'SwitchingThreshold': [self.threshold,
                                           self.threshold,
                                           self.threshold]}

        for p_ix in range(self.num_particles):

            # assert prev[p_ix]['V'] >= 0
            # if prev[p_ix]['V'] < 0:
            #     pdb.set_trace()

            # if prev[p_ix]['V'] != np.round(prev[p_ix]['V']) and prev[p_ix]['V'] < (self.threshold-1):
            #     pdb.set_trace()
            ptcl = prev[p_ix].copy()
            x0 = [ptcl['T'], ptcl['I'], ptcl['V']]
            theta = [
                ptcl['beta'],
                ptcl['p'],
                ptcl['c'],
                ptcl['gamma'],
            ]

            # pdb.set_trace()
            xs, ts = JSF.JumpSwitchFlowSimulator(
                x0,
                lambda x, time: self._rates(x, theta, time),
                self._stoich,
                time_step.dt,
                _my_opts
            )
            
            # print('partile ', p_ix, ' state (', xs[0][-1], xs[1][-1], xs[2][-1], ')')
            print('time = ', time_step.start, ' partile = ', p_ix, end='\r')

            if (xs[0][-1] != np.round(xs[0][-1]) and xs[0][-1] <= (self.threshold-1)) or (xs[1][-1] != np.round(xs[1][-1]) and xs[1][-1] <= (self.threshold-1)) or (xs[2][-1] != np.round(xs[2][-1]) and xs[2][-1] <= (self.threshold-1)):
                pdb.set_trace()

            # if any of the states are negative, go into debug mode
            if (xs[0][-1] < 0 or xs[1][-1] < 0 or xs[2][-1] < 0 ):
                pdb.set_trace()

            ptcl['T'] = xs[0][-1]
            ptcl['I'] = xs[1][-1]
            ptcl['V'] = xs[2][-1]
            curr[p_ix] = ptcl
    
    def can_smooth(self):
        """
        """
        return {'lnV0', 'beta', 'p', 'c', 'gamma'}
    
    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 1
        return np.ceil( (hist['V'] < extinction_thresh)*(hist['I'] < extinction_thresh) )

class TIV_ODE(Model):
    """
    """
    def field_types(self, ctx):
        """
        """
        return [
            ('lnV0', np.float_),
            ('beta', np.float_),
            ('p', np.float_),
            ('c', np.float_),
            ('gamma', np.float_),
            ('T', np.float_),
            ('I', np.float_),
            ('V', np.float_),
        ]

    def init(self, ctx, vec):
        """
        """
        prior = ctx.data['prior']

        vec['beta'] = prior['beta']
        vec['p'] = prior['p']
        vec['lnV0'] = prior['lnV0']
        vec['c'] = prior['c']
        vec['gamma'] = prior['gamma']
        vec['T'] = prior['T0']
        vec['I'] = prior['I0']
        vec['V'] = np.exp(vec['lnV0'])


    def update(self, ctx, time_step, is_forecast, prev, curr):
        """
        """
        curr['beta'] = prev['beta']
        curr['p'] = prev['p']
        curr['lnV0'] = prev['lnV0']
        curr['c'] = prev['c']
        curr['gamma'] = prev['gamma']

        dT_dt = - prev['beta'] * prev['V'] * prev['T']
        dI_dt = prev['beta'] * prev['V'] * prev['T'] - prev['gamma'] * prev['I']
        dV_dt = prev['p'] * prev['I'] - prev['c'] * prev['V']

        curr['T'] = prev['T'] + time_step.dt * dT_dt
        curr['I'] = prev['I'] + time_step.dt * dI_dt
        curr['V'] = prev['V'] + time_step.dt * dV_dt

    def can_smooth(self):
        """
        """
        return {'lnV0', 'beta', 'p', 'c', 'gamma'}
    
    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 10**(-32)
        return np.ceil(hist['V'] < extinction_thresh and hist['I'] < extinction_thresh)


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


class PerfectMeasurement(Univariate):
    """
    Measurement model for perfect measurements.
    """

    def distribution(self, ctx, snapshot):
        expect = np.log10(snapshot.state_vec[self.unit])
        return scipy.stats.norm(loc=expect, scale=0.0)

class Gaussian(Univariate):
    """
    Measurement model for Gaussian measurements.
    """

    def quantiles(self, ctx, snapshot, probs):
        # pdb.set_trace()       
        zero_mask = snapshot.state_vec[self.unit] == 0

        my_y_vals = np.zeros_like(snapshot.state_vec[self.unit])
        # my_y_vals[zero_mask] = -np.inf
        my_y_vals[np.logical_not(zero_mask)] = np.log10(snapshot.state_vec[self.unit][np.logical_not(zero_mask)])
        my_y_vals[my_y_vals <= 0.5 ] = 0.5
        
        my_weighted_sample = smws.DescrStatsW(my_y_vals, weights = snapshot.weights)
        
        # my_y_vals[np.logical_not(zero_mask)] = (snapshot.state_vec[self.unit][np.logical_not(zero_mask)])

        # return (np.quantile(my_y_vals, probs))
        return (my_weighted_sample.quantile(probs, return_pandas=False))

    def distribution(self, ctx, snapshot):
        zero_mask = snapshot.state_vec[self.unit] == 0

        loc = np.zeros_like(snapshot.state_vec[self.unit])
        # loc[zero_mask] = -np.inf
        loc[np.logical_not(zero_mask)] = np.log10(snapshot.state_vec[self.unit][np.logical_not(zero_mask)])
        loc[loc <= 0.5 ] = 0.5
        trunc_mask = loc == 0.5
        scale = np.repeat(ctx.settings['observations'][self.unit]['scale'], len(loc))
        scale[trunc_mask] = 10**(-32)
        # print('loc = ',  loc)
        #loc = np.log10(snapshot.state_vec[self.unit]) 
        scale = ctx.settings['observations'][self.unit]['scale']
        return scipy.stats.norm(loc=loc, scale=scale)
