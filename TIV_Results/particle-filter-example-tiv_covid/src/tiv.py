import scipy.stats              # type: ignore
import numpy as np
import pypfilt                  # type: ignore
from pypfilt.model import Model # type: ignore
from pypfilt.obs import Univariate, Obs # type: ignore
import pdb
import src.JSF_Solver_BasePython as JSF
import statsmodels.stats.weightstats as smws

from pypfilt.io import time_field

class RefractoryCellModel_JSF(Model):
    """
    """
    num_particles = -1
    threshold = 100
    #    T, R, E, I, V
    _nu_reactants = [
        [1, 0, 0, 0, 1],
        [1, 0, 0, 1, 0],
        [0, 1, 0, 0, 0],
        [0, 0, 1, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, 1]
    ]
    _nu_products = [
        [0, 0, 1, 0, 1],
        [0, 1, 0, 1, 0],
        [1, 0, 0, 0, 0],
        [0, 0, 0, 1, 0],
        [0, 0, 0, 0, 0],
        [0, 0, 0, 1, 1],
        [0, 0, 0, 0, 0]
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
        r = x[1]
        e = x[2]
        i = x[3]
        v = x[4]
        m_beta = theta[0]*10**(-9)
        m_phi = theta[1]*10**(-5)
        m_rho = theta[2]
        m_k = theta[3]
        m_delta = theta[4]
        m_pi = theta[5]
        m_c = theta[6]
        return [m_beta*(t*v),
                m_phi*i*t,
                m_rho*r,
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
            ('phi', np.float_),
            ('rho', np.float_),
            ('k', np.float_),
            ('delta', np.float_),
            ('pi', np.float_),
            ('c', np.float_),
            ('T', np.float_),
            ('E', np.float_),
            ('I', np.float_),
            ('R', np.float_),
            ('V', np.float_),
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
        vec['T'] = prior['T0']
        vec['E'] = prior['E0']
        vec['I'] = prior['I0']
        vec['R'] = prior['R0']
        vec['V'] = np.round(np.exp(vec['lnV0']))

    def update(self, ctx, time_step, is_forecast, prev, curr):
        """
        """
        assert time_step.dt==1
        _my_opts = {'EnforceDo': [0, 0, 0, 0, 0],
                    'dt': 0.00005,
                    'SwitchingThreshold': [self.threshold,
                                           self.threshold,
                                           self.threshold,
                                           self.threshold,
                                           self.threshold]}
        for p_ix in range(self.num_particles):

            ptcl = prev[p_ix].copy()
            x0 = [ptcl['T'], ptcl['E'], ptcl['I'], ptcl['R'],  ptcl['V']]
            theta = [
                ptcl['beta'],
                ptcl['phi'],
                ptcl['rho'],
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
            # print clean line
            print('\033[2K', end='\r')

            print('time = ', time_step.start, ' partile = ', p_ix, end='\r')

            # check if any of the xs are nans at the end of the simulation
            if (np.isnan(xs[0][-1]) or np.isnan(xs[1][-1]) or np.isnan(xs[2][-1]) or np.isnan(xs[3][-1]) or np.isnan(xs[4][-1])):
                pdb.set_trace()
            
            # print('time = ', time_step.start, ' partile = ', p_ix)
            

            # if (xs[0][-1] != np.round(xs[0][-1]) and xs[0][-1] <= (self.threshold-1)) or (xs[1][-1] != np.round(xs[1][-1]) and xs[1][-1] <= (self.threshold-1)) or (xs[2][-1] != np.round(xs[2][-1]) and xs[2][-1] <= (self.threshold-1)) or (xs[3][-1] != np.round(xs[3][-1]) and xs[3][-1] <= (self.threshold-1)):
            #     pdb.set_trace()
            # # if any of the states are negative, go into debug mode
            # if (xs[0][-1] < 0 or xs[1][-1] < 0 or xs[2][-1] < 0 or xs[3][-1] < 0):
            #     pdb.set_trace()

            ptcl['T'] = xs[0][-1]
            ptcl['R'] = xs[1][-1]
            ptcl['E'] = xs[2][-1]
            ptcl['I'] = xs[3][-1]
            ptcl['V'] = xs[4][-1]
            # print('T = ' + str(ptcl['T']) + ' R = ' + str(ptcl['R']) + ' E = ' + str(ptcl['E']) + ' I = ' + str(ptcl['I']) + ' V = ' + str(ptcl['V']) )
            # print('Beta = ' + str(ptcl['beta']) + ' phi = ' + str(ptcl['phi']) + ' rho = ' + str(ptcl['rho']) + ' delta = ' + str(ptcl['delta']) + ' pi = ' + str(ptcl['pi']))
            curr[p_ix] = ptcl
    def can_smooth(self):
        """
        """
        return {'lnV0', 'beta', 'phi', 'rho', 'k', 'delta', 'pi', 'c'}
    
    def is_extinct(self, hist):
        """
        identify which of the particles has gone extinct
        """
        extinction_thresh = 1
        return np.ceil( (hist['V'] < extinction_thresh)*(hist['I'] < extinction_thresh)*(hist['E'] < extinction_thresh) )


# class TIV_ODE(Model):
#     """
#     """
#     def field_types(self, ctx):
#         """
#         """
#         return [
#             ('lnV0', np.float_),
#             ('beta', np.float_),
#             ('p', np.float_),
#             ('c', np.float_),
#             ('gamma', np.float_),
#             ('T', np.float_),
#             ('I', np.float_),
#             ('V', np.float_),
#         ]

#     def init(self, ctx, vec):
#         """
#         """
#         prior = ctx.data['prior']

#         vec['beta'] = prior['beta']
#         vec['p'] = prior['p']
#         vec['lnV0'] = prior['lnV0']
#         vec['c'] = prior['c']
#         vec['gamma'] = prior['gamma']
#         vec['T'] = prior['T0']
#         vec['I'] = prior['I0']
#         vec['V'] = np.exp(vec['lnV0'])


#     def update(self, ctx, time_step, is_forecast, prev, curr):
#         """
#         """
#         curr['beta'] = prev['beta']
#         curr['p'] = prev['p']
#         curr['lnV0'] = prev['lnV0']
#         curr['c'] = prev['c']
#         curr['gamma'] = prev['gamma']

#         dT_dt = - prev['beta'] * prev['V'] * prev['T']
#         dI_dt = prev['beta'] * prev['V'] * prev['T'] - prev['gamma'] * prev['I']
#         dV_dt = prev['p'] * prev['I'] - prev['c'] * prev['V']

#         curr['T'] = prev['T'] + time_step.dt * dT_dt
#         curr['I'] = prev['I'] + time_step.dt * dI_dt
#         curr['V'] = prev['V'] + time_step.dt * dV_dt

#     def can_smooth(self):
#         """
#         """
#         return {'lnV0', 'beta', 'p', 'c', 'gamma'}
    
#     def is_extinct(self, hist):
#         """
#         identify which of the particles has gone extinct
#         """
#         extinction_thresh = 10**(-32)
#         return np.ceil(hist['V'] < extinction_thresh and hist['I'] < extinction_thresh)


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
    Measurement model for Gaussian measurements with truncation at the detection limit stored in a class variable.
    """

    limitOfDetection = -0.65

    def quantiles(self, ctx, snapshot, probs):
        """
        The quantiles of the observations adjusting for the limit of detection as described in the `distribution method.
        """
        # pdb.set_trace()       
        zero_mask = snapshot.state_vec[self.unit] == 0

        my_y_vals = np.zeros_like(snapshot.state_vec[self.unit])
        # my_y_vals[zero_mask] = -np.inf
        my_y_vals[np.logical_not(zero_mask)] = np.log10(snapshot.state_vec[self.unit][np.logical_not(zero_mask)])
        my_y_vals[my_y_vals <= self.limitOfDetection ] = self.limitOfDetection
        
        my_weighted_sample = smws.DescrStatsW(my_y_vals, weights = snapshot.weights)
        
        # my_y_vals[np.logical_not(zero_mask)] = (snapshot.state_vec[self.unit][np.logical_not(zero_mask)])

        # return (np.quantile(my_y_vals, probs))
        return (my_weighted_sample.quantile(probs, return_pandas=False))

    def distribution(self, ctx, snapshot):
        """
        For particles where there is virus present, we take the log-base-10 and use that as the mean of the normal distribution with a scale which is read from the TOML.
        For particles where there is no virus present, we set the mean to the limit of detection and the scale to a very small number to indicate that the measurement device should return that value.
        """
        zero_mask = snapshot.state_vec[self.unit] == 0

        loc = np.zeros_like(snapshot.state_vec[self.unit])
        # loc[zero_mask] = -np.inf
        loc[np.logical_not(zero_mask)] = np.log10(snapshot.state_vec[self.unit][np.logical_not(zero_mask)])
        loc[loc <= self.limitOfDetection ] = self.limitOfDetection
        trunc_mask = loc == self.limitOfDetection
        scale = np.repeat(ctx.settings['observations'][self.unit]['scale'], len(loc))
        scale[trunc_mask] = 10**(-32)
        # print('loc = ',  loc)
        #loc = np.log10(snapshot.state_vec[self.unit]) 
        # scale = ctx.settings['observations'][self.unit]['scale']
        return scipy.stats.norm(loc=loc, scale=scale)
