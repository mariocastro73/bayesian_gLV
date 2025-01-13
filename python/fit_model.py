import numpy as np
import scipy.optimize
from scipy.integrate import odeint

def loss_function(params, model, target, y0, time, obs):
    params = list(params)
    y0 = target[0]
    hidden = [i for i in range(len(y0)) if i not in obs]
    hidden.reverse()
    for h in hidden:
        y0[h] = np.exp(params.pop())

    output = odeint(model, y0, time, args=tuple(params))
    return np.sum(np.sum((output[:,obs]-target[:,obs])**2))  

def fit(model, time, target, obs, init_params):
    y0 = target[0]

    hidden = [i for i in range(len(y0)) if i not in obs]
    for h in hidden:
        init_params.append(np.log(1e-16+y0[h]))
    args = [model, target, y0, time, obs]
    minimum = scipy.optimize.fmin(loss_function, init_params, args=tuple(args)).tolist()

    hidden.reverse()
    for h in hidden:
        y0[h] = np.exp(minimum.pop())

    print(minimum)
    return minimum, y0

