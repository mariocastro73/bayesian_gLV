from mbam import mbam_model
from models import *
from plot_model import *
import pickle
import fit_model as fm
#import fit_model_only_params as fm
from scipy.integrate import odeint
from numpy import log, abs, exp, array

# original parameters
init_params_lin = [ 1.00, 0.7, 1.5, 1.3, # r1..r4
                   -1.0 ,-0.5 ,-3.6 ,-1.5,
                   -1.2 ,-0.7 ,-0.1 ,-0.7,
                   -1.5 ,-0.3 ,-1.5 ,-0.4,
                   -0.2 ,-1.0 ,-0.7 ,-1.3
                   ]

init_params = [log(abs(par)) for par in init_params_lin]
str_params = [r"$r_{1}$", r"$r_{2}$", r"$r_{3}$", r"$r_{4}$", r"$b_{11}$", r"$b_{21}$", r"$b_{31}$", r"$b_{41}$", r"$b_{12}$", r"$b_{22}$", r"$b_{32}$", r"$b_{42}$", r"$b_{13}$", r"$b_{23}$",r"$b_{33}$",r"$b_{43}$", r"$b_{14}$", r"$b_{24}$",r"$b_{34}$",r"$b_{44}$"]

# Initial conditions
init_pop = [.05,.05,.05,.05]
str_pop = [r"$x_1$", r"$x_2$", r"$x_3$", r"$x_4$"]
# Observed quantities (I)
observables = [0,1,2,3] # All

taumax = 50 # Geodesic "intrinsic time", it's arbitrary as the code stops earlier
velthresh = 30. # Transtrum et al. usually fix this to 20. 30 is more conservative.
tmax = 30 # Simulation time
dt = tmax # Number of steps
niter = 12

time = np.linspace(0,tmax,dt)
data = odeint(fourpop_0_exp_model,init_pop,time,tuple(init_params))
removed_params = []

#for order in range(niter): # All iterations at once, but with the reductions already known, for figures
for order in [0]:  # put the order you want to explore here
    if order != 0:
        init_params, str_params, removed_params = pickle.load(open(f'fourpop_{order-1}_exp_params.pickle',"rb"))
        print(f"Initial conditions\n{init_pop}")
        print("Fitting....")
        init_params, init_pop = fm.fit(eval(f"fourpop_{order}_exp_model"),time,data,observables,init_params)

    print(f"Initial conditions\n{init_pop}")
    model = mbam_model(0,tmax,dt,taumax,eval(f"fourpop_{order}_exp_model"),init_pop,init_params,
                       dump_name=f"fourpop_{order}_exp_model.pickle",observed = observables,thres=30.) 
    model.integrate()
    
    filename = f'fourpop_{order}_exp_figure.pdf'
    plot_model(model, data, str_params, str_pop, filename, observables,geo=True )
    
    init_params = model.geo.xs[-1].tolist()
    vels = abs(model.geo.vs[-1]).tolist()
    index = vels.index(max(vels))
    removed_params.append(index)
    print("#"*10)
    print(f"Remove variable {str_params[index]}")

    init_params.pop(index)
    str_params.pop(index)
    pickle.dump([init_params, str_params, removed_params],open(f"fourpop_{order}_exp_params.pickle","wb"))

    print(removed_params)
    print("#"*10)
