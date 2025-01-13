from geodesic import geodesic, InitialVelocity
from mpl_toolkits.mplot3d import Axes3D
from scipy.integrate import odeint
import numdifftools as nd
import numpy as np
import pickle
from matplotlib import pyplot as plt,use
#use('tkagg')


class mbam_model:
    def __init__(self,t0,tmax,tlen,taumax,dydt,init_pop,init_params,dump_name="dump.pickle",observed=[],thres = 200.0):
        self.nspecies = len(init_pop)
        self.dydt = dydt
        self.init_pop = init_pop
        self.times = np.linspace(t0,tmax,num=tlen)
        self.taumax = taumax
        self.restart = None
        self.dump_name = dump_name
        self.thres = thres
        
        self.observed = observed 
        if observed == []:
            self.observed = list(range(self.nspecies))

        if type(init_params) == list: 
            self.init_params = init_params
            self.init_vel = InitialVelocity(self.init_params, self.jacobian, self.Avv)
            self.tau0 = []
            self.deltatau = 0
            self.xs0 = []
            print(f"vel:{self.init_vel}")
        elif type(init_params) == str:
            self.restart = self.load_model(init_params)
            self.init_vel = self.restart.geo.vs[-1]
            self.init_params = self.restart.geo.xs[-1]
            self.tau0 = self.restart.geo.ts
            self.deltatau = self.tau0[-1]
            self.xs0 = self.restart.geo.xs
        else:
            print("init_params should be a list or filename of a pickle file")
            return


        self.geo = geodesic(self.ode_model, self.jacobian, self.Avv, len(self.times), len(self.init_params), self.init_params, self.init_vel, atol = 1e-2, rtol = 1e-2, callback = self.callback)  
    
    
    def ode_model_full(self,x):
        thetavals = tuple(x)
        outvals = odeint(self.dydt, self.init_pop, self.times, args=thetavals) # integrate ODE
        return outvals

    def ode_model(self,x):
        thetavals = tuple(x)
        outvals = odeint(self.dydt, self.init_pop, self.times, args=thetavals) # integrate ODE
        aux= outvals[1:,self.observed].T
        return aux.flatten()
        
    # Jacobian
    def jacobian(self,x):
        from FiniteDifference import CD4
        h = 1e-4
        n = len(x)
        vs = np.eye(n)
        aux = []
        for i in range(n):
            aux.append(CD4(self.ode_model, x, vs[:, i], h))
        return np.array(aux).T
    
    # Directional second derivative
    def Avv(self,x,v):
        h = 1e-4
        return (self.ode_model(x + h*v) + self.ode_model(x - h*v) - 2*self.ode_model(x))/h/h
        
    
    def fisher(self,x):
        jj = self.jacobian(x)
        return np.dot(jj.T,jj)
        
    # Callback function used to monitor the geodesic after each step
    def callback(self,geo):
        geo = self.geo
        # Integrate until the norm of the velocity has grown by a factor of 10
        # and print out some diagnotistic along the way
        deltatau = 0
        if self.restart != None:
            deltatau = self.restart.geo.ts[-1]
        print("\x1b[0;32;40mIteration: %i, tau: %f\x1b[0m\n\x1b[0;31;40m|v| = %f\x1b[0m" %(len(self.geo.vs), self.geo.ts[-1]+deltatau, np.linalg.norm(self.geo.vs[-1])))

        print(f"\x1b[0;33;40meigen= {self.geo.eigen}\x1b[0m")
        print(f"\x1b[0;35;40mvelocity= {self.geo.vs[-1]}\x1b[0m")
        print(f"params= {self.geo.xs[-1]}")
        return np.linalg.norm(self.geo.vs[-1]) < self.thres
    
    def integrate(self):
        self.geo.integrate(self.taumax)
        if len(self.xs0) != 0:
            self.geo.ts = np.concatenate((self.tau0,self.geo.ts+self.deltatau))
            self.geo.xs = np.concatenate((self.xs0,self.geo.xs))
        #self.geo.ts +=self.deltatau
        self.dump_model(self.dump_name)
    
    def plot_data(self,sol,sym='-',obs=None):
        plt.interactive(True) 
        plt.show() 
        col = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]*3                                                
        if obs == None:
            obs = range(len(sol[0]))
        if sym == '-':
            for i in range(len(obs)):
                #for i in obs:
                plt.plot(self.times,sol[:,obs[i]],'-',color=col[i])
        else:
            mark = ["o","v","^","<",">","s","p","P","*","X","D","d"]*3
            for i in range(len(obs)):
                #for i in obs:
                plt.plot(self.times,sol[:,obs[i]],sym,color=col[i],marker=mark[i])

        plt.xlabel("Time")
        plt.ylabel("Populations")

    def plot_sol(self,x,sym='-',obs=None):
        plt.interactive(True) 
        plt.show() 
        sol = self.ode_model_full(x)
        col = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]*3                                                
        if obs == None:
            obs = range(len(sol[0]))
        if sym == '-':
            for i in range(len(obs)):
            #for i in obs:
                plt.plot(self.times,sol[:,obs[i]],'-',color=col[i])
        else:
            mark = ["o","v","^","<",">","s","p","P","*","X","D","d"]*3
            for i in range(len(obs)):
            #for i in obs:
                plt.plot(self.times,sol[:,obs[i]],sym,color=col[i],marker=mark[i])                                                                     
        plt.xlabel("Time")
        plt.ylabel("Populations")

    def plot_geo(self,legend=""):
        plt.interactive(True) 
        plt.show() 
        mark = ["o","v","^","<",">","s","p","P","*","X","D","d"]*3
        for i in range(len(self.init_params)): 
            plt.plot(self.geo.ts,self.geo.xs[:,i],marker=mark[i])
            if legend=="":
                plt.legend([str(i) for i in range(len(self.init_params))],loc="lower left")
            else:
                plt.legend(legend,loc="lower left")
            plt.xlabel(r"$\tau$")
            plt.ylabel("Parameter Values")
    
    def plot_all(self,figure=None): 
        plt.subplot(1,2,1) 
        self.plot_geo() 
        plt.subplot(1,2,2) 
        self.plot_sol(self.init_params) 
        self.plot_sol(self.geo.xs[-1],"-")
        plt.xlabel("time")
        plt.ylabel("Species")
        if figure != None:
            plt.savefig(figure)

    def dump_model(self,filename):
        pickle.dump(self, open(filename, "wb")) 
        
    def load_model(self,filename):
        return pickle.load(open(filename, "rb"))

