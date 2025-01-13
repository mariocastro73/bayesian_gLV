import matplotlib.pyplot as plt
from numpy.linalg import norm

def plot_model(model,data,str_params,str_pop,filename,observables,geo=True):
    plt.close()
    plt.rcParams["figure.figsize"] = (16,9)
    plt.interactive(True)
    if geo is True:
        plt.subplot(231)
        model.plot_geo(str_params)
        plt.title("Geodesics")

    if geo is True:
        plt.subplot(232)
    else:
        plt.subplot(221)
    model.plot_data(data,'o')
    model.plot_sol(model.geo.xs[-1],'-')
    plt.legend(str_pop)#,loc="center left")
    plt.title("Data vs Model solution at the MB")
    if geo is True:
        plt.subplot(233)
    else:
        plt.subplot(222)
    model.plot_data(data,'o',obs=observables)
    model.plot_sol(model.geo.xs[-1],'-',obs=observables)
    plt.legend([str_pop[i] for i in observables],loc="upper left")
    plt.title("Data vs Model solution at the MB (observables)")
    if geo is True:
        plt.subplot(234)
        frame1 = plt.gca() 
        frame1.axes.xaxis.set_ticklabels([])
        l0 = [model.geo.eigen0,model.geo.eigen0]
        l = [model.geo.eigen[:-1],model.geo.eigen[:-1]] 
        plt.semilogy([0,1],l0,'b-')
        plt.semilogy([2,3],l,'b-')
        plt.subplot(235)
    else:
        plt.subplot(223)
    g = model.geo
    xaxis = [i for i in range(len(g.xs[-1]))]
    plt.bar(xaxis,g.vs[0]/norm(g.vs[0]),width=0.4)
    plt.xlabel("Parameters")
    plt.ylabel("Normalized velocity")
    plt.xticks([i for i in range(len(str_params))],str_params)
    plt.title(r"Geodesic velocity at initial $\tau$")
    if geo is True:
        plt.subplot(236)
    else:
        plt.subplot(224)
    plt.bar(xaxis,g.vs[-1]/norm(g.vs[-1]),width=0.4)
    plt.xlabel("Parameters")
    plt.ylabel("Normalized velocity")
    plt.xticks([i for i in range(len(str_params))],str_params)
    plt.title(r"Geodesic velocity at final $\tau$")

    plt.tight_layout()  
    plt.show()
    plt.savefig(filename)

