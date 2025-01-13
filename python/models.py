from matplotlib import pyplot as plt, use
#use('tkagg')
import numpy as np
exp = np.exp


def lotka3(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    x1, x2, x3 = y
    dx1 = x1*(r1+b11*x1+b12*x2+b13*x3)
    dx2 = x2*(r2+b21*x1+b22*x2+b23*x3)
    dx3 = x3*(r3+b31*x1+b32*x2+b33*x3)
    return [dx1, dx2, dx3]

def lotka4(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b32, b42, b13, b23, b33, b43, b14, b24, b34, b44):
    x1, x2, x3, x4 = y
    dx1 = x1*(r1+b11*x1+b12*x2+b13*x3+b14*x3)
    dx2 = x2*(r2+b21*x1+b22*x2+b23*x3+b24*x3)
    dx3 = x3*(r3+b31*x1+b32*x2+b33*x3+b34*x3)
    dx4 = x4*(r4+b41*x1+b42*x2+b43*x3+b44*x4)
    return [dx1, dx2, dx3, dx4]

def simple_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5,0.019,0.9,-0.02,-0.001,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_0_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5,0.019,0.9,-0.02,-0.001,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_1_exp_model(y, t, r1, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5, 0.9,-0.02,-0.001,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(        -exp(b21)*x1-exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_2_exp_model(y, t, r1, r3, b11, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5, 0.9,-0.02,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(                    -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_3_exp_model(y, t, r1, r3, b11, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5, 0.9,-0.02,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(                    -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)            +exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_4_exp_model(y, t, r1, r3, b11, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5, 0.9,-0.02,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1            +exp(b13)*x3)
    dx2 = x2*(                    -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)            +exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_5_exp_model(y, t, r1, r3, b11, b22, b32, b23, b33):
    """
    init_params_lin = [0.5, 0.9,-0.02,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1                        )
    dx2 = x2*(                    -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)            +exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_6_exp_model(y, t, r1, r3, b11, b22, b23, b33):
    """
    init_params_lin = [0.5, 0.9,-0.02,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1                        )
    dx2 = x2*(                    -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)                        -exp(b33)*x3)

    return [dx1, dx2, dx3]

def simple_sign_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.5,0.019,0.9,-0.02,+0.001,0.0015,0.0089,-0.05,0.025,0.0049,0.0057,-0.02]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(+exp(r2)+exp(b21)*x1-exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]


def remien_0_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [6,4,2,-0.05,0.15,-0.20,-0.01,-0.026,0.05,0.10,-0.10,-0.0148]
    init_params_lin = [6,4,2,-0.05,-0.01,0.10,0.15,-0.027,-0.10,-0.20,0.05,-0.0148]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1-exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]


def remien_1_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23):
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1-exp(b32)*x2            )

    return [dx1, dx2, dx3]


def remien_2_exp_model(y, t, r1, r2, r3, b11, b31, b12, b22, b32, b13, b23):
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(+exp(r2)            -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(+exp(r3)+exp(b31)*x1-exp(b32)*x2            )

    return [dx1, dx2, dx3]

def remien_3_exp_model(y, t, r1, r2, b11, b31, b12, b22, b32, b13, b23):
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(+exp(r2)            -exp(b22)*x2+exp(b23)*x3)
    dx3 = x3*(        +exp(b31)*x1-exp(b32)*x2            )

    return [dx1, dx2, dx3]

def remien_4_exp_model(y, t, r1, r2, b11, b31, b12, b22, b32, b13):
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1+exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(+exp(r2)            -exp(b22)*x2            )
    dx3 = x3*(        +exp(b31)*x1-exp(b32)*x2            )

    return [dx1, dx2, dx3]


def decay_osc_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.44,-0.07,0.27,0.07,0.14,-0.18,-0.19,-0.18,-0.25,0.17,-0.65,0.20]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)+exp(b11)*x1-exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(+exp(r3)-exp(b31)*x1-exp(b32)*x2+exp(b33)*x3)

    return [dx1, dx2, dx3]

def decayosc_0_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.44,-0.07,0.27,0.07,0.14,-0.18,-0.19,-0.18,-0.25,0.17,-0.65,0.20]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)+exp(b11)*x1-exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(+exp(r3)-exp(b31)*x1-exp(b32)*x2+exp(b33)*x3)

    return [dx1, dx2, dx3]

def decayosc_1_exp_model(y, t, r1, r2, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.44,-0.07,0.27,0.07,0.14,-0.18,-0.19,-0.18,-0.25,0.17,-0.65,0.20]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)+exp(b11)*x1-exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(        -exp(b31)*x1-exp(b32)*x2+exp(b33)*x3)

    return [dx1, dx2, dx3]

def decayosc_2_exp_model(y, t, r1, r2, b11, b21, b31, b12, b22, b32, b13, b23):
    """
    init_params_lin = [0.44,-0.07,0.27,0.07,0.14,-0.18,-0.19,-0.18,-0.25,0.17,-0.65,0.20]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)+exp(b11)*x1-exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(        -exp(b31)*x1-exp(b32)*x2            )

    return [dx1, dx2, dx3]

def decayosc_3_exp_model(y, t, r1, r2, b11, b21, b31, b12, b22, b13, b23):
    """
    init_params_lin = [0.44,-0.07,0.27,0.07,0.14,-0.18,-0.19,-0.18,-0.25,0.17,-0.65,0.20]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)+exp(b11)*x1-exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(        -exp(b31)*x1                        )

    return [dx1, dx2, dx3]

def decayosc_4_exp_model(y, t, r1, b11, b21, b31, b12, b22, b13, b23):
    """
    init_params_lin = [0.44,-0.07,0.27,0.07,0.14,-0.18,-0.19,-0.18,-0.25,0.17,-0.65,0.20]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)+exp(b11)*x1-exp(b12)*x2+exp(b13)*x3)
    dx2 = x2*(        +exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(        -exp(b31)*x1                        )

    return [dx1, dx2, dx3]
   
def cycle_0_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def cycle_1_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2            )

    return [dx1, dx2, dx3]


def cycle_2_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b32, b13, b23):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1            -exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2            )

    return [dx1, dx2, dx3]


def cycle_3_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b32, b13):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1                        )
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2            )

    return [dx1, dx2, dx3]


def cycle_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def cycle_obs_0_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23, b33):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2-exp(b33)*x3)

    return [dx1, dx2, dx3]

def cycle_obs_1_exp_model(y, t, r1, r2, r3, b11, b21, b31, b12, b22, b32, b13, b23):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2            )

    return [dx1, dx2, dx3]

def cycle_obs_2_exp_model(y, t, r1, r2, r3, b21, b31, b12, b22, b32, b13, b23):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)            -exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1+exp(b32)*x2            )

    return [dx1, dx2, dx3]

def cycle_obs_3_exp_model(y, t, r1, r2, r3, b21, b31, b12, b22, b13, b23):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)            -exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2-exp(b23)*x3)
    dx3 = x3*(-exp(r3)+exp(b31)*x1                        )

    return [dx1, dx2, dx3]

def cycle_obs_4_exp_model(y, t, r1, r2, r3, b21, b31, b12, b22, b13):
    """
    init_params_lin = [0.25, -0.5, -0.5, -0.001, 0.04, 0.02, -0.04, -0.002, 0.04, -0.04, -0.02, -0.003]
    """
    x1, x2, x3 = y
    dx1 = x1*(+exp(r1)            -exp(b12)*x2-exp(b13)*x3)
    dx2 = x2*(-exp(r2)+exp(b21)*x1-exp(b22)*x2            )
    dx3 = x3*(-exp(r3)+exp(b31)*x1                        )

    return [dx1, dx2, dx3]



def fourpop_0_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b32, b42, b13, b23, b33, b43, b14, b24, b34, b44):
    """
    init_params_lin = [ 1.00, 0.7, 1.5, 1.3, # r1..r4
                   -1.0 ,-0.5 ,-3.6 ,-1.5,
                   -1.2 ,-0.7 ,-0.1 ,-0.7,
                   -1.5 ,-0.3 ,-1.5 ,-0.4,
                   -0.2 ,-1.0 ,-0.7 ,-1.3
                   ]
    """
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3-exp(b14)*x4)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2-exp(b23)*x3-exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1-exp(b32)*x2-exp(b33)*x3-exp(b34)*x4)
    dx4 = x4*(+exp(r4)-exp(b41)*x1-exp(b42)*x2-exp(b43)*x3-exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_1_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b42, b13, b23, b33, b43, b14, b24, b34, b44):
    # Remove b32
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3-exp(b14)*x4)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2-exp(b23)*x3-exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3-exp(b34)*x4)
    dx4 = x4*(+exp(r4)-exp(b41)*x1-exp(b42)*x2-exp(b43)*x3-exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_2_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b42, b13, b33, b43, b14, b24, b34, b44):
    # Remove b23
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3-exp(b14)*x4)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2            -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3-exp(b34)*x4)
    dx4 = x4*(+exp(r4)-exp(b41)*x1-exp(b42)*x2-exp(b43)*x3-exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_3_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b42, b13, b33, b43, b14, b24, b44):
    # Remove b34
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3-exp(b14)*x4)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2            -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3            )
    dx4 = x4*(+exp(r4)-exp(b41)*x1-exp(b42)*x2-exp(b43)*x3-exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_4_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b42, b13, b33, b14, b24, b44):
    # Remove b43
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3-exp(b14)*x4)
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2            -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3            )
    dx4 = x4*(+exp(r4)-exp(b41)*x1-exp(b42)*x2            -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_5_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b42, b13, b33, b24, b44):
    # Remove b14
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2            -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3            )
    dx4 = x4*(+exp(r4)-exp(b41)*x1-exp(b42)*x2            -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_6_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b22, b13, b33, b24, b44):
    # Remove b42
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
    dx2 = x2*(+exp(r2)-exp(b21)*x1-exp(b22)*x2            -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3            )
    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_7_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b13, b33, b24, b44):
    # Remove b22
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
    dx2 = x2*(+exp(r2)-exp(b21)*x1                        -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3            )
    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_8_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b13, b33, b24, b44):
    # Remove b21
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
    dx2 = x2*(+exp(r2)                                    -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1            -exp(b33)*x3            )
    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_9_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b13, b24, b44):
    # Remove b33
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
    dx2 = x2*(+exp(r2)                                    -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_10_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b24, b44):
    # Remove b13
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2                        )
    dx2 = x2*(+exp(r2)                                    -exp(b24)*x4)
    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

def fourpop_11_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b44):
    # Remove b24
    x1, x2, x3, x4 = y
    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2                        )
    dx2 = x2*(+exp(r2)                                                )
    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)

    return [dx1, dx2, dx3, dx4]

## Eliminar
#def fourpop_8_exp_model(y, t, r1, r2, r3, r4, b11, b21, b31, b41, b12, b13, b24, b44):
#    # Remove b33
#    x1, x2, x3, x4 = y
#    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
#    dx2 = x2*(+exp(r2)-exp(b21)*x1                        -exp(b24)*x4)
#    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
#    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)
#
#    return [dx1, dx2, dx3, dx4]
#
#
#def fourpop_9_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b13, b24, b44):
#    # Remove b21
#    x1, x2, x3, x4 = y
#    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2-exp(b13)*x3            )
#    dx2 = x2*(+exp(r2)                                    -exp(b24)*x4)
#    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
#    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)
#
#    return [dx1, dx2, dx3, dx4]
#
#def fourpop_10_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b24, b44):
#    # Remove b13
#    x1, x2, x3, x4 = y
#    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2                        )
#    dx2 = x2*(+exp(r2)                                    -exp(b24)*x4)
#    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
#    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)
#
#    return [dx1, dx2, dx3, dx4]
#
#def fourpop_11_exp_model(y, t, r1, r2, r3, r4, b11, b31, b41, b12, b44):
#    # Remove b24
#    x1, x2, x3, x4 = y
#    dx1 = x1*(+exp(r1)-exp(b11)*x1-exp(b12)*x2                        )
#    dx2 = x2*(+exp(r2)                                                )
#    dx3 = x3*(+exp(r3)-exp(b31)*x1                                    )
#    dx4 = x4*(+exp(r4)-exp(b41)*x1                        -exp(b44)*x4)
#
#    return [dx1, dx2, dx3, dx4]
#
#
# Auxiliary functions
def save_solution(filename,time,output,header=None):
    out = np.array(time)
    for i in range(len(output[0])):
        out = np.vstack((out,output[:,i]))

    if header == None:
        np.savetxt(filename,out.T)
    else:
        np.savetxt(filename,out.T,header=header,comments="")


def plot_solution(time,output,legend=None,sym='-',col = ["tab:blue","tab:orange","tab:green","tab:red","tab:purple","tab:brown","tab:pink","tab:gray","tab:olive","tab:cyan"]):
    mark = ["o","v","^","<",">","s","p","P","*","X","D","d"]
    
    for i in range(len(output[0])):
        if sym=='-':
            plt.plot(time,output[:,i],'-',color=col[i])
        else:
            plt.plot(time,output[:,i],'o',color=col[i],marker=mark[i])
    plt.interactive(True)
    if legend!=None:
        plt.legend(legend)
    plt.show()

