import numpy as np
from scipy.integrate import odeint

g = 9.8
dt = 0.05
Kx = [0.5,0,0] # Kpx, Kix, Kdx (=0 because that's just velocity)
Kv = [5,0,0.1] # Kpv, Kiv (=0 because that's just diplacement), Kdv

# The following function gives the ordinary differential
# equation that our plant follows. Do not meddle with this.
def f(x, t, theta):
    return (x[1], (-5 * g / 7) * np.radians(theta))



# Write your function here.

def solve(theta_prev,x_prev,v_prev,ex_prev,ev_prev,Ix_prev,Iv_prev,X,Y):
    temp = odeint(f, [x_prev,v_prev], [0,dt], args=(theta_prev,))
    x,v = temp[1]

    ex = X - 400 - x
    ev = -v

    Ix = Ix_prev + Kx[1]*ex*dt
    Iv = Iv_prev + Kv[1]*ev*dt
    
    theta = -( (Kx[0]*ex + Ix - Kx[2]*(x-x_prev)/dt) + (Kv[0]*ev + Iv - Kv[2]*(v-v_prev)/dt) )

    if abs(x) > 300:
        x = 300*np.sign(x)
        Ix = Ix_prev
    if abs(theta - theta_prev) > 1:
        theta = theta_prev + 1*np.sign(theta - theta_prev)
    if abs(theta) > 15:
        theta = 15*np.sign(theta)
    
    return theta, x, v, ex, ev, Ix, Iv

# print(solve(0.1, 1, 1, 0, 0, 0, 0, 1, 0))
