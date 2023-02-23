import numpy as np
import rebound
from .constants import *
from .log import *

def orbital_charcteristics(m1, m2, d, e=0, phase=0, Omega=0, inc=0):
    r1 = m2/(m1+m2)*d
    r2 = m1/(m1+m2)*d
    
    a1 = get_semi_major(r1, e, phase)
    a2 = get_semi_major(r2, e, phase)

    T1 = np.sqrt(4*np.pi**2/(G_AU_YR_MSUN * (m1+m2)) * (a1+a2)**3)
    T2 = T1
    
    
    v1 = get_velocity(r1, a1, e, T1)
    v2 = get_velocity(r2, a2, e, T2)
    

    
    xy1 = np.array([r1*np.cos(phase), r1*np.sin(phase), 0])
    xy2 = np.array([r2*np.cos(phase+np.pi), r2*np.sin(phase+np.pi), 0])
    
    

    rot = np.array([[np.cos(Omega), -np.sin(Omega)*np.cos(inc), np.sin(Omega)*np.sin(inc)],
                    [np.sin(Omega), np.cos(Omega)*np.cos(inc), -np.cos(Omega)*np.sin(inc)],
                    [0,             np.sin(inc),                np.cos(inc)]])
    
    x_vec1 = np.matmul(rot, xy1)
    x_vec2 = np.matmul(rot, xy2)
    
    v_dir1 = np.cross(x_vec1, np.matmul(rot, [0,0,1]))
    v_dir2 = np.cross(x_vec2, np.matmul(rot, [0,0,1]))
    
    
    v_hat1 = v_dir1/np.linalg.norm(v_dir1)
    v_hat2 = v_dir2/np.linalg.norm(v_dir2)
    
    v_vec1 = v_hat1 * v1
    v_vec2 = v_hat2 * v2

    return [x_vec1, v_vec1], [x_vec2, v_vec2]
    
    
def get_period(m1, m2, d, e=0, phase=0):
    m1=m1*MASS_E/MASS_SUN
    m2=m2*MASS_E/MASS_SUN
    
    r1 = m2/(m1+m2)*d
    r2 = m1/(m1+m2)*d
    
    a1 = get_semi_major(r1, e, phase)
    a2 = get_semi_major(r2, e, phase)

    return np.sqrt(4*np.pi**2/(G_AU_YR_MSUN * (m1+m2)) * (a1+a2)**3)
    
    
    
def get_semi_major(r, e, theta):
    return r*(1 + e*np.cos(theta))/(1 - e**2)

def get_velocity(r, a , e, T):
    return 2* np.pi * (a**2) * np.sqrt(1-e**2)/(r*T)



def init_binary_planet(m1, m2, d, e=0, phase=0, Omega=0, inc=0, n_log=1000, 
                       integrator="whfast", dt=1e-4):

    sim=rebound.Simulation()
    print(sim.status())
    sim.units = ('yr', 'AU', 'Msun')
    sim.integrator = integrator
    sim.dt = dt
    sim.ri_whfast.safe_mode = 0
    sim.ri_whfast.corrector = 11
    sim.add(m=1)
    sim.move_to_hel()

    m1=m1*MASS_E/MASS_SUN
    m2=m2*MASS_E/MASS_SUN

    [x1, v1], [x2, v2] = orbital_charcteristics(m1, m2, d, inc=inc, e=e, 
                                                phase=phase, Omega=Omega)


    sim.add(m=m1, 
            x=1+x1[0], y=x1[1], z=x1[2], 
            vx=v1[0],vy=2*np.pi+v1[1], vz=v1[2])
    
    sim.add(m=m2,
            x=1+x2[0], y=x2[1], z=x2[2], 
            vx=v2[0], vy=2*np.pi+v2[1], vz=v2[2])

    sim.move_to_com()
    log = init_log(n_log, 3)

    
    return sim, log

def simulate(sim, log, t_end):
    n_log = get_log_len(log)
    
    for t in np.linspace(0, t_end, n_log):
        sim.integrate(t, exact_finish_time=0)
        log = log_elements(sim, log)
    
    return sim, log

