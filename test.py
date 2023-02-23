"""test.py
A test system with two Earth-like planets orbiting eachother
"""

import rebound
import numpy as np
from binary_planets.constants import *
import imageio
from matplotlib import pyplot as plt
from tqdm import tqdm
import time
from scipy.stats import describe

def animate_separation(name, sim, t_start, t_end, n_frames=100):
    n_particles = sim.N
    x = np.zeros((2, n_frames))
    y = np.zeros((2, n_frames))

    dt = (t_end - t_start)/n_frames

    for i in tqdm(range(n_frames)):
        sim.integrate(t_start + i*dt, exact_finish_time=0)
        p1 = sim.particles[1]
        p2 = sim.particles[2]
        x_com = (p1.x*p1.m + p2.x*p2.m)/(p1.m+p2.m)
        y_com = (p1.y*p1.m + p2.y*p2.m)/(p1.m+p2.m)

        
        x[0, i] = p1.x - x_com
        y[0, i] = p1.y - y_com
        
        x[1, i] = p2.x - x_com
        y[1, i] = p2.y - y_com
        
        plt.plot(x[0, :i+1], y[0, :i+1], color="r")
        plt.plot(x[0, i], y[0, i], marker = "o", color="r")
        
        plt.plot(x[1, :i+1], y[1, :i+1], color="k")
        plt.plot(x[1, i], y[1, i], marker = "o", color="k")
        plt.text(-.0075,-.0075, f"{(i*dt) % 1:.2f}")

        plt.axis("square")
        plt.xlim(-.01, .01)
        plt.ylim(-0.01, .01)
        plt.tight_layout()
        plt.savefig(f"anim/{i}.png")
        plt.close()

    frames = [imageio.imread(f"anim/{i}.png") for i in range(n_frames)]
    imageio.mimsave(f'{name}.gif', frames, fps=10)

def animate_simulation(name, sim, t_start, t_end, n_frames=100):
    n_particles = sim.N
    x = np.zeros((n_particles,n_frames))
    y = np.zeros((n_particles, n_frames))

    dt = (t_end - t_start)/n_frames

    for i in range(n_frames):
        sim.integrate(t_start + i*dt, exact_finish_time=0)
        for n_p in range(n_particles):
            p = sim.particles[n_p]
            x[n_p, i] = p.x
            y[n_p, i] = p.y
            plt.plot(x[n_p, :i+1], y[n_p, :i+1], color="k")
            plt.plot(x[n_p, i], y[n_p, i], marker = "o", color="k")
        plt.axis("square")
        plt.xlim(.95, 1.05)
        plt.ylim(-.05, .05)
        plt.savefig(f"anim/{i}.png")
        plt.close()

    frames = [imageio.imread(f"anim/{i}.png") for i in range(n_frames)]
    imageio.mimsave(f'{name}.gif', frames, fps=10)


def init_log(n_log, n_particles):
    return np.zeros((n_log, n_particles-1, 5)), 0

def log_elements(sim, log):
    particles = sim.particles
    elements, t_step = log
    n_log, n_particles, _ = elements.shape
    for i in range(n_particles):
        p = particles[i+1]
        o = p.calculate_orbit(primary=particles[0])
        elements[t_step, i] = [o.a, o.e, o.inc, o.Omega, o.omega]
    t_step += 1
    return elements, t_step

def save_log(log):
    elements, _ = log
    np.save("elements.npy", elements)

def calc_moments(log, file="stats.npy"):
    elements, _ = log
    n_log, n_particles, n_elements = elements.shape
    summary_stats = np.zeros((n_particles, n_elements, 4))
    for i in range(n_particles):
        for j in range(n_elements):
            summary_stats[i,j,:] = [i for i in describe(elements[:,i,j])][-4:]
    if file:
        np.save(file, summary_stats)
    return summary_stats
    
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
    
def get_semi_major(r, e, theta):
    return r*(1 + e*np.cos(theta))/(1 - e**2)

def get_velocity(r, a , e, T):
    return 2* np.pi * (a**2) * np.sqrt(1-e**2)/(r*T)
    
t_end = 10000

n_log = 1000    

sim=rebound.Simulation()
sim.units = ('yr', 'AU', 'Msun')
sim.integrator = "ias15"
sim.dt = 1e-5
sim.ri_whfast.safe_mode = 0
sim.ri_whfast.corrector = 11
sim.add(m=1)
sim.move_to_hel()
# Earth - Moon
# sim.add(m=MASS_E/MASS_SUN, x=1., vy=2*np.pi*(1))
# sim.add(m=MASS_E/MASS_SUN/6, x=1+.00257, vy=2*np.pi + .215)

# 2 Earths

m1=MASS_E/MASS_SUN
m2=MASS_E/MASS_SUN

# [x1, v1], [x2, v2] = orbital_charcteristics(MASS_E/MASS_SUN, MASS_E/MASS_SUN/6, .00257)
[x1, v1], [x2, v2] = orbital_charcteristics(m1, m2, .005, inc=-np.pi)



sim.add(m=m1, x=1+x1[0], y=x1[1], z=x1[2], vx=v1[0], vy=2*np.pi+v1[1], vz=v1[2])
sim.add(m=m2, x=1+x2[0], y=x2[1], z=x2[2], vx=v2[0], vy=2*np.pi+v2[1], vz=v2[2])

sim.move_to_com()
op = rebound.OrbitPlot(sim)
op.fig.savefig("initial_orbit.png")

print("Initial Orbital Elements")
for p in sim.particles[2:]:
    o = p.calculate_orbit(primary=sim.particles[1])
    print(f"a:{o.a}, e:{o.e}")

# animate_separation("start", sim, 0, 3) 
log = init_log(n_log, 3)
t0 = time.time()
for t in np.linspace(0, t_end, n_log):
    sim.integrate(t, exact_finish_time=0)
    log = log_elements(sim, log)

t1 = time.time()
print(f"Time elapesed: {t1-t0}")
save_log(log)
moments = calc_moments(log)

print(moments)


op = rebound.OrbitPlot(sim)
op.fig.savefig("final_orbit.png")

print("\nFinal Orbital Elements")
for p in sim.particles[2:]:
    o = p.calculate_orbit(primary=sim.particles[1])
    print(f"a:{o.a}, e:{o.e}")
animate_separation("end", sim, t_end, t_end+.4)




    
