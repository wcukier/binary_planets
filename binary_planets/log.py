import numpy as np
from scipy.stats import describe



def init_log(n_log, n_particles):
    return np.zeros((n_log, n_particles-1, 5)), 0

def get_log_len(log):
    return log[0].shape[0]

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

def save_log(log, file="output/elements.npy"):
    elements, _ = log
    np.save(file, elements)
    
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