import numpy as np
from scipy.stats import describe
import corner
from matplotlib import pyplot as plt


def init_log(n_log, n_particles):
    return np.zeros((n_log, n_particles-1, 6)), 0

def get_log_len(log):
    return log[0].shape[0]

def log_elements(sim, log):
    sim.status()
    particles = sim.particles
    elements, t_step = log
    n_log, n_particles, _ = elements.shape
    halt = 0
    for i in range(n_particles)[1:]:
        p = particles[i+1]
        
        if i == 0:
            primary = particles[2]
        # elif i==1:
        #     primary = particles[1]
        else:
            primary = particles[0]
            
        o = p.calculate_orbit(primary=primary)
        elements[t_step, i] = [o.a, o.e, o.inc, o.Omega, o.omega, sim.t]
        if o.a<0: halt += 1
    t_step += 1
    # halt = o.a < 0
    print("WARNING, binary bound check off")
    return [elements, t_step], halt>0

def save_log(log, file="output/elements.npy"):
    elements, _ = log
    np.save(file, elements)
    
def calc_moments(log, file=False):
    elements, _ = log
    n_log, n_particles, n_elements = elements.shape
    summary_stats = np.zeros((n_particles, n_elements-1, 4))
    for i in range(n_particles):
        for j in range(n_elements-1):
            summary_stats[i,j,:] = [i for i in describe(elements[:,i,j])][-4:]
    if file:
        np.save(file, summary_stats)
    return summary_stats

def plot_corner(log, file):
    elements, _ = log
    corner.corner(np.hstack((elements[:,0,:5], elements[:,1,:5])), 
                  range=[.999]*10, 
                  labels=["a1", "e1", "inc1", "Omega1", "omega1", 
                          "a2", "e2", "inc2", "Omega2", "omega2"])
    plt.savefig(file)
    plt.close()
    
def get_derivatives(log):
    elements, _ = log
    decile = int(round(elements.shape[0] / 10))
    first_i = np.mean((elements[:decile-2, :, :5] - elements[2:decile, :, :5]) 
                      / np.mean(elements[:decile-2, :, 5] - elements[2:decile, :, 5]), 
                      axis=0)
    first_f = np.mean((elements[-decile:-2, :, :5] - elements[-decile+2:, :, :5]) 
                      / np.mean(elements[-decile:-2, :, 5] - elements[-decile+2:, :, 5]), 
                      axis=0)
    
    second_i = np.mean((elements[:decile-2, :, :5] -2*elements[1:decile-1, :, :5]
                        +elements[2:decile, :, :5]) 
                      / np.mean(elements[:decile-2, :, 5] - elements[2:decile, :, 5])**2/4, 
                      axis=0)
    
    second_f = np.mean((elements[-decile:-2, :, :5] - 2*elements[-decile+1:-1, :, :5]
                       + elements[-decile+2:, :, :5]) 
                    / np.mean(elements[-decile:-2, :, 5] - elements[-decile+2:, :, 5])**2/4, 
                    axis=0)
    
    return first_i, first_f, second_i, second_f