from binary_planets.run_model import run_model
from binary_planets.sim import get_hill_radius
from binary_planets.utils import *

from os import sys
import json
import os
import numpy as np

if __name__ == "__main__":
    try:
        print(f"Base Config file: {sys.argv[1]}", file=sys.stderr)
        with open(sys.argv[1]) as f:
            config = json.load(f)
        if len(sys.argv) > 2:
            config["name"] = sys.argv[2]
            
        run_num = sys.argv[3]
        sys_num = int(int(run_num) % 9)
        batch_num = int(int(run_num) / 9) + 1
        
        compact_sys = np.load("data/compact_systems.npy", allow_pickle=True)
        system = compact_sys[sys_num]
        config["name"] = f"{config['name']}/{system['name']}"
        print(f"Run Name: {config['name']}. Batch #: {batch_num}", 
              file=sys.stderr)

        mode = sys.argv[4] # 0 for no inj, 1 for single inj, 2 for binary inj

    except Exception as e:
        print("Error in loading config file")
        print(f"Error was {e}")
        raise()

    n_secondary = len(system["a"])
    

    
    config["n_secondary"] = n_secondary 
    
    if mode == 1:
        config["n_secondary"] += 1
    
    name = config["name"]
    for i in range(400):
        config["name"] = name
        config["m_star"] = float(get_stellar_mass(system))
        print(f"m_star: {config['m_star']}", file=sys.stderr)
        semi_majors = get_semimajor(system)
        es = get_eccen(system)
        incs = get_inc(system)
        masses = get_mass(system)
        for j in range(n_secondary):
            config[f"secondary_{j}"] = {}
            config[f"secondary_{j}"]["m"] = masses[j]
            config[f"secondary_{j}"]["a"] = semi_majors[j]
            config[f"secondary_{j}"]["e"] = es[j]
            config[f"secondary_{j}"]["inc"] = incs[j]
            config[f"secondary_{j}"]["omega"] = np.random.uniform(-np.pi, np.pi)
            config[f"secondary_{j}"]["Omega"] = np.random.uniform(-np.pi, np.pi)        
            
        separation = np.random.uniform(.15, .3)

        if mode == 1:
            j = n_secondary
            config[f"secondary_{j}"] = {}
            config[f"secondary_{j}"]["m"] = np.random.normal(np.mean(masses),
                                                             np.std(masses))
            config[f"secondary_{j}"]["a"] = np.random.uniform(system["gap"][0], 
                                                              system["gap"][1])
            config[f"secondary_{j}"]["e"] = np.random.uniform(0, 1) 
            config[f"secondary_{j}"]["inc"] = np.random.uniform(-np.pi, np.pi)
            config[f"secondary_{j}"]["omega"] = np.random.uniform(-np.pi, np.pi)
            config[f"secondary_{j}"]["Omega"] = np.random.uniform(-np.pi, np.pi)
        
        if mode == 2:
            mass_total = np.random.uniform(0.38, 72)
            q = np.random.uniform(0.25, 0.5)
            config["binary"]["m1"] = q * mass_total
            config["binary"]["m2"] = (1 - q) * mass_total
            config["binary"]["e"] = np.random.uniform(0,1)
            config["binary"]["e_sys"] = np.random.uniform(0,1)
            config["binary"]["phase"] = np.random.uniform(-np.pi, np.pi)
            config["binary"]["Omega"] = np.random.uniform(-np.pi, np.pi)
            config["binary"]["inc"] = np.random.uniform(-np.pi, np.pi)
            config["binary"]["a"] = np.random.uniform(system["gap"][0], 
                                                      system["gap"][1])
            r_hill = get_hill_radius(config["binary"]["a"], 
                                     config["binary"]["e"],
                                     config["binary"]["m2"],
                                     config["m_star"])
            config["binary"]["d"] = r_hill * np.random.uniform(0, 1)
        
        else:
            config["binary"]["m1"] = 1e-20
            config["binary"]["m2"] = 1e-20
            
        
        

        i = 1
        try:
            os.mkdir(f"output/{config['name']}")
        except:
            pass
        while(True):
            try:
                os.mkdir(f"output/{config['name']}/{i}")
                break
            except:
                i += 1
        config["name"] = f"{config['name']}/{i}"
        
        print(config, file=sys.stderr)
        run_model(config, mode)