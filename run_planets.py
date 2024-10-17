from binary_planets.run_model import run_model
from binary_planets.sim import get_hill_radius
from binary_planets.utils import *

from os import sys
import time
import json
import os
import numpy as np
from multiprocessing import Pool


global config
global mode
global pl_num
compact_sys = np.load("data/compact_systems_run_composite.npy", allow_pickle=True)
# compact_sys = np.load("data/TOI-178.npy", allow_pickle=True)
# compact_sys = np.load("data/Kepler-11.npy", allow_pickle=True)
n_sys = len(compact_sys)
n_runs = 10000
seq = np.random.SeedSequence().generate_state(n_runs)



# def one_run(run_num, config, mode, pl_num, debug=0):
def one_run(run_num):
    seed = seq[run_num]
    print(f"Seed: {seed}")
    np.random.seed(seed)
    cfg = dict(config)
    sys_num = pl_num#np.random.randint(n_sys)
    batch_num = int(int(run_num) / n_sys) + 1

    system = compact_sys[sys_num]
    cfg["name"] = f"{cfg['name']}/{system['name']}"
    print(f"Sys_Num: {sys_num} Run Number: {run_num}. Run Name: {cfg['name']}. Batch #: {batch_num}",
            file=sys.stderr)

    n_secondary = len(system["a"])


    cfg["n_secondary"] = n_secondary

    if mode == 1:
        cfg["n_secondary"] += 1

    # cfg["n_secondary"] = 0 #DEBUG

    name = cfg["name"]
    for i in range(1):
        cfg["name"] = name
        cfg["m_star"] = float(get_stellar_mass(system))
        print(f"m_star: {cfg['m_star']}", file=sys.stderr)
        semi_majors = get_semimajor(system)
        es = get_eccen(system)
        incs = get_inc(system)
        masses = get_mass(system)
        for j in range(n_secondary):
            cfg[f"secondary_{j}"] = {}
            cfg[f"secondary_{j}"]["m"] = masses[j]
            cfg[f"secondary_{j}"]["a"] = semi_majors[j]
            cfg[f"secondary_{j}"]["e"] = es[j]
            cfg[f"secondary_{j}"]["inc"] = incs[j] - np.pi/2
            cfg[f"secondary_{j}"]["omega"] = np.random.uniform(-np.pi, np.pi)
            cfg[f"secondary_{j}"]["Omega"] = np.random.uniform(-np.pi, np.pi)


        if mode == 1:
            j = n_secondary
            cfg[f"secondary_{j}"] = {}
            cfg[f"secondary_{j}"]["m"] = np.random.normal(np.mean(masses),
                                                            np.std(masses))

            diff = np.abs(system["gap"][1] - system["gap"][0])

            cfg[f"secondary_{j}"]["a"] = np.random.uniform(system["gap"][0] + 0.2*diff,
                                                            system["gap"][1] - 0.2*diff)

            cfg[f"secondary_{j}"]["e"] = np.random.uniform(0, .3)
            cfg[f"secondary_{j}"]["inc"] = 0 #np.random.uniform(-np.pi, np.pi)
            cfg[f"secondary_{j}"]["omega"] = np.random.uniform(-np.pi, np.pi)
            cfg[f"secondary_{j}"]["Omega"] = np.random.uniform(-np.pi, np.pi)

        if mode == 2:
            mass_total = np.random.uniform(0.38, 72)
            q = np.random.uniform(0.25, 0.5)
            cfg["binary"]["m1"] = q * mass_total
            cfg["binary"]["m2"] = (1 - q) * mass_total
            cfg["binary"]["e"] = np.random.uniform(0, 0.4)
            cfg["binary"]["e_sys"] = np.random.uniform(0, 0.4)
            cfg["binary"]["phase"] = np.random.uniform(-np.pi, np.pi)

            cfg["binary"]["Omega"] = np.random.uniform(-np.pi, np.pi)
            cfg["binary"]["inc"] = 0 #np.random.uniform(-np.pi, np.pi) #TODO
            cfg["binary"]["bin_inc"] =0 # np.random.uniform(-np.pi, np.pi) #TODO

            diff = np.abs(np.log10(system["gap"][1]) - np.log10(system["gap"][0]))

            cfg["binary"]["a"] = 10 ** np.random.uniform(np.log10(system["gap"][0]),
                                                   np.log10(system["gap"][1]))
            r_hill = get_hill_radius(cfg["binary"]["a"],
                                    cfg["binary"]["e"],
                                    mass_total,
                                    cfg["m_star"])
            cfg["binary"]["d"] = r_hill * np.random.uniform(0, 1.5)

        else:
            cfg["binary"]["m1"] = 1e-20
            cfg["binary"]["m2"] = 1e-20




        i = 1
        try:
            os.mkdir(f"output/{cfg['name']}")
        except:
            pass
        while(True):
            try:
                os.mkdir(f"output/{cfg['name']}/{i}")
                break
            except:
                i += 1
        cfg["name"] = f"{cfg['name']}/{i}"

        print(f"system configuration: {cfg}", file=sys.stderr)
        run_model(cfg, mode)
        return

if __name__ == "__main__":
    try:
        print(f"Base Config file: {sys.argv[1]}", file=sys.stderr)
        with open(sys.argv[1]) as f:
            config = json.load(f)
        if len(sys.argv) > 2:
            config["name"] = sys.argv[2]




        mode = int(sys.argv[3]) # 0 for no inj, 1 for single inj, 2 for binary inj
        pl_num = int(sys.argv[4])
        print(f"Simulating {compact_sys[pl_num]['name']}, mode = {mode}", file=sys.stderr)

    except Exception as e:
        print("Error in loading config file")
        print(f"Error was {e}")
        raise()

    with Pool() as p:
        p.map(one_run, range(0, 1000))


