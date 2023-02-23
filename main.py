from binary_planets.anim import *
from binary_planets.constants import *
from binary_planets.log import *
from binary_planets.sim import *

import sys
import json
import time
import os
import datetime

if __name__ == "__main__":
    try:
        print(f"Config file: {sys.argv[1]}", file=sys.stderr)
        with open(sys.argv[1]) as f:
            config = json.load(f)
        if len(sys.argv) > 2:
            config["name"] = sys.argv[2]
        print(f"Run Name: {config['name']}", file=sys.stderr)

    except Exception as e:
        print("Error in loading config file")
        print(f"Error was {e}")
        raise()


    i = 1
    if not os.path.exists(f"output/{config['name']}"):
        os.mkdir(f"output/{config['name']}")
    while(True):
        if not os.path.exists(f"output/{config['name']}/{i}"):
            config["name"] = f"{config['name']}/{i}"
            os.mkdir(f"output/{config['name']}")
            break
        i += 1

    with open(f"output/{config['name']}/run_notes.out", 'w+') as run_notes:
        run_notes.write(f"{sys.argv}\n")
        run_notes.write(f"Run start time: {datetime.datetime.now()}\n")
        for key in config.keys():
            run_notes.write(f"{key}: {config[key]}\n")
        run_notes.write("*************************************************\n\n")

    sim, log = init_binary_planet(config["m1"], config["m2"], config["d"], 
                                config["e"], config["phase"], 
                                config["Omega"], config["inc"], 
                                config["n_log"], config["integrator"], 
                                config["dt"])
    
    t0 = time.time()
    simulate(sim, log, config["t_end"])
    t1 = time.time()

    save_log(log, f"output/{config['name']}/elements.npy")

    with open(f"output/{config['name']}/summary.out", 'a+') as summary:
        summary.write(f"Time elapsed: {t1-t0}\n")
        summary.write("\nStatistical Moments:\n")
        summary.write(f"{calc_moments(log)}")
        
    T = get_period(config["m1"], config["m2"], config["d"], 
                   config["e"], config["phase"])
    print(T)
    animate_separation(f"output/{config['name']}/end", sim, 
                       config['t_end'], config['t_end']+.2)