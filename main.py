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
    
    E0 = sim.energy()
    Lx0, Ly0, Lz0, = sim.angular_momentum()
    
    t0 = time.time()
    simulate(sim, log, config["t_end"])
    t1 = time.time()

    save_log(log, f"output/{config['name']}/elements.npy")

    E1 = sim.energy()
    Lx1, Ly1, Lz1, = sim.angular_momentum()

    E_err = (E0 - E1)/E0
    L_err = np.sqrt(((Lx0 - Lx1)**2 + (Ly0 - Ly1)**2 + (Lz0 - Lz1)**2)
                    /(Lx0**2+ Ly0**2 + Lx0**2+1e-20))

    with open(f"output/{config['name']}/run_notes.out", 'a+') as summary:
        summary.write(f"Time elapsed: {t1-t0}\n")
        
        summary.write(f"Energy Error: {E_err:.2e}\tE_i: {E0:.2e}\
            \tE_f: {E1:.2e}\n")
        
        summary.write(f"Angular Mom Err: {L_err:.2e}\t\
                      L_i: ({Lx0:.2e}, {Ly0:.2e}, {Lz0:.2e})\t\
                          L_f: ({Lx1:.2e}, {Ly1:.2e}, {Lz1:.2e})\n")
        
        summary.write("\nStatistical Moments:\n")
        summary.write(f"{calc_moments(log)}")
        
    T = get_period(config["m1"], config["m2"], config["d"], 
                   config["e"], config["phase"])
    print(T)
    # animate_separation(f"output/{config['name']}/end", sim, 
    #                    config['t_end'], config['t_end']+.2)