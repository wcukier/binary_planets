from binary_planets.run_model import run_model
from binary_planets.sim import get_hill_radius


from os import sys
import json
import os
import numpy as np

if __name__ == "__main__":
    try:
        print(f"Config file: {sys.argv[1]}", file=sys.stderr)
        with open(sys.argv[1]) as f:
            config = json.load(f)
        if len(sys.argv) > 2:
            config["name"] = sys.argv[2]
            
        run_num = sys.argv[3]
        print(f"Run Name: {config['name']}. Run #: {run_num}", file=sys.stderr)

    except Exception as e:
        print("Error in loading config file")
        print(f"Error was {e}")
        raise()

    r_hill = get_hill_radius(config["binary"]["a"], 
                             config["binary"]["e"], 
                             config["binary"]["m1"])
   
    # jupiter_sep = np.linspace(.85, 0.05, 96)[int(run_num)]
    
    # config["secondary_0"]["a"] = 1 - jupiter_sep
    # config["secondary_1"]["a"] = 1 + jupiter_sep
    
    
    # ds = np.linspace(0.05*r_hill, 1.2*r_hill, 20)
    # q = np.linspace(0, 1, 5)[1:][int(int(run_num)%5)]
    # e0 = np.linspace(0, 1, 20)[int(int(run_num)/5)]
    # e1 = np.linspace(0, 1, 20)[int(int(run_num)/5)+1] 
    # es = np.linspace(e0, e1, 5)
    # np.random.shuffle(ds)
    # np.random.shuffle(es)


    # ds = np.linspace(.25, .35, 100)  * r_hill
    # q = np.linspace(0, 0.05, 11)[1:][int(run_num)]

    # config["m2"] = config["m1"] * q
    # config["dt"] = 0.01
    # config["t_end"] = 10000000


    if len(sys.argv) > 5:
        key = sys.argv[4]
        val = sys.argv[5]
        # print(f"{key}: {val}", file=sys.stderr)
        # print(f"{config[key]}", file=sys.stderr)
        print(type(val), file=sys.stderr)
        config[key] = float(val)
        
        
    name = config['name']
    ds = [1]
    es = [1]
    for i in range(400):
        # d = ds[int(i % len(ds))]
        # e = es[int(i / len(ds))]
        # # config["d"] = d
        # # config["e"] = e

        config["binary"]["d"] = np.random.uniform(0.05, 1.0) * r_hill
        config["binary"]["e"] = np.random.uniform(0.0, 0.2)
        config["binary"]["inc"] = np.random.uniform(-0.25, 0.25)
        config["binary"]["e_sys"] = np.random.uniform(0.0, 0.2)

        inc = np.random.uniform(-0.25, 0.25)
        
        config["secondary_0"]["omega"] = np.random.uniform(-np.pi, np.pi)
        config["secondary_0"]["Omega"] = np.random.uniform(-np.pi, np.pi)
        config["secondary_0"]["inc"] = inc + np.random.uniform(-0.25, 0.25)
        config["secondary_1"]["omega"] = np.random.uniform(-np.pi, np.pi)
        config["secondary_1"]["Omega"] = np.random.uniform(-np.pi, np.pi)
        config["secondary_1"]["inc"] = inc + np.random.uniform(-0.25, 0.25)
        
        
        # separation = np.random.uniform(.15, .3)
        # contrast = 3 ** np.random.uniform(-1, 1)
        # config["secondary_0"]["a"] = 1 +(separation * contrast)
        # config["secondary_1"]["a"] = 1 - (separation)

        config["secondary_0"]["a"] = np.random.uniform(1, 1.5)
        config["secondary_1"]["a"] = np.random.uniform(0.5, 1)


        mass = (2 ** np.random.uniform(np.log2(0.1), np.log2(20))) * 2
        m1 = mass * np.random.uniform(0.2, 0.8)
        m2 = mass - m1

        config["secondary_0"]["m"] = m1
        config["secondary_1"]["m"] = m2
        
        m_p = 2 ** np.random.uniform(np.log2(0.05), np.log2(10))
        q = np.random.uniform(0.5, 1)
        m2 = m1 * (1 - q)/q
        
        config["binary"]["m1"] = m_p
        config["binary"]["m2"] = m2
        
        config['name'] = name 

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
        
        run_model(config)