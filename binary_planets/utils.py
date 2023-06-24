from scipy.stats import skewnorm
from scipy.optimize import curve_fit
from scipy.interpolate import interp1d
import numpy as np
from .forecaster import mr_forecast as mr
import sys as system
def _cdf(x,a,b,c):
    return skewnorm.cdf(x,a,b,c)

mass, cdf = np.load("data/st_mass_dist.npy")
mass_func = interp1d(cdf, mass)
MASS_CDF_MIN = np.min(cdf)
MASS_CDF_MAX = np.max(cdf)

e, cdf = np.load("data/e_dist.npy")
e_func = interp1d(cdf, e)
E_CDF_MIN = np.min(cdf)
E_CDF_MAX = np.max(cdf)

def draw_rv(mid, lo, hi, num=1):
    if (np.isnan(lo) * np.isnan(hi)* (mid==0)): return 0
    if np.isnan(lo): lo = -.2*mid
    if np.isnan(hi): hi = .2*mid
    try:
        [a,b,c], _ = curve_fit(_cdf, [mid+lo, mid, mid+hi], [.16, .5, .84]) 
        print(a,b,c)
        return skewnorm.rvs(a,b,c, size=num)   
    except:
        return np.random.normal(mid, np.max([lo, hi]))
    
def get_stellar_mass(sys):
    mass = sys["st_mass"]
    if (np.isnan(mass) + (mass<.1)):
        x = np.random.uniform(MASS_CDF_MIN, MASS_CDF_MAX)
        return mass_func(x)
    return draw_rv(mass, sys["st_lower"], sys["st_upper"])
    
def get_eccen(sys):
    mids = sys["e"]
    los = sys["e_lower"]
    his = sys["e_upper"]
    es = np.zeros(len(mids))
    for i in range(len(es)):
        if np.isnan(mids[i]):
            x = np.random.uniform(E_CDF_MIN, E_CDF_MAX)
            es[i] = e_func(x)
        else:
            while True:
                es[i] = draw_rv(mids[i], los[i], his[i])
                if es[i]>0: break

    return es

def get_inc(sys):
    mids = np.array(sys["inc"])
    
    if np.all(np.isnan(mids)):
        return np.ones(len(mids)) * 90
    
    incs = np.zeros(len(mids))
    los = sys["inc_lower"]
    his = sys["inc_upper"]
    
    mean = np.mean(mids[~np.isnan(mids)])
    std = np.std(mids[~np.isnan(mids)])
    
    for i in range(len(mids)):
        if np.isnan(mids[i]):
            incs[i] = np.random.normal(mean, std)
        else:
            incs[i] = draw_rv(mids[i], los[i], his[i])
            
    return incs

def get_mass(sys):
    mass_mid = sys["mass"]
    mass_upper = sys["mass_upper"]
    mass_lower = sys["mass_lower"]
    r_mid = sys["r"]
    r_upper = sys["r_upper"]
    r_lower = sys["r_lower"]
    mass = np.zeros(len(mass_mid))
    for i in range(len(mass)):
        while True:
            if np.isnan(mass_mid[i]):
                r = draw_rv(r_mid[i], r_lower[i], r_upper[i])
                mass[i] = mr.Rpost2M([r], "Earth", 1e5, classify="No")
            else:
                mass[i] = draw_rv(mass_mid[i], mass_lower[i], mass_upper[i])
            if mass[i] > 0: break
    return mass
            
def get_semimajor(sys):
    semi_majors = np.zeros(len(sys["a"]))
    for i in range(len(sys["a"])):
        semi_majors[i] = draw_rv(sys["a"][i], sys["a_lower"][i], sys["a_upper"][i])
    return semi_majors