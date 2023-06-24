# utilities
import os, time, datetime
import pickle
from tqdm import tqdm

## numerics
import numpy as np
import scipy as sp
from scipy.special import erfi
import scipy.fftpack
import scipy.stats
import pandas as pd

# plotting
import matplotlib
matplotlib.rcParams['figure.dpi']= 200

import matplotlib.pyplot as plt
plt.rc('text', usetex=True)
plt.rc('text.latex', preamble=r'\usepackage{amsmath}')

# paralel processing
import multiprocessing

# astropy
import astropy.coordinates, astropy.units
import astropy.io

import sklearn


class gdatstrt(object):

    def __init__(self):
        #self.boollockmodi = False
        pass
    
    def __setattr__(self, attr, valu):
        super(gdatstrt, self).__setattr__(attr, valu)


class gdatstrtpcat(object):
    
    def __init__(self):
        self.boollockmodi = False
        pass

    def lockmodi(self):
        self.boollockmodi = True

    def unlkmodi(self):
        self.boollockmodi = False

    def __setattr__(self, attr, valu):
        
        if hasattr(self, attr) and self.boollockmodi and attr != 'boollockmodi':
            raise KeyError('{} has already been set'.format(attr))
        
        if len(attr) % 4 != 0 and not attr.startswith('path'):
            raise Exception('')
        
        #if attr == 'thislliktotl' and hasattr(self, attr) and getattr(self, attr) - 100. > valu:
        #    raise Exception('Trying to decrease lliktotl too much...')
       
        # temp
        #if attr == 'indxsampmodi':
        #    if not isinstance(valu, int64) and len(valu) > 6:
        #        raise Exception('Setting gdatmodi.indxsampmodi to an array!')
        
        super(gdatstrt, self).__setattr__(attr, valu)


def retr_dictstrg():
    
    dicttdpy = dict()
    dicttdpy['raww'] = 'Raw'
    dicttdpy['lcur'] = 'LightCurve'
    dicttdpy['cosc'] = 'CompactObjectWithStellarCompanion'
    dicttdpy['psys'] = 'SystemOfPlanets'
    dicttdpy['psyspcur'] = 'SystemOfPlanetsWithPhaseCurve'
    dicttdpy['psysdiskedgehori'] = 'SystemOfPlanetsWithEdgeOnHorizontalDisks'

    return dicttdpy


class datapara(object):

    def __init__(self, numbpara):
        
        self.numbpara = numbpara
        
        self.indx = dict()
        self.name = np.empty(numbpara, dtype=object)
        self.minm = np.zeros(numbpara)
        self.maxm = np.zeros(numbpara)
        self.scal = np.empty(numbpara, dtype=object)
        self.labl = np.empty(numbpara, dtype=object)
        self.unit = np.empty(numbpara, dtype=object)
        self.vari = np.zeros(numbpara)
        self.true = np.zeros(numbpara)
        self.strg = np.empty(numbpara, dtype=object)
        
        self.cntr = 0
    
    def defn_para(self, name, minm, maxm, scal, labl, unit, vari, true):
        
        self.indx[name] = self.cntr
        self.name[self.cntr] = name
        self.minm[self.cntr] = minm
        self.maxm[self.cntr] = maxm
        self.scal[self.cntr] = scal
        self.labl[self.cntr] = labl
        self.unit[self.cntr] = unit
        self.vari[self.cntr] = vari
        self.true[self.cntr] = true
        self.strg[self.cntr] = labl + ' ' + unit
        self.cntr += 1


def time_func(func, *args):
    '''
    Measure the execution time of a function.
    '''    
    
    numbiter = 100
    timediff = np.empty(numbiter)
    for k in range(numbiter):
        timeinit = time.time()
        func(*args)
        timediff[k] = time.time() - timeinit
    
    return mean(timediff), std(timediff)


def time_func_verb(func, *args):
    
    meantimediff, stdvtimediff = time_func(func, *args)


# astrophysics
def retr_specbbod(tmpt, wlen):
    '''
    Calculate the spectrum of a blackbody
    '''
    
    #0.0143877735e6 # [um K]
    spec = 3.742e11 / wlen**5 / (np.exp(0.0143877735e6 / (wlen * tmpt)) - 1.)
    
    return spec


def calc_visitarg(rasctarg, decltarg, latiobvt, longobvt, strgtimeobvtyear, listdelttimeobvtyear, heigobvt=None):
    '''
    Calculate the visibility of a celestial target.
    '''
    
    print('rasctarg')
    print(rasctarg)
    print('decltarg')
    print(decltarg)
    print('latiobvt')
    print(latiobvt)
    print('longobvt')
    print(longobvt)
    print('strgtimeobvtyear')
    print(strgtimeobvtyear)
    print('listdelttimeobvtyear')
    summgene(listdelttimeobvtyear)
    print('heigobvt')
    print(heigobvt)
    if heigobvt is None:
        heigobvt = 0.

    # location object for the observatory
    objtlocaobvt = astropy.coordinates.EarthLocation(lat=latiobvt*astropy.units.deg, lon=longobvt*astropy.units.deg, height=heigobvt*astropy.units.m)
    
    # time object for the year
    timeyear = astropy.time.Time(strgtimeobvtyear).jd + listdelttimeobvtyear
    objttimeyear = astropy.time.Time(timeyear, format='jd', location=objtlocaobvt)
    timeside = objttimeyear.sidereal_time('mean')
    
    # delt time arry for night
    timedeltscal = 0.5
    timedelt = np.linspace(-12., 12. - timedeltscal, int(24. / timedeltscal))
    
    # frame object for the observatory during the year
    objtframobvtyear = astropy.coordinates.AltAz(obstime=objttimeyear, location=objtlocaobvt)
    
    # alt-az coordinate object for the target
    objtcoorplanalazyear = astropy.coordinates.SkyCoord(ra=rasctarg, dec=decltarg, frame='icrs', unit='deg').transform_to(objtframobvtyear)
    
    # air mass of the target during the time interval
    massairr = objtcoorplanalazyear.secz
    
    return massairr


# statistics
def cdfn_self(parascal, minmpara, maxmpara):
    '''
    Return the CDF transform of a variable drawrn from the uniform distribution
    '''
    paraunit = (parascal - minmpara) / (maxmpara - minmpara)
    
    return paraunit


def cdfn_powr(parascal, minm, maxm, slop):
        
    paraunit = (parascal**(1. - slop) - minm**(1. - slop)) / (maxm**(1. - slop) - minm**(1. - slop))
    
    return paraunit


def cdfn_dpow(parascal, minm, maxm, brek, sloplowr, slopuppr):
    
    if np.isscalar(parascal):
        parascal = np.array([parascal])
    
    faca = 1. / (brek**(sloplowr - slopuppr) * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) / (1. - sloplowr) + \
                                                            (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    facb = faca * brek**(sloplowr - slopuppr) / (1. - sloplowr)

    paraunit = np.empty_like(parascal)
    indxlowr = np.where(parascal <= brek)[0]
    indxuppr = np.where(parascal > brek)[0]
    
    if indxlowr.size > 0:
        paraunit[indxlowr] = facb * (parascal[indxlowr]**(1. - sloplowr) - minm**(1. - sloplowr))
    if indxuppr.size > 0:
        paraunitbrek = facb * (brek**(1. - sloplowr) - minm**(1. - sloplowr))
        paraunit[indxuppr] = paraunitbrek + faca / (1. - slopuppr) * (parascal[indxuppr]**(1. - slopuppr) - brek**(1. - slopuppr))
    
    return paraunit


def cdfn_igam(parascal, slop, cutf):
    
    paraunit = sp.stats.invgamma.cdf(parascal, slop - 1., scale=cutf)
    
    return paraunit


def cdfn_expo(parascal, maxm, scal):
    '''
    Return the CDF transform of a variable drawrn from the exponential distribution
    '''

    paraunit = (1. - np.exp(-parascal / maxm)) / (1. - np.exp(-maxm / scal))

    return paraunit


def cdfn_dexp(parascal, maxm, scal):
    
    if parascal < 0.:
        paraunit = cdfn_expo(-parascal, maxm, scal)
    else:
        paraunit = cdfn_expo(parascal, maxm, scal)
    
    return paraunit


def cdfn_lnor(parascal, meanpara, stdvpara):
   
    paraunit = cdfn_gaus(np.log(parascal), np.log(meanpara), stdvpara)
    
    return paraunit


def cdfn_gaus(parascal, meanpara, stdvpara):
   
    paraunit = 0.5  * (1. + sp.special.erf((parascal - meanpara) / np.sqrt(2) / stdvpara))
    
    return paraunit


def cdfn_lgau(parascal, mean, stdv):
    
    paraunit = cdfn_gaus(np.log(parascal), np.log(mean), stdv)

    return paraunit


def cdfn_eerr(parascal, meanpara, stdvpara, paraunitnormminm, paraunitnormdiff):
    
    tranpara = (parascal - meanpara) / stdvpara
    paraunitnormpara = 0.5 * (sp.special.erf(tranpara / np.sqrt(2.)) + 1.)
    paraunit = (paraunitnormpara - paraunitnormminm) / paraunitnormdiff

    return paraunit


def cdfn_logt(parascal, minmpara, maxmpara):
    
    paraunit = (np.log(parascal) - np.log(minmpara)) / (np.log(maxmpara) - np.log(minmpara))
    
    return paraunit


def cdfn_atan(parascal, minmpara, maxmpara):
    
    paraunit = (np.arctan(parascal) - np.arctan(minmpara)) / (np.arctan(maxmpara) - np.arctan(minmpara))
    
    return paraunit


def icdf_self(paraunit, minm, maxm):
    
    para = paraunit * (maxm - minm) + minm
    
    return para


def icdf_powr(paraunit, minm, maxm, slop):
    '''
    Inverse CDF for P(para) = (1 - slop) / (maxmpara^(1 - slop) - minmpara^(1 - slop)) * para^(-slop)
    '''
    
    para = (paraunit * (maxm**(1. - slop) - minm**(1. - slop)) + minm**(1. - slop))**(1. / (1. - slop))
    
    return para


def samp_powr(numbsamp, minm, maxm, slop):

    paraunit = np.random.rand(numbsamp)
    
    para = icdf_powr(paraunit, minm, maxm, slop)
    
    return para


def icdf_dpow(paraunit, minm, maxm, brek, sloplowr, slopuppr):
    
    if np.isscalar(paraunit):
        paraunit = np.array([paraunit])
    
    faca = 1. / (brek**(sloplowr - slopuppr) * (brek**(1. - sloplowr) - minm**(1. - sloplowr)) \
                                / (1. - sloplowr) + (maxm**(1. - slopuppr) - brek**(1. - slopuppr)) / (1. - slopuppr))
    facb = faca * brek**(sloplowr - slopuppr) / (1. - sloplowr)

    para = np.empty_like(paraunit)
    paraunitbrek = facb * (brek**(1. - sloplowr) - minm**(1. - sloplowr))
    indxlowr = np.where(paraunit <= paraunitbrek)[0]
    indxuppr = np.where(paraunit > paraunitbrek)[0]
    if indxlowr.size > 0:
        para[indxlowr] = (paraunit[indxlowr] / facb + minm**(1. - sloplowr))**(1. / (1. - sloplowr))
    if indxuppr.size > 0:
        para[indxuppr] = ((1. - slopuppr) * (paraunit[indxuppr] - paraunitbrek) / faca + brek**(1. - slopuppr))**(1. / (1. - slopuppr))
    
    return para


def icdf_dexp(paraunit, maxm, scal):
    
    if paraunit < 0.5:
        icdf = -icdf_expo(2. * paraunit, maxm, scal)
    else:
        icdf = icdf_expo(2. * (paraunit - 0.5), maxm, scal)
    
    return icdf


def icdf_igam(xdat, slop, cutf):
    
    icdf = sp.stats.invgamma.ppf(xdat, slop - 1., scale=cutf)
    
    return icdf


def icdf_lgau(paraunit, mean, stdv):
    
    icdf = np.exp(icdf_gaus(paraunit, np.log(mean), stdv))

    return icdf


def icdf_lnor(paraunit, meanpara, stdvpara):
    
    para = np.exp(icdf_gaus(paraunit, np.log(meanpara), stdvpara))

    return para


def icdf_gaus(paraunit, meanpara, stdvpara):
    
    para = meanpara + stdvpara * np.sqrt(2) * sp.special.erfinv(2. * paraunit - 1.)

    return para


def icdf_expo(paraunit, maxm, scal):

    para = -scal * np.log(1. - paraunit * (1. - np.exp(-maxm / scal)))

    return para


def icdf_eerr(paraunit, meanpara, stdvpara, paraunitnormminm, paraunitnormdiff):
    
    paraunitnormpara = paraunit * paraunitnormdiff + paraunitnormminm
    tranpara = sp.special.erfinv(2. * paraunitnormpara - 1.) * np.sqrt(2)
    para = tranpara * stdvpara + meanpara
   
    return para


def icdf_atan(paraunit, minmpara, maxmpara):

    para = np.tan((np.arctan(maxmpara) - np.arctan(minmpara)) * paraunit + np.arctan(minmpara))
    
    return para


def retr_llikgaustrun(para, gdat):
    
    meantria = para[0]
    stdvtria = para[1]
    
    meanresu, stdvresu = scipy.stats.truncnorm.stats(gdat.a, gdat.b, loc=meantria, scale=stdvtria, moments='mv')
    resimean = (meanresu - gdat.meanpara) / gdat.meanpara
    resistdv = (stdvresu - gdat.stdvpara) / gdat.stdvpara
    cost = -10. * (resimean**2 + resistdv**2)
    
    return cost


def samp_gaustrun(numbsamp, meanpara, stdvpara, minmpara, maxmpara):
    '''
    Sample from a truncated Gaussian
    '''

    if not np.isfinite(meanpara).any():
        raise Exception('')

    if not np.isfinite(stdvpara).any():
        raise Exception('')

    # scaled minimum and maximum of the truncated Gaussian
    #gdat = gdatstrt()
    #gdat.a = (minmpara - meanpara) / stdvpara
    #gdat.b = (maxmpara - meanpara) / stdvpara
    
    #gdat.meanpara = meanpara
    #gdat.stdvpara = stdvpara

    # solve for the Gaussian mean and standard deviation that will produce the desired mean and standard deviation shen truncated
    #from scipy.optimize import minimize
    #obtjmini = minimize(cost_gaustrun, (meanpara, stdvpara), args=(a, b, meanpara, stdvpara), tol=1e-6, method='Nelder-Mead')
    #meanresu, stdvresu = obtjmini.x
   
    ## list of parameter names
    #listnamepara = ['meantria', 'stdvtria']
    ## list of parameter labels and units
    #listlablpara = [['$\mu$', ''], ['$\sigma$', '']]
    ## list of parameter scalings
    #listscalpara = ['self', 'self']
    ## list of parameter minima
    #listminmpara = [0., 0.]
    ## list of parameter maxima
    #listmaxmpara = [3. * meanpara, 3. * stdvpara]
        
    ##gdat.typemodl = typemodl
    #path = '/Users/tdaylan/Desktop/'
    #path = None
    #dictparafitt, dictvarbderi = samp(gdat, path, 100, retr_llikgaustrun, \
    #                                            listnamepara, listlablpara, listscalpara, listminmpara, listmaxmpara, \
    #                                            numbsampburnwalk=80)#, strgextn=gdat.strgextn)
    #meantria = np.median(dictparafitt['meantria'])
    #stdvtria = np.median(dictparafitt['stdvtria'])
    
    a = (minmpara - meanpara) / stdvpara
    b = (maxmpara - meanpara) / stdvpara
    meantria = meanpara
    stdvtria = stdvpara
    # sample from truncated Gaussian
    para = scipy.stats.truncnorm.rvs(a, b, loc=meantria, scale=stdvtria, size=numbsamp)
    
    return para



def retr_pctlvarb(listpara):

    shap = np.zeros(len(listpara.shape), dtype=int)
    shap[0] = 3
    shap[1:] = listpara.shape[1:]
    shap = list(shap)
    postvarb = np.zeros(shap)
    
    postvarb[0, ...] = np.percentile(listpara, 50., axis=0)
    postvarb[1, ...] = np.percentile(listpara, 16., axis=0)
    postvarb[2, ...] = np.percentile(listpara, 84., axis=0)

    return postvarb


def prnt_list(listpara, strg):
    
    medi = np.median(listpara)
    lowr = np.percentile(listpara, 16.)
    uppr = np.percentile(listpara, 84.)
    print('%s: %.3g +%.3g -%.3g' % (strg, medi, medi - lowr, uppr - medi))


def retr_errrvarb(inpt, samp=False):

    if samp:
        postvarb = retr_pctlvarb(inpt)
    else:
        postvarb = inpt

    errr = np.abs(postvarb[0, ...] - postvarb[1:3, ...])

    return errr


def retr_kdegpdfn(listsamp, binsvarb, stdv):
    
    meanvarb = (binsvarb[1:] + binsvarb[:-1]) / 2.
    deltvarb = (binsvarb[1:] - binsvarb[:-1]) / 2.
    kdeg = retr_kdeg(listsamp, meanvarb, stdv)
    kdegpdfn = kdeg / np.sum(kdeg)
    kdegpdfn /= deltvarb
    
    return kdeg


def retr_kdeg(listsamp, varb, stdv):
    
    if np.isscalar(varb):
        varb = np.array([varb])
    kdeg = np.sum(np.exp(-0.5 * (varb[None, :] - listsamp[:, None])**2 / stdv**2), axis=0)
    
    return kdeg


def plot_recaprec( \
                  # base path for placing the plots
                  pathimag, \
                  
                  # a string describing the run
                  strgextn, \
                  
                  # list of parameters for relevant samples
                  listpararele=None, \
                  
                  # list of parameters for analyzed samples
                  listparaanls=None, \
                  
                  # list of parameter names for relevant samples
                  listnamepararele=None, \
                  
                  # list of parameter names for analyzed samples
                  listnameparaanls=None, \

                  # list of parameter names for relevant samples
                  listlablvarbrele=None, \

                  # list of parameter names for analyzed samples
                  listlablvarbanls=None, \
                  
                  # list of Booleans for all relevant samples indicating whether they are positive
                  boolposirele=None, \

                  # list of Booleans for all positive samples indicating whether they are relevant
                  boolreleposi=None, \

                  strgplotextn='pdf', \
                  typeverb=1, \
                  numbbins=10, \
                  strgreca='Recall', \
                  listparadete=None, \
                 ):
    
    if isinstance(boolreleposi, list):
        boolreleposi = np.array(boolreleposi)

    if isinstance(boolposirele, list):
        boolposirele = np.array(boolposirele)
    
    # sanity checks
                                #listpararele is not None and boolposirele.size != listpararele.shape[0] or \
    if np.sum(boolposirele) != np.sum(boolreleposi) or \
                                listparaanls is not None and boolreleposi.size != listparaanls.shape[0] or \
                                isinstance(listpararele, list) or isinstance(listparaanls, list) or \
                                listlablvarbrele is not None and len(listlablvarbrele) != listpararele.shape[1] or \
                                listlablvarbanls is not None and len(listlablvarbanls) != listparaanls.shape[1]:
        print('listlablvarbrele')
        print(listlablvarbrele)
        print('listpararele')
        summgene(listpararele)
        print('listlablvarbanls')
        print(listlablvarbanls)
        print('listparaanls')
        summgene(listparaanls)
        print('boolposirele')
        summgene(boolposirele)
        print('boolreleposi')
        summgene(boolreleposi)
        print('np.sum(boolreleposi)')
        print(np.sum(boolreleposi))
        print('np.sum(boolposirele)')
        print(np.sum(boolposirele))
        #raise Exception('')

    # get the labels and scalings if not already defined
    if listlablvarbrele is None:
        listlablvarbrele, _, _, _, _ = retr_listlablscalpara(listnamepararele)
    
    if listlablvarbanls is None:
        listlablvarbanls, _, _, _, _ = retr_listlablscalpara(listnameparaanls)
    
    numbvarbreca = listpararele.shape[1]
    indxvarbreca = np.arange(numbvarbreca)
    listlablvarbreletotl = retr_labltotl(listlablvarbrele)
    
    if listparaanls is not None:
        numbvarbprec = listparaanls.shape[1]
        indxvarbprec = np.arange(numbvarbprec)
        listlablvarbanlstotl = retr_labltotl(listlablvarbanls)
    
    if typeverb > 1:
        print('numbbins')
        print(numbbins)
    
    if isinstance(boolposirele, list):
        raise Exception('')

    indxbins = np.arange(numbbins)
    indxiter = range(2)
    
    metr = np.zeros((3, numbbins)) + np.nan
    for c in indxiter:
        
        if listparaanls is None and c == 1:
            continue
        
        if c == 0:
            if boolposirele.size == 0:
                continue
            listpara = listpararele
            listnamepara = listnamepararele
            listlablvarbtemp = listlablvarbreletotl
            strgmetr = 'reca'
            lablyaxi = strgreca + ' [\%]'
        if c == 1:
            if boolreleposi.size == 0:
                continue
            listpara = listparaanls
            listnamepara = listnameparaanls
            listlablvarbtemp = listlablvarbanlstotl
            strgmetr = 'prec'
            lablyaxi = 'Precision [\%]'
        if c == 2:
            listpara = listpararele
            listnamepara = listnamepararele
            listlablvarbtemp = listlablvarbreletotl
            strgmetr = 'occu'
            lablyaxi = 'Occurence rate'
            
        for k, namepara in enumerate(listnamepara):
        
            bins = np.linspace(np.nanmin(listpara[:, k]), np.nanmax(listpara[:, k]), numbbins + 1)
            meanvarb = (bins[1:] + bins[:-1]) / 2.
            
            if not np.isfinite(bins).all():
                print('np.nanmin(listpara[:, k])')
                print(np.nanmin(listpara[:, k]))
                raise Exception('')

            if c == 2 and listparadete is not None:
                boolupprlimt = np.zeros(numbbins, dtype=bool)
            
            for a in indxbins:
                
                if c == 0 or c == 1:
                    indx = np.where((bins[a] < listpara[:, k]) & (listpara[:, k] < bins[a+1]))[0]
                    numb = indx.size
                    if numb > 0:
                        # recall
                        if c == 0:
                            metr[c, a] = float(np.sum(boolposirele[indx].astype(float))) / numb
                        # precision
                        if c == 1:
                            metr[c, a] = float(np.sum(boolreleposi[indx].astype(float))) / numb
                    
                if c == 2 and listparadete is not None:
                    if len(listparadete) == 0:
                        numb = 0
                    else:
                        indx = np.where((bins[a] < listparadete[:, k]) & (listparadete[:, k] < bins[a+1]))[0]
                        numb = indx.size
                    if numb > 0:
                        boolupprlimt[a] = False
                        metr[c, a] = numb / listpararele[k].size * metr[1, a] / metr[0, a]
                    else:
                        boolupprlimt[a] = True
                        metr[c, a] = 1. / listpararele[k].size * metr[1, a] / metr[0, a]
            
            if (c == 0 or c == 1) and listparadete is None:
                figr, axis = plt.subplots(figsize=(6, 4))
                axis.plot(meanvarb, metr[c, :] * 100.)
                axis.set_xlabel(listlablvarbtemp[k])
                axis.set_ylabel(lablyaxi)
                path = pathimag + '%s_%s_%s.%s' % (strgmetr, namepara, strgextn, strgplotextn) 
                print('Writing to %s...' % path)
                plt.savefig(path)
                plt.close()
            
            if c == 2 and listparadete is not None:
                figr, axis = plt.subplots(figsize=(6, 4))
                axis.errorbar(meanvarb, metr[c, :], marker='o', uplims=boolupprlimt)
                axis.set_ylabel('Occurence rate')
                axis.set_xlabel(listlablvarbreletotl[k])
                path = pathimag + 'occu_%s_%s.%s' % (namepara, strgextn, strgplotextn) 
                print('Writing to %s...' % path)
                plt.savefig(path)
                plt.close()


def prep_mask(data, epoc=None, peri=None, duramask=None, limttime=None):
    '''
    Read the data, mask out the transits
    '''
    
    time = data[:, 0]
    
    if epoc is not None:
        listindxtimethis = []
        for n in range(-2000, 2000):
            timeinit = epoc + n * peri - duramask / 24.
            timefinl = epoc + n * peri + duramask / 24.
            indxtimethis = np.where((time > timeinit) & (time < timefinl))[0]
            listindxtimethis.append(indxtimethis)
        indxtimethis = np.concatenate(listindxtimethis)
    else:
        indxtimethis = np.where((time < limttime[1]) & (time > limttime[0]))[0]

    numbtime = time.size
    indxtime = np.arange(numbtime)
    
    indxtimegood = np.setdiff1d(indxtime, indxtimethis)
    
    dataoutp = data[indxtimegood, :]
    
    return dataoutp


def retr_nfwp(nfwg, numbside, norm=None):
    
    edenlocl = 0.3 # [GeV/cm^3]
    radilocl = 8.5 # [kpc]
    rscl = 23.1 # [kpc]
    
    nradi = 100
    minmradi = 1e-2
    maxmradi = 1e2
    radi = np.logspace(np.log10(minmradi), np.log10(maxmradi), nradi)
    
    nsadi = 100
    minmsadi = 0.
    maxmsadi = 2. * radilocl
    sadi = np.linspace(minmsadi, maxmsadi, nsadi)
    
    lghp, bghp, numbpixl, apix = retr_healgrid(numbside)
    
    cosigahp = cos(np.deg2rad(lghp)) * cos(np.deg2rad(bghp))
    gahp = np.rad2deg(arccos(cosigahp))
    
    eden = 1. / (radi / rscl)**nfwg / (1. + radi / rscl)**(3. - nfwg)
    eden *= edenlocl / interp1d(radi, eden)(radilocl)
    
    edengrid = np.zeros((nsadi, numbpixl))
    for i in range(nsadi):
        radigrid = np.sqrt(radilocl**2 + sadi[i]**2 - 2 * radilocl * sadi[i] * cosigahp)
        edengrid[i, :] = interp1d(radi, eden)(radigrid)

    edengridtotl = sum(edengrid**2, axis=0)

    if norm != None:
        jgahp = np.argsort(gahp)
        edengridtotl /= interp1d(gahp[jgahp], edengridtotl[jgahp])(5.)
        
    return edengridtotl


def retr_lablmexp(numb):
    '''
    Return a string for a number including the mantissa and exponent
    '''
    
    if numb == 0.:
        strg = '$0$'
    else:
        logn = np.log10(np.fabs(numb))
        expo = np.floor(logn)
        expo = int(expo)
        mant = 10**(logn - expo) * numb / np.fabs(numb)
        
        if np.fabs(numb) >= 1e4 or np.fabs(numb) < 1e-3:
            if mant == 1. or mant == -1.:
                strg = r'$10^{%d}$' % expo
            else:
                strg = r'$%.3g \times 10^{%d}$' % (mant, expo)
        else:
            strg = r'$%g$' % numb

    return strg


class varb(object):
    
    def __init__(self, numb=None):
        
        self.name = []
        self.strg = []
        self.para = []
        self.scal = []
        self.strg = []
        if numb != None:
            self.numb = numb

    
    def defn_para(self, name, minm, maxm, numb=None, strg='', scal='self'):
        
        if numb != None:
            numbtemp = numb
        else:
            numbtemp = self.numb

        if scal == 'logt':
            arry = np.logspace(np.log10(minm), np.log10(maxm), numbtemp) 
        if scal == 'self':
            arry = np.linspace(minm, maxm, numbtemp) 
        
        self.name.append(name)
        self.para.append(arry)
        self.scal.append(scal)
        self.strg.append(strg)
        self.size = len(self.para)

    
def retr_listscalpara(listnamepara):
    
    numbpara = len(listnamepara)
    indxpara = np.arange(numbpara)
    listscalpara = ['self' for k in indxpara]
    
    return listscalpara


def retr_labltotl(listlablpara):
    
    listlablparatotl = [[] for k in range(len(listlablpara))]
    for k in range(len(listlablpara)):
        listlablparatotl[k] = retr_labltotlsing(listlablpara[k][0], listlablpara[k][1])
   
    return listlablparatotl


def retr_labltotlsing(labl, lablunit):
    
    if lablunit == '' or lablunit is None:
        labltotl = '%s' % labl
    else:
        labltotl = '%s [%s]' % (labl, lablunit)

    return labltotl


def retr_factconv():
    '''
    Make a dictionary of relevant conversion factors.
    '''
    
    dictfact = dict()

    dictfact['rsrj'] = 9.731
    dictfact['rjre'] = 11.21
    dictfact['rsre'] = dictfact['rsrj'] * dictfact['rjre']
    dictfact['msmj'] = 1048.
    dictfact['mjme'] = 317.8
    dictfact['msme'] = dictfact['msmj'] * dictfact['mjme']
    
    dictfact['pcau'] = 206265.
    
    # Astronomical Unit in Solar radius
    dictfact['aurs'] = 215.
    
    # Astronomical Unit in kilometers
    dictfact['aukm'] = 149597870700

    return dictfact


def retr_dictturk():
    '''
    Dictionary for English-Turkish translation.
    '''

    dictturk = dict()
    
    dictturk['Planet candidate'] = 'Gezegen adayı'
    dictturk['False positive'] = 'Yanlış pozitif'
    dictturk['Known planet'] = 'Bilinen gezegen'
    dictturk['Confirmed planet'] = 'Doğrulanmış gezegen'
    
    dictturk['Number of TOIs'] = 'TOI sayısı'
    dictturk['Number of exoplanets'] = 'Ötegezegen sayısı'

    dictturk['Faint-star search during PM'] = 'Sönük-yıldız taraması, ana görev'
    dictturk['Faint-star search during EM1'] = 'Sönük-yıldız taraması, uzatma görev 1'
    dictturk['All'] = 'Tüm'
    dictturk['Detected'] = 'Keşfedilen'
    dictturk['Exoplanets'] = 'Ötegezegenler'
    dictturk['Stellar age'] = 'Yıldızın yaşı'
    dictturk['Age of the host star'] = 'Barınak yıldızının yaşı'
    dictturk['Gyr'] = 'Milyar yıl'
    dictturk['Exoplanet Detections'] = 'Ötegezegen Keşifleri'
    dictturk['Transit'] = 'Geçiş'
    dictturk['Radial Velocity'] = 'Dikine hız'
    dictturk['Other'] = 'Diğer'
    dictturk['Exoplanet detections other than those via transits and radial velocity'] = 'Geçiş ve dikine hız dışındaki ötegezegen keşifleri'
    dictturk['Imaging'] = 'Görüntüleme'
    dictturk['Microlensing'] = 'Mikromerceklenme'
    dictturk['Transiting'] = 'Geçiş Yapan'
    dictturk['exoplanet'] = 'ötegezegen'
    dictturk['Exoplanet'] = 'Ötegezegen'
    dictturk['Exoplanets with precise density'] = 'Yoğunluk ölçümü hassas olan ötegezegenler'
    dictturk['Exoplanets with precise mass'] = 'Kütle ölçümü hassas olan ötegezegenler'
    dictturk['Exoplanets with weak mass'] = 'Kütle ölçümü zayıf olan ötegezegenler'
    dictturk['Young'] = 'Genç'
    dictturk['Old'] = 'Yaşlı'
    dictturk['Discovery Year'] = 'Keşif Yılı'
    dictturk['Discovered by TESS'] = 'TESS tarafından keşfedilen'
    dictturk['Discovered by Kepler'] = 'Kepler tarafından keşfedilen'
    dictturk['With precise mass and radius'] = 'Hassas kütle ve yarıçap ölçümüne sahip'
                
    dictturk['Sun-like host'] = 'Güneş-benzeri barınak yıldızı'
    dictturk['Habitable zone'] = 'Yaşanabilir bölge'
    dictturk['Hot'] = 'Sıcak'
    dictturk['Low irradiation'] = 'Az derecede aydınlatılmış'
    dictturk['Medium irradiation'] = 'Orta derecede aydınlatılmış'
    dictturk['High irradiation'] = 'Yüksek derecede aydınlatılmış'
    
    dictturk['Faint-star search'] = 'Sönük-yıldız taraması'
    dictturk['faint-star search'] = 'sönük-yıldız taraması'
    dictturk['Other TOIs'] = 'Diğerleri'
    dictturk['Planetary radius'] = 'Gezegen yarıçapı'
    dictturk['Orbital period'] = 'Yörünge periyodu'
    dictturk['Distance'] = 'Uzaklık'
    dictturk['Transit depth'] = 'Geçiş derinliği'
    
    dictturk['Irradiation'] = 'Aydınlatma'
    
    dictturk['TOI'] = 'TOI'
    dictturk["TOIs"] = 'TESS İlginç Nesneleri (TOI)'
    
    dictturk['Transiting exoplanets discovered by TESS'] = 'TESS tarafından keşfedilen geçiş yapan ötegezegenler'
    dictturk['Bright'] = 'Parlak'
    dictturk['Bright and small'] = 'Parlak ve küçük'
    dictturk['Bright, small, and visible from LCO'] = "Parlak, küçük ve LCO'dan gözlemlenebilir"
    dictturk['Bright, small, visible from LCO, and favorable for AC'] = "Parlak, küçük ve LCO'dan gözlemlenebilir ve atmosfer nitelendirmesine elverişli"
    dictturk['High TSM/ESM'] = 'Yüksek TSM/ESM'
    dictturk['High TSM/ESM \& Precise mass'] = 'Yüksek TSM/ESM ve hassas kütle ölçümü'
    dictturk['High TSM/ESM \& Weak mass'] = 'Yüksek TSM/ESM ve zayıf kütle ölçümü'
    dictturk['Precise mass'] = 'Hassas Kütle Ölçümü'
    dictturk['Weak mass'] = 'Zayıf Kütle Ölçümü'
    dictturk['Time from midtransit'] = 'Geçiş ortasına göre zaman'
    dictturk['hour'] = 'saat'
    dictturk['days'] = 'gün'
    dictturk['degree'] = 'derece'
    dictturk['Right Ascension'] = 'Bahar açısı'
    dictturk['Declination'] = 'Yükselim'
    dictturk['Astrometry'] = 'Astrometri'
    dictturk['Orbital Brightness Modulation'] = 'Yörünge parlaklık değişimi'
    dictturk['Transit Timing Variations'] = 'Geçiş zamanlama değişimleri'
    dictturk['Eclipse Timing Variations'] = 'Örtme zamanlama değişimleri'
    
    return dictturk


def sign_code(axis, typesigncode):
    '''
    Sign the axis to indicate the generating code.
    ''' 

    bbox = dict(boxstyle='round', ec='k', fc='white')
    axis.text(0.97, 0.05, r'github.com/tdaylan/\textbf{%s}' % (typesigncode), bbox=bbox, \
                                                                        transform=axis.transAxes, color='firebrick', ha='right', size='small')


def retr_listlablscalpara(listnamepara, listlablpara=None, dictdefa=None, booldiag=True, typelang='English', boolmath=False):
    
    if dictdefa is not None:
        if not isinstance(dictdefa, dict):
            raise Exception('')
    
    if booldiag:
        if listlablpara is not None and len(listnamepara) != len(listlablpara):
            print('')
            print('')
            print('listlablpara is not None and len(listnamepara) != len(listlablpara)')
            print('listnamepara')
            print(listnamepara)
            print('listlablpara')
            print(listlablpara)
            raise Exception('')
                
    numbpara = len(listnamepara)
    indxpara = np.arange(numbpara)
    
    if listlablpara is None:
        listlablpara = [[] for k in indxpara]
    listlabltotlpara = [[] for k in indxpara]
    listlablrootpara = [[] for k in indxpara]
    listlablunitpara = [[] for k in indxpara]
    listscalpara = [[] for k in indxpara]
    for k in indxpara:
        if listlablpara[k] is not None and len(listlablpara[k]) > 0:
            continue
        listlablpara[k] = [None, None]
        if dictdefa is not None and listnamepara[k] in dictdefa:
            listlablpara[k] = dictdefa[listnamepara[k]]['labl']
            listscalpara[k] = dictdefa[listnamepara[k]]['scal']
        elif listnamepara[k] == 'toii':
            listlablpara[k] = ['TOI ID', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'HostRedshift':
            listlablpara[k] = ['HostRedshift', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'Redshift':
            listlablpara[k] = ['Redshift', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'ID':
            listlablpara[k] = ['ID', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'Sector':
            listlablpara[k] = ['Sector', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'nameplan':
            listlablpara[k] = ['Planet', '']
            listscalpara[k] = ''
        # TOI
        elif listnamepara[k] == 'numbobsvtime':
            listlablpara[k] = ['Number of time-series observations', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'numbobsvspec':
            listlablpara[k] = ['Number of spectroscopic observations', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'numbobsvimag':
            listlablpara[k] = ['Number of imaging observations', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'yearaler':
            listlablpara[k] = ['Time of Alert', '']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'tolerrat':
            listlablpara[k] = ['$\delta_{rr}$', '']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'namestar':
            listlablpara[k] = ['Star', '']
            listscalpara[k] = 'self'

        # discovery magnitude
        elif listnamepara[k] == 'magtdisc':
            listlablpara[k] = [r'$m_{disc}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'laecstar':
            listlablpara[k] = [r'$\beta$', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'loecstar':
            listlablpara[k] = [r'$\lambda$', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'rasc' or listnamepara[k] == 'rascstar':
            listlablpara[k] = ['Right Ascension', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'decl' or listnamepara[k] == 'declstar':
            listlablpara[k] = ['Declination', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'lgal' or listnamepara[k] == 'lgalstar':
            listlablpara[k] = ['$l$', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'bgal' or listnamepara[k] == 'bgalstar':
            listlablpara[k] = ['$b$', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'yeardisc':
            listlablpara[k] = ['Discovery Year', '']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'duraobsv':
            listlablpara[k] = ['Duration of Observation', 'days']
            listscalpara[k] = 'logt'
        
        elif listnamepara[k] == 'amplflar':
            listlablpara[k] = ['Flare amplitude', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'rateflar':
            listlablpara[k] = ['Flare rate', 'day$^{-1}$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'fwhmflar':
            listlablpara[k] = ['FWHM of flare', 'seconds']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'enerflar':
            listlablpara[k] = ['Energy of flare', 'erg']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'amplstar':
            listlablpara[k] = ['$A_{f}$', '']
            listscalpara[k] = 'logt'
        
        elif listnamepara[k] == 'areasrch':
            listlablpara[k] = ['Area to be searched', 'degree$^2$']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'ratimass':
            listlablpara[k] = ['$q$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'radiplan':
            if boolmath:
                listlablpara[k][0] = '$R_p$'
            else:
                listlablpara[k][0] = 'Planetary radius'
            listlablpara[k][1] = '$R_\oplus$'
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'stdvradiplan':
            listlablpara[k] = ['$\sigma_{R_p}$', '$R_\oplus$']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'radistar':
            listlablpara[k] = ['$R_\star$', '$R_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'tagestar':
            listlablpara[k] = ['Age of the host star', 'Gyr']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'massplan':
            listlablpara[k] = ['$M_p$', '$M_\oplus$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'stdvmassplan':
            listlablpara[k] = ['$\sigma_{M_p}$', '$M_\oplus$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'masssyst':
            listlablpara[k] = ['$M_{sys}$', '$M_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'massstar':
            listlablpara[k] = ['$M_\star$', '$M_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'stdvmassstar':
            listlablpara[k] = ['$\sigma_{M_\star}$', '$M_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'massbhol':
            listlablpara[k] = ['$M_{BH}$', '$M_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'densstar':
            listlablpara[k] = ['$d_\star$', 'g cm$^{-3}$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'tmptplan':
            listlablpara[k] = ['$T_p$', 'K']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'lumistar':
            listlablpara[k] = ['$L_\star$', '$L_{\odot}$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'tmptstar':
            listlablpara[k] = ['$T_\star$', 'K']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'loggstar':
            listlablpara[k] = ['$\log g$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'metastar':
            listlablpara[k] = ['[M/H]', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'esmm':
            listlablpara[k] = ['Emission Spectroscopy Metric (ESM)', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'tsmm':
            listlablpara[k] = ['Transmission Spectroscopy Metric (TSM)', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'esmmacwg':
            listlablpara[k] = ['ESM$_{ACWG}$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'tsmmacwg':
            listlablpara[k] = ['TSM$_{ACWG}$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'densplan' or listnamepara[k] == 'denscomp':
            listlablpara[k] = ['d', 'g cm$^{-3}$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'epocmtratess':
            listlablpara[k] = ['$T_0$', 'BJD-2457000']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'offsphascomp':
            listlablpara[k] = ['$\phi_{off}$', '$^\circ$']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'epocmtra':
            listlablpara[k] = ['$T_0$', 'BJD']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'timesupn':
            listlablpara[k] = ['$T_0$', 'BJD']
            listscalpara[k] = 'self'
        
        # predicted radial velocity semi-amplitude
        elif listnamepara[k] == 'rvelsemapred':
            listlablpara[k] = ['Predicted RV semi-amplitude', 'm s$^{-1}$']
            listscalpara[k] = 'logt'
        
        # time-series constants
        elif listnamepara[k].startswith('cons') or listnamepara[k].startswith('timestep') or listnamepara[k].startswith('scalstep') or \
                                                listnamepara[k].startswith('coeflmdk') or listnamepara[k].startswith('rratcom'):
            ## white
            if listnamepara[k].endswith('whit'):
                strg = ',w'
            ## energy-dependent
            elif listnamepara[k][-4:-2] == 'en' and listnamepara[k][-2:].isnumeric():
                strg = ',%d' % int(listnamepara[k][-2:])
            ## single-wavelength
            else:
                strg = ''
            
            if listnamepara[k].startswith('consfrst'):
                listlablpara[k] = ['$C_{1%s}$' % strg, '']
            elif listnamepara[k].startswith('consseco'):
                listlablpara[k] = ['$C_{2%s}$' % strg, '']
            elif listnamepara[k].startswith('cons'):
                listlablpara[k] = ['$C_{%s}$' % strg, '']
            # the time at which a time-series baseline has a step (discontinuity)
            elif listnamepara[k].startswith('timestep'):
                listlablpara[k] = ['$t_{b%s}$' % strg, '']
            # the slope of the step
            elif listnamepara[k].startswith('scalstep'):
                listlablpara[k] = ['$A_{b%s}$' % strg, '']
            elif listnamepara[k].startswith('rratcomp'):
                listlablpara[k] = ['$R_{p%s}/R_\star$' % strg, '']
            elif listnamepara[k].startswith('rratcom'):
                listlablpara[k] = ['$R_{%d%s}/R_\star$' % (int(listnamepara[k][7]), strg), '']
            # linear limb-darkening coefficient
            elif listnamepara[k].startswith('coeflmdklinr'):
                listlablpara[k] = ['$u_{l%s}$' % strg, '']
            # quadratic limb-darkening coefficient
            elif listnamepara[k].startswith('coeflmdkquad'):
                listlablpara[k] = ['$u_{q%s}$' % strg, '']
            
            if listnamepara[k].startswith('rratcom'):
                listscalpara[k] = 'logt'
            else:
                listscalpara[k] = 'self'
            
        elif listnamepara[k] == 'limtvmag':
            listlablpara[k] = ['Limiting V Magnitude', '']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'timeexpo':
            listlablpara[k] = ['Exposure Time', 'seconds']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'tmagsyst' or listnamepara[k] == 'tmag' or listnamepara[k] == 'tmagstar':
            listlablpara[k] = ['TESS Magnitude', '']
            listscalpara[k] = 'self'
        # number of companions per star
        elif listnamepara[k] == 'numbcompstar':
            listlablpara[k] = ['$N_{c}$', '']
            listscalpara[k] = 'self'
        # mean number of companions per star
        elif listnamepara[k] == 'numbcompstarmean':
            listlablpara[k] = ['$\mu_{c}$', '']
            listscalpara[k] = 'self'
        # number of transiting companions per star
        elif listnamepara[k] == 'numbcomptranstar':
            listlablpara[k] = ['$N_{ct}$', '']
            listscalpara[k] = 'self'
        # number of planets per star
        elif listnamepara[k] == 'numbplanstar':
            listlablpara[k] = ['$N_{p}$', '']
            listscalpara[k] = 'self'
        # number of transiting planets per star
        elif listnamepara[k] == 'numbplantranstar':
            listlablpara[k] = ['$N_{pt}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'hmagsyst' or listnamepara[k] == 'hmag':
            listlablpara[k] = ['H', 'mag']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'kmagsyst' or listnamepara[k] == 'kmag':
            listlablpara[k] = ['K', 'mag']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'jmagsyst' or listnamepara[k] == 'jmag':
            listlablpara[k] = ['J', 'mag']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'vmagsyst' or listnamepara[k] == 'vmag':
            listlablpara[k] = ['V', 'mag']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'rmagsyst' or listnamepara[k] == 'rmag':
            listlablpara[k] = ['r', 'mag']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'dist' or listnamepara[k] == 'distsyst':
            if boolmath:
                listlablpara[k][0] = '$d$'
            else:
                listlablpara[k][0] = 'Distance'
            #listlablpara[k][1] = 'pc'
            listlablpara[k][1] = 'Mpc'
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'noisphot':
            listlablpara[k] = ['Photometric Noise', 'ppt']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'nois':
            listlablpara[k] = ['$\sigma$', 'ppt']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'proboccu':
            listlablpara[k] = ['$P_{occ}$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'boolfpos':
            listlablpara[k] = ['$B_{FP}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'booloccu':
            listlablpara[k] = ['$B_{occ}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'numbtsec':
            listlablpara[k] = ['$N_{sector}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'numbtran':
            listlablpara[k] = ['$N_{tr}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'incl' or listnamepara[k] == 'inclcomp':
            listlablpara[k] = ['$i$', 'degree']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'cosi':
            listlablpara[k] = ['$\cos i$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'booltran':
            listlablpara[k] = ['Whether transiting', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'irra':
            if boolmath:
                listlablpara[k][0] = '$I_{irr}$'
            else:
                listlablpara[k][0] = 'Irradiation'
            listlablpara[k][1] = '$I_{\oplus}$'
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'sdee' or listnamepara[k] == 'sdeepboxprim':
            listlablpara[k] = ['SDE', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'powrlspeprim':
            listlablpara[k] = ['$\eta_{LS,max}$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'peripboxprim':
            listlablpara[k] = ['$P_{BLS,max}$', 'days']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'perilspeprim':
            listlablpara[k] = ['$P_{LS,max}$', 'days']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'peri' or listnamepara[k] == 'pericomp' or listnamepara[k] == 'periplan':
            if boolmath:
                listlablpara[k][0] = '$P$'
            else:
                listlablpara[k][0] = 'Orbital period'
            listlablpara[k][1] = 'days'
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'masscomp':
            listlablpara[k] = ['$M_{comp}$', '$M_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'masstotl':
            listlablpara[k] = ['$M_{tot}$', '$M_\odot$']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'smax' or listnamepara[k] == 'smaxcomp':
            listlablpara[k] = ['$a$', 'AU']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'rsma':
            listlablpara[k] = ['$(R_\star+R_p)/a$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'rs2a':
            listlablpara[k] = ['$R_\star/a$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'dcyc':
            listlablpara[k] = ['$\Delta \phi_{tr}$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'amplslen':
            listlablpara[k] = ['$A_{SL}$', 'ppt']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'duratrancomp' or listnamepara[k] == 'duratrancomptotl' or listnamepara[k] == 'duratrantotl':
            listlablpara[k] = ['Total Transit Duration', 'hours']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'typebrgtcomp':
            listlablpara[k] = ['Brightness Type of the Companion', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'duratrancompfull':
            listlablpara[k] = ['Full Transit Duration', 'hours']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'radicomp':
            if boolmath:
                listlablpara[k][0] = '$R_c$'
            else:
                listlablpara[k][0] = 'Companion radius'
            listlablpara[k][1] = '$R_\oplus$'
            listscalpara[k] = 'logt'
        
        # orbital parameters for each component
        elif listnamepara[k][:7] == 'radicom' and len(listnamepara[k]) == 8 and (listnamepara[k][7].isnumeric() or listnamepara[k][7] == 'p'):
            listlablpara[k] = ['$R_{%s}$' % listnamepara[k][7], '']
            listscalpara[k] = 'logt'
        elif listnamepara[k][:7] == 'rratcom' and len(listnamepara[k]) == 8 and (listnamepara[k][7].isnumeric() or listnamepara[k][7] == 'p'):
            listlablpara[k] = ['$R_{%s}/R_\star$' % listnamepara[k][7], '']
            listscalpara[k] = 'self'
        elif listnamepara[k][:7] == 'rsmacom' and len(listnamepara[k]) == 8 and (listnamepara[k][7].isnumeric() or listnamepara[k][7] == 'p'):
            listlablpara[k] = ['$(R_\star+R_{%s})/a$' % listnamepara[k][7], '']
            listscalpara[k] = 'self'
        elif listnamepara[k][:7] == 'cosicom' and len(listnamepara[k]) == 8 and (listnamepara[k][7].isnumeric() or listnamepara[k][7] == 'p'):
            listlablpara[k] = ['$\cos i_{%s}$' % listnamepara[k][7], '']
            listscalpara[k] = 'self'
        elif listnamepara[k][:11] == 'epocmtracom' and len(listnamepara[k]) == 12 and (listnamepara[k][11].isnumeric() or listnamepara[k][11] == 'p'):
            listlablpara[k] = ['$T_{0;%s}$' % listnamepara[k][11], 'BJD']
            listscalpara[k] = 'self'
        elif listnamepara[k][:7] == 'pericom' and len(listnamepara[k]) == 8 and (listnamepara[k][7].isnumeric() or listnamepara[k][7] == 'p'):
            listlablpara[k] = ['$P_{%s}$' % listnamepara[k][7], 'days']
            listscalpara[k] = 'self'
        elif listnamepara[k][:11] == 'massplancom' and len(listnamepara[k]) == 12 and (listnamepara[k][7].isnumeric() or listnamepara[k][7] == 'p'):
            listlablpara[k] = ['$M_{%s}$' % listnamepara[k][11], '$M_\oplus$']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'masscomp':
            listlablpara[k] = ['$M_{C}$', '$M_\odot$']
            listscalpara[k] = 'self'
        
        elif listnamepara[k] == 'depttran' or listnamepara[k] == 'depttrancomp':
            if boolmath:
                listlablpara[k][0] = '$\delta_{tr}$'
            else:
                listlablpara[k][0] = 'Transit depth'
            listlablpara[k][1] = 'ppt'
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 's2nr':
            listlablpara[k] = ['S/N', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 's2nrblss':
            listlablpara[k] = ['S/N$_{BLS}$', '']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'booldete':
            listlablpara[k] = ['$B_{det}$', '']
            listscalpara[k] = 'self'
        elif listnamepara[k] == 'duraslen':
            listlablpara[k] = ['$D_{sl}$', 'hours']
            listscalpara[k] = 'logt'
        elif listnamepara[k] == 'tici':
            listlablpara[k] = ['TIC ID', '']
            listscalpara[k] = 'self'
        else:
            listlablpara[k] = [''.join(listnamepara[k].split('_')), '']
            listscalpara[k] = 'self'
            print('Warning! Unrecognized parameter name: %s' % listnamepara[k])
        
        if len(listscalpara[k]) == 0:
            print('')
            print('')
            print('')
            print('Scale (listscalpara[k]) is not defined for the parameter.')
            print('listnamepara[k]')
            print(listnamepara[k])
            raise Exception('')

        if len(listlablpara[k]) == 0:
            raise Exception('')

    if typelang == 'Turkish':
        dictturk = retr_dictturk()
        for k in indxpara:
            if listlablpara[k][0] in dictturk:
                listlablpara[k][0] = dictturk[listlablpara[k][0]]
            if listlablpara[k][1] in dictturk:
                listlablpara[k][1] = dictturk[listlablpara[k][1]]

    for k in indxpara:
        listlablrootpara[k] = listlablpara[k][0]
        listlablunitpara[k] = listlablpara[k][1]
        if listlablunitpara[k] == '':
            listlabltotlpara[k] = listlablpara[k][0]
        else:
            listlabltotlpara[k] = '%s [%s]' % (listlablpara[k][0], listlablpara[k][1])
            
    return listlablpara, listscalpara, listlablrootpara, listlablunitpara, listlabltotlpara


def summgene(varb, boolslin=False, namevarb=None, varbcomp=None, boolshowuniq=False):
    '''
    Output to the command line the content of a varible including its type and summary statistics
    '''
    def show_uniq(varb):
        valuuniq, cntsuniq = np.unique(np.array(varb), return_counts=True)
        print('Unique values:')
        print(valuuniq)
        print('Counts:')
        print(cntsuniq)
    
    if isinstance(varb, dict):
        print('Type is dict with keys:')
        print(list(varb.keys()))
    elif isinstance(varb, list):
        print('Type is list with length %d.' % len(varb))
        #if isinstance(varb[0], str):
        #    print('First element is a string.')
        if boolshowuniq:
            show_uniq(varb)
    
    elif isinstance(varb, tuple):
        print('Type is tuple with length %d.' % len(varb))
    else:
        if boolshowuniq:
            show_uniq(varb)
        if boolslin:
            if namevarb != None:
                namevarb = '%s: ' % namevarb
            else:
                namevarb = ''
            if len(varb) == 0:
                print('Empty variable with size 0.')
                return
            medi = np.nanpercentile(varb, 50)
            minm = np.nanmin(varb)
            maxm = np.nanmax(varb)
            lowr = np.nanpercentile(varb, 16)
            uppr = np.nanpercentile(varb, 84)
            
            if not np.isfinite(varb).all():
                indx = np.where(~np.isfinite(varb))[0]
                strginfi = ', %d infinite samples!' % indx.size
            else:
                strginfi = ''
            
            if np.isnan(varb).any():
                indx = np.where(np.isnan(varb))[0]
                strgnann = ', %d NaN samples!' % indx.size
            else:
                strgnann = ''
            
            if varbcomp is None:
                strgcomp = ''
            else:
                medicomp = np.nanpercentile(varbcomp, 50)
                lowrcomp = np.nanpercentile(varbcomp, 16)
                upprcomp = np.nanpercentile(varbcomp, 84)
                stdvcomp = (upprcomp - lowrcomp) / 2.
                sigm = (medi - medicomp) / np.sqrt(stdvcomp**2 + ((uppr - lowr) / 2.)**2)
                strgcomp = ', %.3g sigma diff' % sigm
            print('%s%.3g +%.3g -%.3g, limt: %.3g %.g, %d samples%s%s%s' % (namevarb, medi, uppr - medi, medi - lowr, \
                                                                                minm, maxm, varb.size, strginfi, strgnann, strgcomp))
        else:
            
            #print('varb')
            #print(varb)
            #print('type(varb)')
            #print(type(varb))
            #print('varb.dtype')
            #print(varb.dtype)
            
            if isinstance(varb, np.ndarray):
                if varb.size == 0:
                    print('An empty number array with size 0.')
                    print('')
                elif varb.dtype.type is np.string_ or varb.dtype.type is np.str_:
                    print('string type numpy array')
                elif varb.dtype.type is np.object_:
                    print('Object type numpy array')
                else:
                    print('varb.dtype')
                    print(varb.dtype)
                    if varb.dtype.names is None:
                        if not np.isfinite(varb).all():
                            boolfini = np.isfinite(varb)
                            indxfini = np.where(boolfini)
                            indxinfi = np.where(~boolfini)[0]
                            print('%d elements in total.' % varb.size)
                            print('%d elements are not finite!' % indxinfi.size)
                            if indxfini[0].size > 0:
                                minmtemp = np.amin(varb[indxfini])
                                maxmtemp = np.nanmax(varb[indxfini])
                                meantemp = np.nanmean(varb[indxfini])
                                print('Among %d finite samples, 0p: %g, 0.3p: %g, 5p: %g, 50p: %g, 95p: %g, 99.7p: %g, 100p: %g, mean: %g' % \
                                                                                                                (indxfini[0].size, minmtemp, \
                                                                                                                np.percentile(varb[indxfini], 0.3), \
                                                                                                                np.percentile(varb[indxfini], 5.), \
                                                                                                                np.percentile(varb[indxfini], 50.), \
                                                                                                                np.percentile(varb[indxfini], 95.), \
                                                                                                                np.percentile(varb[indxfini], 99.7), \
                                                                                                                maxmtemp, meantemp, \
                                                                                                                ))
                        if np.isfinite(varb).any():
                            print(np.nanmin(varb))
                            print(np.nanmax(varb))
                            print(np.nanmean(varb))
                            print(varb.shape)
                            print('%d unique elements.' % np.unique(varb).size)
                    else:
                        print('This is a numpy record!')
                        for name in varb.dtype.names:
                            print(name)
                            summgene(varb[name])

                    print('')
            elif isinstance(varb, np.floating):
                print('Numpy float with type %s' % type(varb))
                print(varb)
            else:
                print('Type is %s' % type(varb))


def retr_p4dm_spec(anch, part='el'):
    
    pathimag, pathdata = retr_path('tdpy')
    if part == 'el':
        strg = 'AtProduction_positrons'
    if part == 'ph':
        strg = 'AtProduction_gammas'
    name = pathdata + 'p4dm/' + strg + '.dat'
    p4dm = loadtxt(name)
    
    p4dm[:, 0] *= 1e3 # [MeV]
    
    mass = unique(p4dm[:, 0])
    nmass = mass.size
    numbener = p4dm.shape[0] / nmass
    
    mult = np.zeros((numbener, nmass))
    for k in range(nmass):
        jp4dm = np.where(abs(p4dm[:, 0] - mass[k]) == 0)[0]

        if anch == 'e':
            mult[:, k] = p4dm[jp4dm, 4]
        if anch == 'mu':
            mult[:, k] = p4dm[jp4dm, 7]
        if anch == 'tau':
            mult[:, k] = p4dm[jp4dm, 10]
        if anch == 'b':
            mult[:, k] = p4dm[jp4dm, 13]
        
    enerscal = 10**p4dm[jp4dm, 1]

    return mult, enerscal, mass


def show_prog(cntr, maxmcntr, thiscntr, nprog=20, indxprocwork=None, showmemo=False, accp=None, accpprio=None):

    nextcntr = int(nprog * float(cntr + 1) / maxmcntr) * 100 / nprog
    if nextcntr > thiscntr:
        if indxprocwork != None:
            print('Process %d is %3d%% completed.' % (indxprocwork, nextcntr))
        else:
            print('%3d%% completed.' % nextcntr)
        if accp != None:
            print('Acceptance ratio: %.3g%%' % accp)
            print('Acceptance through prior boundaries: %.3g%%' % accpprio)
        thiscntr = nextcntr
        if showmemo:
            show_memo_simp()
        
    return thiscntr            


def retr_galcfromequc(rasc, decl):

    icrs = astropy.coordinates.SkyCoord(ra=rasc*astropy.units.degree, dec=decl*astropy.units.degree)

    lgal = icrs.galactic.l.degree
    bgal = icrs.galactic.b.degree
    
    return lgal, bgal


def regr(xdat, ydat, ordr):
    
    coef = polyfit(xdat, ydat, ordr)
    func = poly1d(coef)
    strg = '$y = '
    if ordr == 0:
        strg += '%.5g$'
    if ordr == 1:
        strg += '%.5g x + %.5g$' % (coef[0], coef[1])
    if ordr == 2:
        strg += '%.5g x^2 + %.5g x + %.5g$' % (coef[0], coef[1], coef[2])

    return coef, func, strg


def corr_catl(lgalseco, bgalseco, lgalfrst, bgalfrst, anglassc=np.pi/180., typeverb=1):

    numbfrst = lgalfrst.size

    indxsecoassc = np.zeros(numbfrst, dtype=int) - 1
    numbassc = np.zeros(numbfrst, dtype=int)
    distassc = np.zeros(numbfrst) + 1000.
    lgalbgalfrst = np.array([lgalfrst, bgalfrst])
    thisfraccomp = -1
    numbseco = lgalseco.size
    for k in range(numbseco):
        lgalbgalseco = np.array([lgalseco[k], bgalseco[k]])
        dist = angdist(lgalbgalfrst, lgalbgalseco, lonlat=True)
        thisindxfrst = np.where(dist < anglassc)[0]
        
        if thisindxfrst.size > 0:
            
            # if there are multiple associated true PS, sort them
            indx = np.argsort(dist[thisindxfrst])
            dist = dist[thisindxfrst][indx]
            thisindxfrst = thisindxfrst[indx]
                
            # store the index of the model PS
            numbassc[thisindxfrst[0]] += 1
            if dist[0] < distassc[thisindxfrst[0]]:
                distassc[thisindxfrst[0]] = dist[0]
                indxsecoassc[thisindxfrst[0]] = k

        nextfraccomp = int(100 * float(k) / numbseco)
        if typeverb > 1 and nextfraccomp > thisfraccomp:
            thisfraccomp = nextfraccomp
            print('%02d%% completed.' % thisfraccomp)

    return indxsecoassc


def show_memo_simp():
    
    memoresi, memoresiperc = retr_memoresi()

    strgmemo = retr_strgmemo(memoresi)

    print('Resident memory: %s, %4.3g%%' % (strgmemo, memoresiperc))


def retr_strgmemo(memo):

    if memo >= float(2**30):
        memonorm = memo / float(2**30)
        strg = 'GB'
    elif memo >= float(2**20):
        memonorm = memo / float(2**20)
        strg = 'MB'
    elif memo >= float(2**10):
        memonorm = memo / float(2**10)
        strg = 'KB'
    else:
        memonorm = memo
        strg = 'B'
    strgmemo = '%d %s' % (memonorm, strg)
    return strgmemo


def retr_axis(minm=None, maxm=None, numb=None, bins=None, scal='self'):
    
    if bins == None:
        if scal == 'self':
            bins = np.linspace(minm, maxm, numb + 1)
            mean = (bins[1:] + bins[:-1]) / 2.
        else:
            bins = np.logspace(np.log10(minm), np.log10(maxm), numb + 1)
            mean = np.sqrt(bins[1:] * bins[:-1])
    else:
        if scal == 'self':
            mean = (bins[1:] + bins[:-1]) / 2.
        else:
            mean = np.sqrt(bins[1:] * bins[:-1])
        numb = mean.size
        
    indx = np.arange(numb)
   
    return bins, mean, diff(bins), numb, indx


def retr_psfngausnorm(angl):

    norm = np.sqrt(2. / np.pi**3) / angl / exp(-0.5 * angl**2) / \
                    real(-erfi((angl**2 - np.pi * 1j) / np.sqrt(2) / angl) - erfi((angl**2 + np.pi * 1j) / np.sqrt(2) / angl) + 2. * erfi(angl / np.sqrt(2.)))

    return norm


def retr_mapspnts(lgal, bgal, stdv, flux, numbside=256, typeverb=0):
    
    # lgal, bgal and stdv are in degrees
    numbpnts = lgal.size
    lgalheal, bgalheal, numbpixl, apix = retr_healgrid(numbside)
    gridheal = np.array([lgalheal, bgalheal])
    stdvradi = np.deg2rad(stdv)
    mapspnts = np.zeros(numbpixl)
    for n in range(numbpnts):
        gridpnts = np.array([lgal[n], bgal[n]])
        angl = angdist(gridheal, gridpnts, lonlat=True)
        norm = retr_psfngausnorm(stdvradi)
        mapspnts += apix * norm * flux[n] * exp(-0.5 * angl**2 / stdvradi**2)

    return mapspnts


def retr_mapsplnkfreq(indxpixloutprofi=None, numbsideoutp=256, indxfreqrofi=None):

    numbside = 2048
    numbpixl = 12 * numbside**2
    meanfreq = np.array([30, 44, 70, 100, 143, 217, 353, 545, 857])
    numbfreq = meanfreq.size
    indxfreq = np.arange(numbfreq)
    strgfreq = ['%04d' % meanfreq[k] for k in indxfreq]
    
    indxpixloutp = np.arange(numbsideoutp)

    if indxfreqrofi == None:
        indxfreqrofi = indxfreq

    if indxpixloutprofi == None:
        indxpixloutprofi = indxpixloutp

    rtag = '_%04d' % (numbsideoutp)

    path = retr_path('tdpy', onlydata=True)
    pathsbrt = path + 'plnksbrt_%s.fits' % rtag
    pathsbrtstdv = path + 'plnksbrtstdv_%s.fits' % rtag
    if os.path.isfile(pathsbrt) and os.path.isfile(pathsbrtstdv):
        print('Reading %s...' % pathsbrt)
        mapsplnkfreq = pf.getdata(pathsbrt)
        print('Reading %s...' % pathsbrtstdv)
        mapsplnkfreqstdv= pf.getdata(pathsbrtstdv)
    else:
        mapsplnkfreq = np.zeros((numbfreq, numbpixl))
        mapsplnkfreqstdv = np.zeros((numbfreq, numbpixl))
        for k in indxfreq:
            print('Processing Planck Map at %d GHz...' % (meanfreq[k]))
            # read sky maps
            if strgfreq[k][1] == '0':
                strgfrst = 'plnk/LFI_SkyMap_' 
                strgseco = '-BPassCorrected-field-IQU_0256_R2.01_full.fits'
                strgcols = 'TEMPERATURE'
            elif strgfreq[k][1] == '1' or strgfreq[k][1] == '2' or strgfreq[k][1] == '3':
                strgfrst = 'plnk/HFI_SkyMap_'
                strgseco = '-field-IQU_2048_R2.02_full.fits'
                strgcols = 'I_STOKES'
            else:
                strgfrst = 'plnk/HFI_SkyMap_'
                strgseco = '-field-Int_2048_R2.02_full.fits'
                strgcols = 'I_STOKES'
            strg = strgfrst + '%s' % strgfreq[k][1:] + strgseco
        
            mapsinpt = pf.getdata(path + strg, 1)[strgcols]
            numbpixlinpt = mapsinpt.size
            numbsideinpt = int(np.sqrt(numbpixlinpt / 12))
            mapsplnkfreq[k, :] = pf.getdata(path + strg, 1)[strgcols]
            mapsplnkfreq[k, :] = hp.reorder(mapsplnkfreq[k, :], n2r=True)
        
            # change units of the sky maps to Jy/sr
            if strgfreq[k] != '0545' and strgfreq[k] != '0857':
                ## from Kcmb
                if calcfactconv:
                    # read Planck band transmission data
                    if strgfreq[k][1] == '0':
                        strg = 'LFI_RIMO_R2.50'
                        strgextn = 'BANDPASS_%s' % strgfreq[k][1:]
                        freqband = pf.open(path + '%s.fits' % strg)[strgextn].data['WAVENUMBER'][1:] * 1e9
                    else:
                        strg = 'plnk/HFI_RIMO_R2.00'
                        strgextn = 'BANDPASS_F%s' % strgfreq[k][1:]
                        freqband = 1e2 * velolght * pf.open(path + '%s.fits' % strg)[strgextn].data['WAVENUMBER'][1:]
                    tranband = pf.open(path + '%s.fits' % strg)[strgextn].data['TRANSMISSION'][1:]
                    indxfreqbandgood = np.where(tranband > 1e-6)[0]
                    indxfreqbandgood = np.arange(amin(indxfreqbandgood), amax(indxfreqbandgood) + 1)
        
                    # calculate the unit conversion factor
                    freqscal = consplnk * freqband[indxfreqbandgood] / consbolt / tempcmbr
                    freqcntr = float(strgfreq[k]) * 1e9
                    specdipo = 1e26 * 2. * (consplnk * freqband[indxfreqbandgood]**2 / velolght / tempcmbr)**2 / consbolt / (exp(freqscal) - 1.)**2 * exp(freqscal) # [Jy/sr]
                    factconv = trapz(specdipo * tranband[indxfreqbandgood], freqband[indxfreqbandgood]) / \
                                                    trapz(freqcntr * tranband[indxfreqbandgood] / freqband[indxfreqbandgood], freqband[indxfreqbandgood]) # [Jy/sr/Kcmb]
                else:
                    # read the unit conversion factors provided by Planck
                    factconv = factconvplnk[k, 1] * 1e6
            else:
                ## from MJy/sr
                factconv = 1e6
            mapsplnk[k, :] *= factconv
        
        pf.writeto(pathsbrt, mapsplnkfreq, clobber=True)
        pf.writeto(pathsbrtstdv, mapsplnkfreqstdv, clobber=True)
        
    return mapsplnkfreq, mapsplnkfreqstdv


def rebn_arry(arry, shap):
    '''
    Rebin a numpy array into a new shape by averaging
    '''

    shaptemp = shap[0], arry.shape[0] // shap[0], shap[1], arry.shape[1] // shap[1]
    arryrebn = arry.reshape(shaptemp).mean(-1).mean(1)

    return arryrebn


def retr_indximagmaxm(data):

    sizeneig = 10
    cntpthrs = 10
    maxmdata = sp.ndimage.filters.maximum_filter(data, sizeneig)
    
    boolmaxm = (data == maxmdata)
    minmdata = sp.ndimage.filters.minimum_filter(data, sizeneig)
    diff = ((maxmdata - minmdata) > cntpthrs)
    boolmaxm[diff == 0] = 0
    mapslabl, numbobjt = sp.ndimage.label(boolmaxm)
    mapslablones = np.zeros_like(mapslabl)
    mapslablones[np.where(mapslabl > 0)] = 1.
    indxmaxm = np.array(sp.ndimage.center_of_mass(data, mapslabl, range(1, numbobjt+1))).astype(int)
    if len(indxmaxm) == 0:
        indxydatmaxm = np.array([0])
        indxxdatmaxm = np.array([0])
    else:
        indxydatmaxm = indxmaxm[:, 1]
        indxxdatmaxm = indxmaxm[:, 0]
    return indxxdatmaxm, indxydatmaxm


def plot_timeline(
                  path, \
                  dictrows, \
                  # separation in months between successive ticks on the horizontal axis
                  delttime=None, \
                  # number of tick marks along the horizontal axis
                  numbtime=None, \
                  # tuple indicating the size of the figure
                  sizefigr=None, \
                  # type of the times to be highlighted
                  typetimeshow=None, \
                  # a list of times and labels to highligh with vertical lines
                  listjdatlablhigh=None, \
                  ):
    '''
    Make a timeline (Gantt) chart
    '''
    
    if sizefigr is None:
        sizefigr = (8., 4.)

    listnamerows = list(dictrows.keys())
    numbrows = len(dictrows)
    indxrows = np.arange(numbrows)
    minmjdat = 1e100
    maxmjdat = -1e100
    listlablrows = [[] for k in indxrows]
    listsizerows = [[] for k in indxrows]
    listcolrrows = [[] for k in indxrows]
    for k in indxrows:
        
        if 'size' in dictrows[listnamerows[k]]:
            listsizerows[k] = dictrows[listnamerows[k]]['size']
        else:
            listsizerows[k] = 1.
        
        if 'colr' in dictrows[listnamerows[k]]:
            listcolrrows[k] = dictrows[listnamerows[k]]['colr']
        else:
            listcolrrows[k] = 'gray'
        
        if 'labl' in dictrows[listnamerows[k]]:
            listlablrows[k] = dictrows[listnamerows[k]]['labl']
        else:
            listlablrows[k] = listnamerows[k]
            
        for l in range(len(dictrows[listnamerows[k]]['listelem'])):
            minmjdat = min(minmjdat, astropy.time.Time(dictrows[listnamerows[k]]['listelem'][l][0], format='iso').jd)
            maxmjdat = max(maxmjdat, astropy.time.Time(dictrows[listnamerows[k]]['listelem'][l][1], format='iso').jd)
    
    #dicttemp = dict()
    #dicttemp['init'] = dict()
    #dicttemp['finl'] = dict()
    #for limt in ['minm', 'maxm']:
    #    if limt == 'minm':
    #        limt = minmjdat
    #    if limt == 'maxm':
    #        limt = maxmjdat
    #    for strg in ['year', 'month', 'day']:
    #        dicttemp[limt][strg] = limtjdatdict
    #    dicttemp['month'] = 1
    
    #objttimeminm = astropy.time.Time(minmjdat, format='jd')
    #objttimemaxm = astropy.time.Time(maxmjdat, format='jd')
    #
    #objtdatetimedelt = datetime.datetime(0, 1, 0, 0, 0, 0)
    #
    #dicttemp = dict()
    #dicttemp['month'] = 1
    #objttimedelt = astropy.time.TimeDelta(objtdatetimedelt, format='ymdhms')
    #listobjttime = []
    #objttimetemp = objttimeminm
    #while True:
    #    objttimetemp = objttimetemp + objttimedelt
    #    if objttimetemp > objttimemaxm:
    #        break
    #    listobjttime.append(objttimetemp)
    
    #dicttemp['init']['year'] = 
    limtjdat = [minmjdat, maxmjdat]
    if delttime is not None:
        listjdat = np.arange(minmjdat, maxmjdat, delttime / 12. * 365.25)
    elif numbtime is not None:
        listjdat = np.linspace(minmjdat, maxmjdat, numbtime)
    else:
        print('delttime')
        print(delttime)
        print('numbtime')
        print(numbtime)
        raise Exception('')
    listlabltick = [[] for jdat in listjdat]
    
    listjdatvert = []
    if typetimeshow == 'summer':
        jdatinitintg = int(minmjdat)
        jdatinit = minmjdat
        jdatfinl = minmjdat
        cntr = 0
        
        while True:
            strgtimethis = astropy.time.Time(jdatinitintg + cntr * 365.25, format='jd').to_value('iso', subfmt='date')
            print('strgtimethis')
            print(strgtimethis)
            
            #strgtimeinit = str(strgtimethis)
            #strgtimefinl = str(strgtimethis)
            strgtimeinit = strgtimethis[:4] + '-06-01'
            strgtimefinl = strgtimethis[:4] + '-08-31'
            jdatinit = astropy.time.Time(strgtimeinit, format='iso').jd
            jdatfinl = astropy.time.Time(strgtimefinl, format='iso').jd
            if jdatinit < maxmjdat:
                listjdatvert.append(jdatinit)
            if jdatfinl < maxmjdat:
                listjdatvert.append(jdatfinl)
            if jdatfinl > maxmjdat:
                break
            cntr += 1

    #listlabltick[kk] = astropy.time.Time(jdat, format='TimeYMDHMS').to_value('iso', subfmt='date')
    
    for kk, jdat in enumerate(listjdat):
        print('temp: Trucating the date labels...')
        listlabltick[kk] = astropy.time.Time(jdat, format='jd').to_value('iso', subfmt='date')[:-3]
    
    figr, axis = plt.subplots(1, figsize=sizefigr)
    
    ydattext = 0.
    listtick = [[] for k in indxrows]
    for k in indxrows:
        
        for l in range(len(dictrows[listnamerows[k]]['listelem'])):
            
            jdatfrst = astropy.time.Time(dictrows[listnamerows[k]]['listelem'][l][0], format='iso').jd
            jdatseco = astropy.time.Time(dictrows[listnamerows[k]]['listelem'][l][1], format='iso').jd
            
            xdattext = jdatfrst + 0.5 * (jdatseco - jdatfrst)
            
            axis.barh(ydattext, jdatseco - jdatfrst, left=jdatfrst, color=listcolrrows[k], height=listsizerows[k], edgecolor='k', alpha=0.5)
            
            listlablxaxi = axis.set_xticks(listjdat, rotation=45)
            
            axis.set_xticklabels(listlabltick)
            
            axis.text(xdattext, ydattext, dictrows[listnamerows[k]]['listelem'][l][2], color='k', ha='center', va='center')
        
        if k == 0:
            minmydat = ydattext - 0.5 * listsizerows[k]
        if k == numbrows - 1:
            maxmydat = ydattext + 0.5 * listsizerows[k]
        
        listtick[k] = ydattext

        if k < numbrows - 1:
            ydattext += 0.5 * (listsizerows[k] + listsizerows[k+1])
    
    for jdatvert in listjdatvert:
        axis.axvline(jdatvert, ls='--', color='gray', alpha=0.4)
    
    if listjdatlablhigh is not None:
        for jdatlabl in listjdatlablhigh:
            if jdatlabl == 'now':
                jdat = astropy.time.Time.now().jd
            else:
                jdat = astropy.time.Time(jdatlabl[0], format='iso').jd
            axis.axvline(jdat, ls='--')
            axis.text(jdat, -1, jdatlabl[1], ha='center', va='center')
    
    limtydat = [minmydat, maxmydat]
    axis.set_yticks(listtick)
    axis.set_yticklabels(listlablrows)
    axis.set_xlim(limtjdat)
    axis.set_ylim(limtydat)
    axis.set_xlabel('Time')
    print('Writing to %s...' % path)
    plt.tight_layout()
    plt.savefig(path, dpi=300)
    plt.close()


def plot_gene(path, xdat, ydat, yerr=None, scalxdat=None, scalydat=None, \
                                            lablxdat='', lablydat='', plottype=None, limtxdat=None, limtydat=None, colr=None, listlinestyl=None, \
                                            alph=None, listlegd=None, listvlinfrst=None, listvlinseco=None, listhlin=None, drawdiag=False):
    
    if not isinstance(ydat, list):
        listydat = [ydat]
    else:
        listydat = ydat
  
    if yerr != None and not isinstance(ydat, list):
        listyerr = [yerr]
    else:
        listyerr = yerr

    numbelem = len(listydat)
    
    if listlinestyl == None:
        listlinestyl = [None for k in range(numbelem)]

    if listlegd == None:
        listlegd = [None for k in range(numbelem)]

    if not isinstance(xdat, list):
        listxdat = [xdat for k in range(numbelem)]
    else:
        listxdat = xdat
    
    if plottype == None:
        listplottype = ['line' for k in range(numbelem)]
    elif not isinstance(plottype, list):
        listplottype = [plottype]
    else:
        listplottype = plottype
    
    figr, axis = plt.subplots(figsize=(6, 6))
    
    for k in range(numbelem):
        xdat = listxdat[k]
        ydat = listydat[k]
        plottype = listplottype[k]
        legd = listlegd[k]
        if plottype == 'scat':
            axis.scatter(xdat, ydat, color=colr, alpha=alph, label=legd, s=2)
        elif plottype == 'hist':
            deltxdat = xdat[1] - xdat[0]
            axis.bar(xdat - deltxdat / 2., ydat, deltxdat, color=colr, alpha=alph, label=legd)
        else:
            if listyerr != None:
                axis.errorbar(xdat, ydat, yerr=listyerr[k], color=colr, lw=2, alpha=alph, label=legd, ls=listlinestyl[k])
            else:
                axis.plot(xdat, ydat, color=colr, lw=2, alpha=alph, label=legd, ls=listlinestyl[k])
    
    if listlegd != None:
        axis.legend(framealpha=1.)

    if scalxdat == 'logt':
        axis.set_xscale('log')
    if scalydat == 'logt':
        axis.set_yscale('log')
    
    if limtxdat == None:
        limtxdat = [np.amin(np.concatenate(listxdat)), np.amax(np.concatenate(listxdat))]
    if limtydat == None:
        limtydat = [np.amin(np.concatenate(listydat)), np.amax(np.concatenate(listydat))]
    
    if drawdiag:
        axis.plot(limtxdat, limtxdat, ls='--', alpha=0.3, color='black')
    
    axis.set_xlim(limtxdat)
    axis.set_ylim(limtydat)
    
    if listhlin != None:
        if isscalar(listhlin):
            listhlin = [listhlin]
        for k in range(len(listhlin)):
            axis.axhline(listhlin[k], ls='--', alpha=0.2, color=colr)
    
    if listvlinfrst != None:
        if isscalar(listvlinfrst):
            listvlinfrst = [listvlinfrst]
        for k in range(len(listvlinfrst)):
            axis.axvline(listvlinfrst[k], ls='--', alpha=0.2, color=colr)
    
    if listvlinseco != None:
        if isscalar(listvlinseco):
            listvlinseco = [listvlinseco]
        for k in range(len(listvlinseco)):
            axis.axvline(listvlinseco[k], ls='-.', alpha=0.2, color=colr)
    
    axis.set_xlabel(lablxdat)
    axis.set_ylabel(lablydat)

    figr.tight_layout()
    figr.savefig(path)
    plt.close(figr)


def cart_heal(cart, minmlgal=-180., maxmlgal=180., minmbgal=-90., maxmbgal=90., nest=False, numbside=256):
    
    numbbgcr = cart.shape[0]
    numblgcr = cart.shape[1]
    lghp, bghp, numbpixl, apix = retr_healgrid(numbside)
    
    indxlgcr = (numblgcr * (lghp - minmlgal) / (maxmlgal - minmlgal)).astype(int)
    indxbgcr = (numbbgcr * (bghp - minmbgal) / (maxmbgal - minmbgal)).astype(int)
    
    indxpixlrofi = np.where((minmlgal <= lghp) & (lghp <= maxmlgal) & (minmbgal <= bghp) & (bghp <= maxmbgal))[0]
    
    heal = np.zeros(numbpixl)
    heal[indxpixlrofi] = np.fliplr(cart)[indxbgcr[indxpixlrofi], indxlgcr[indxpixlrofi]]
    
    return heal


class cntr():
    def gets(self):
        return self.cntr

    def incr(self, valu=1):
        temp = self.cntr
        self.cntr += valu
        return temp
    def __init__(self):
        self.cntr = 0


def retr_evttferm(recotype):
    
    if recotype == 'rec7':
        evtt = np.array([1, 2])
    if recotype == 'rec8':
        evtt = np.array([16, 32])
    if recotype == 'manu':
        evtt = np.array([0, 1])
    
    numbevtt = evtt.size
    indxevtt = np.arange(numbevtt)

    return evtt, numbevtt, indxevtt


def writ_sbrtfdfm(numbside=256, regitype='igal', binsenertype='full', recotype='rec7'):
    
    pathdata = os.environ["PCAT_DATA_PATH"] + '/data/'
    
    evtt, numbevtt, indxevtt = retr_evttferm(recotype)
    
    binsener = np.array([0.1, 0.3, 1., 3., 10., 100.])
    
    meanener = np.sqrt(binsener[1:] * binsener[:-1])
    numbpixl = 12 * numbside**2
    numbener = binsener.size - 1

    # get the Fermi-LAT diffuse model
    sbrtfdfm = retr_sbrtfdfm(binsener, numbside)
    
    # rotate if necessary
    for m in indxevtt:
        if regitype == 'ngal':
            for i in range(numbener):
                almc = hp.map2alm(sbrtfdfm[i, :])
                hp.rotate_alm(almc, 0., 0.5 * np.pi, 0.)
                sbrtfdfm[i, :] = hp.alm2map(almc, numbside)

    # smooth the model
    sbrtfdfm = smth_ferm(sbrtfdfm[:, :, None], meanener, recotype)

    path = pathdata + 'sbrtfdfm%04d%s%s%s.fits' % (numbside, recotype, binsenertype, regitype)
    pf.writeto(path, sbrtfdfm, clobber=True)


def retr_strgtimestmp():

    strgtimestmp = datetime.datetime.now().strftime('%Y%m%d_%H%M%S')

    return strgtimestmp


def read_fits(path, pathimag=None, full=False, typeverb=0):
    '''
    Read FITS file
    '''
    
    if pathimag != None:
        os.system('mkdir -p ' + pathimag)
    
    print('Reading from %s...' % path)
    hdun = astropy.io.fits.open(path)
    numbhead = len(hdun)
    listdata = []
    for k in range(numbhead):
        head = hdun[k].header
        data = hdun[k].data

        print('Extension %d...' % k)
        listdata.append(data)

        print('Header:')
        print(repr(head))
        print('')
        print('')
        print('')
        print('Data:')
        if data is not None:
            if isinstance(data, np.ndarray):
                print('It is a numpy array')
                summgene(data)
            else:
                print('data.names')
                print(data.names)
                for name in data.names:
                    print(name)
                    summgene(data[name])
        print('')
        print('')
        print('')
        print('')
        print('')
        print('')
        continue

        if data == None:
            continue

        numb = len(list(head.keys()))
        listkeys = list(head.keys())
        listvalu = list(head.values())
        #print('listkeys')
        #print(listkeys)
        #print('listvalu')
        #print(listvalu)
        
        for m, keys in enumerate(listkeys):
            print(keys)
            print(listvalu[m])
            print('')
            
        listtype = []
        listform = []
        listunit = []
       
        for n in range(numb):
            if listkeys[n].startswith('TTYPE'):
                listtype.append(listvalu[n])
            if listkeys[n].startswith('TFORM'):
                listform.append(listvalu[n])
            if listkeys[n].startswith('TUNIT'):
                listunit.append(listvalu[n])

        if pathimag != None:
            cmnd = 'convert -density 300'
    
            for n in range(len(listtype)):
                if not listform[n].endswith('A') and isfinite(data[listtype[n]]).all():
                    figr, axis = plt.subplots()
                    bins = np.linspace(amin(data[listtype[n]]), amax(data[listtype[n]]), 100)
                    axis.hist(data[listtype[n]], bins=bins)
                    #axis.set_xlabel('%s [%s]' % (listtype[n], listunit[n]))
                    axis.set_yscale('log')
                    
                    axis.set_xlabel('%s' % (listtype[n]))
                    plt.tight_layout()
                    path = pathimag + 'readfits_%s.pdf' % listtype[n]
                    cmnd += ' ' + path
                    figr.savefig(path)
                    plt.close(figr)
        
            cmnd += ' ' + pathimag + 'merg.pdf'
            os.system(cmnd)

    return listdata


def plot_maps(path, maps, pixltype='heal', scat=None, indxpixlrofi=None, numbpixl=None, titl='', minmlgal=None, maxmlgal=None, minmbgal=None, maxmbgal=None, \
                                                                                        resi=False, satu=False, numbsidelgal=None, numbsidebgal=None, igal=False):
   
    if minmlgal == None:
        if not igal:
            minmlgal = -180.
            minmbgal = -90.
            maxmlgal = 180.
            maxmbgal = 90.
        else:
            minmlgal = -20.
            minmbgal = -20.
            maxmlgal = 20.
            maxmbgal = 20.
            
    asperati = (maxmbgal - minmbgal) / (maxmlgal - minmlgal)
    
    if indxpixlrofi != None:
        mapstemp = np.zeros(numbpixl)
        mapstemp[indxpixlrofi] = maps
        maps = mapstemp
    else:
        numbpixl = maps.size
    
    if numbsidelgal == None:
        numbsidelgal = min(4 * int((maxmlgal - minmlgal) / np.rad2deg(np.sqrt(4. * np.pi / numbpixl))), 2000)
    if numbsidebgal == None:
        numbsidebgal = min(4 * int((maxmbgal - minmbgal) / np.rad2deg(np.sqrt(4. * np.pi / numbpixl))), 2000)
    
    # saturate the map
    if satu:
        mapstemp = copy(maps)
        maps = mapstemp
        if not resi:
            satu = 0.1 * amax(maps)
        else:
            satu = 0.1 * min(np.fabs(amin(maps)), amax(maps))
            maps[np.where(maps < -satu)] = -satu
        maps[np.where(maps > satu)] = satu

    exttrofi = [minmlgal, maxmlgal, minmbgal, maxmbgal]

    if pixltype == 'heal':
        cart = retr_cart(maps, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal, numbsidelgal=numbsidelgal, numbsidebgal=numbsidebgal)
    else:
        numbsidetemp = int(np.sqrt(maps.size))
        cart = maps.reshape((numbsidetemp, numbsidetemp)).T

    sizefigr = 8
    if resi:
        cmap = 'RdBu'
    else:
        cmap = 'Reds'

    if asperati < 1.:   
        factsrnk = 1.5 * asperati
    else:
        factsrnk = 0.8
    figr, axis = plt.subplots(figsize=(sizefigr, asperati * sizefigr))
    if resi:
        valu = max(np.fabs(amin(cart)), np.fabs(amax(cart)))
        imag = plt.imshow(cart, origin='lower', cmap=cmap, extent=exttrofi, interpolation='none', vmin=-valu, vmax=valu)
    else:
        imag = plt.imshow(cart, origin='lower', cmap=cmap, extent=exttrofi, interpolation='none')
    
    if scat != None:
        numbscat = len(scat)
        for k in range(numbscat):
            axis.scatter(scat[k][0], scat[k][1])

    cbar = plt.colorbar(imag, shrink=factsrnk) 
    
    plt.title(titl, y=1.08)

    figr.savefig(path)
    plt.close(figr)
    

def rebn(arry, shapoutp, totl=False):
    
    shaptemp = shapoutp[0], arry.shape[0] // shapoutp[0], shapoutp[1], arry.shape[1] // shapoutp[1]
    
    if totl:
        arryoutp = arry.reshape(shaptemp).sum(-1).sum(1)
    else:
        arryoutp = arry.reshape(shaptemp).mean(-1).mean(1)

    return arryoutp


def cart_heal_depr(cart, minmlgal=-180., maxmlgal=180., minmbgal=-90., maxmbgal=90., nest=False, numbside=256):
    
    nbgcr = cart.shape[0]
    nlgcr = cart.shape[1]
    lghp, bghp, numbpixl, apix = retr_healgrid(numbside)
    heal = np.zeros(numbpixl)
    jpixl = np.where((minmlgal < lghp) & (lghp < maxmlgal) & (minmbgal < bghp) & (bghp < maxmbgal))[0]
    jlgcr = (nlgcr * (lghp[jpixl] - minmlgal) / (maxmlgal - minmlgal)).astype(int)
    jbgcr = (nbgcr * (bghp[jpixl] - minmbgal) / (maxmbgal - minmbgal)).astype(int)
    
    heal[jpixl] = np.fliplr(cart)[jbgcr, jlgcr]
    
    return heal


def retr_healgrid(numbside):
    
    numbpixl = 12 * numbside**2
    apix = 4. * np.pi / numbpixl # [sr]
    thhp, phhp = hp.pixelfunc.pix2ang(numbside, np.arange(numbpixl), nest=False) # [rad]
    lghp = np.rad2deg(phhp)
    lghp = 180. - ((lghp - 180.) % 360.)# - 180. # [deg]
    bghp = 90. - np.rad2deg(thhp)

    return lghp, bghp, numbpixl, apix


def retr_isot(binsener, numbside=256):
    
    diffener = binsener[1:] - binsener[:-1]
    numbpixl = 12 * numbside**2
    numbener = binsener.size - 1
    numbsamp = 10

    # get the best-fit isotropic surface brightness given by the Fermi-LAT collaboration

    pathdata = retr_path('tdpy', onlydata=True)
    path = pathdata + 'iso_P8R2_SOURCE_V6_v06.txt'
    isotdata = loadtxt(path)
    enerisot = isotdata[:, 0] * 1e-3 # [GeV]
    sbrtisottemp = isotdata[:, 1] * 1e3 # [1/cm^2/s/sr/GeV]
    
    # sampling energy grid
    binsenersamp = np.logspace(np.log10(amin(binsener)), np.log10(amax(binsener)), numbsamp * numbener)
    
    # interpolate the surface brightness over the sampling energy grid
    sbrtisottemp = interp(binsenersamp, enerisot, sbrtisottemp)
    
    # take the mean surface brightness in the desired energy bins
    sbrtisot = np.empty((numbener, numbpixl))
    for i in range(numbener):
        sbrtisot[i, :] = trapz(sbrtisottemp[i*numbsamp:(i+1)*numbsamp], binsenersamp[i*numbsamp:(i+1)*numbsamp]) / diffener[i]
        
    return sbrtisot


def retr_cart(hmap, indxpixlrofi=None, numbsideinpt=None, minmlgal=-180., maxmlgal=180., minmbgal=-90., maxmbgal=90., nest=False, \
                                                                                                            numbsidelgal=None, numbsidebgal=None):
   
    shap = hmap.shape
    if indxpixlrofi == None:
        numbpixlinpt = shap[0]
        numbsideinpt = int(np.sqrt(numbpixlinpt / 12.))
    else:
        numbpixlinpt = numbsideinpt**2 * 12
    
    if numbsidelgal == None:
        numbsidelgal = 4 * int((maxmlgal - minmlgal) / np.rad2deg(np.sqrt(4. * np.pi / numbpixlinpt)))
    if numbsidebgal == None:
        numbsidebgal = 4 * int((maxmbgal - minmbgal) / np.rad2deg(np.sqrt(4. * np.pi / numbpixlinpt)))
    
    lgcr = np.linspace(minmlgal, maxmlgal, numbsidelgal)
    indxlgcr = np.arange(numbsidelgal)
    
    bgcr = np.linspace(minmbgal, maxmbgal, numbsidebgal)
    indxbgcr = np.arange(numbsidebgal)
    
    lghp, bghp, numbpixl, apix = retr_healgrid(numbsideinpt)

    bgcrmesh, lgcrmesh = np.meshgrid(bgcr, lgcr, indexing='ij')
    
    indxpixlmesh = hp.ang2pix(numbsideinpt, np.pi / 2. - np.deg2rad(bgcrmesh), np.deg2rad(lgcrmesh))
    
    if indxpixlrofi == None:
        indxpixltemp = indxpixlmesh
    else:
        pixlcnvt = np.zeros(numbpixlinpt, dtype=int) - 1
        for k in range(indxpixlrofi.size):
            pixlcnvt[indxpixlrofi[k]] = k
        indxpixltemp = pixlcnvt[indxpixlmesh]
    
    indxbgcrgrid, indxlgcrgrid = np.meshgrid(indxbgcr, indxlgcr, indexing='ij')

    if hmap.ndim == 2:
        hmapcart = np.empty((numbsidebgal, numbsidelgal, shap[1]))
        for k in range(shap[1]):
            hmapcart[np.meshgrid(indxbgcr, indxlgcr, k, indexing='ij')] = hmap[indxpixltemp, k][:, :, None]
    else:
        hmapcart = np.empty((numbsidebgal, numbsidelgal))
        hmapcart[np.meshgrid(indxbgcr, indxlgcr, indexing='ij')] = hmap[indxpixltemp]

    return np.fliplr(hmapcart).T


def retr_sbrtfdfm(binsener, numbside=256, vfdm=7):                    
    
    diffener = binsener[1:] - binsener[0:-1]
    numbener = diffener.size

    numbpixl = numbside**2 * 12
    
    path = os.environ["TDPY_DATA_PATH"] + '/data/'
    if vfdm == 2:
        path += 'gll_iem_v02.fit'
    if vfdm == 3:
        path += 'gll_iem_v02_P6_V11_DIFFUSE.fit'
    if vfdm == 4:
        path += 'gal_2yearp7v6_v0.fits'
    if vfdm == 5:
        path += 'gll_iem_v05.fit'
    if vfdm == 6:
        path += 'gll_iem_v05_rev1.fit'
    if vfdm == 7:
        path += 'gll_iem_v06.fits'
    
    sbrtcart = pf.getdata(path, 0) * 1e3 # [1/cm^2/s/sr/GeV]
    enerfdfm = np.array(pf.getdata(path, 1).tolist()).flatten() * 1e-3 # [GeV]
    sbrtfdfmheal = np.zeros((enerfdfm.size, numbpixl))
    for i in range(enerfdfm.size):
        sbrtfdfmheal[i, :] = cart_heal(np.fliplr(sbrtcart[i, :, :]), numbside=numbside)
    
    sbrtfdfm = np.empty((numbener, numbpixl))
    numbsampbins = 10
    enersamp = np.logspace(np.log10(amin(binsener)), np.log10(amax(binsener)), numbsampbins * numbener)
    
    sbrtfdfmheal = interpolate.interp1d(enerfdfm, sbrtfdfmheal, axis=0)(enersamp)
    for i in range(numbener):
        sbrtfdfm[i, :] = trapz(sbrtfdfmheal[i*numbsampbins:(i+1)*numbsampbins, :], enersamp[i*numbsampbins:(i+1)*numbsampbins], axis=0) / diffener[i]

    return sbrtfdfm


def plot_matr(axis, xdat, ydat, labl, loc=1):
    
    listlinestyl = [':', '--', '-']
    listcolr = ['b', 'r', 'g']
    
    for i in range(3):
        for  j in range(3):
            if len(xdat.shape) == 3:
                axis.plot(xdat[i, j, :], ydat[i, j, :], color=listcolr[j], ls=listlinestyl[i])
            else:
                axis.plot(xdat, ydat[i, j, :], color=c[j], ls=ls[i])

    line = []
    for k in np.arange(3):
        line.append(plt.Line2D((0,1),(0,0), color='k', ls=listlinestyl[k]))
    for l in range(3):
        line.append(plt.Line2D((0,1),(0,0), color=listcolr[l]))
    axis.legend(line, labl, loc=loc, ncol=2, framealpha=1.)


def plot_braz(axis, xdat, ydat, yerr=None, numbsampdraw=0, lcol='yellow', dcol='green', mcol='black', labl=None, alpha=None):
    
    if numbsampdraw > 0:
        jsampdraw = choice(np.arange(ydat.shape[0]), size=numbsampdraw)
        axis.plot(xdat, ydat[jsampdraw[0], :], alpha=0.1, color='b', label='Samples')
        for k in range(1, numbsampdraw):
            axis.plot(xdat, ydat[jsampdraw[k], :], alpha=0.1, color='b')

    if yerr != None:
        axis.plot(xdat, ydat - yerr[0, :], color=dcol, alpha=alpha)
        line, = axis.plot(xdat, ydat, color=mcol, label=labl, alpha=alpha)
        axis.plot(xdat, ydat + yerr[1, :], color=dcol, alpha=alpha)
        ptch = mpl.patches.Patch(facecolor=dcol, alpha=alpha, linewidth=0)
        #ax.legend([(p1,p2)],['points'],scatterpoints=2, framealpha=1.)
        #plt.legend([(ptch, line)], ["Theory"], handler_map = {line : mpl.legend_handler.HandlerLine2D(marker_pad = 0)} )
        
        axis.fill_between(xdat, ydat - yerr[0, :], ydat + yerr[1, :], color=dcol, alpha=alpha)
        
        hand, labl = axis.get_legend_handles_labels()
        hand[0] = [hand[0], ptch]

        return ptch, line
    else:
        axis.plot(xdat, np.percentile(ydat, 2.5, 0), color=lcol, alpha=alpha)
        axis.plot(xdat, np.percentile(ydat, 16., 0), color=dcol, alpha=alpha)
        axis.plot(xdat, np.percentile(ydat, 50., 0), color=mcol, label=labl, alpha=alpha)
        axis.plot(xdat, np.percentile(ydat, 84., 0), color=dcol, alpha=alpha)
        axis.plot(xdat, np.percentile(ydat, 97.5, 0), color=lcol, alpha=alpha)
        axis.fill_between(xdat, np.percentile(ydat, 2.5, 0), np.percentile(ydat, 97.5, 0), color=lcol, alpha=alpha)#, label='95% C.L.')
        axis.fill_between(xdat, np.percentile(ydat, 16., 0), np.percentile(ydat, 84., 0), color=dcol, alpha=alpha)#, label='68% C.L.')


def retr_psfn(gdat, psfp, indxenertemp, thisangl, psfntype, binsoaxi=None, oaxitype=None, strgmodl='fitt'):

    numbpsfpform = getattr(gdat, strgmodl + 'numbpsfpform')
    numbpsfptotl = getattr(gdat, strgmodl + 'numbpsfptotl')
    
    indxpsfpinit = numbpsfptotl * (indxenertemp[:, None] + gdat.numbener * gdat.indxevtt[None, :])
    if oaxitype:
        indxpsfponor = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp]
        indxpsfpoind = numbpsfpform + numbpsfptotl * gdat.indxener[indxenertemp] + 1
    
    if gdat.exprtype == 'ferm':
        scalangl = 2. * arcsin(np.sqrt(2. - 2. * cos(thisangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
        scalanglnorm = 2. * arcsin(np.sqrt(2. - 2. * cos(gdat.binsangl)) / 2.)[None, :, None] / gdat.fermscalfact[:, None, :]
    else:
        if oaxitype:
            scalangl = thisangl[None, :, None, None]
        else:
            scalangl = thisangl[None, :, None]
    
    if oaxitype:
        factoaxi = retr_factoaxi(gdat, binsoaxi, psfp[indxpsfponor], psfp[indxpsfpoind])
   
    if psfntype == 'singgaus':
        sigc = psfp[indxpsfpinit]
        if oaxitype:
            sigc = sigc[:, None, :, None] * factoaxi[:, None, :, :]
        else:
            sigc = sigc[:, None, :]
        psfn = retr_singgaus(scalangl, sigc)
    
    elif psfntype == 'singking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        
        psfn = retr_singking(scalangl, sigc, gamc)
    
    elif psfntype == 'doubking':
        sigc = psfp[indxpsfpinit]
        gamc = psfp[indxpsfpinit+1]
        sigt = psfp[indxpsfpinit+2]
        gamt = psfp[indxpsfpinit+3]
        frac = psfp[indxpsfpinit+4]
        sigc = sigc[:, None, :]
        gamc = gamc[:, None, :]
        sigt = sigt[:, None, :]
        gamt = gamt[:, None, :]
        frac = frac[:, None, :]
        
        psfn = retr_doubking(scalangl, frac, sigc, gamc, sigt, gamt)
        if gdat.exprtype == 'ferm':
            psfnnorm = retr_doubking(scalanglnorm, frac, sigc, gamc, sigt, gamt)
    
    # normalize the PSF
    if gdat.exprtype == 'ferm':
        fact = 2. * np.pi * trapz(psfnnorm * sin(gdat.binsangl[None, :, None]), gdat.binsangl, axis=1)[:, None, :]
        psfn /= fact

    return psfn


def retr_psfpferm(gdat):
   
    gdat.exproaxitype = False
    
    if gdat.recotype == 'rec8':
        path = gdat.pathdata + 'expr/irfn/psf_P8R2_SOURCE_V6_PSF.fits'
    else:
        path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
    irfn = pf.getdata(path, 1)
    minmener = irfn['energ_lo'].squeeze() * 1e-3 # [GeV]
    maxmener = irfn['energ_hi'].squeeze() * 1e-3 # [GeV]
    enerirfn = np.sqrt(minmener * maxmener)

    numbpsfpscal = 3
    numbpsfpform = 5
    
    fermscal = np.zeros((gdat.numbevtt, numbpsfpscal))
    fermform = np.zeros((gdat.numbener, gdat.numbevtt, numbpsfpform))
    
    parastrg = ['score', 'gcore', 'stail', 'gtail', 'ntail']
    for m in gdat.indxevtt:
        if gdat.recotype == 'rec7':
            if m == 0:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_front.fits'
            elif m == 1:
                path = gdat.pathdata + 'expr/irfn/psf_P7REP_SOURCE_V15_back.fits'
            irfn = pf.getdata(path, 1)
            fermscal[m, :] = pf.getdata(path, 2)['PSFSCALE']
        if gdat.recotype == 'rec8':
            irfn = pf.getdata(path, 1 + 3 * gdat.indxevtt[m])
            fermscal[m, :] = pf.getdata(path, 2 + 3 * gdat.indxevtt[m])['PSFSCALE']
        for k in range(numbpsfpform):
            fermform[:, m, k] = interp1d(enerirfn, mean(irfn[parastrg[k]].squeeze(), axis=0))(gdat.meanener)
    # convert N_tail to f_core
    for m in gdat.indxevtt:
        for i in gdat.indxener:
            fermform[i, m, 4] = 1. / (1. + fermform[i, m, 4] * fermform[i, m, 2]**2 / fermform[i, m, 0]**2)

    # calculate the scale factor
    gdat.fermscalfact = np.sqrt((fermscal[None, :, 0] * (10. * gdat.meanener[:, None])**fermscal[None, :, 2])**2 + fermscal[None, :, 1]**2)
    
    # store the fermi PSF parameters
    gdat.psfpexpr = np.zeros(gdat.numbener * gdat.numbevtt * numbpsfpform)
    for m in gdat.indxevtt:
        for k in range(numbpsfpform):
            indxfermpsfptemp = m * numbpsfpform * gdat.numbener + gdat.indxener * numbpsfpform + k
            #if k == 0 or k == 2:
            #    gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k] * gdat.fermscalfact[:, m]
            #else:
            #    gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k]
            gdat.psfpexpr[indxfermpsfptemp] = fermform[:, m, k]
    

def retr_fwhm(psfn, binsangl):

    if psfn.ndim == 1:
        numbener = 1
        numbevtt = 1
        psfn = psfn[None, :, None]
    else:
        numbener = psfn.shape[0]
        numbevtt = psfn.shape[2]
    indxener = np.arange(numbener)
    indxevtt = np.arange(numbevtt)
    wdth = np.zeros((numbener, numbevtt))
    for i in indxener:
        for m in indxevtt:
            indxanglgood = np.argsort(psfn[i, :, m])
            intpwdth = max(0.5 * amax(psfn[i, :, m]), amin(psfn[i, :, m]))
            if intpwdth > amin(psfn[i, indxanglgood, m]) and intpwdth < amax(psfn[i, indxanglgood, m]):
                wdth[i, m] = interp1d(psfn[i, indxanglgood, m], binsangl[indxanglgood])(intpwdth)
    return wdth


def retr_doubking(scaldevi, frac, sigc, gamc, sigt, gamt):

    psfn = frac / 2. / np.pi / sigc**2 * (1. - 1. / gamc) * (1. + scaldevi**2 / 2. / gamc / sigc**2)**(-gamc) + \
        (1. - frac) / 2. / np.pi / sigt**2 * (1. - 1. / gamt) * (1. + scaldevi**2 / 2. / gamt / sigt**2)**(-gamt)
    
    return psfn


def retr_path(strg, pathextndata=None, pathextnimag=None, rtag=None, onlyimag=False, onlydata=False):
    
    pathbase = os.environ[strg.upper() + '_DATA_PATH'] + '/'

    if not onlyimag:
        pathdata = pathbase
        if pathextndata != None:
            pathdata += pathextndata
        pathdata += 'data/'
        os.system('mkdir -p %s' % pathdata)
    if not onlydata:        
        pathimag = pathbase
        if pathextnimag != None:
            pathimag += pathextnimag
        pathimag += 'imag/'
        if rtag != None:
            pathimag += rtag + '/'
        os.system('mkdir -p %s' % pathimag)

    if not onlyimag and not onlydata:
        return pathimag, pathdata
    elif onlyimag:
        return pathimag
    else:
        return pathdata


def conv_rascdecl(args):

    rasc = args[0] * 8 + args[1] / 60. + args[2] / 3600.
    decl = args[3] + args[4] / 60. + args[5] / 3600.

    return rasc, decl


def smth_heal(maps, scalsmth, mpol=False, retrfull=False, numbsideoutp=None, indxpixlmask=None):

    if mpol:
        mpolsmth = scalsmth
    else:
        mpolsmth = 180. / scalsmth

    numbpixl = maps.size
    numbside = int(np.sqrt(numbpixl / 12))
    numbmpol = 3 * numbside
    maxmmpol = 3. * numbside - 1.
    mpolgrid, temp = hp.Alm.getlm(lmax=maxmmpol)
    mpol = np.arange(maxmmpol + 1.)
    
    if numbsideoutp == None:
        numbsideoutp = numbside
    
    if indxpixlmask != None:
        mapsoutp = copy(maps)
        mapsoutp[indxpixlmask] = hp.UNSEEN
        mapsoutp = hp.ma(mapsoutp)
        almctemp = hp.map2alm(maps)
    else:
        mapsoutp = maps
    
    almc = hp.map2alm(mapsoutp)

    wght = exp(-0.5 * (mpolgrid / mpolsmth)**2)
    almc *= wght
    
    mapsoutp = hp.alm2map(almc[np.where(mpolgrid < 3 * numbsideoutp)], numbsideoutp, verbose=False)

    if retrfull:
        return mapsoutp, almc, mpol, exp(-0.5 * (mpol / mpolsmth)**2)
    else:
        return mapsoutp


def smth_ferm(mapsinpt, meanener, recotype, maxmmpol=None, makeplot=False, kerntype='ferm'):

    numbpixl = mapsinpt.shape[1]
    numbside = int(np.sqrt(numbpixl / 12))
    if maxmmpol == None:
        maxmmpol = 3 * numbside - 1

    numbener = meanener.size
    indxener = np.arange(numbener)
    
    evtt, numbevtt, indxevtt = retr_evttferm(recotype)
    numbalmc = (maxmmpol + 1) * (maxmmpol + 2) / 2
    
    mapsoutp = np.empty_like(mapsinpt)
    if kerntype == 'gaus':
        angl = np.pi * np.linspace(0., 10., 100) / 180.
        
        if recotype != 'manu':
            gdat = gdatstrt()
            gdat.pathdata = os.environ["PCAT_DATA_PATH"] + '/data/'
            gdat.numbener = numbener
            gdat.indxener = indxener
            gdat.numbevtt = numbevtt
            gdat.indxevtt = indxevtt
            gdat.meanener = meanener
            gdat.recotype = recotype
            retr_psfpferm(gdat)
            gdat.exprtype = 'ferm'
            gdat.fittnumbpsfpform = 5
            gdat.fittnumbpsfptotl = 5
            gdat.binsangl = angl
            psfn = retr_psfn(gdat, gdat.psfpexpr, gdat.indxener, angl, 'doubking', None, False)
            fwhm = retr_fwhm(psfn, angl) 
        for i in indxener:
            for m in indxevtt:
                if recotype == 'manu':
                    if i == 0:
                        if m == 0:
                            sigm = 100.
                        if m == 1:
                            sigm = 100.
                    if i == 1:
                        if m == 0:
                            sigm = 1.05
                        if m == 1:
                            sigm = 0.7
                    if i == 2:
                        if m == 0:
                            sigm = 0.47
                        if m == 1:
                            sigm = 0.35
                    if i == 3:
                        if m == 0:
                            sigm = 0.06
                        if m == 1:
                            sigm = 0.04
                    if i == 4:
                        if m == 0:
                            sigm = 0.04
                        if m == 1:
                            sigm = 0.03
                    fwhmtemp = 2.355 * sigm * np.pi / 180.
                else:
                    fwhmtemp = fwhm[i, m]
                mapsoutp[i, :, m] = hp.smoothing(mapsinpt[i, :, m], fwhm=fwhmtemp)
    
    if kerntype == 'ferm':
        # get the beam
        beam = retr_beam(meanener, evtt, numbside, maxmmpol, recotype)
        # construct the transfer function
        tranfunc = ones((numbener, numbalmc, numbevtt))
        cntr = 0
        for n in np.arange(maxmmpol+1)[::-1]:
            tranfunc[:, cntr:cntr+n+1, :] = beam[:, maxmmpol-n:, :]
            cntr += n + 1


        indxener = np.arange(numbener)
        indxevtt = np.arange(numbevtt)
        for i in indxener:
            for m in indxevtt:
                
                # temp
                if i != 0 or m != 0:
                    continue
                # temp
                #mapsoutp[i, :, m] = hp.smoothing(mapsinpt[i, :, m], fwhm=radians(1.))
                #mapsoutp[i, :, m] = mapsinpt[i, :, m]
                almc = hp.map2alm(mapsinpt[i, :, m], lmax=maxmmpol)
                almc *= tranfunc[i, :, m]
                mapsoutp[i, :, m] = hp.alm2map(almc, numbside, lmax=maxmmpol)
    
    return mapsoutp


def retr_beam(meanener, evtt, numbside, maxmmpol, recotype, fulloutp=False, evaltype='invt'):
    
    numbpixl = 12 * numbside**2
    apix = 4. * np.pi / numbpixl

    numbener = meanener.size
    numbevtt = evtt.size
    indxevtt = np.arange(numbevtt)
    numbmpol = int(maxmmpol) + 1

    # alm of the delta function at the North Pole
    mapsinpt = np.zeros(numbpixl)
    mapsinpt[:4] = 1.
    mapsinpt /= sum(mapsinpt) * apix
    almcinpt = real(hp.map2alm(mapsinpt, lmax=maxmmpol)[:maxmmpol+1])
    
    # alm of the point source at the North Pole
    if evaltype != 'invt':
        lgalgrid, bgalgrid, numbpixl, apix = retr_healgrid(numbside)
        dir1 = np.array([lgalgrid, bgalgrid])
        dir2 = np.array([0., 90.])
        angl = hp.rotator.angdist(dir1, dir2, lonlat=True)
    else:
        angl = np.pi / np.linspace(1., maxmmpol, maxmmpol + 1)
    mapsoutp = retr_psfnferm(meanener, angl)
    if evaltype != 'invt':
        almcoutp = np.empty((numbener, maxmmpol+1, numbevtt))
        for i in range(numbener):
            for m in indxevtt:
                almcoutp[i, :, m] = real(hp.map2alm(mapsoutp[i, :, m], lmax=maxmmpol)[:maxmmpol+1])
        tranfunc = almcoutp / almcinpt[None, :, None]
        # temp
        tranfunc /= tranfunc[:, 0, :][:, None, :]
    else:    
        numbangl = angl.size
        matrdesi = np.empty((numbmpol, numbener, numbangl, numbevtt))
        tranfunc = np.empty((numbener, numbangl, numbevtt))
        for n in range(numbmpol):
            temp = 1. / np.sqrt(2. * n + 1.) * np.sqrt(4. * np.pi) / sp.special.lpmv(0, n, cos(angl))
            matrdesi[n, :, :, :] = temp[None, :, None]
        for i in range(numbener):
            for m in indxevtt:
                # temp
                if i != 0 or m != 0:
                    continue
                tranfunc[i, :, m] = matmul(linalg.inv(matrdesi[:, i, :, m]), mapsoutp[i, :, m])
        tranfunc = tranfunc**2
        for i in range(numbener):
            for m in indxevtt:
                tranfunc[i, :, m] /= tranfunc[i, 0, m]

    if fulloutp:
        return tranfunc, almcinpt, almcoutp
    else:
        return tranfunc


def plot_fermsmth():

    numbside = 256
    numbpixl = 12 * numbside**2
    maxmmpol = 3 * numbside - 1
    mpol = np.arange(maxmmpol + 1)
    
    listrecotype = ['rec7', 'rec8']
    for recotype in listrecotype:
        evtt, numbevtt, indxevtt = retr_evttferm(recotype)
        
        binsenerplot = np.array([0.3, 1., 3., 10.])
        meanenerplot = np.sqrt(binsenerplot[1:] * binsenerplot[:-1])
        numbenerplot = meanenerplot.size
        tranfunc, almcinpt, almcoutp = retr_beam(meanenerplot, evtt, numbside, maxmmpol, recotype, fulloutp=True)

        figr, axis = plt.subplots()

        plt.loglog(mpol, almcinpt, label='HealPix')
        plt.loglog(mpol, np.sqrt((2. * mpol + 1.) / 4. / np.pi), label='Analytic')
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/almcinpt.pdf'
        plt.legend(framealpha=1.)
        figr.savefig(path)
        plt.close(figr)
        
        figr, axis = plt.subplots()
        for i in np.arange(meanenerplot.size):
            for m in indxevtt:
                plt.loglog(mpol, almcoutp[i, :, m], label='$E=%.3g$, PSF%d' % (meanenerplot[i], indxevtt[m]))
        plt.legend(loc=3, framealpha=1.)
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/almcoutp.pdf'
        figr.savefig(path)
        plt.close(figr)
            
        figr, axis = plt.subplots()
        for i in np.arange(meanenerplot.size):
            for m in indxevtt:
                plt.loglog(mpol, tranfunc[i, :, m], label='$E=%.3g$, PSF%d' % (meanenerplot[i], indxevtt[m]))
        plt.legend(loc=3, framealpha=1.)
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/tranfunc.pdf'
        figr.savefig(path)
        plt.close(figr)
        
        maxmgang = 20.
        minmlgal = -maxmgang
        maxmlgal = maxmgang
        minmbgal = -maxmgang
        maxmbgal = maxmgang
            
        # get the Planck radiance map
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/HFI_CompMap_ThermalDustModel_2048_R1.20.fits'
        maps = pf.getdata(path, 1)['RADIANCE']
        mapstemp = hp.ud_grade(maps, numbside, order_in='NESTED', order_out='RING')
        maps = np.empty((numbenerplot, numbpixl, numbevtt))
        for i in np.arange(numbenerplot):
            for m in indxevtt:
                 maps[i, :, m] = mapstemp

        # smooth the map with the Fermi-LAT kernel
        mapssmthferm = smth_ferm(maps, meanenerplot, recotype)
        
        # smooth the map with the Gaussian kernel
        mapssmthgaus =  hp.sphtfunc.smoothing(mapstemp, sigma=np.deg2rad(0.5))

        # plot the maps
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/maps.pdf'
        plot_maps(path, mapstemp, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)

        for i in np.arange(meanenerplot.size):
            for m in indxevtt:
                path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/mapssmthferm%d%d.pdf' % (i, m)
                plot_maps(path, mapssmthferm[i, :, m], minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)
                
        path = os.environ["FERM_IGAL_DATA_PATH"] + '/imag/mapssmthgaus.pdf'
        plot_maps(path, mapssmthgaus, minmlgal=minmlgal, maxmlgal=maxmlgal, minmbgal=minmbgal, maxmbgal=maxmbgal)


def icdf_self(paraunit, minmpara, maxmpara):
    para = (maxmpara - minmpara) * paraunit + minmpara
    return para


def icdf_logt(paraunit, minmpara, maxmpara):
    
    para = minmpara * np.exp(paraunit * np.log(maxmpara / minmpara))
    
    return para


def icdf_atan(paraunit, minmpara, maxmpara):
    para = tan((arctan(maxmpara) - arctan(minmpara)) * paraunit + arctan(minmpara))
    return para


def icdf_gaus(cdfn, meanpara, stdvpara):
    
    para = meanpara + stdvpara * np.sqrt(2) * sp.special.erfinv(2. * cdfn - 1.)

    return para


def cdfn_self(para, minmpara, maxmpara):
    paraunit = (para - minmpara) / (maxmpara - minmpara)
    return paraunit


def cdfn_logt(para, minmpara, maxmpara):
    paraunit = np.log(para / minmpara) / np.log(maxmpara / minmpara)
    return paraunit


def cdfn_atan(para, minmpara, maxmpara):
    paraunit = (arctan(para) - arctan(minmpara)) / (arctan(maxmpara) - arctan(minmpara))
    return paraunit


def cdfn_samp(sampvarb, datapara, k=None):
    
    if k is None:
        samp = empty_like(sampvarb)
        for k in range(sampvarb.size):
            samp[k] = cdfn_samp_sing(sampvarb[k], k, datapara)
    else:
        samp = cdfn_samp_sing(sampvarb[k], k, datapara)
    return samp


def cdfn_samp_sing(sampvarb, k, datapara):
    
    if datapara.scal[k] == 'self':
        samp = cdfn_self(sampvarb, datapara.minm[k], datapara.maxm[k])
    if datapara.scal[k] == 'logt':
        samp = cdfn_logt(sampvarb, datapara.minm[k], datapara.maxm[k])
    if datapara.scal[k] == 'atan':
        samp = cdfn_atan(sampvarb, datapara.minm[k], datapara.maxm[k])
        
    return samp


def icdf_samp(samp, datapara, k=None):
    
    if k is None:
        sampvarb = empty_like(samp)
        for k in range(sampvarb.size):
            sampvarb[k] = icdf_samp_sing(samp[k], k, datapara)
    else:
        sampvarb = icdf_samp_sing(samp[k], k, datapara)
    return sampvarb


def icdf_samp_sing(samp, k, datapara):

    if datapara.scal[k] == 'self':
        sampvarb = icdf_self(samp, datapara.minm[k], datapara.maxm[k])
    if datapara.scal[k] == 'logt':
        sampvarb = icdf_logt(samp, datapara.minm[k], datapara.maxm[k])
    if datapara.scal[k] == 'atan':
        sampvarb = icdf_atan(samp, datapara.minm[k], datapara.maxm[k])
        
    return sampvarb


def gmrb_test(griddata):
    
    withvari = np.mean(var(griddata, 0))
    btwnvari = griddata.shape[0] * var(np.mean(griddata, 0))
    wgthvari = (1. - 1. / griddata.shape[0]) * withvari + btwnvari / griddata.shape[0]
    psrf = sqrt(wgthvari / withvari)

    return psrf


def retr_atcr_neww(listpara):

    numbsamp = listpara.shape[0]
    four = sp.fftpack.fft(listpara - np.mean(listpara, axis=0), axis=0)
    atcr = sp.fftpack.ifft(four * np.conjugate(four), axis=0).real
    atcr /= np.amax(atcr, 0)
    
    return atcr[:int(numbsamp/2), ...]


def retr_timeatcr(listpara, typeverb=1, atcrtype='maxm'):

    numbsamp = listpara.shape[0]
    listpara = listpara.reshape((numbsamp, -1))
    numbpara = listpara.shape[1]

    boolfail = False
    if listpara.shape[0] == 1:
        boolfail = True

    atcr = retr_atcr_neww(listpara)
    indxatcr = np.where(atcr > 0.2)
     
    if indxatcr[0].size == 0:
        boolfail = True
        timeatcr = 0
    else:
        if atcrtype == 'nomi':
            timeatcr = np.argmax(indxatcr[0], axis=0)
        if atcrtype == 'maxm':
            indx = np.argmax(indxatcr[0])
            indxtimemaxm = indxatcr[0][indx]
            indxparamaxm = indxatcr[1][indx]
            atcr = atcr[:, indxparamaxm]
            timeatcr = indxtimemaxm
   
    if boolfail:
        if atcrtype == 'maxm':
            return np.zeros((1, 1)), 0.
        else:
            return np.zeros((1, numbpara)), 0.
    else:
        return atcr, timeatcr


def retr_numbsamp(numbswep, numbburn, factthin):
    
    numbsamp = int((numbswep - numbburn) / factthin)
    
    return numbsamp


def plot_gmrb(path, gmrbstat):

    numbbinsplot = 40
    bins = np.linspace(1., np.amax(gmrbstat), numbbinsplot + 1)
    figr, axis = plt.subplots()
    axis.hist(gmrbstat, bins=bins)
    axis.set_title('Gelman-Rubin Convergence Test')
    axis.set_xlabel('PSRF')
    axis.set_ylabel('$N_p$')
    figr.savefig(path + 'gmrb.pdf')
    plt.close(figr)


def plot_atcr(path, atcr, timeatcr, strgextn=''):

    numbsampatcr = atcr.size
    
    figr, axis = plt.subplots(figsize=(6, 4))
    axis.plot(np.arange(numbsampatcr), atcr)
    axis.set_xlabel(r'$\tau$')
    axis.set_ylabel(r'$\xi(\tau)$')
    axis.text(0.8, 0.8, r'$\tau_{exp} = %.3g$' % timeatcr, ha='center', va='center', transform=axis.transAxes)
    axis.axhline(0., ls='--', alpha=0.5)
    plt.tight_layout()
    pathplot = path + 'atcr%s.pdf' % strgextn
    figr.savefig(pathplot)
    plt.close(figr)
    
        
def plot_propeffi(path, numbswep, numbpara, listaccp, listindxparamodi, namepara):

    indxlistaccp = np.where(listaccp == True)[0]
    binstime = np.linspace(0., numbswep - 1., 10)
    
    numbcols = 2
    numbrows = (numbpara + 1) / 2
    figr, axgr = plt.subplots(numbrows, numbcols, figsize=(16, 4 * (numbpara + 1)))
    if numbrows == 1:
        axgr = [axgr]
    for a, axrw in enumerate(axgr):
        for b, axis in enumerate(axrw):
            k = 2 * a + b
            if k == numbpara:
                axis.axis('off')
                break
            indxlistpara = np.where(listindxparamodi == k)[0]
            indxlistintc = intersect1d(indxlistaccp, indxlistpara, assume_unique=True)
            histotl = axis.hist(indxlistpara, binstime, color='b')
            histaccp = axis.hist(indxlistintc, binstime, color='g')
            axis.set_title(namepara[k])
    figr.subplots_adjust(hspace=0.3)
    figr.savefig(path + 'propeffi.pdf')
    plt.close(figr)


def plot_trac(path, listpara, labl, truepara=None, scalpara='self', titl=None, \
                        
                        # Boolean flag to overplot quantiles
                        boolplotquan=False, \
                        
                        listparadraw=None, listlabldraw=None, numbbinsplot=20, logthist=False, listcolrdraw=None):
    
    if not np.isfinite(listpara).all():
        return
    
    if not np.isfinite(listpara).all():
        raise Exception('')
    
    if listpara.size == 0:
        return

    maxmpara = np.amax(listpara)
    if scalpara == 'logt':
        minmpara = np.amin(listpara[np.where(listpara > 0.)])
        bins = icdf_logt(np.linspace(0., 1., numbbinsplot + 1), minmpara, maxmpara)
    else:
        minmpara = np.amin(listpara)
        bins = icdf_self(np.linspace(0., 1., numbbinsplot + 1), minmpara, maxmpara)
    limspara = np.array([minmpara, maxmpara])
        
    if boolplotquan:
        quanarry = sp.stats.mstats.mquantiles(listpara, prob=[0.025, 0.16, 0.84, 0.975])

    if scalpara == 'logt':
        numbtick = 5
        listtick = np.logspace(np.log10(minmpara), np.log10(maxmpara), numbtick)
        listlabltick = ['%.3g' % tick for tick in listtick]
    
    figr, axrw = plt.subplots(1, 2, figsize=(14, 7))
    if titl is not None:
        figr.suptitle(titl, fontsize=18)
    for n, axis in enumerate(axrw):
        if n == 0:
            axis.plot(listpara, lw=0.5)
            axis.set_xlabel('$i_{samp}$')
            axis.set_ylabel(labl)
            if truepara is not None and not np.isnan(truepara):
                axis.axhline(y=truepara, color='g', lw=4)
            if scalpara == 'logt':
                axis.set_yscale('log')
                axis.set_yticks(listtick)
                axis.set_yticklabels(listlabltick)
            axis.set_ylim(limspara)
            if listparadraw is not None:
                for k in range(len(listparadraw)):
                    axis.axhline(listparadraw[k], label=listlabldraw[k], color=listcolrdraw[k], lw=3)
            if boolplotquan:
                axis.axhline(quanarry[0], color='b', ls='--', lw=2)
                axis.axhline(quanarry[1], color='b', ls='-.', lw=2)
                axis.axhline(quanarry[2], color='b', ls='-.', lw=2)
                axis.axhline(quanarry[3], color='b', ls='--', lw=2)
        else:
            axis.hist(listpara, bins=bins)
            axis.set_xlabel(labl)
            if logthist:
                axis.set_yscale('log')
            axis.set_ylabel('$N_{samp}$')
            if truepara is not None and not np.isnan(truepara):
                axis.axvline(truepara, color='g', lw=4)
            if scalpara == 'logt':
                axis.set_xscale('log')
                axis.set_xticks(listtick)
                axis.set_xticklabels(listlabltick)
            axis.set_xlim(limspara)
            if listparadraw is not None:
                for k in range(len(listparadraw)):
                    axis.axvline(listparadraw[k], label=listlabldraw[k], color=listcolrdraw[k], lw=3)
            if boolplotquan:
                axis.axvline(quanarry[0], color='b', ls='--', lw=2)
                axis.axvline(quanarry[1], color='b', ls='-.', lw=2)
                axis.axvline(quanarry[2], color='b', ls='-.', lw=2)
                axis.axvline(quanarry[3], color='b', ls='--', lw=2)
                
    figr.subplots_adjust()#top=0.9, wspace=0.4, bottom=0.2)

    figr.savefig(path + '_trac.pdf')
    plt.close(figr)


def plot_plot(path, xdat, ydat, lablxdat, lablydat, scalxaxi, titl=None, linestyl=[None], colr=[None], legd=[None], **args):
    
    if not isinstance(ydat, list):
        listydat = [ydat]
    else:
        listydat = ydat

    figr, axis = plt.subplots(figsize=(6, 6))
    for k, ydat in enumerate(listydat):
        if k == 0:
            linestyl = '-'
        else:
            linestyl = '--'
        axis.plot(xdat, ydat, ls=linestyl, color='k', **args)
        # temp
        #axis.plot(xdat, ydat, ls=linestyl[k], color=colr[k], label=legd[k], **args)
    axis.set_ylabel(lablydat)
    axis.set_xlabel(lablxdat)
    if scalxaxi == 'logt':
        axis.set_xscale('log')
    if titl is not None:
        axis.set_title(titl)
    plt.tight_layout()
    figr.savefig(path + '.pdf')
    plt.close(figr)


def plot_hist( \
              path, listpara, strg, titl=None, numbbins=20, truepara=None, \
              
              # Boolean flag to overplot quantiles
              boolplotquan=False, 
              
              # type of the file for plots
              typefileplot='pdf', \
                                            scalpara='self', listparadraw=None, listlabldraw=None, listcolrdraw=None):

    minmvarb = np.amin(listpara)
    maxmvarb = np.amax(listpara)
    if scalpara == 'logt':
        bins = icdf_logt(np.linspace(0., 1., numbbins + 1), minmvarb, maxmvarb)
    else:
        bins = icdf_self(np.linspace(0., 1., numbbins + 1), minmvarb, maxmvarb)
    figr, axis = plt.subplots(figsize=(6, 6))
    axis.hist(listpara, bins=bins)
    axis.set_ylabel(r'$N_{samp}$')
    axis.set_xlabel(strg)
    if truepara is not None:
        axis.axvline(truepara, color='g', lw=4)
    if listparadraw is not None:
        for k in range(len(listparadraw)):
            axis.axvline(listparadraw[k], label=listlabldraw[k], color=listcolrdraw[k], lw=3)
    if boolplotquan:
        quanarry = sp.stats.mstats.mquantiles(listpara, prob=[0.025, 0.16, 0.84, 0.975])
        axis.axvline(quanarry[0], color='b', ls='--', lw=2)
        axis.axvline(quanarry[1], color='b', ls='-.', lw=2)
        axis.axvline(quanarry[2], color='b', ls='-.', lw=2)
        axis.axvline(quanarry[3], color='b', ls='--', lw=2)
    if titl is not None:
        axis.set_title(titl)
    plt.tight_layout()
    figr.savefig(path + '_hist.%s' % typefileplot)
    plt.close(figr)


def retr_limtpara(scalpara, minmpara, maxmpara, meanpara, stdvpara):
    
    numbpara = len(scalpara)
    limtpara = np.empty((2, numbpara))
    indxpara = np.arange(numbpara)
    for n in indxpara:
        if scalpara[n] == 'self' or scalpara[n] == 'logt':
            limtpara[0, n] = minmpara[n]
            limtpara[1, n] = maxmpara[n]
        if scalpara[n] == 'gaus':
            limtpara[0, n] = meanpara[n] - 10 * stdvpara[n]
            limtpara[1, n] = meanpara[n] + 10 * stdvpara[n]
    
    return limtpara


def retr_lpos(para, *dictlpos):
     
    gdat, indxpara, scalpara, minmpara, maxmpara, meangauspara, stdvgauspara, retr_llik, retr_lpri = dictlpos
    
    boolreje = False
    for k in indxpara:
        if scalpara[k] != 'gaus':
            if para[k] < minmpara[k] or para[k] > maxmpara[k]:
                lpos = -np.inf
                boolreje = True
    
    if not boolreje:
        llik = retr_llik(para, gdat)
        lpri = 0.
        if retr_lpri is None:
            for k in indxpara:
                if scalpara[k] == 'gaus':
                    lpri += (para[k] - meangauspara[k]) / stdvgauspara[k]**2
        else:
            lpri = retr_lpri(para, gdat)
        lpos = llik + lpri
    
    #print('lpos')
    #print(lpos)
    #print('')
    
    return lpos


def retr_icdfunif(cdfn, minm, maxm):

    icdf = minm + cdfn * (maxm - minm)
    
    return icdf


def retr_icdf(cdfn, scalpara, minm, maxm):

    numbpara = len(scalpara)
    indxpara = np.arange(numbpara)
    icdf = np.empty(numbpara)
    for k in indxpara:
        if scalpara[k] == 'self':
            icdf[k] = minm[k] + cdfn[k] * (maxm[k] - minm[k])
    
    return icdf


def opti(pathimag, retr_llik, minmpara, maxmpara, numbtopp=3, numbiter=5):

    numbsamp = 4
    indxsamp = np.arange(numbsamp)
    numbpara = minmpara.size
    indxiter = np.arange(numbiter)
    indxtopp = np.arange(numbtopp)

    # seeds
    listfact = []
    listparacent = []
    listopen = []
    listllikmaxmseed = []
    # all samples
    #listllik = np.empty(0)
    #listpara = np.empty((0, numbpara))
    
    for i in indxiter:
        print('i')
        print(i)
        #print('listllik')
        #print(listllik)
        #print('listpara')
        #print(listpara)
        
        print('listfact')
        print(listfact)
        print('listparacent')
        print(listparacent)
        print('listllikmaxmseed')
        print(listllikmaxmseed)
        print('listopen')
        print(listopen)
        if i == 0:
            minmparatemp = minmpara
            maxmparatemp = maxmpara
            thisindxseed = 0
            paramidi = (maxmpara + minmpara) / 2.
            listfact.append(1.)
            listparacent.append(paramidi)
            listopen.append(True)
            listllikmaxmseed.append([])
        else:
            indxopen = np.where(listopen)[0]
            thisindxseed = np.random.choice(indxopen)
            print('thisindxseed')
            print(thisindxseed)
            maxmparatemp = listparacent[thisindxseed] + listfact[thisindxseed] * (maxmpara - listparacent[thisindxseed])
            minmparatemp = listparacent[thisindxseed] - listfact[thisindxseed] * (listparacent[thisindxseed] - minmpara)
        print('thisindxseed')
        print(thisindxseed)
        para = np.random.rand(numbpara * numbsamp).reshape((numbsamp, numbpara)) * (maxmpara[None, :] - minmpara[None, :]) + minmpara[None, :]
        print('para')
        print(para)
        
        print('Evaluating samples...')
        llik = np.empty(numbsamp)
        for k in indxsamp:
            llik[k] = retr_llik(para[k, :])
        print('llik')
        print(llik)
        listllikmaxmseed[thisindxseed] = np.amax(llik)
        
        # add new seeds
        if i == 0:
            for k in indxtopp:
                # factor
                listfact.append(0.5 * listfact[thisindxseed])
                
                # parameters of the seeds
                indxsampsort = np.argsort(llik)[::-1]
                listparacent.append(para[indxsampsort, :])
                
                # determine if still open
                boolopen = np.amax(llik) >= listllikmaxmseed[thisindxseed]
                listopen.append(boolopen)
                
                listllikmaxmseed.append([])

        #listllik = np.concatenate((listllik, llik))
        #listpara = np.concatenate((listpara, para), 0)
        
        #print('listllik')
        #print(listllik)
        #print('listpara')
        #print(listpara)
        print('listfact')
        print(listfact)
        print('listparacent')
        print(listparacent)
        print('listllikmaxmseed')
        print(listllikmaxmseed)
        print('listopen')
        print(listopen)
        print('')
        print('')
        print('')
        if not np.array(listopen).any():
            break
    
    return listparatopp


def retr_modl_corr(gdat, feat, inpt):
    '''
    Return a linear model
    '''    
    angl = feat[0]
    gamm = feat[1]

    slop = -1. / np.tan(angl)
    intc = gamm / np.sin(angl)
    
    line = inpt * slop + intc
    
    return line, []


def retr_llik_corr(feat, gdat):
    '''
    Return the likelihood for a linear model
    '''
    angl = feat[0]
    gamm = feat[1]

    dist = np.cos(angl) * gdat.tempfrst + np.sin(angl) * gdat.tempseco - gamm
    vari = np.cos(angl)**2 * gdat.tempfrststdv**2 + np.sin(angl)**2 * gdat.tempsecostdv**2

    llik = -0.5 * np.sum(dist**2 / vari)

    return llik


def samp( \
         gdat, \

         numbsampwalk, \

         retr_llik, \
              
         # model parameters
         ## list of names of parameters
         listnamepara, \
         ## list of labels of parameters
         listlablpara, \
         ## list of scalings of parameters
         scalpara, \
         ## list of minima of parameters
         minmpara, \
         ## list of maxima of parameters
         maxmpara, \
         
         # base path for placing plot and data files
         pathbase=None, \
         
         # Boolean flag to enforce a reprocess and overwrite
         boolforcrepr=False, \

         # Boolean flag to make plots
         boolplot=True, \
         
         numbsamppostwalk=None, \

         meangauspara=None, \
         stdvgauspara=None, \
         retr_lpri=None, \
         # Boolean flag to turn on multiprocessing
         boolmult=False, \
         # burn-in
         ## number of samples in a precursor run whose final state will be used as the initial state of the actual sampler
         numbsampburnwalkinit=0, \
         ## number of initial samples to be burned
         numbsampburnwalk=0, \
         # function to return derived variables from the parameter vector
         retr_dictderi=None, \
         # dictionary of labels and scalings for derived parameters
         dictlablscalparaderi=None, \
         # Boolean flag to diagnose the code
         booldiag=True, \
         # a string used to save and retrieve results
         strgextn='', \
         typesamp='mcmc', \

         # type of the file for plots
         typefileplot='png', \
         
         # verbosity level
         typeverb=1, \
        ):
    '''
    Sample from a posterior using MCMC
    '''
    
    numbpara = len(listlablpara)
   
    if numbsampwalk <= numbsampburnwalk:
        raise Exception('Burn-in samples cannot outnumber samples.')
    
    if isinstance(minmpara, list):
        minmpara = np.array(minmpara)

    if isinstance(maxmpara, list):
        maxmpara = np.array(maxmpara)

    if numbpara != minmpara.size:
        raise Exception('')
    if numbpara != maxmpara.size:
        raise Exception('')
    
    if dictlablscalparaderi is not None and retr_dictderi is None: 
        raise Exception('')

    indxpara = np.arange(numbpara)
    
    if typesamp == 'mcmc':
        numbwalk = max(20, 2 * numbpara)
        indxwalk = np.arange(numbwalk)
        numbsamptotl = numbsampwalk * numbwalk
        
        if typeverb > 0:
            print('numbsampwalk')
            print(numbsampwalk)
            print('numbwalk')
            print(numbwalk)
            print('numbsamptotl')
            print(numbsamptotl)

    if pathbase is not None:
        pathbasesamp = pathbase + '%s/' % typesamp
        pathimag = pathbasesamp + 'imag/'
        pathdata = pathbasesamp + 'data/'
        os.system('mkdir -p %s' % pathimag)
        os.system('mkdir -p %s' % pathdata)

    # plotting
    ## plot limits 
    limtpara = retr_limtpara(scalpara, minmpara, maxmpara, meangauspara, stdvgauspara)

    ## plot bins
    numbbins = 20
    indxbins = np.arange(numbbins)
    binspara = np.empty((numbbins + 1, numbpara)) 
    for k in indxpara:
        if limtpara[0, k] == limtpara[1, k]:
            raise Exception('')

        binspara[:, k] = np.linspace(limtpara[0, k], limtpara[1, k], numbbins + 1)
    meanpara = (binspara[1:, :] + binspara[:-1, :]) / 2.
    
    for k in indxpara:
        if minmpara[k] >= maxmpara[k]:
            raise Exception('')
    
    dictlpos = [gdat, indxpara, scalpara, minmpara, maxmpara, meangauspara, stdvgauspara, retr_llik, retr_lpri]
    
    if numbsamppostwalk is None:
        numbsamppostwalk = 100

    if typeverb > 0:
        print('listnamepara')
        print(listnamepara)
        print('listlablpara')
        print(listlablpara)
        print('scalpara')
        print(scalpara)
        if minmpara is not None:
            print('minmpara')
            print(minmpara)
            print('maxmpara')
            print(maxmpara)
        if meangauspara is not None:
            print('meangauspara')
            print(meangauspara)
            print('stdvgauspara')
            print(stdvgauspara)
    
    if strgextn is not '':
        strgextn = '_%s' % strgextn
    
    # path of the posterior
    if pathbase is not None:
        pathpost = pathdata + 'postpara%s.csv' % strgextn
        pathdict = pathdata + 'samppostpara%s.csv' % strgextn
        pathpickderi = pathdata + 'postderi%s.pickle' % strgextn
    
    if boolforcrepr or not boolforcrepr and (pathbase is None or not os.path.exists(pathdict)):
        
        # initialize
        if pathbase is not None and os.path.exists(pathpost) and np.loadtxt(pathpost, delimiter=',').shape[0] == limtpara.shape[1]:
            print('Reading the initial state from %s...' % pathpost)
            parainitcent = np.loadtxt(pathpost, delimiter=',')[:, 0]
            paraunitinitcent = np.empty_like(parainitcent)
            for m in indxpara:
                if scalpara[m] == 'self':
                    paraunitinitcent[m] = cdfn_self(parainitcent[m], limtpara[0, m], limtpara[1, m])
                if scalpara[m] == 'logt':
                    paraunitinitcent[m] = cdfn_logt(parainitcent[m], limtpara[0, m], limtpara[1, m])
                if scalpara[m] == 'gaus':
                    paraunitinitcent[m] = cdfn_logt(parainitcent[m], meanpara[m], gauspara[m])
        else:
            paraunitinitcent = np.full(numbpara, 0.5)
        
        parainit = [np.empty(numbpara) for k in indxwalk]
        for k in indxwalk:
            for m in indxpara:
                paraunit = paraunitinitcent[m] + 0.05 * scipy.stats.norm.rvs()
                paraunit = paraunit % 1.
                if scalpara[m] == 'self':
                    parainit[k][m] = icdf_self(paraunit, limtpara[0, m], limtpara[1, m])
                if scalpara[m] == 'logt':
                    if booldiag:
                        if limtpara[0, m] <= 0.:
                            raise Exception('')
                        if limtpara[1, m] <= 0.:
                            raise Exception('')
                    parainit[k][m] = icdf_logt(paraunit, limtpara[0, m], limtpara[1, m])
                if scalpara[m] == 'gaus':
                    parainit[k][m] = icdf_logt(paraunit, meanpara[m], gauspara[m])

        if typesamp == 'mcmc':
            if typeverb > 0:
                progress = True
            else:
                progress = False
        
            import emcee

            numbsamp = numbwalk * numbsampwalk
            indxsampwalk = np.arange(numbsampwalk)
            indxsamp = np.arange(numbsamp)
            numbsampburn = numbsampburnwalkinit * numbwalk
            if booldiag:
                if numbsampwalk == 0:
                    raise Exception('')
        
            if boolmult:
                pool = multiprocessing.Pool(processes=multiprocessing.cpu_count()-1)
            else:
                pool = None
            objtsamp = emcee.EnsembleSampler(numbwalk, numbpara, retr_lpos, args=dictlpos, pool=pool)
            if numbsampburnwalkinit > 0:
                parainitburn, prob, state = objtsamp.run_mcmc(parainit, numbsampburnwalkinit, progress=progress)
                if typeverb == 1:
                    print('Parameter states from the burn-in:')
                    print('parainitburn')
                    print(parainitburn)
                parainit = np.array(parainitburn)
                indxwalkmpos = np.argmax(objtsamp.lnprobability[:, -1], 0)
                parainittemp = parainit[indxwalkmpos, :]
                parainitburn = [[[] for m in indxpara] for k in indxwalk]
                for m in indxpara:
                    for k in indxwalk:
                        parainitburn[k][m] = parainittemp[m] * (1. + 1e-5 * np.random.randn())
                objtsamp.reset()
            else:
                parainitburn = parainit
            objtsamp.run_mcmc(parainitburn, numbsampwalk, progress=progress)
            listlposwalk = objtsamp.lnprobability
            listparafittwalk = objtsamp.chain
            
            # get rid of burn-in and thin
            indxsampwalkkeep = np.linspace(numbsampburnwalk, numbsampwalk - 1, numbsamppostwalk).astype(int)
            listparafitt = listparafittwalk[:, indxsampwalkkeep, :].reshape((-1, numbpara))
            
            numbsampkeep = listparafitt.shape[0]
            indxsampkeep = np.arange(numbsampkeep)
            
            listparaderi = None
            dictvarbderi = dict()
            if retr_dictderi is not None:
                numbsampderi = 10
                indxsampderi = np.random.choice(indxsampkeep, size=numbsampderi)
                #numbsampderi = indxsampderi.size

                print('Evaluating derived variables...')
                listdictvarbderi = [[] for n in indxsampderi]
                for n in tqdm(range(numbsampderi)):
                    listdictvarbderi[n] = retr_dictderi(listparafitt[indxsampderi[n], :], gdat)

                print('Placing the evaluated derived variables in the output dictionary...')
                for strg, valu in listdictvarbderi[0].items():
                    if booldiag and isinstance(valu, list):
                        print('strg')
                        print(strg)
                        print('valu')
                        print(valu)
                        raise Exception('')

                    if valu is None:
                        print('Encountered a None derived variable (%s). Skipping...' % strg)
                        continue
                    if np.isscalar(valu):
                        dictvarbderi[strg] = np.empty((numbsampkeep, 1))
                    else:
                        dictvarbderi[strg] = np.empty([numbsampderi] + list(valu.shape))
                    
                    for n in range(numbsampderi):
                        dictvarbderi[strg][n, ...] = listdictvarbderi[n][strg]
                
                listnamevarbderi = listdictvarbderi[0].keys()
                listnameparaderi = []
                for name, valu in listdictvarbderi[0].items():
                    if np.isscalar(valu) and valu is not None:
                        listnameparaderi.append(name)
                numbparaderi = len(listnameparaderi)
                print('Placing the evaluated derived parameters in the output dictionary...')
                listparaderi = np.empty((numbsampkeep, numbparaderi)) 
                k = 0
                for strg, valu in listdictvarbderi[0].items():
                    if valu is None:
                        continue
                    for n in range(numbsampderi):
                        dictvarbderi[strg][n, ...] = listdictvarbderi[n][strg]
                    
                    #if np.isscalar(listdictvarbderi[n][strg]) and numbsamp == numbsampderi:
                    #    for n in indxsampkeep:
                    #        print('')
                    #        print('')
                    #        print('')
                    #        print('n, k')
                    #        print(n, k)
                    #        print('strg')
                    #        print(strg)
                    #        print('listparaderi')
                    #        summgene(listparaderi)
                    #        print('listdictvarbderi')
                    #        summgene(listdictvarbderi)

                    #        listparaderi[n, k] = listdictvarbderi[n][strg]
                    #        if n == numbsampkeep - 1:
                    #            k += 1
                
            indxsampwalk = np.arange(numbsampwalk)
            
            if boolplot and pathimag is not None:
                # plot the posterior
                ### trace
                figr, axis = plt.subplots(numbpara + 1, 1, figsize=(12, (numbpara + 1) * 4))
                for i in indxwalk:
                    axis[0].plot(indxsampwalk, listlposwalk[i, :])
                axis[0].axvline(numbsampburnwalk, color='k')
                axis[0].set_ylabel('log P')
                for k in indxpara:
                    for i in indxwalk:
                        axis[k+1].plot(indxsampwalk, listparafittwalk[i, :, k])
                    labl = listlablpara[k][0]
                    if listlablpara[k][1] != '':
                        labl += ' [%s]' % listlablpara[k][1]
                    axis[k+1].axvline(numbsampburnwalk, color='k')
                    axis[k+1].set_ylabel(labl)
                path = pathimag + 'trac%s.%s' % (strgextn, typefileplot)
                if typeverb == 1:
                    print('Writing to %s...' % path)
                plt.savefig(path)
                plt.close()
                
                # plot the posterior
                ### trace
                if numbsampburnwalk > 0:
                    figr, axis = plt.subplots(numbpara + 1, 1, figsize=(12, (numbpara + 1) * 4))
                    for i in indxwalk:
                        axis[0].plot(indxsampwalk[numbsampburnwalk:], listlposwalk[i, numbsampburnwalk:])
                    axis[0].set_ylabel('log P')
                    for k in indxpara:
                        for i in indxwalk:
                            axis[k+1].plot(indxsampwalk[numbsampburnwalk:], listparafittwalk[i, numbsampburnwalk:, k])
                        labl = listlablpara[k][0]
                        if listlablpara[k][1] != '':
                            labl += ' [%s]' % listlablpara[k][1]
                        axis[k+1].set_ylabel(labl)
                    path = pathimag + 'tracgood%s.%s' % (strgextn, typefileplot)
                    if typeverb == 1:
                        print('Writing to %s...' % path)
                    plt.savefig(path)
                    plt.close()
        
        if boolplot and pathimag is not None:
            ## joint PDF
            strgplot = 'postparafitt' + strgextn
            plot_grid(pathimag, strgplot, listparafitt, listlablpara, numbbinsplot=numbbins)
            
            # derived
            if dictlablscalparaderi is not None:
                listnameparaderi = dictlablscalparaderi.keys()
                listlablparaderi = []
                for name in listnameparaderi:
                    listlablparaderi.append(dictlablscalparaderi[name][0])
                listlablparatotl = listlablpara + listlablparaderi
                listparatotl = np.concatenate([listparafitt, listparaderi], 1)
                strgplot = 'postparaderi' + strgextn
                plot_grid(pathimag, strgplot, listparaderi, listlablparaderi, numbbinsplot=numbbins)
                strgplot = 'postparatotl' + strgextn
                plot_grid(pathimag, strgplot, listparatotl, listlablparatotl, numbbinsplot=numbbins)
            else:
                listparatotl = listparafitt
        
        if pathbase is not None:
            if dictlablscalparaderi is not None:
                numbparapost = listparatotl.shape[1]
            else:
                numbparapost = numbpara
            print('Writing the posterior to %s...' % pathpost)
            arry = np.empty((numbparapost, 3))
            arry[:, 0] = np.median(listparatotl, 0)
            arry[:, 1] = np.percentile(listparatotl, 84.) - arry[:, 0]
            arry[:, 2] = arry[:, 0] - np.percentile(listparatotl, 16.)
            np.savetxt(pathpost, arry, delimiter=',')
            
            print('Writing the posterior derived variables to %s...' % pathpickderi)
            with open(pathpickderi, 'wb') as objtfile:
                pickle.dump(dictvarbderi, objtfile, protocol=pickle.HIGHEST_PROTOCOL)


        dictsamp = dict()
        for k, name in enumerate(listnamepara):
            dictsamp[name] = listparafitt[:, k]
        dictsamp['lpos'] = listlposwalk[:, indxsampwalkkeep].flatten()
        pd.DataFrame.from_dict(dictsamp).to_csv(pathdict, index=False)
        for name, valu in dictvarbderi.items():
            dictsamp[name] = valu

    else:
        if typeverb > 0:
            print('Reading from %s...' % pathdict)
        dictsamp = pd.read_csv(pathdict).to_dict(orient='list')
        for name in dictsamp.keys():
            dictsamp[name] = np.array(dictsamp[name])
        
        if retr_dictderi is not None:
            objtfile = open(pathpickderi, "rb")
            if typeverb > 0:
                print('Reading from %s...' % pathpickderi)
            dictvarbderi = pickle.load(objtfile)
            for strgvarbderi, valuvarbderi in dictvarbderi.items():
                dictsamp[strgvarbderi] = valuvarbderi

    return dictsamp


def setp_para_defa(gdat, strgmodl, strgvarb, valuvarb):
    
    if strgmodl == 'fitt' or strgmodl == 'true':
        setp_para_defa_wrap(gdat, strgmodl, strgvarb, valuvarb)
    elif strgmodl == 'full':
        for strgmodl in gdat.liststrgmodl:
            setp_para_defa_wrap(gdat, strgmodl, strgvarb, valuvarb)


def retr_matrrota(thet, \
                  # Boolean indicating that the input is in radians
                  boolradn=False, \
                 ):
    
    matrrota = np.array([[np.cos(thet), np.sin(thet)], [np.sin(thet), np.cos(thet)]])
    
    return matrrota


def retr_sperfromcart(xaxi, yaxi, zaxi):
    '''
    Convert Cartesian coordinates to spherical coordinates
    '''
    radi = np.sqrt(xaxi**2 + yaxi**2 + zaxi**2)
    phii = (360. + np.arctan2(yaxi, xaxi) * 180. / np.pi) % 360.
    thet = 90. - np.arccos(zaxi / radi) * 180. / np.pi
    
    return phii, thet, radi


def setp_para_defa_wrap(gdat, strgmodl, strgvarb, valuvarb):
    
    gmod = getattr(gdat, strgmodl)
    
    if strgmodl == 'fitt':
        dicttemp = gdat.dictfitt
    if strgmodl == 'true':
        dicttemp = gdat.dicttrue

    if strgvarb in dicttemp:
        valuvarbdefa = dicttemp[strgvarb]
    else:
        valuvarbdefa = valuvarb
    
    setattr(gmod, strgvarb, valuvarbdefa)


def plot_grid_diag(k, axis, listpara, truepara, listparadraw, boolplotquan, \
                            listlablpara, listtypeplottdim, indxpopl, listcolrpopl, listmrkrpopl, listlablpopl, boolmakelegd, listsizepopl, bins=None):
                    
    for u in indxpopl:
        if indxpopl.size > 1:
            labl = listlablpopl[u]
        else:
            labl = None
        axis.hist(listpara[u][:, k], bins=bins[k], label=labl, \
                                                               edgecolor=matplotlib.colors.to_rgba(listcolrpopl[u], 1.0), lw=2, ls='-', \
                                                               facecolor=matplotlib.colors.to_rgba(listcolrpopl[u], 0.2))
    if boolmakelegd and indxpopl.size > 1:
        axis.legend(framealpha=1.)
    
    if truepara is not None and truepara[k] is not None and not np.isnan(truepara[k]):
        axis.axvline(truepara[k], color='g', lw=4)
    
    # draw the provided reference values
    if listparadraw is not None:
        for m in indxdraw:
            axis.axvline(listparadraw[m][k], color='r', lw=3)
    
    if boolplotquan:
        quan = np.empty(4)
        quan[0] = np.nanpercentile(listpara[0][:, k], 2.5)
        quan[1] = np.nanpercentile(listpara[0][:, k], 16.)
        quan[2] = np.nanpercentile(listpara[0][:, k], 84.)
        quan[3] = np.nanpercentile(listpara[0][:, k], 97.5)
        axis.axvline(quan[0], color='r', ls='--', lw=2)
        axis.axvline(quan[1], color='r', ls='-.', lw=2)
        axis.axvline(quan[2], color='r', ls='-.', lw=2)
        axis.axvline(quan[3], color='r', ls='--', lw=2)
        medivarb = np.nanmedian(listpara[0][:, k])
    
        if listlablpara[k][1] != '':
            strgunit = ' ' + listlablpara[k][1]
        else:
            strgunit = ''
        axis.set_title(r'%s = %.3g $\substack{+%.2g \\\\ -%.2g}$ %s' % (listlablpara[k][0], medivarb, quan[2] - medivarb, medivarb - quan[1], strgunit))
    
    axis.set_yscale('log')
    
                

def plot_grid_pair(k, l, axis, limt, listmantlabl, listpara, truepara, listparadraw, boolplotquan, listlablpara, \
                                     listscalpara, boolsqua, listvectplot, listtypeplottdim, indxpopl, listcolrpopl, \
                                     listmrkrpopl, listcolrpopltdim, listlablpopl, boolmakelegd, bins=None, \
                                     listlablsamp=None, boolcbar=True, \
                                    ):
    
    if not np.isfinite(limt).all():
        print('limt')
        print(limt)
        raise Exception('')
    
    if (listtypeplottdim == 'hist').any():
        binstemp = [bins[l], bins[k]]

    for u in indxpopl:
         
        if listpara[u][:, l].size == 0:
            continue
            
        if indxpopl.size > 1:
            labl = listlablpopl[u]
            alph = 0.5
        else:
            labl = None
            alph = 1.

        if listtypeplottdim[u] == 'scat':
        
            if listpara[u][:, l].size >= 1e6:
                print('Warning! Skipped the scatter plot because there are too many points to plot!')
                return

            axis.scatter(listpara[u][:, l], listpara[u][:, k], s=1, color=listcolrpopl[u], label=labl, alpha=alph, marker=listmrkrpopl[u])
        
            # add text labels on outliers
            if listlablsamp is not None and indxpopl.size == 1:
                
                # remove infinite samples
                indx = np.where(np.isfinite(listpara[u][:, l]) & np.isfinite(listpara[u][:, k]))[0]
                listparapair = listpara[u][indx, :]
                print('listlablsamp[u]')
                print(listlablsamp[u])
                print('indx')
                summgene(indx)
                listlablsamppair = listlablsamp[u][indx]
                listparapair = listparapair[:, np.array([l, k])]
                
                # transform from data to axis coordinate positions
                ## display coordinate positions
                posidisp = axis.transData.transform(listparapair)
                ## axis coordinate positions
                listparapairoutl = axis.transAxes.inverted().transform(posidisp)
                
                from sklearn.neighbors import LocalOutlierFactor
                numbsamp = listparapairoutl.shape[0]
                
                numboutf = min(0, numbsamp)
                if numboutf > 0:
                    if numbsamp > numboutf:
                        n_neighbors = min(numbsamp, 100)
                        objtfore = LocalOutlierFactor(n_neighbors=n_neighbors)
                        objtfore.fit(listparapairoutl)
                        louf = objtfore.negative_outlier_factor_
                    else:
                        louf = np.zeros(numbsamp)
                    listindxsamplouf = np.argsort(louf)[:numboutf]
                else:
                    listindxsamplouf = np.array([], dtype=int)

                indxtemp = np.array([np.argmin(listparapair[:, 0]), np.argmin(listparapair[:, 1]), np.argmax(listparapair[:, 0]), np.argmax(listparapair[:, 1])])
                listindxsamplouf = np.concatenate([listindxsamplouf, indxtemp])
                listindxsamplouf = np.unique(listindxsamplouf)
                numboutf = listindxsamplouf.size
                listparapair = listparapair[listindxsamplouf, :]
                listlablsamppair = listlablsamppair[listindxsamplouf]
            
                ## coordinates of the samples to be annotated
                if listscalpara[l] == 'self' or listscalpara[l] == 'gaus':
                    xpossamp = cdfn_self(listparapair[:, 0], limt[0, l], limt[1, l])
                elif listscalpara[l] == 'logt':
                    xpossamp = cdfn_logt(listparapair[:, 0], limt[0, l], limt[1, l])
                elif listscalpara[l] == 'atan':
                    xpossamp = cdfn_atan(listparapair[:, 0], limt[0, l], limt[1, l])
                else:
                    raise Exception('Unrecognized scaling: %s' % listscalpara[l])
                
                if listscalpara[k] == 'self' or listscalpara[k] == 'gaus':
                    ypossamp = cdfn_self(listparapair[:, 1], limt[0, k], limt[1, k])
                elif listscalpara[k] == 'logt':
                    ypossamp = cdfn_logt(listparapair[:, 1], limt[0, k], limt[1, k])
                elif listscalpara[k] == 'atan':
                    ypossamp = cdfn_atan(listparapair[:, 1], limt[0, k], limt[1, k])
                else:
                    raise Exception('Unrecognized scaling: %s' % listscalpara[k])
                
                # place larger markers at the positions of the outliers
                axis.scatter(listparapair[:, 0], listparapair[:, 1], s=3, color=listcolrpopl[u], marker=listmrkrpopl[u])
                
                print('Automatically positioning the annotations...')
                # automatically position the annotations

                # tunable parameters
                sizelablxpos = 0.3 
                sizelablypos = 0.1
                
                # minimum horizontal position of the label
                minmlablxpos = 0.#0.5 * sizelablxpos
                # minimum vertical position of the label
                minmlablypos = 0.#0.5 * sizelablypos

                # maximum horizontal position of the label
                maxmlablxpos = 1.# - 0.5 * sizelablxpos
                # maximum vertical position of the label
                maxmlablypos = 1.# - 0.5 * sizelablypos
                
                # list of trial maximum horizontal distances between the sample and label
                listdistxpos = np.array([1.5, 3., 100.]) * sizelablxpos
                # list of trial maximum vertical distances between the sample and label
                listdistypos = np.array([1.5, 3., 100.]) * sizelablypos
                for ll in range(len(listdistxpos)):
                    distxpos = listdistxpos[ll]
                    distypos = listdistypos[ll]
                    
                    numbtria = 0
                    while True:
                        
                        # trial coordinates of the annotations
                        listxposlabl = xpossamp + distxpos * (2. * np.random.rand(numboutf) - 1.)
                        listyposlabl = ypossamp + distypos * (2. * np.random.rand(numboutf) - 1.)
                        
                        indxxpos = np.where((listxposlabl > maxmlablxpos) | (listxposlabl < minmlablxpos))[0]
                        listxposlabl[indxxpos] = minmlablxpos + listxposlabl[indxxpos] % (maxmlablxpos - minmlablxpos)

                        indxypos = np.where((listyposlabl > maxmlablypos) | (listyposlabl < minmlablypos))[0]
                        listyposlabl[indxypos] = minmlablypos + listyposlabl[indxypos] % (maxmlablypos - minmlablypos)

                        # concatenated list of coordinates of samples and trial annotations
                        xpos = np.concatenate((listxposlabl, xpossamp))
                        ypos = np.concatenate((listyposlabl, ypossamp))
                        
                        # check if each trial annotation is sufficiantly distant from all samples and other trial annotations
                        boolgood = True
                        for n in range(numboutf):
                            xposdiff = abs(listxposlabl[n] - xpos)
                            yposdiff = abs(listyposlabl[n] - ypos)
                            if ((xposdiff < sizelablxpos) & (xposdiff > 0) & (yposdiff < sizelablypos) & (yposdiff > 0)).any():
                                boolgood = False
                        if boolgood or numbtria > 99999:
                            break
                        numbtria += 1
                    
                    if boolgood:
                        print('Automatically positioned the annotations in %d trials, ll=%d...' % (numbtria, ll))
                        break
                if not boolgood:
                    print('Failed to automatically positioned the annotations...')
                
                #indxoutf = np.arange(numboutf)
                #indxsampassi = []
                #indxlablassi = []
                #for n in indxoutf:
                #    indxlabllive = np.setdiff1d(indxoutf, np.array(indxlablassi))
                #    indxsamplive = np.setdiff1d(indxoutf, np.array(indxsampassi))
                #    distlabl = np.sqrt((listxposlabl[indxlabllive, None] - xpossamp[None, indxsamplive])**2 + (listyposlabl[indxlabllive, None] - ypossamp[None, indxsamplive])**2)
                #    
                #    indxlablminm, indxsampminm = np.unravel_index(distlabl.argmin(), distlabl.shape)
                #    indxlablassi.append(indxlabllive[indxlablminm])
                #    indxsampassi.append(indxsamplive[indxsampminm])
                #listxposlabl = listxposlabl[indxlablassi]
                #listyposlabl = listyposlabl[indxlablassi]
                
                for n in range(numboutf):
                    axis.annotate(listlablsamppair[n], xy=(xpossamp[n], ypossamp[n]), \
                                  xytext=(listxposlabl[n], listyposlabl[n]),
                                  textcoords=axis.transAxes, \
                                  xycoords=axis.transAxes, \
                                  ha='center', va='center', \
                                  color=listcolrpopl[u], \
                                  bbox=dict(boxstyle='round,pad=0.2', fc='yellow', ec='darkblue', alpha=1.), \
                                  arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.1', color='red'), \
                                 )

        else:
            hist = np.histogram2d(listpara[u][:, l], listpara[u][:, k], bins=binstemp)[0]
            objtaxispcol = axis.pcolor(bins[l], bins[k], hist.T, cmap=listcolrpopltdim[u], label=labl)
    
    if boolcbar and (listtypeplottdim == 'hist').any():
        cbar = plt.colorbar(objtaxispcol)
        
    if boolmakelegd and indxpopl.size > 1:
        axis.legend(framealpha=1.)
    
    if truepara is not None and truepara[l] is not None and not np.isnan(truepara[l]) and truepara[k] is not None and not np.isnan(truepara[k]):
        axis.scatter(truepara[l], truepara[k], color='g', marker='x', s=500)
    
    # draw the provided reference values
    if listparadraw is not None:
        for m in indxdraw:
            axis.scatter(listparadraw[m][l], listparadraw[m][k], color='r', marker='x', s=350)
    
    if listvectplot is not None:
        for vectplot in listvectplot:
            axis.arrow(vectplot[0], vectplot[1], vectplot[2], vectplot[3])
    
    if boolsqua:
        axis.set_aspect('equal')
    
    if listscalpara[l] == 'logt':
        setp_axislogt(axis, limt[:, l], 'x', listmantlabl)
    
    if listscalpara[k] == 'logt':
        setp_axislogt(axis, limt[:, k], 'y', listmantlabl)
    

def retr_listvalutickmajr(minmlogt, maxmlogt):
    
    if np.ceil(maxmlogt) - np.floor(minmlogt) < 8:
        diff = 1
    else:
        diff = 2

    listvalutickmajr = 10**(np.arange(np.floor(minmlogt), np.ceil(maxmlogt) + 1, diff))
    
    return listvalutickmajr


def retr_valulabltick( \
                      # minimum value of the axis
                      minm, \
                      # maximum values of the axis
                      maxm, \
                      # a string indicating the scaling
                      ## 'powr': power-law
                      ## 'logt': logarithmic
                      ## 'asnh': arcsinh
                      scal, \
                      # optional list of mantissa to be used
                      listmantlabl=None, \
                     ):
    
    if scal == 'logt' or scal == 'powr':
        minm = np.log10(minm)
        maxm = np.log10(maxm)
    elif scal == 'asnh':
        minm = np.arcsinh(minm)
        maxm = np.arcsinh(maxm)
    elif scal == 'self':
        raise Exception('This function is not supposed to be used for linear scaling.')
    else:
        raise Exception('Scaling is not recognized.')
    
    # determine major ticks and labels
    ## list of values for the major ticks
    listvalutickmajr = retr_listvalutickmajr(minm, maxm)
    ## list of labels for the major ticks
    listlabltickmajr = [retr_lablmexp(listvalutickmajr[a]) for a in range(len(listvalutickmajr))]
    
    # determine minor ticks and labels
    ## minor tick mantissa
    if listmantlabl is None:
        if len(listvalutickmajr) > 2:
            listmantlabl = []
        elif len(listvalutickmajr) == 2:
            listmantlabl = [3.]
        else:
            listmantlabl = [2., 5.]
    listmantminr = np.arange(2., 10.)
    listexpominr = np.arange(np.floor(minm), np.ceil(maxm) + 1)
    
    # list of values for the minor ticks
    listvalutickminr = []
    listvalutickminrlabl = []
    indxlablmant = []
    indxlablexpo = []
    for kk in range(len(listmantminr)):
        for mm in range(len(listexpominr)):
            if scal == 'logt' or scal == 'powr':
                valu = listmantminr[kk] * 10**(listexpominr[mm])
            if scal == 'asnh':
                valu = listmantminr[kk] * np.sinh(listexpominr[mm])
            listvalutickminr.append(valu)
            if listmantminr[kk] in listmantlabl:
                listvalutickminrlabl.append(valu)
                indxlablmant.append(kk)
                indxlablexpo.append(mm)
    indxlabl = np.array(indxlablexpo, dtype=int) * len(listmantminr) + np.array(indxlablmant, dtype=int)
    listvalutickminr = np.array(listvalutickminr)
    listvalutickminr = np.sort(listvalutickminr)
    
    # list of labels for the minor ticks
    listlabltickminr = np.empty_like(listvalutickminr, dtype=object)
    listlabltickminr[:] = ''
    for a in range(len(listvalutickminrlabl)):
        listlabltickminr[indxlabl[a]] = retr_lablmexp(listvalutickminrlabl[a])
    
    return listvalutickmajr, listlabltickmajr, listvalutickminr, listlabltickminr


def setp_axislogt(axis, limt, typeaxis, listmantlabl):
    
    minm = limt[0]
    maxm = limt[1]
    
    listvalutickmajr, listlabltickmajr, listvalutickminr, listlabltickminr = retr_valulabltick(minm, maxm, 'logt', listmantlabl=listmantlabl)

    if typeaxis == 'x':
        axis.set_xscale('log', base=10)
        axis.set_xticks(listvalutickmajr)
        axis.set_xticklabels(listlabltickmajr)
        axis.set_xticks(listvalutickminr, minor=True)
        axis.set_xticklabels(listlabltickminr, minor=True)
    else:
        axis.set_yscale('log', base=10)
        axis.set_yticks(listvalutickmajr)
        axis.set_yticklabels(listlabltickmajr)
        axis.set_yticks(listvalutickminr, minor=True)
        axis.set_yticklabels(listlabltickminr, minor=True)


def plot_grid_histodim(listmantlabl, listpara, k, listlablparatotl, indxpopl, listlablpopl, bins, \
                                    listcolrpopl, listparadraw, lablnumbsamp, lablsampgene, boolinte, \
                                    boolmakelegd, listscalpara, factulimyaxihist, titl, plotsize, limt, limtrims, boolcumu=False, path=None):
    
    figr, axis = plt.subplots(figsize=(plotsize, plotsize))
    for u in indxpopl:
        
        if listpara[u][:, k].size == 0:
            continue

        if indxpopl.size > 1:
            labl = listlablpopl[u]
            alph = 0.7
        else:
            labl = None
            alph = 1.
        axis.hist(listpara[u][:, k], bins=bins[k], label=labl, cumulative=boolcumu, \
                                                               edgecolor=matplotlib.colors.to_rgba(listcolrpopl[u], 1.0), lw=2, ls='-', \
                                                               facecolor=matplotlib.colors.to_rgba(listcolrpopl[u], 0.2))
    if listparadraw is not None:
        for m in indxdraw:
            axis.axvline(listparadraw[m][k], color='orange', lw=3)
    axis.set_xlabel(listlablparatotl[k])
    if lablnumbsamp is not None:
        axis.set_ylabel(lablnumbsamp)
    elif lablsampgene is not None:
        #axis.set_ylabel(r'$N_{\rm{%s}}$' % lablsampgene)
        axis.set_ylabel('Number of %s' % lablsampgene)
    
    if boolinte[k]:
        axis.xaxis.get_major_locator().set_params(integer=True)

    if boolmakelegd and indxpopl.size > 1:
        axis.legend(framealpha=1.)
    if listscalpara[k] == 'logt':
        setp_axislogt(axis, limtrims[:, k], 'x', listmantlabl)
    limtyaxi = axis.get_ylim()
    
    boolyaxilogt = limtyaxi[1] - limtyaxi[0] > 30.

    # rescale the upper limit of the vertical axis
    if factulimyaxihist != 1.:
        limtyaxiprim = np.array(limtyaxi)
        if boolyaxilogt:
            # find the minimum value of the histogram for all populations
            llimyaxihist = 1e100
            for u in indxpopl:
                hist = np.histogram(listpara[u][:, k], bins=bins[k])[0]
                if np.amax(hist) > 0:
                    llimyaxihist = min(np.amin(hist[np.where(hist > 0)[0]]), llimyaxihist)
            llimyaxihist *= 0.5
        else:
            llimyaxihist = None
        limtyaxiprim = [llimyaxihist, limtyaxi[1] * factulimyaxihist]
        axis.set_ylim(limtyaxiprim)

    if boolyaxilogt:
        axis.set_yscale('log')
    
    if titl is not None:
        axis.set_title(titl)
    axis.set_xlim(limt[:, k]) 
    plt.tight_layout()
    if path is not None:
        print('Writing to %s...' % path)
        plt.savefig(path)
        plt.close()
    else:
        plt.show()


def plot_grid(
              # two dimensional numpy array of samples, where the first dimension is the sample and second dimension is the parameter
              listpara, \
              
              # a list with length equal to the number of parameters, 
              # Each element of the list should itself be list of two strings, where
              # the first string is the label for the parameter and the second string is the unit
              listlablpara, \
              
              # an optional string indicating the path of the folder in which to write the plot
              pathbase=None, \
              
              # an optional base string to include in the file name
              strgplot=None, \
              
              # the limits for the parameters
              limt=None, \

              # list of scalings for the parameters
              listscalpara=None, \
              
              # size of each subplot
              plotsize=3.5, \
              
              # type of the file for plots
              typefileplot='png', \
              
              # Boolean flag to generate the lower-triangle plot
              boolplottria=None, \

              # Boolean flag to generate individual histograms
              boolplothistodim=None, \
              
              # Boolean flag to generate individual pair-wise plots
              boolplotpair=False, \
              
              # list of base file names for the individual histograms
              listnamepara=None, \
                  
              # Boolean flag to make the two-dimensional plots square
              boolsqua=False, \
              
              # label to be used to denote the number of samples (takes priority over lablsampgene)
              lablnumbsamp=None, \

              # label to be used as subscript to denote the number of samples
              lablsampgene=None, \

              # list of markers for populations
              listmrkrpopl=None, \

              # list of colors for populations
              listcolrpopl=None, \

              # list of colors for populations to be used in two histograms
              listcolrpopltdim=None, \

              # list of vectors to overplot
              listvectplot=None, \
              
              # list of labels for each sample
              listlablsamp=None, \
              
              # list of feature names for which a cumulative histogram will be made
              listnamefeatcumu=None, \

              # list of pairs of feature names to be skipped
              listnamefeatskip=None, \
              
              # list of labels for populations
              listlablpopl=None, \
              
              # Boolean flag to include a legend
              boolmakelegd=True, \
             
              # Boolean flag to indicate that the populations are mutually-exclusive
              boolpoplexcl=False, \
             
              # factor by which to multiply the upper limit of the y-axis
              factulimyaxihist=None, \

              # title for the plots
              titl=None, \
              
              # list of tick mantices (other than 1) to show in the label when the axis is log-streched
              listmantlabl=None, \
              
              # number of bins
              numbbinsplot=40, \

              # Boolean flag to overplot quantiles
              boolplotquan=False, \
              
              # list of parameters to overplot
              listparadraw=None, \
              
              # true parameters to highlight
              truepara=None, \
              
              # a list of pairs of feature names, which enforces the given order (first item y-axis, second item y-axis)
              listnameordrpair=None, \

              # verbosity level
              typeverb=1, \
              
              # type of the two-dimensional plots
              ## 'scat': all populations are scatter plots.
              ## 'hist': all populations are histogram plots.
              ## 'best': the largest population is histogram if there are too many samples, scatter otherwise. All other populations are scatter.
              typeplottdim='best', \
              
              # Boolean flag to diagnose the code
              booldiag=True, \
              
             ):
    
    '''
    Make a corner plot of a multivariate distribution.
    '''
    
    if typeverb > 1:
        print('tdpy.plot_grid() initialized...')
    
    if pathbase is not None:
        if not os.path.exists(pathbase):
            os.system('mkdir -p %s' % pathbase)
    else:
        path = None

    # check whether there is a single population or multiple populations
    if isinstance(listpara, list):
        boolmpop = True
    else:
        boolmpop = False
        listpara = [listpara]
        
    # preclude quantile lines if there are multiple populations
    if boolplotquan and boolmpop:
        raise Exception('')
    
    if listpara[0].ndim == 1:
        listpara[0] = listpara[0][:, None]
        listlablpara = [listlablpara]

    # temp: number of parameters should be able to be different for different populations
    numbpara = listpara[0].shape[1]
    indxpara = np.arange(numbpara)
    
    if boolplothistodim is None:
        if boolmpop:
            boolplothistodim = False
        else:
            if numbpara == 1:
                boolplothistodim = True
            else:
                boolplothistodim = False
    
    if boolplottria is None:
        boolplottria = numbpara > 1
    
    if lablnumbsamp is None:
        lablnumbsamp = 'Number of samples'

    # check inputs
    if (boolplotpair or boolplothistodim) and listnamepara is None and pathbase is not None:
        raise Exception('You should define listnamepara for individual and pairwise plots.')
    
    listlablparatotl = retr_labltotl(listlablpara)
    
    if listmrkrpopl is None:
        listmrkrpopl = np.array(['o', 'x', '+', 'D', '^', '*', '<', '>', 's', 'p'])
    
    if listcolrpopl is None:
        listcolrpopl = np.array(['g', 'r', 'b', 'purple', 'orange', 'pink', 'magenta', 'olive', 'cyan', 'teal'])
    
    if listcolrpopltdim is None:
        listcolrpopltdim = np.array(['Greens', 'Blues', 'Purples', 'Oranges'])
    
    numbpopl = len(listpara)
    
    if numbpopl > 10:
        print('')
        print('')
        print('')
        print('Number of populations must be less than or equal to 10.')
        raise Exception('')

    if numbpopl > 1 and listlablpopl is None:
        print('listlablpopl should be defined when there are more than one populations.')
        raise Exception('')

    indxpopl = np.arange(numbpopl)
    listcolrpopl = listcolrpopl[:numbpopl]
    for u in indxpopl:
        for k in indxpara:
            if not np.isfinite(listpara[u][:, k]).all():
                print('tdpy.plot_grid(): Not all samples are finite for parameter %d!' % k)
    
    if listscalpara is None:
        listscalpara = ['self'] * numbpara

    if len(listscalpara) != len(listlablpara):
        print('listscalpara')
        print(listscalpara)
        print('listlablpara')
        print(listlablpara)
        raise Exception('')
    if len(listscalpara) != numbpara:
        print('listscalpara')
        print(listscalpara)
        print('len(listscalpara)')
        print(len(listscalpara))
        for u in indxpopl:
            print('listpara[u]')
            summgene(listpara[u])
        print('numbpara')
        print(numbpara)
        raise Exception('len(listscalpara) != numbpara')
    
    # number of samples in each population
    numbsamp = np.empty(numbpopl, dtype=int)
    for u in indxpopl:
        numbsamp[u] = listpara[u][:, 0].size

    # determine the type of 2-dimensional plots
    if typeplottdim == 'scat' or numbpopl > 1:
        listtypeplottdim = np.array(['scat' for u in indxpopl])
    elif typeplottdim == 'hist':
        listtypeplottdim = np.array(['hist' for u in indxpopl])
    elif typeplottdim == 'best':
        listtypeplottdim = np.array(listtypeplottdim)
        indxpoplmaxm = np.argmax(numbsamp)
        if numbsamp[indxpoplmaxm] >= 1e3:
            listtypeplottdim[indxpoplmaxm] = 'hist'
    
    # sort the populations in decreasing order of size
    indxpopl = np.argsort(numbsamp)[::-1]

    listindxgood = [[[] for k in indxpara] for u in indxpopl]
    if limt is None:
        limt = np.empty((2, numbpara))
        limt[0, :] = 1e100
        limt[1, :] = -1e100
        
        boolinte = [[] for k in indxpara]
        for k in indxpara:
            for u in indxpopl:
                if listscalpara[k] == 'logt':
                    boolsampfini = np.isfinite(listpara[u][:, k])
                    boolsampposi = listpara[u][:, k] > 0
                    listindxgood[u][k] = np.where(boolsampposi & boolsampfini)[0]
                    if (listpara[u][boolsampfini, k] <= 0).any():
                        print('')
                        print('')
                        print('')
                        print('')
                        print('tdpy.plot_grid(): Parameter %d (%s) has a log scaling but also nonpositive elements!' % (k, listlablpara[k]))
                        print('listpara[u][:, k]')
                        summgene(listpara[u][:, k])
                        print('listindxgood[u][k]')
                        summgene(listindxgood[u][k])
                        print('np.where(boolsampfini)[0]')
                        summgene(np.where(boolsampfini)[0])
                        print('')
                        print('')
                        print('')
                        print('')
                        #raise Exception('')
                else:
                    listindxgood[u][k] = np.arange(listpara[u][:, k].size)
            
                if listindxgood[u][k].size > 0:
                    limt[0, k] = min(limt[0, k], np.nanmin(listpara[u][listindxgood[u][k], k], 0))
                    limt[1, k] = max(limt[1, k], np.nanmax(listpara[u][listindxgood[u][k], k], 0))
            
            if booldiag:
                if limt[0, k] == limt[1, k]:
                    print('')
                    print('')
                    print('')
                    print('')
                    print('Warning! The lower and upper limits for parameters %s are the same: %g.' % (listnamepara[k], limt[0, k]))
                    print('listnamepara[k]')
                    print(listnamepara[k])
                    print('listpara[u][:, k]')
                    summgene(listpara[u][:, k])
                    print('listpara[u][listindxgood[u][k], k]')
                    summgene(listpara[u][listindxgood[u][k], k])
                    print('')
                    #raise Exception('')

            # list of Booleans for each parameter indicating whether it is a list of integers
            boolinte[k] = True
            for u in indxpopl:
                if ((listpara[u][:, k] - listpara[u][:, k].astype(int)) != 0).any():
                    boolinte[k] = False
        
        # sanity checks
        if booldiag:
            for k in indxpara:
                #for u in indxpopl:
                #    if listpara[u][:, k].size == 0:
                #        print('')
                #        print('')
                #        print('')
                #        print('')
                #        print('')
                #        print('')
                #        print('')
                #        print('tdpy.plot_grid(): size is 0 for parameter %d!' % k)
                #        print('k')
                #        print(k)
                #        print('listindxgood[u][k]')
                #        summgene(listindxgood[u][k])
                #        print('listpara[u][:, k]')
                #        summgene(listpara[u][:, k])
                #        print('numbbinsplot')
                #        print(numbbinsplot)
                #        print('listnamepara[k]')
                #        print(listnamepara[k])
                #        print('listlablpara[k]')
                #        print(listlablpara[k])
                #        print('listscalpara[k]')
                #        print(listscalpara[k])
                #        print('limt[:, k]')
                #        print(limt[:, k])
                #        raise Exception('')
                
                for u in indxpopl:
                    if listpara[u][:, k].size > 0:
                        minmtemp = np.amin(np.abs(listpara[u][:, k]))
                        if not np.isfinite(limt[:, k]).all() or minmtemp < 1e-100 and minmtemp > 0.:
                            print('')
                            print('')
                            print('')
                            print('')
                            print('')
                            print('')
                            print('')
                            print('tdpy.plot_grid(): limit for parameter %d is infinite!' % k)
                            print('k')
                            print(k)
                            print('listindxgood[u][k]')
                            summgene(listindxgood[u][k])
                            print('listpara[u][:, k]')
                            summgene(listpara[u][:, k])
                            print('numbbinsplot')
                            print(numbbinsplot)
                            print('listlablpara[k]')
                            print(listlablpara[k])
                            print('listscalpara[k]')
                            print(listscalpara[k])
                            print('limt[:, k]')
                            print(limt[:, k])
                            raise Exception('')

        # sanity checks
        if booldiag:
            for k in indxpara:
                if not np.isfinite(limt[:, k]).all():
                    print('')
                    print('')
                    print('')
                    print('')
                    print('')
                    print('')
                    print('')
                    print('tdpy.plot_grid(): limit for parameter %d is infinite!' % k)
                    print('k')
                    print(k)
                    print('numbpopl')
                    print(numbpopl)
                    for u in indxpopl:
                        print('listindxgood[u][k]')
                        summgene(listindxgood[u][k])
                        print('listpara[u][:, k]')
                        summgene(listpara[u][:, k])
                    print('numbbinsplot')
                    print(numbbinsplot)
                    print('listlablpara[k]')
                    print(listlablpara[k])
                    print('listscalpara[k]')
                    print(listscalpara[k])
                    print('limt[:, k]')
                    print(limt[:, k])
                    raise Exception('')

        if (listtypeplottdim == 'scat').any():
            # update limits to leave white space in the rims which is good for scatter plots
            for k in indxpara:
                
                if limt[0, k] == 1e100 and limt[1, k] == -1e100:
                    continue

                if listscalpara[k] == 'self':
                    if boolinte[k]:
                        delt = 0.5
                    else:
                        delt = 0.05 * (limt[1, k] - limt[0, k])
                    limt[0, k] -= delt
                    limt[1, k] += delt
                if listscalpara[k] == 'logt':
                    fact = np.exp(0.05 * np.log(limt[1, k] / limt[0, k]))
                    limt[0, k] /= fact
                    limt[1, k] *= fact
        
        ## if the limits are finite
        #for k in indxpara:
        #    if not np.isfinite(limt[:, k]).all():
        #        print('The limits for parameter %s were not finite. Fixing them to [0, 1]...' % listnamepara[k])
        #        limt[0, k] = 0.
        #        limt[1, k] = 1.
            
        # limits that do not leave any room for white space which is good for histograms
        limtrims = np.copy(limt)

        # sanity checks
        if booldiag:
            for k in indxpara:
                if not np.isfinite(limt[:, k]).all():
                    print('')
                    print('')
                    print('')
                    print('')
                    print('')
                    print('')
                    print('')
                    print('tdpy.plot_grid(): limit for parameter %d is infinite!' % k)
                    print('k')
                    print(k)
                    print('numbpopl')
                    print(numbpopl)
                    for u in indxpopl:
                        print('listindxgood[u][k]')
                        summgene(listindxgood[u][k])
                        print('listpara[u][:, k]')
                        summgene(listpara[u][:, k])
                    print('numbbinsplot')
                    print(numbbinsplot)
                    print('listlablpara[k]')
                    print(listlablpara[k])
                    print('listscalpara[k]')
                    print(listscalpara[k])
                    print('limt[:, k]')
                    print(limt[:, k])
                    raise Exception('')

        for k in indxpara:
            if listlablpopl is not None:
                for u in indxpopl:
                    if not np.isfinite(listpara[u][:, k]).all():
                        print('Warning! Population %d (%s), parameter %d (%s) has nonfinite samples.' % (u, listlablpopl[u], k, listlablpara[k][0]))
                        summgene(listpara[u][:, k])
    
    for k in indxpara:
        if limt[0, k] == limt[1, k]:
            print('')
            print('')
            print('')
            print('')
            print('tdpy.plot_grid(): WARNING! Lower and upper limits are the same for the following parameter.')
            print('k')
            print(k)
            print('limt[:, k]')
            print(limt[:, k])
            print('listlablpara[k]')
            print(listlablpara[k])
            for u in indxpopl:
                print('u')
                print(u)
                print('listpara[u]')
                summgene(listpara[u])
            print('')
    
    if truepara is not None:
        for k in indxpara:
            if truepara[k] is not None:
                if truepara[k] < limt[0, k]:
                    limt[0, k] = truepara[k] - 0.1 * (limt[1, k] - truepara[k]) 
                if truepara[k] > limt[1, k]:
                    limt[1, k] = truepara[k] + 0.1 * (truepara[k] - limt[0, k])
    
    if listparadraw is not None:
        numbdraw = len(listparadraw)
        indxdraw = np.arange(numbdraw)

    # list of Booleans indicating whether a parameter is good to plot
    boolparagood = np.ones(numbpara, dtype=bool)
    for k in indxpara:
        booltemp = False
        for u in indxpopl:
            if np.isfinite(listpara[u][:, k]).any():
                booltemp = True
        boolparagood[k] = booltemp
    boolparagood[np.where(limt[0, :] == limt[1, :])] = False
        
    if boolplottria or boolplothistodim or boolplotpair and (listtypeplottdim == 'hist').any():
        bins = [[] for k in indxpara]
        for k in indxpara:
            
            if not boolparagood[k]:
                continue
            
            if boolinte[k]:
                bins[k] = np.linspace(limt[0, k], limt[1, k], int(limt[1, k] - limt[0, k] + 1))
            
            else:
                if listscalpara[k] == 'self' or listscalpara[k] == 'gaus':
                    bins[k] = icdf_self(np.linspace(0., 1., numbbinsplot + 1), limt[0, k], limt[1, k])
                elif listscalpara[k] == 'logt':
                    bins[k] = icdf_logt(np.linspace(0., 1., numbbinsplot + 1), limt[0, k], limt[1, k])
                elif listscalpara[k] == 'atan':
                    bins[k] = icdf_atan(np.linspace(0., 1., numbbinsplot + 1), limt[0, k], limt[1, k])
                else:
                    raise Exception('Unrecognized scaling: %s' % listscalpara[k])
            
            if booldiag:
                boolgood = True
                if bins[k].size > 1e6:
                    boolgood = False
                if not np.isfinite(bins[k]).all() or bins[k][0] > bins[k][-1]:
                    boolgood = False
                    
                if not boolgood:
                    print('')
                    print('')
                    print('')
                    print('k')
                    print(k)
                    print('numbbinsplot')
                    print(numbbinsplot)
                    print('listnamepara[k]')
                    print(listnamepara[k])
                    print('listlablpara[k]')
                    print(listlablpara[k])
                    print('listscalpara[k]')
                    print(listscalpara[k])
                    print('limt[:, k]')
                    print(limt[:, k])
                    print('bins[k]')
                    summgene(bins[k])
                    for u in indxpopl:
                        print('listpara[u][:, k]')
                        summgene(listpara[u][:, k])
                    raise Exception('')
            if np.amin(bins[k]) == 0 and np.amax(bins[k]) == 0:
                print('Lower and upper limits of the bins are the same for %s. Grid plot failed.' % listlablpara[k][0])
                print('k')
                print(k)
                print('bins[k]')
                print(bins[k])
                print('limt[:, k]')
                print(limt[:, k])
                print('listscalpara[k]')
                print(listscalpara[k])
                print('listlablpara[k]')
                print(listlablpara[k])
                #return
    else:
        bins = None

    if boolplothistodim:
        
        if numbpopl == 1:
            factulimyaxihist = 1.
        else:
            factulimyaxihist = 1.
            # temp
            #factulimyaxihist = 200.

        # one dimensional histograms
        for k in indxpara:
            if not boolparagood[k]:
                continue
            
            if pathbase is not None:
                path = pathbase + 'hist_%s_%s.%s' % (listnamepara[k], strgplot, typefileplot)
            if not (pathbase is not None and os.path.exists(path)):
                plot_grid_histodim(listmantlabl, listpara, k, listlablparatotl, indxpopl, listlablpopl, \
                                                bins, listcolrpopl, listparadraw, lablnumbsamp, lablsampgene, boolinte, \
                                                boolmakelegd, listscalpara, factulimyaxihist, titl, plotsize, limt, limtrims, boolcumu=False, path=path)
            if listnamefeatcumu is not None:
                if listnamepara[k] in listnamefeatcumu:
                    if pathbase is not None:
                        path = pathbase + 'histcumu_%s_%s.%s' % (listnamepara[k], strgplot, typefileplot)
                    if not (pathbase is not None and os.path.exists(path)):
                        plot_grid_histodim(listmantlabl, listpara, k, listlablparatotl, indxpopl, listlablpopl, bins, \
                                                    listcolrpopl, listparadraw, lablnumbsamp, lablsampgene, boolinte, \
                                                    boolmakelegd, listscalpara, factulimyaxihist, titl, plotsize, limt, limtrims, boolcumu=True, path=path)
                    
    listsizepopl = []
    for u in indxpopl:
        listsizepopl.append(listpara[u][:, 0].size)
    
    if numbpopl > 1:
        print('Number of samples in the populations:')
        for u in indxpopl:
            if listlablpopl is not None:
                print(listlablpopl[u])
            print(listsizepopl[u])
    
    # make pie-chart of the populations if populations are mutually-exclusive
    if boolpoplexcl:
    
        #def make_autopct(listsizepopl):
        #    
        #    def my_autopct(pct):
        #        total = sum(listsizepopl)
        #        val = int(round(pct*total/100.0))
        #        return '{p:.1f}% hey  ({v:d})'.format(p=pct,v=val)
        #    
        #    return my_autopct
        
        def make_autopct(listsizepopl):
            def my_autopct(pct):
                total = sum(listsizepopl)
                val = int(round(pct*total/100.0))
                return '%.3g%%' % val
                #return '{p:.2f}%%  ({v:d})'.format(p=pct,v=val)
            return my_autopct

        figr, axis = plt.subplots(figsize=(plotsize, plotsize))
        objt = axis.pie(listsizepopl, labels=listlablpopl, \
                                                                #autopct=make_autopct(listsizepopl), \
                                                                colors=listcolrpopl, \
                                                                autopct='%.3g%%', \
                                                                )
        for kk in range(len(objt[0])):
            objt[0][kk].set_alpha(0.5)

        axis.axis('equal')  # Equal aspect ratio ensures that pie is drawn as a circle.
        if titl is not None:
            axis.set_title(titl)
        path = pathbase + 'pies_%s.%s' % (strgplot, typefileplot)
        print('Writing to %s...' % path)
        
        figr.savefig(path, bbox_inches='tight')
        plt.close(figr)
    
    if boolplotpair:
        if listlablsamp is None or numbpopl > 1:
            liststrgtext = ['']
        else:
            liststrgtext = ['', '_anno']
        for strgtext in liststrgtext:
            if strgtext == '':
                listlablsamptemp = None
            else:
                listlablsamptemp = listlablsamp
            
            for k in indxpara:
                for l in indxpara:
                    if not boolparagood[k] or not boolparagood[l]:
                        continue
                    if k <= l:
                        continue
                    
                    # skip the feature pair if specified by the user
                    if listnamefeatskip is not None:
                        boolskip = False
                        for gg in range(len(listnamefeatskip)):
                            if listnamepara[k] in listnamefeatskip[gg] and listnamepara[l] in listnamefeatskip[gg]:
                                boolskip = True
                        if boolskip:
                            continue
                    
                    if not (listnamepara[l] == 'rascstar' and listnamepara[k] == 'declstar'):
                        numbiter = 1
                    else:
                        numbiter = 2
                    
                    for e in range(numbiter):
                        
                        if e == 0:
                            projection = None
                            strgiter = ''
                        if e == 1:
                            projection = 'aitoff'
                            strgiter = '_aito'
                        path = pathbase + 'pmar_%s_%s_%s_%s%s%s.%s' % (typeplottdim, listnamepara[k], \
                                                    listnamepara[l], strgplot, strgtext, strgiter, typefileplot)
                        
                        if not os.path.exists(path):
                        
                            figr = plt.figure(figsize=(plotsize, plotsize))
                            axis = figr.add_subplot(111, projection=projection)
                            
                            plot_grid_pair(k, l, axis, limt, listmantlabl, listpara, truepara, listparadraw, boolplotquan, listlablpara, \
                                                             listscalpara, boolsqua, listvectplot, listtypeplottdim, indxpopl, listcolrpopl, \
                                                             listmrkrpopl, listcolrpopltdim, listlablpopl, boolmakelegd, bins=bins, \
                                                             listlablsamp=listlablsamptemp, boolcbar=True)
                            
                            if e == 0:
                                axis.set_xlim(limt[:, l])
                                axis.set_ylim(limt[:, k])
                
                            axis.set_xlabel(listlablparatotl[l])
                            axis.set_ylabel(listlablparatotl[k])
                            
                            if titl is not None:
                                axis.set_title(titl)
                            
                            print('Writing to %s...' % path)
                            figr.savefig(path, bbox_inches='tight')
                            plt.close(figr)
                        
    if boolplottria:
        figr, axgr = plt.subplots(numbpara, numbpara, figsize=(0.6*plotsize*numbpara, 0.6*plotsize*numbpara))
        if numbpara == 1:
            axgr = [[axgr]]
        for k, axrw in enumerate(axgr):
            for l, axis in enumerate(axrw):
                if not boolparagood[k] or not boolparagood[l]:
                    continue
                if k < l:
                    axis.axis('off')
                    continue

                if k == l:
                    plot_grid_diag(k, axis, listpara, truepara, listparadraw, boolplotquan, listlablpara, listtypeplottdim, indxpopl, \
                                                                              listcolrpopl, listmrkrpopl, listlablpopl, boolmakelegd, listsizepopl, bins=bins)
                else:
                    plot_grid_pair(k, l, axis, limt, listmantlabl, listpara, truepara, listparadraw, boolplotquan, listlablpara, \
                                                     listscalpara, boolsqua, listvectplot, listtypeplottdim, indxpopl, listcolrpopl, \
                                                     listmrkrpopl, listcolrpopltdim, listlablpopl, boolmakelegd, bins=bins, \
                                                     boolcbar=False)
                    
                    axis.set_xlim(limt[:, l])
                    axis.set_ylim(limt[:, k])
                
                    
                if k == numbpara - 1:
                    axis.set_xlabel(listlablparatotl[l])
                else:
                    axis.set_xticklabels([])
                
                if (l == 0 and k != 0):
                    axis.set_ylabel(listlablparatotl[k])
                else:
                    if k != 0:
                        axis.set_yticklabels([])
        
        figr.tight_layout()
        plt.subplots_adjust(wspace=0.05, hspace=0.05)
        
        if pathbase is not None:
            path = pathbase + 'pmar_%s_%s.%s' % (typeplottdim, strgplot, typefileplot)
            print('Writing to %s...' % path)
            figr.savefig(path, dpi=300)
            plt.close(figr)
        else:
            plt.show()
        

