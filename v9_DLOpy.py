import builtins
#!/usr/bin/env python
# coding: utf-8

import obspy
from obspy import read as r
from obspy.core import UTCDateTime
from obspy.clients.fdsn import Client

from numpy import *
import numpy as np

from scipy import signal
import scipy.signal as sig
from scipy.stats import circmean as cmean
from scipy.stats import circstd as cstd
from scipy.stats import hmean as hm

from geographiclib.geodesic import *
import subprocess
import matplotlib.pyplot as plt
from matplotlib.ticker import FuncFormatter, MultipleLocator

import os
import sys

import yaml
import argparse

import csv
import re

SELF_CHECK= True
BOOTNUM_DEFAULT = 10000

# temporary
import time, builtins
from contextlib import contextmanager

_profile_totals = {}

@contextmanager
def timed(label: str):
    t0=time.perf_counter()
    try:
        yield
    finally:
        dt=time.perf_counter() -t0
        _profile_totals[label]=_profile_totals.get(label, 0.0)+float(dt)

def print_profile(header="Station"):
    tot=float(builtins.sum((float(v) for v in _profile_totals.values()), 0.0))
    rows=sorted(_profile_totals.items(),key=lambda x: x[1],reverse=True)
    print(f"\n[{header}]")
    print(f"Total Time: {tot}")
    for k,v in rows:
        print(f"{k:<20} {float(v):8.3f} s")


#############################

def getmsd1(t1):
    dir1='/Volumes/AKD1/AKD-data/'
    t1=UTCDateTime(t1)

    d=r(dir1+'AKD.CH?.%i.%03d.*.msd'%(t1.year,t1.julday))

	# include next day of data if close to edge
    if t1.hour>=20:
        d2=r(dir1+'AKD.CH?.%i.%03d.*.msd'%(t1.year,t1.julday+1))
#       d3=d+a2 
        d3=d+d2
        d=d3.merge()

	# rename channels
    C1=d.select(channel='CH1')
    C1[0].stats.channel='BH1'
    
    C2=d.select(channel='CH2')
    C2[0].stats.channel='BH2'
    
    C3=d.select(channel='CH3')
    C3[0].stats.channel='BHZ'
    
    st1=C1+C2+C3

    # set sampling rate
    sps=20
    # define four hours of data
    t2=60*60*4

    return st1.slice(starttime=t1,endtime=t1+t2)

#############################
# Functions for driver program
#############################

# First pass at detrending
def detr(T):
    T.detrend()
    T.detrend('linear')
    return T

# Organize channels by coordinate system
def org(T,nameconv):
    # Need to add to this in case other type (eg LHZ) is used
    if nameconv==1:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??N")
        bh2=T.select(channel="??E")    
    if nameconv==2:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??1")
        bh2=T.select(channel="??2")
    if nameconv==3:
        bhz=T.select(channel="??Z")
        bh1=T.select(channel="??2")
        bh2=T.select(channel="??1")
    
    s=bh1+bh2+bhz

    # decimate just to make things easier
    # for now assume all channels are same sampling rate
    while s[0].stats.sampling_rate >= 10.0 and s[0].stats.sampling_rate / 10.0 == int(s[0].stats.sampling_rate / 10.0):
        s.decimate(10)

    while s[0].stats.sampling_rate / 5.0 == int(s[0].stats.sampling_rate / 5.0):
        s.decimate(5)

    return s
  
# used by catclean
def close(x1,x2):
    if abs(x1-x2)<0.8:
        return True
    else:
        return False

# look for repeat events
def catclean(stime,lat,lon,mag):
    rep=array([],dtype=int)
    for k in arange((len(stime))):
        for kk in arange((len(stime))):
            if stime[kk]-stime[k]<60*15 and close(lat[kk],lat[k]) and close(lon[kk],lon[k]) and mag[kk]<mag[k]-0.3:
                rep=append(rep,kk)
    return rep

#############################
# PHASE VELOCITY CALCULATIONS
#############################

def getf(freq,A):
    for i in arange((len(A))):
        if A[i,0]==freq:
            v=A[i,1]
    return v            

###
# find nearest value in an array
def nv(x,v):
    # x is array
    # v is value
    idx = (abs(x-v)).argmin()
    return x[idx]

# Overall function to get path-averaged group velocity
def pathvels(lat1,lon1,lat2,lon2,map10,map15,map20,map25,map30,map35,map40):

    Rearth=6371.25;
    circE=2*np.pi*Rearth;
    
    # Get distance and azimuth
    p=Geodesic.WGS84.Inverse(lat1,lon1,lat2,lon2)
    
    minor=float(p['s12']) / 1000.00
    major=circE-minor
            
    l1=Geodesic.WGS84.Line(lat1,lon1,p['azi1'])         
    l2=Geodesic.WGS84.Line(lat1,lon1,p['azi1']-180)     
    
    deg=111.17
    D1=zeros([1,2])
    for i in arange((361)):
        b=l1.Position(deg*1000*i)
        p2=array([b['lat2'],b['lon2']])
    
        if i==0:
            D1[0,0]=p2[0]
            D1[0,1]=p2[1]
        else:
            D1=vstack((D1,p2))
        
        bb=Geodesic.WGS84.Inverse(lat2,lon2,b['lat2'],b['lon2'])
        if bb['s12'] <= deg*1000.0 :
            break
    
    
    D2=zeros([1,2])
    for i in arange((361)):
        b=l2.Position(deg*1000*i)
        p2=array([b['lat2'],b['lon2']])
    
        if i==0:
            D2[0,0]=p2[0]
            D2[0,1]=p2[1]
        else:
            D2=vstack((D2,p2))
        
        bb=Geodesic.WGS84.Inverse(lat2,lon2,b['lat2'],b['lon2'])
        if bb['s12'] <= deg*1000.0 :
            break
    
    """
    We now have lat and lon points along the major and minor great circles.
    We calcaulte the group velocity of at each point, and then
    find the average velocities. 
    """
    
    for k in arange((len(D1))):
        if D1[k,1]<0:
            D1[k,1]+=360
    for k in arange((len(D2))):
        if D2[k,1]<0:
            D2[k,1]+=360

    def Ray(D):
        U1=zeros([len(D),7])
        for k in arange((len(D))):
            # do latitude first
            
            ## get correct precision
            ## designed to match results of Ma et al codes
            if abs(D[k,1]) < 10:
                D[k,1]=round(D[k,1],5)
            if abs(D[k,1]) >= 10 and abs(D[k,1]) < 100:
                D[k,1]=round(D[k,1],4)
            if abs(D[k,1]) > 100:
                D[k,1]=round(D[k,1],3)
            if abs(D[k,0]) < 10:
                D[k,0]=round(D[k,0],5)
            if abs(D[k,0]) >= 10:
                D[k,0]=round(D[k,0],4)
            #    
            # find right index
            q=where( map10[:,1] == nv(map10[:,1], (D[k,0])  ))[0]
            qq=where( map10[q,0]==nv(map10[q,0],  (D[k,1])  ))[0]
            idx=q[qq]
            
            # update path
            idx=int(np.atleast_1d(idx)[0])
            U1[k,0]= map10[ idx, 2]
            U1[k,1]= map15[ idx, 2]
            U1[k,2]= map20[ idx, 2]
            U1[k,3]= map25[ idx, 2]
            U1[k,4]= map30[ idx, 2]
            U1[k,5]= map35[ idx, 2]
            U1[k,6]= map40[ idx, 2]
        mhz=array([10,15,20,25,30,35,40])
        return array((mhz, hm(U1,axis=0))).T

    return Ray(D1),Ray(D2)

#####################
## ANGLE CALCULATION FUNCTIONS
#####################

# keep eqs above certain cc limit
# also keep which earthquakes were kept
# Different version of C1
def C1_2(phi,cc,clim):
    PHI=array([]); C=array([]); ix=array([])
    for i in arange((len(phi))):
        if cc[i]>clim:
            PHI=append(PHI,phi[i])
            C=append(C,cc[i])
            ix=append(ix,i)
    return PHI, C, ix

# C2 culling - keep values within 95% circ conf of circ mean
# also keep which earthquakes were kept

# Get unique events used in final calculation
def uniqueevents(phis,ccs,n,R1cc,R2cc):
    L=len(phis)/2
    ii=zeros((len(n)))
    for i in arange((len(n))):
        if n[i]<L:
            ii[i]=where(R1cc==ccs[int(n[i])])[0][0]
        else:
            ii[i]=where(R2cc==ccs[int(n[i])])[0][0]

    return unique(ii)


# Plotting function
def centerat(phi,m=0):
    phinew=copy(phi)
    if len(shape(phi))==1:
        for i in arange((len(phi))):
            if phi[i]>=m+180:
                phinew[i]-=360
            if phi[i]<=m-180:
                phinew[i]+=360
    else:
        for k in arange((shape(phi)[1])):
            for i in arange((shape(phi)[0])):
                if phi[i,k]>=m+180:
                    phinew[i,k]-=360
                if phi[i,k]<=m-180:
                    phinew[i,k]+=360
    return phinew


# Function to flatten result arrays
def flatten(X):
    return reshape(X,[X.shape[0]*X.shape[1],1])

# median absolute deviation
def mad(x):
    return median(abs(x-median(x)))

# remove outliars
def outlier1(Tphi,ix,lim=5.0):
    devs=abs(Tphi-median(Tphi))/mad(Tphi)
    ixs=where(devs<lim)
    return Tphi[ixs],ix[ixs]

# bootstrap mean
def boot1(phi, bootnum):
    import numpy as np
    phi = np.asarray(phi, dtype=float) % 360.0
    L = len(phi)
    if L == 0:
        return np.array([], dtype=float)
    rng = np.random.default_rng()
    out = np.empty(bootnum, dtype=float)

    chunk = 1000  # tune
    rad = np.deg2rad(phi)
    for start in range(0, bootnum, chunk):
        b = chunk if chunk < (bootnum - start) else (bootnum - start)
        idx = rng.integers(0, L, size=(b, L))
        rr  = rad[idx]  
        C   = np.cos(rr).mean(axis=1)
        S   = np.sin(rr).mean(axis=1)
        out[start:start+b] = (np.rad2deg(np.arctan2(S, C)) % 360.0)
    return out


# reorganize results
def resort(phi,col2):
    phi2=centerat(phi,m=cmean(phi,high=360))
    t=zeros((len(phi2),2))
    t[:,0]=phi2; t[:,1]=col2
    t = t[t[:,0].argsort()]
    return t[:,0], t[:,1]

# final Doran-Laske calculation
def fcalc1(phi, cc, lim, R1cc, R2cc):
    import numpy as np

    # Keep cc over limit
    Tphi, Tcc, ii = C1_2(phi, cc, lim)

    # -------- Defaults so every return path has 8 values --------
    ci_low_pct  = np.nan
    ci_high_pct = np.nan
    m           = np.array([])   # bootstrap samples (may be empty)
    Tphi2       = np.array([])   # filtered angles (may be empty)

    # Empty input after thresholding
    if len(Tphi) == 0:
        ue = np.array([0])   # no unique events
        n_used = 0
        return 0.0, 180.0, ue, n_used, ci_low_pct, ci_high_pct, m, Tphi2

    # Degenerate: zero MAD (all angles effectively identical)
    if mad(Tphi) == 0:
        m0 = mean(Tphi)      # matches your original semantics
        ue = np.array([1])   # single "event" marker to mimic prior behavior
        n_used = len(Tphi)
        # Use the point value as both CI bounds; keep your original 90 error
        ci_low_pct  = m0
        ci_high_pct = m0
        m     = np.asarray(Tphi)     # treat the observed set as the bootstrap sample
        Tphi2 = np.asarray(Tphi)
        err_full = 90.0
        return m0, err_full, ue, n_used, ci_low_pct, ci_high_pct, m, Tphi2

    # Remove outliers using MAD
    Tphi, ii  = resort(Tphi, ii)
    Tphi2, ii2 = outlier1(Tphi, ii)

    # If outlier removal zeroes out the data, return safe defaults
    if len(Tphi2) == 0:
        # circular mean of the pre-outlier set is a reasonable center
        m0 = cmean(Tphi, high=360)
        ue = np.array([0])
        n_used = 0
        err_full = 180.0
        return m0, err_full, ue, n_used, ci_low_pct, ci_high_pct, m, Tphi2

    # Bootstrap; ensure array type even if a function would return a scalar
    m = np.atleast_1d(boot1(Tphi2, 10000))
    m0 = cmean(m, high=360)

    # CIs on wrapped deltas
    if m.size > 0:
        delta = ((m - m0 + 180.0) % 360.0) - 180.0
        q_lo, q_hi = np.percentile(delta, [2.5, 97.5])
        ci_low_pct  = (m0 + q_lo) % 360.0
        ci_high_pct = (m0 + q_hi) % 360.0
        err_full = 2 * 1.96 * cstd(m, high=360)
    else:
        err_full = 180.0

    # Unique events & counts
    ue = uniqueevents(phi, cc, ii2, R1cc, R2cc)
    n_used = len(Tphi2)

    return m0, err_full, ue, n_used, ci_low_pct, ci_high_pct, m, Tphi2



# create histogram of results
def fhist1(phi,cc,lim):
    Tphi,Tcc,ii=C1_2(phi,cc,lim)
    Tphi,ii=resort(Tphi,ii)
    Tphi2,ii2=outlier1(Tphi,ii)
#    plt.figure()
    plt.hist(Tphi2,bins=25)
    plt.xlabel('Orientation Estimate')
    plt.ylabel('Counts')
    return 

#############################
# A few other random necessary ones
#############################

# Define rotation function
#   -rotates horizontal components CW from N by alpha (in degrees)
def rot2d(N,E,alpha):
    a=deg2rad(alpha)  # convert from degrees to radians
    r=cos(a)*N - sin(a)*E
    t=sin(a)*N + cos(a)*E
    return(r,t)

# root mean square
def rms(x):
    return sqrt(mean(abs(x)**2))

def find_nearest(array,value):
    return (abs(array-value)).argmin()

def checklen(st,hrs):
    # checks to see if there is enough downloaded data to run program
    L=len(st)
    for i in arange((L)):
        if (UTCDateTime(st[i].stats.endtime)-UTCDateTime(st[i].stats.starttime))+100 < hrs:
            return True        
    if var(st[0].data)<1 or var(st[1].data)<1 or var(st[2].data)<1:
        return True
    return False

# save resutls
def saved(R1cc,R2cc,R1phi,R2phi,loc='temp.dir'):
    if not os.path.exists(loc):
        os.mkdir(loc)
    savetxt(loc+'/'+'R1cc',R1cc)
    savetxt(loc+'/'+'R2cc',R2cc)
    savetxt(loc+'/'+'R1phi',R1phi)
    savetxt(loc+'/'+'R2phi',R2phi)
    return

"""
Functions from file formerly called Phases
Mostly deal with computing correlations
and rotating data
"""

# Resize arrays to all identical shapes
def resiz(x1,x2,x3):
    a1=len(x1); a2=len(x2); a3=len(x3)
    L=min(array([a1,a2,a3]))
    return x1[0:L], x2[0:L], x3[0:L]
    
# preprocess segments of data
# taper, zerophase filter, detrend
def sw1proc(T,LPF,HPF,corn=4):
    T.taper(type='hann',max_percentage=0.05)
    T.filter("lowpass",freq=LPF,corners=corn,zerophase=True)
    T.filter("highpass",freq=HPF,corners=corn,zerophase=True)
    T.detrend()
    
    return T

# DORAN-LASKE calculation for one freq, one orbit of surface wave
def SW1(TT,Rf,LPF,HPF,daz1,A,nameconv,winlen=10.0,ptype=0):
    # event info
    daz=daz1[0]/1000.0 # convert to KM
    baz=daz1[1] # angle from station to event        

    Rvel=getf(Rf,A) # Group velocity at Rf
    R1window=(1.0/(Rf/1000.0))*winlen

    # Process
    T=sw1proc(TT.copy(),LPF,HPF)

    # Window info
    arv=1.0/Rvel * daz
    r1=arv-R1window/2.0
    r2=arv+R1window/2.0

    dt=T[0].stats.starttime
    P=T.slice(starttime=dt+r1,endtime=dt+r2)
    
    rdat=P[0].data
    rdat2=P[1].data
    vdat=imag(sig.hilbert(P[2].data))
    
    # Ensure all data vectors are same length
    rdat,rdat2,vdat=resiz(rdat,rdat2,vdat)
    
    # rotate through and find max cc *****
    # one of the biggest clogs >> 1440 iterations each data point
    degs=360*4
    ang=arange((degs))
    cc=zeros((degs)); cc2=zeros((degs)); # cc3=zeros((degs))
    for k in ang:
        n,e=rot2d(rdat,rdat2,k/4.0)
        covmat=corrcoef(n,vdat)
        cc[k]=covmat[0,1]
        cstar=cov(vdat,n)/cov(vdat)
        cc2[k]=cstar[0,1]
        # cstar=cov(vdat,e)/cov(vdat)
        # cc3[k]=cstar[0,1]

    # Keep angle determined by cstar, but use rating from corrcoef
    #   Formulas in Stachnik paper
    ANG=cc2.argmax(); #CC[j]=cc[ANG]
    # correct for angles above 360
    or_ang= (baz- (360 - ANG/4.0) )

    # ADJUST FOR NAMING CONVENTION
    if nameconv==3: or_ang+=180

    if or_ang<0: or_ang+=360
    if or_ang>=360: or_ang-=360
    # Can plot xc

    # plotting:
    # ptype=0, no plot
    # ptype=1, Rayleigh plot
    # ptype=2, Love plot

    if ptype==1:
        import matplotlib.dates as dat
        X=P[0].times()
        T=zeros((len(X)))
        for q in arange((len(T))):
            T[q]=dt+r1+X[q]
        ZZ=dat.epoch2num(T)
        Z=dat.num2date(ZZ)
        n,e=rot2d(rdat,rdat2,ANG/4.0)
#        plt.figure()
        plt.plot(Z,vdat,label='Vertical')
#        savetxt('/Users/adoran/Desktop/T.txt',T)
#        savetxt('/Users/adoran/Desktop/vdat.txt',vdat)
#        savetxt('/Users/adoran/Desktop/n.txt',n)
        plt.hold("on")
        plt.plot(Z,n,label='BH1')    
        plt.legend(loc=4)
        plt.xlabel('Time')
        plt.ylabel('Counts')
        plt.title('D-L Results (%i mHz)' %(Rf))
    elif ptype==2:
        import matplotlib.dates as dat
        X=P[0].times()
        T=zeros((len(X)))
        for q in arange((len(T))):
            T[q]=dt+r1+X[q]
        ZZ=dat.epoch2num(T)
        Z=dat.num2date(ZZ)
        n,e=rot2d(rdat,rdat2,ANG/4.0)
        plt.figure()
        plt.subplot(121)
        plt.plot(Z,vdat,label='Vertical')
        plt.hold("on")
        plt.plot(Z,n,label='BH1')    
        plt.legend(loc=4)
        plt.xlabel('Time')
        plt.suptitle('D-L Results (%i mHz)' %(Rf))
        plt.subplot(122)
        plt.plot(Z,e,label='BH2')
        plt.xlabel('Time')
        plt.ylabel('Counts')
        plt.legend(loc=4)
    elif ptype==3:
        import matplotlib.dates as dat
        X=P[0].times()
        T=zeros((len(X)))
        for q in arange((len(T))):
            T[q]=dt+r1+X[q]
        n,e=rot2d(rdat,rdat2,ANG/4.0)
        plt.figure()
        plt.plot(T,vdat,label='Vertical')
    return or_ang, cc[ANG]

#############################
# Main function - import parameters
#############################

file_path = os.path.realpath(__file__)
directory = os.path.dirname(file_path)

flags = {"config": ("-c", str),
         "net": ("-n", str),
         "sta": ("-s", str),
         "sta_file": ("-sf", str), # ***
         "cha": ("-ch", str),
         "com": ("-co", str),
         "loc": ("-l", str),
         "Rdir": ("-r", str),
         "outdir": ("-o", str),
         "nameconv": ("-nc", int),
         "cat_client": ("-cc", str),
         "wf_client": ("-wf", str),    ## wc_client -> wf_client
         "time1": ("-t1", str),
         "time2": ("-t2", str),
         "localdata": ("-ldopt", bool),
         "localcoords": ("-lcopt", bool),
         "inputlat": ("-lat", float),
         "inputlon": ("-lon", float),
         "minmag": ("-mm", float),
         "mindeg_sw": ("-mndgsw", float),
         "maxdeg_sw": ("-mxdgsw", float),
         "maxdep_sw": ("-mxdpsw", float),
         "verb": ("-v", int),
         "constsave": ("-cs", bool),
         "finplot": ("-fp", bool),
         "savecat": ("-sc", bool),
         "multiple": ("-m", bool), # ***
        }

parser = argparse.ArgumentParser(description="Flags for cal_PSD.py. By default, these variables are established in \
                                 the cal_PSD.yml config file, and the flags serve as overrides for these values.")

for key in flags.keys():
    curr_var = key
    flag, vartype = flags[key]
    parser.add_argument(flag, f"--{curr_var}", metavar='', type=vartype)

args, unknown = parser.parse_known_args()
args_dict = args.__dict__

config_name = "DLOPy.yml"

if args_dict["config"]:
    config_name = args_dict["config"]

with open(f"{directory}/{config_name}") as f:
    config = yaml.safe_load(f)

for curr_var in flags.keys():
    if args_dict[curr_var]:
        config[curr_var] = args_dict[curr_var]

net = config["net"]
sta = config["sta"]
cha = config["cha"]
com = config["com"]
loc = config["loc"]

Rdir = config["Rdir"]

sta_file = config.get("sta_file", None)             # ***
multiple  = bool(config.get("multiple", False))

#****
MAP10 = np.loadtxt(os.path.join(Rdir, 'R.gv.10.txt'), skiprows=1)
MAP15 = np.loadtxt(os.path.join(Rdir, 'R.gv.15.txt'), skiprows=1)
MAP20 = np.loadtxt(os.path.join(Rdir, 'R.gv.20.txt'), skiprows=1)
MAP25 = np.loadtxt(os.path.join(Rdir, 'R.gv.25.txt'), skiprows=1)
MAP30 = np.loadtxt(os.path.join(Rdir, 'R.gv.30.txt'), skiprows=1)
MAP35 = np.loadtxt(os.path.join(Rdir, 'R.gv.35.txt'), skiprows=1)
MAP40 = np.loadtxt(os.path.join(Rdir, 'R.gv.40.txt'), skiprows=1)

SMALL_SIZE, MEDIUM_SIZE, BIGGER_SIZE = 16, 18, 20
plt.rc('font', size=SMALL_SIZE)
plt.rc('axes', titlesize=SMALL_SIZE, labelsize=MEDIUM_SIZE)
plt.rc('xtick', labelsize=SMALL_SIZE)
plt.rc('ytick', labelsize=SMALL_SIZE)
plt.rc('legend', fontsize=SMALL_SIZE)
plt.rc('figure', titlesize=BIGGER_SIZE)
plt.rcParams['figure.figsize'] = (10, 4)


if sta in ["BRIB", "TRIN", "SCOT", "SUTB", "DCMP", "SHEP"]:   # special cases
    loc = "01"

print("# sta = {0} net = {1} com = {2} loc = {3}".format(sta, net, com, loc))

nameconv = config["nameconv"]

cat_client = config["cat_client"]
wf_client = config["wf_client"]

time1 = config["time1"]
time2 = config["time2"]

print("# time1 = {0} time2 = {1}".format(time1, time2))

localdata = config["localdata"]
if localdata:
    from readlocal import *
    LF = getmsd1

localcoords = config["localcoords"]
if localcoords:
    inputlat = config["inputlat"]
    inputlon = config["inputlon"]

minmag = config["minmag"]

mindeg_sw = config["mindeg_sw"]
maxdeg_sw = config["maxdeg_sw"]
maxdep_sw = config["maxdep_sw"]

verb = config["verb"]

constsave = config["constsave"]
saveloc = sta

finplot = config["finplot"]

savecat = config["savecat"]
catname = sta + 'cat.txt'





"""
IMARY ORIENTATION PROGRAM
ADRIAN. K. DORAN
GABI LASKE
VERSION 1.0
RELEASED APRIL 2017
"""



#####################
# LOOP OVER MULTIPLE STATIONS ***
#####################
def loopstations(sta_value):

    # OLD  variables and special cases ***
    global sta
    global net
    global loc
    global catname
    global saveloc

    old_sta = sta
    old_net = net
    old_loc = loc

    sta = sta_value

    if sta in ["BRIB", "TRIN", "SCOT", "SUTB", "DCMP", "SHEP"]:
        loc = "01"

    saveloc = sta
    catname = sta + 'cat.txt'

    ########################

    # Define variables using param file inputs1
    client1=Client(cat_client)      # Catalog Client
    client2=Client(wf_client)     # Waveform Client
    t1=UTCDateTime(time1)           # start date
    t2=UTCDateTime(time2)           # end date




    with timed("get_events"):
        cat=client1.get_events(starttime=t1,endtime=t2,minmagnitude=minmag)#,maxmagnitude=maxmag)

    print("get_events done...")



    # get station info
    if localdata or localcoords:
        sta_lat=inputlat
        sta_lon=inputlon
    else:
        # Get station and event data
        inv=client1.get_stations(network=net,station=sta)
        sta_lat=inv[0][0].latitude
        sta_lon=inv[0][0].longitude

    # Organize station and event data
    # organize data from catalog
    L=len(cat.events)

    print("total EQ L = {0}".format(L))

    lat=zeros(L); lon=zeros(L); mag=zeros(L); stime=zeros(L); dep=zeros(L); deg=zeros(L); baz=zeros(L)
    for i in arange(L):
        lat[i]=cat.events[i].origins[0].latitude                            # latitude
        lon[i]=cat.events[i].origins[0].longitude                            # longitude
        dep[i]=cat.events[i].origins[0].depth                            # depth
        stime[i]=UTCDateTime(cat.events[i].origins[0].time)                # event start time
        mag[i]=cat.events[i].magnitudes[0].mag                            # magnitude
        #daz1=obspy.core.util.gps2DistAzimuth(sta_lat,sta_lon,lat[i],lon[i])   # depricated version
        daz1=obspy.geodetics.gps2dist_azimuth(sta_lat,sta_lon,lat[i],lon[i])   # distance b/t station and event
        #deg[i]=obspy.core.util.kilometer2degrees(daz1[0]/1000)                # depricated version
        deg[i]=obspy.geodetics.kilometer2degrees(daz1[0]/1000)               
        baz[i]=daz1[1]                                                     # angle from station to event
    # get index of repeat events, save for later
    reps=unique(catclean(stime,lat,lon,mag))

    # to save catalog:
    if savecat:
        ts=array([],dtype=object)
        h1=array(['Time','Lon','Lat','Dep(km)','Mag'],dtype=object)
        for i in arange(L):
            ts=append(ts,UTCDateTime(stime[i]))
        catprint2=array((ts,lon,lat,dep/1000,mag),dtype=object).T
        catprint=vstack((h1,catprint2))
        savetxt(catname,catprint,fmt="%s")

    print("save catalog done...")
    #   INITIALIZE INITIALIZE
    # Initialize surface wave arrays
    numsurfcalcs=7
    R1phi=zeros([L,numsurfcalcs]); R1cc=zeros([L,numsurfcalcs])
    R2phi=zeros([L,numsurfcalcs]); R2cc=zeros([L,numsurfcalcs])
    # Initialize Stachnik arrays
    R4phi=zeros((L)); R4cc=zeros((L))

    hrs=4*60*60     # Length of data to download



    # map10=loadtxt(Rdir+'/R.gv.10.txt'); map15=loadtxt(Rdir+'/R.gv.15.txt')
    # map20=loadtxt(Rdir+'/R.gv.20.txt'); map25=loadtxt(Rdir+'/R.gv.25.txt')
    # map30=loadtxt(Rdir+'/R.gv.30.txt'); map35=loadtxt(Rdir+'/R.gv.35.txt')
    # map40=loadtxt(Rdir+'/R.gv.40.txt')

    # LOOP OVER ALL EVENTS
    for j in arange((L)):
        print("# {0}-the event of total {1}".format(j, L))
        #GET WAVEFORMS, including protections
        try:
            if not localdata:
                # download data from client
                with timed("get_waveforms"):
                    s=client2.get_waveforms(net,sta,loc,cha,UTCDateTime(stime[j]),UTCDateTime(stime[j]+hrs))
            else:
                # access local data
                s=getmsd1(UTCDateTime(stime[j]))
                
            s.merge()
            
            # don't want any masked data or data with nans
            for q in arange((len(s))):
                if ma.count_masked(s[q].data)>0:
                    continue
                
            # remove mean and trend
            s.detrend()
            s.detrend('linear')
            if len(s)<3:
                continue
            
            # [0] bh1, [1] bh2, [2] bhz
            st=org(s.copy(),nameconv) 
                # organizes data by coordinate system
                # also downsamples to <10 Hz

            # check data length, data quality        
            if checklen(st,hrs):
                continue

        except Exception as e:
            continue
            
        # get some additional parameters
        #daz1=obspy.core.util.gps2DistAzimuth(sta_lat,sta_lon,lat[j],lon[j])
        daz1=obspy.geodetics.gps2dist_azimuth(sta_lat,sta_lon,lat[j],lon[j])
        daz2=copy(daz1)
        Rearth=6371.25*1000; circE=2*np.pi*Rearth;
        daz2[0]=circE-daz2[0]; daz2[1]=daz2[1]+180  # major & minor arc calculation
        if daz2[1]>=360: daz2[1]-=360
        
        # SURFACE WAVE CALCULATIONS
        
        # conditions
        # minimum distance, maximum distance, and maximum depth
        if deg[j]<mindeg_sw or deg[j]>maxdeg_sw or dep[j]>=maxdep_sw*1000:
            continue
        # clean catalog of repeats (repeat conditions set in catclean)
        if j in reps:
            continue
        
        with timed("pathvels"):
            # get path-averaged group velocities
            Ray1,Ray2 = pathvels(sta_lat, sta_lon, lat[j], lon[j], MAP10, MAP15, MAP20, MAP25, MAP30, MAP35, MAP40)

        with timed("SW1"):
            # FOR EACH FREQUENCY AND ORBIT, calculate arrival angle

            # freq 1 (40 mHz)
            Rf=40.0; HPF=0.035; LPF=0.045
            
            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=20.0,ptype=0)
            R1phi[j,0]=ANG; R1cc[j,0]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=24.0,ptype=0)
            R2phi[j,0]=ANG; R2cc[j,0]=cc

            # freq 2 (35 mHz)
            Rf=35.0; HPF=0.030; LPF=0.040

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=17.0,ptype=0)
            R1phi[j,1]=ANG; R1cc[j,1]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=20.0,ptype=0)
            R2phi[j,1]=ANG; R2cc[j,1]=cc


            # freq 3 (30 mHz)
            Rf=30.0; HPF=0.025; LPF=0.035

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=14.0,ptype=0)
            R1phi[j,2]=ANG; R1cc[j,2]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=16.0,ptype=0)
            R2phi[j,2]=ANG; R2cc[j,2]=cc


            # freq 4 (25 mHz)
            Rf=25.0; HPF=0.020; LPF=0.030

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=12.0,ptype=0)
            R1phi[j,3]=ANG; R1cc[j,3]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=13.0,ptype=0)
            R2phi[j,3]=ANG; R2cc[j,3]=cc

            # freq 5 (20 mHz)
            Rf=20.0; HPF=0.015; LPF=0.025

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=10.0,ptype=0)
            R1phi[j,4]=ANG; R1cc[j,4]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=10.0,ptype=0)
            R2phi[j,4]=ANG; R2cc[j,4]=cc

            # freq 6 (15 mHz)
            Rf=15.0; HPF=0.020; LPF=0.010

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=10.0,ptype=0)
            R1phi[j,5]=ANG; R1cc[j,5]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=10.0,ptype=0)
            R2phi[j,5]=ANG; R2cc[j,5]=cc

            # freq 7 (10 mHz)
            Rf=10.0; HPF=0.005; LPF=0.015

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz1,Ray1,nameconv,winlen=7.0,ptype=0)
            R1phi[j,6]=ANG; R1cc[j,6]=cc

            ANG,cc=SW1(st.copy(),Rf,LPF,HPF,daz2,Ray2,nameconv,winlen=7.0,ptype=0)
            R2phi[j,6]=ANG; R2cc[j,6]=cc

        # save up that Data
        if constsave:
            saved(R1cc,R2cc,R1phi,R2phi,loc=str(saveloc))          
            
        if verb >= 2 and (j % 10 == 0 or j == L - 1):
            print(f"{sta}: {j+1} / {L}")

    # PLOT ALL RESULTS
    if finplot==1:
        plt.figure()
        plt.subplot(1,1,1)
        plt.title('Surf Waves')
        plt.plot(R1cc,R1phi,'x',R2cc,R2phi,'x')
        plt.ylim([0,360]); plt.xlim([0,1])

    # SAVE DATA
    saved(R1cc,R2cc,R1phi,R2phi,loc=str(saveloc))


    ### from plotonly

    # font size
    SMALL_SIZE = 16
    MEDIUM_SIZE = 18
    BIGGER_SIZE = 20

    plt.rc('font', size=SMALL_SIZE)          # controls default text sizes
    plt.rc('axes', titlesize=SMALL_SIZE)     # fontsize of the axes title
    plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
    plt.rc('xtick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('ytick', labelsize=SMALL_SIZE)    # fontsize of the tick labels
    plt.rc('legend', fontsize=SMALL_SIZE)    # legend fontsize
    plt.rc('figure', titlesize=BIGGER_SIZE)  # fontsize of the figure title

    plt.rcParams['figure.figsize'] = 10, 4

    # com = "HHN"








    """
    Final Orientation Calculation File
    A. Doran and G. Laske
    """

    #####################
    ## ANGLE CALCULATION PARAMETERS
    #####################

    # location of result files
    loc1= sta+'/'

    # create output directory
    if not os.path.exists(loc1):
        os.makedirs(loc1)

    # PARAM
    cc_min = 0.95
    # PARAM
    LIM=0.8         # CC limit for Surface wave calculations
    LIM=0.95         # CC limit for Surface wave calculations
    LIM = cc_min
    ### Specify phases to use
    R1use=1
    R1_40=1; R1_35=1; R1_30=1; R1_25=1; R1_20=1; R1_15=1; R1_10=1
    R2use=1
    R2_40=1; R2_35=1; R2_30=1; R2_25=1; R2_20=1; R2_15=1; R2_10=1


    ## Load files
    R1phi=loadtxt(loc1+'R1phi')
    R1cc=loadtxt(loc1+'R1cc')

    R2phi=loadtxt(loc1+'R2phi')
    R2cc=loadtxt(loc1+'R2cc')

    ######################
    ### FINAL ANGLE CALCULATIONS
    ######################

    # Initialize arrays
    L=len(R1phi)
    phis=array([])
    ccs=array([])
    finval=array([]); finerr=array([])
    N=array([])

    N=full((L,L),-1.0)
    LN=zeros((L))
    phases=array([])

    startL=0
    endL=0

    A=array([L])
    if endL!=0:
        A=array([endL])
        

    # If not all calculations are desired, adjust accordingly
    sha=shape(R1phi)
    if R1use==0: R1cc=zeros(sha)
    if R1use==1 and R1_40==0: R1cc[:,0]=zeros((sha[0]))
    if R1use==1 and R1_35==0: R1cc[:,1]=zeros((sha[0]))
    if R1use==1 and R1_30==0: R1cc[:,2]=zeros((sha[0]))
    if R1use==1 and R1_25==0: R1cc[:,3]=zeros((sha[0]))
    if R1use==1 and R1_20==0: R1cc[:,4]=zeros((sha[0]))
    if R1use==1 and R1_15==0: R1cc[:,5]=zeros((sha[0]))
    if R1use==1 and R1_10==0: R1cc[:,6]=zeros((sha[0]))

    if R2use==0: R2cc=zeros(sha)    
    if R2use==1 and R2_40==0: R2cc[:,0]=zeros((sha[0]))
    if R2use==1 and R2_35==0: R2cc[:,1]=zeros((sha[0]))
    if R2use==1 and R2_30==0: R2cc[:,2]=zeros((sha[0]))
    if R2use==1 and R2_25==0: R2cc[:,3]=zeros((sha[0]))
    if R2use==1 and R2_20==0: R2cc[:,4]=zeros((sha[0]))
    if R2use==1 and R2_15==0: R2cc[:,5]=zeros((sha[0]))
    if R2use==1 and R2_10==0: R2cc[:,6]=zeros((sha[0]))

    for i in A:
        # create one massive list with necessary angles and cc values
        phis=concatenate((flatten(R1phi[startL:i,:]),flatten(R2phi[startL:i,:])))
        ccs=concatenate((flatten(R1cc[startL:i,:]),flatten(R2cc[startL:i,:])))
        
        # Doran-Laske calculation
        val, err, n, ph, ci_low_pct, ci_high_pct, m, Tphi2 = fcalc1(phis, ccs, LIM, R1cc, R2cc)




        finval=append(finval,val)
        finerr=append(finerr,err)
        phases=append(phases,ph)
        for k in arange((len(n))):
            N[k,i-1]=n[k]
        LN[i-1]=len(n)

    
        
    # output results to termianl
    print ("Station = ", loc1)
    #print ("D-L mean, error, data included, unique events: %.2f, %.2f, %i, %i" %(finval[-1],finerr[-1],phases[-1],max(LN))
    #print ("D-L CC level: %f" %(LIM))

    print ("D-L CC level = ",LIM)
    print ("D-L mean = ",finval[-1])
    print ("D-L error = ",finerr[-1])
    print ("D-L phase = ",phases[-1])
    print ("D-L max = ",max(LN))
    print ("# Rsta = ",sta, finval[-1], finerr[-1], phases[-1], max(LN))

    az_mean = finval[-1]
    az_error = finerr[-1]

    # ****

    if Tphi2 is not None and len(Tphi2) > 0:
        sample_delta = ((Tphi2 - val + 180.0) % 360.0) - 180.0

        q_lo = ((ci_low_pct  - val + 180.0) % 360.0) - 180.0
        q_hi = ((ci_high_pct - val + 180.0) % 360.0) - 180.0

        binw = 1
        
        W = builtins.max(abs(q_lo), abs(q_hi), float(np.percentile(np.abs(sample_delta), 95)))
        W = float(np.ceil(W / binw) * binw) + binw
        bins = np.arange(-W, W + binw, binw)

        fig_hist, axh = plt.subplots()

        axh.hist(sample_delta, bins=bins, color='0.85', edgecolor='k', linewidth=0.4)
        axh.axvline(0.0, color='r', linewidth=0.8)
        axh.axvline(q_lo, color='r', linestyle='--', linewidth=0.6)
        axh.axvline(q_hi, color='r', linestyle='--', linewidth=0.6)
        axh.set_xlim(-W*1.5, W*1.5)

        def abs_deg_tick(x, pos):
            return f"{(val + x) % 360:.0f}°"
        
        axh.xaxis.set_major_formatter(FuncFormatter(abs_deg_tick))
        axh.xaxis.set_major_locator(MultipleLocator(5))
        axh.xaxis.set_minor_locator(MultipleLocator(1))
    
        axh.set_xlabel("Spread of Azimuth", fontsize=8)
        axh.set_ylabel("Num of Data Points Used", fontsize=8)
        axh.set_title(f"Bootstrap Mean and Confidence Intervals: {net}.{sta_line}", fontsize=10)
        axh.tick_params(axis='both', labelsize=8)
        axh.grid(alpha=0.2)

        anno = (f"Mean: {val:.2f}°\n"
            f"CI2.5: {ci_low_pct:.2f}°\n"
            f"CI97.5: {ci_high_pct:.2f}°\n"
            f"N: {len(Tphi2)}")
        axh.text(0.98, 0.98, anno,
            transform=axh.transAxes, ha='right', va='top',
            fontsize=8, bbox=dict(fc="white", ec="0.5", alpha=0.9))

        out_hist = os.path.join(loc1+f"{sta_line}.bootstrap_hist.pdf")
        fig_hist.savefig(out_hist, bbox_inches="tight")
        plt.close(fig_hist)

        if m is None or len(m) == 0:
            # BOOTNUM_DEFAULT is optional; fall back to 10000 if not defined
            bootnum_plot = BOOTNUM_DEFAULT if 'BOOTNUM_DEFAULT' in globals() else 10000
            m_plot = boot1(Tphi2, bootnum_plot)
        else:
            m_plot = m

        m_delta   = ((m_plot - val + 180.0) % 360.0) - 180.0

        q_lo_mean = ((ci_low_pct  - val + 180.0) % 360.0) - 180.0
        q_hi_mean = ((ci_high_pct - val + 180.0) % 360.0) - 180.0

        ci_lo_abs = (val + q_lo_mean) % 360.0
        ci_hi_abs = (val + q_hi_mean) % 360.0

        binw = 0.25
        Wm = builtins.max(abs(q_lo_mean), abs(q_hi_mean), float(np.percentile(np.abs(m_delta), 95)))
        Wm   = float(np.ceil(Wm / binw) * binw) + binw
        bins_boot = np.arange(val - Wm, val + Wm + binw, binw)
        m_abs = m_delta + val

        fig_boot, axb = plt.subplots()


        axb.hist(m_abs, bins=bins_boot, color='0.85', edgecolor='k', linewidth=0.4)
        axb.axvline(val,      color='r', linewidth=0.8, label='Mean')
        axb.axvline(ci_lo_abs, color='r', linestyle='--', linewidth=0.6, label='CI 2.5%')
        axb.axvline(ci_hi_abs, color='r', linestyle='--', linewidth=0.6, label='CI 97.5%')
        axb.set_xlim(val - Wm, val + Wm)

        def abs_deg_tick(x, pos):
            return f"{(val + x) % 360:.0f}°"
        
        axh.xaxis.set_major_formatter(FuncFormatter(abs_deg_tick))
        axh.xaxis.set_major_locator(MultipleLocator(5))
        axh.xaxis.set_minor_locator(MultipleLocator(1))

        axb.set_xlabel(f"Azimuth (bootstrap means)", fontsize=8)
        axb.set_ylabel("Num of Bootstrap Resamples", fontsize=8)
        axb.set_title(f"Bootstrap Means/95% CI: {net}.{sta}", fontsize=10)
        axb.tick_params(axis='both', labelsize=8)
        axb.grid(alpha=0.2)

        anno_boot = (f"Mean: {val:.2f}°\n"
                    f"CI2.5: {ci_low_pct:.2f}°\n"
                    f"CI97.5: {ci_high_pct:.2f}°\n"
                    f"Boot Resamples: {len(m_plot)}")
        axb.text(0.98, 0.98, anno_boot, transform=axb.transAxes, ha='right', va='top',
                fontsize=8, bbox=dict(fc='white', ec='0.5', alpha=0.9))

        out_boot = os.path.join(loc1, f"{sta}.bootstrap_means_hist.pdf")
        fig_boot.savefig(out_boot, bbox_inches="tight")
        plt.close(fig_boot)


    ci_low  = (az_mean - (az_error*2.0)) % 360.0  
    ci_high = (az_mean + (az_error*2.0)) % 360.0 

    n_used=int(phases[-1])


    xcmin = 0

    xcmax = 1

    amin = az_mean - 180
    amax = az_mean + 180

    az_min = az_mean - (az_error*2)
    az_max= az_mean + (az_error*2)

    print("# az_mean = ", az_mean, " az_min = ", az_min, " az_max = ", az_max)

    # create figure

    CEN=finval[-1]
    YLIM1=[-10+CEN,10+CEN]

    plt.figure(figsize=(15, 5))
    plt.subplot(1,1,1)
    plt.title("DLOPy result "+net+"."+sta)

    plt.plot(R1cc,centerat(R1phi,m=CEN),'x', color="black")
    plt.plot(R2cc,centerat(R2phi,m=CEN),'x', color="black")

    plt.plot([xcmin, xcmax], [az_mean, az_mean], color="red", linewidth=1.55, linestyle="solid")
    plt.plot([xcmin, xcmax], [az_min, az_min], color="red", linewidth=1.15, linestyle="dashed")
    plt.plot([xcmin, xcmax], [az_max, az_max], color="red", linewidth=1.15, linestyle="dashed")

    plt.plot([cc_min, cc_min], [amin, amax], color="black", linewidth=1.15, linestyle="dashed")

    #plt.ylabel(com+" orientation ($^\circ$)")
    plt.ylabel(com+" orientation (deg.)")

    plt.xlabel("Cross-correlation value")

    plt.ylim(amin, amax) 

    plt.xlim(xcmin,xcmax)

    az_mean_out = "{:.2f}".format(az_mean)
    az_min_out = "{:.2f}".format(az_min)
    az_max_out = "{:.2f}".format(az_max)


    t=plt.text(xcmin+0.01, amin+20, "Az.: "+str(az_mean_out)+" ("+az_min_out+"-"+az_max_out+")"  )
    t.set_bbox(dict(facecolor='white', alpha=0.7, edgecolor='white'))


    plt.grid()
    plt.savefig(loc1+"/"+sta+".v2.cluster.pdf")
    
    az_mean  = finval[-1]
    az_error = finerr[-1]

    # restore old globals ***
    sta = old_sta
    net = old_net
    loc = old_loc

    print_profile(header=f"Station: {net}.{sta_value}")
    _profile_totals.clear()

    return az_mean, az_error, ci_low, ci_high, az_min, az_max, n_used






'''
Reads R1phi, R1cc, R2phi, R2cc files before my code can even create them, immediately
running when the code runs: persistent "FileNotFoundError: YUBA/R1phi not found."


'''

# ***
# Final Orientation Calculation File
# A. Doran and G. Laske


# #####################
# ## ANGLE CALCULATION PARAMETERS
# #####################

# loc1= sta+'/'

# LIM=0.8         # CC limit for Surface wave calculations
# LIM=0.95         # CC limit for Surface wave calculations

# ### Specify phases to use
# R1use=1
# R1_40=1; R1_35=1; R1_30=1; R1_25=1; R1_20=1; R1_15=1; R1_10=1
# R2use=1
# R2_40=1; R2_35=1; R2_30=1; R2_25=1; R2_20=1; R2_15=1; R2_10=1

# ## Load files
# R1phi=loadtxt(loc1+'R1phi')
# R1cc=loadtxt(loc1+'R1cc')

# R2phi=loadtxt(loc1+'R2phi')
# R2cc=loadtxt(loc1+'R2cc')



# ######################
# ### FINAL ANGLE CALCULATIONS
# ######################

# # Initialize arrays
# L=len(R1phi)
# phis=array([])
# ccs=array([])
# finval=array([]); finerr=array([])
# N=array([])

# N=full((L,L),-1.0)
# LN=zeros((L))
# phases=array([])

# startL=0
# endL=0

# A=array([L])
# if endL!=0:
#     A=array([endL])
    

# # If not all calculations are desired, adjust accordingly
# sha=shape(R1phi)
# if R1use==0: R1cc=zeros(sha)
# if R1use==1 and R1_40==0: R1cc[:,0]=zeros((sha[0]))
# if R1use==1 and R1_35==0: R1cc[:,1]=zeros((sha[0]))
# if R1use==1 and R1_30==0: R1cc[:,2]=zeros((sha[0]))
# if R1use==1 and R1_25==0: R1cc[:,3]=zeros((sha[0]))
# if R1use==1 and R1_20==0: R1cc[:,4]=zeros((sha[0]))
# if R1use==1 and R1_15==0: R1cc[:,5]=zeros((sha[0]))
# if R1use==1 and R1_10==0: R1cc[:,6]=zeros((sha[0]))

# if R2use==0: R2cc=zeros(sha)    
# if R2use==1 and R2_40==0: R2cc[:,0]=zeros((sha[0]))
# if R2use==1 and R2_35==0: R2cc[:,1]=zeros((sha[0]))
# if R2use==1 and R2_30==0: R2cc[:,2]=zeros((sha[0]))
# if R2use==1 and R2_25==0: R2cc[:,3]=zeros((sha[0]))
# if R2use==1 and R2_20==0: R2cc[:,4]=zeros((sha[0]))
# if R2use==1 and R2_15==0: R2cc[:,5]=zeros((sha[0]))
# if R2use==1 and R2_10==0: R2cc[:,6]=zeros((sha[0]))

# for i in A:
#     # create one massive list with necessary angles and cc values
#     phis=concatenate((flatten(R1phi[startL:i,:]),flatten(R2phi[startL:i,:])))
#     ccs=concatenate((flatten(R1cc[startL:i,:]),flatten(R2cc[startL:i,:])))
    
#     # Doran-Laske calculation
#     val,err,n,ph=fcalc1(phis,ccs,LIM,R1cc,R2cc)
#     finval=append(finval,val)
#     finerr=append(finerr,err)
#     phases=append(phases,ph)
#     for k in arange((len(n))):
#         N[k,i-1]=n[k]
#     LN[i-1]=len(n)
    
# # output results to termianl
# print ("Station = ", loc1)
# #print ("D-L mean, error, data included, unique events: %.2f, %.2f, %i, %i" %(finval[-1],finerr[-1],phases[-1],max(LN))
# #print ("D-L CC level: %f" %(LIM))

# print ("D-L CC level = ",LIM)
# print ("D-L mean = ",finval[-1])
# print ("D-L error = ",finerr[-1])
# print ("D-L phase = ",phases[-1])
# print ("D-L max = ",max(LN))
# print ("# Rsta = ",sta, finval[-1], finerr[-1], phases[-1], max(LN))

# # create figure

# CEN=finval[-1]
# YLIM1=[-10+CEN,10+CEN]

# plt.figure()
# plt.subplot(1,1,1)
# plt.title('DLOPy results '+sta, fontsize=16)

# plt.plot(R1cc,centerat(R1phi,m=CEN),'x',R2cc,centerat(R2phi,m=CEN),'x')
# plt.ylabel('BH1 Orientation \n Angle ($^\\circ$)',fontsize=16)         # compiler is reading "\c"ircs as escape; \circs -> \\circs
# plt.ylim([CEN-180,CEN+180]); plt.xlim([0,1])
# plt.yticks(fontsize=16)
# plt.xticks(fontsize=16)




############################ CONTROLLER ############################

def parse(line, default_net, default_loc):
    s = line.strip()
    if not s or s.startswith("#"):
        return None
    return default_net, s, default_loc

def summarize(rows, path):
    with open(path, "w", newline="") as f:
        
        w = csv.writer(f)
        w.writerow(["Network", "Station", "Location", "Mean Azimuth", "CI95 Low Interval", "CI95 High Interval", "Data Points Used", "Within5?", "AbsOffset"])
        for r in rows:
            
            w.writerow([r["net"], r["sta"], r["loc"], f"{r['az_mean']:.4f}", f"{r['ci_low']:.4f}", f"{r['ci_high']:.4f}", r["n_used"],
                        "YES" if r["within_5deg"] else "NO", f"{r['abs_offset']:.4f}",])

def azwrap(x):
    y = ((x + 180.0) % 360.0) - 180.0
    return y


if __name__ == "__main__":
    if multiple and sta_file:
        summary_rows = []
        with open(sta_file, "r") as fh:
            lines = fh.readlines()

        for line in lines:
            parsed = parse(line, net, loc)
            if not parsed:
                continue
            net_line, sta_line, loc_line = parsed

            print(f"station: {net_line}.{sta_line}.{loc_line}")

            # temporarily override globals net/loc for this run
            old_net = net
            old_loc = loc
            net = net_line
            loc = loc_line

            try:
                az_mean, az_error, ci_low, ci_high, az_min, az_max, n_used = loopstations(sta_line)

            finally:
                net=old_net
                loc = old_loc

            offset = azwrap(az_mean)
            summary_rows.append({
            "net": net_line,
            "sta": sta_line,
            "loc": loc_line,
            "az_mean": az_mean,
            "ci_low": ci_low,
            "ci_high": ci_high,
            # "az_min": az_min,
            # "az_max": az_max,
            "n_used": n_used,
            "abs_offset": abs(offset),
            "within_5deg": abs(offset) <= 5.0,
        })

        summarize(summary_rows, "summary.csv")

        print("done.")

    else:
        # FOR SINGLE STATION RUNS
        loopstations(sta)


