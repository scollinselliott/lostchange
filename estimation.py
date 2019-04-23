# -*- coding: utf-8 -*-

# estimation.py
# this script contains the method to simulate U(t) as well as produce graphics
# in the paper "Quantifying Artifacts over Time: Interval Estimation of a
# Poisson Distribution using the Jeffreys Prior" submitted to Archaeometry

import csv
import numpy
import random
import matplotlib.pyplot as plt
import collections
import scipy.misc
import scipy.stats

from rpy2.robjects import r
import rpy2.robjects as robjects
from rpy2.robjects import *

r = robjects.r

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter

##################################################################################
##################################################################################
##################################################################################




# load data



#test data: comment out to use RPP data
################################################################################
#rawtable = numpy.array([['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101'],
#       ['test', 'coin', '1', '1', '101']],
#      dtype='<U4')
#contextresult = [['test']]
#assemblagesize = [len(rawtable)]
#sitearea = [1]





#RPP data
###################################################################################
reader = csv.reader(open("input-rpp.csv", "r"), delimiter=",")
rawtable = list(reader)
rawtable = numpy.array(rawtable[1:])
reader = csv.reader(open("outputcontext-rpp.csv", "r"), delimiter=",")
z = list(reader)
contextresult = z[1:]
assemblagesize = [36, 33, 5, 9, 83]
sitearea = [550.49, 764.63, 636.65, 995.29, 2947.06]













# estimation algorithm
################################################################################
################################################################################

alldates = numpy.array(rawtable[:,3:], dtype = int)

startdate = min(alldates.flatten())
enddate = max(alldates.flatten()) + 2
daterange = range(startdate,enddate)

#geometric rate
gamma = 0.1

k = 1000
betaminusalpha = 0

for site in contextresult:
        
    simuls = numpy.zeros((len(daterange),0))   

    results = numpy.zeros((len(daterange),3)) #date | musub | varsub 
    results[:,0] = daterange
    
    n = assemblagesize[contextresult.index(site)]
    area = sitearea[contextresult.index(site)]

    ktotal = numpy.zeros((len(daterange),2)) #keeps track of the sums of artifacts over the course of all k simulation runs
    ktotal[:,0] = daterange
    
    dates = rawtable[rawtable[:,0] == site][:,3:].astype(int)
    quant = rawtable[rawtable[:,0] == site][:,2].astype(int)
           
    tau = numpy.sum(dates[:,1] - dates[:,0])
    
    for j in range(k):
        ktotal[:,1] = 0
        
        #generate a random value of gamma according to a probability distribution
        #gamma = random.uniform(0.01,0.05)

        #this is to generate the random interval of dates from the subsample: each coin will have its own random interval
        dates2 = numpy.empty([0,2]) #the simulated dates 

        rowcount = 0
        for row in dates: 
            a = row[0]
            b = row[1]
            if a < 5000:
                if b < 5000:
                    xquant = quant[rowcount]
                    while xquant > 0:
                        alpha = random.randint(a,b)
                        beta =  alpha + numpy.random.geometric(gamma,1)
                        dates2 = numpy.vstack((dates2,numpy.array([alpha,beta])))
                        xquant = xquant - 1
            rowcount = rowcount + 1
        
        #to calculate the number of artifacts per year according to those random simulated dates
        for date in dates2: 
            alpha = date[0]
            beta = date[1]
            start = numpy.where(ktotal[:,0] == alpha)[0][0]
            if beta >= enddate: #truncate the value of beta if it surpasses the enddate b
                beta = enddate - 1
            end = numpy.where(ktotal[:,0] == beta)[0][0]
            ktotal[start:end,1] = ktotal[start:end,1] + 1 #this adds one to each row, adding one throughout the entire interval

        simul = ktotal[:,1][numpy.newaxis].transpose()
        simuls = numpy.array(numpy.hstack((simuls,simul)))

    #simuls = simuls / area ############### moved to graph

    for j in range(simuls.shape[0]):  
        tempdata = simuls[j,:]            
        if sum(tempdata) != 0: 
            mu = numpy.mean(tempdata)
            sigma2 = numpy.var(tempdata, ddof = 1)
            results[j,1] = mu
            results[j,2] = sigma2

    #results = results[:-1,:]
    mu = results[:,1]
    sigma2 = results[:,2]
        
    if site == ['MZ']:
        MZarea = area
        MZint = numpy.zeros((results.shape[0],results.shape[1])) #lower 1 # upper 2
        MZint[:,0] = results[:,0]
        
        for row in range(simuls.shape[0]):
            sumw = sum(simuls[row,:])
            r.assign('remotek',k)
            r.assign('remotew',sumw)
            lower = r('qgamma(0.05,remotew + 0.5,remotek)')
            upper = r('qgamma(0.95,remotew + 0.5,remotek)')
        
            MZint[row,1] = numpy.asarray(lower)
            MZint[row,2] = numpy.asarray(upper)
    
    if site == ['CN']:
        CNarea = area
        CNint = numpy.zeros((results.shape[0],results.shape[1])) #lower 1 # upper 2
        CNint[:,0] = results[:,0]
        
        for row in range(simuls.shape[0]):
            sumw = sum(simuls[row,:])
            r.assign('remotek',k)
            r.assign('remotew',sumw)
            lower = r('qgamma(0.05,remotew + 0.5,remotek)')
            upper = r('qgamma(0.95,remotew + 0.5,remotek)')
        
            CNint[row,1] = numpy.asarray(lower)
            CNint[row,2] = numpy.asarray(upper)

    if site == ['PI']:
        PIarea = area
        PIint = numpy.zeros((results.shape[0],results.shape[1])) #lower 1 # upper 2
        PIint[:,0] = results[:,0]
        
        for row in range(simuls.shape[0]):
            sumw = sum(simuls[row,:])
            r.assign('remotek',k)
            r.assign('remotew',sumw)
            lower = r('qgamma(0.05,remotew + 0.5,remotek)')
            upper = r('qgamma(0.95,remotew + 0.5,remotek)')
        
            PIint[row,1] = numpy.asarray(lower)
            PIint[row,2] = numpy.asarray(upper)
            
    if site == ['PT']:
        PTarea = area
        PTint = numpy.zeros((results.shape[0],results.shape[1])) #lower 1 # upper 2
        PTint[:,0] = results[:,0]
        
        for row in range(simuls.shape[0]):
            sumw = sum(simuls[row,:])
            r.assign('remotek',k)
            r.assign('remotew',sumw)
            lower = r('qgamma(0.05,remotew + 0.5,remotek)')
            upper = r('qgamma(0.95,remotew + 0.5,remotek)')
        
            PTint[row,1] = numpy.asarray(lower)
            PTint[row,2] = numpy.asarray(upper)





#   Figures
################################################################################
################################################################################




#   Figure 1: Coin Counts by Issue Date
################################################################################
################################################################################

periods = numpy.array(['3rd-2nd c. BCE', '1st c. BCE', 'Triumviral', '27 BCE - 14',
       '14-41', '41-54', '54-69', '69-96', '96-117', '117-138', '138-161',
       '161-180', '180-192', '193-222', '222-238', '238-259', '259-275',
       '275-294', '294-317', '317-330', '330-348', '348-364', '364-378',
       '378-388', '388-402'],
      dtype='<U15')

issuecounts = numpy.array([[0, 1, 0, 1],
       [0, 1, 0, 2],
       [0, 0, 0, 1],
       [2, 0, 1, 4],
       [2, 1, 0, 3],
       [1, 1, 0, 4],
       [0, 0, 0, 1],
       [0, 0, 0, 1],
       [0, 0, 0, 0],
       [0, 0, 0, 1],
       [0, 1, 0, 1],
       [0, 1, 0, 0],
       [0, 0, 0, 0],
       [0, 0, 0, 0],
       [0, 0, 0, 3],
       [0, 0, 0, 1],
       [0, 1, 0, 3],
       [0, 1, 0, 1],
       [0, 1, 0, 4],
       [0, 0, 0, 0],
       [0, 3, 0, 0],
       [0, 3, 0, 0],
       [0, 0, 0, 0],
       [0, 1, 1, 0],
       [0, 1, 0, 0]])


ind = numpy.arange(periods.shape[0]) 

fig = plt.figure(figsize=(13,8))
plt.bar(ind, issuecounts[:,2], bottom=issuecounts[:,1] + issuecounts[:,0] + issuecounts[:,3], color = 'red', alpha = 0.5, label = 'Case Nuove')
plt.bar(ind, issuecounts[:,1], bottom=issuecounts[:,0] + issuecounts[:,3], color = 'blue', alpha = 0.5, label = 'Pievina')
plt.bar(ind, issuecounts[:,3], bottom=issuecounts[:,0], color = 'black', alpha = 0.4, label = 'Podere Marzuolo')
plt.bar(ind, issuecounts[:,0], color = 'orange', alpha = 0.5, label = 'Podere Terrato')
plt.xticks(ind, periods)
plt.xticks(rotation=90)
plt.ylabel('Count')
plt.title(u'Coin Count by Issue Date (Roman Peasant Project, 2009-2014), Periodization after Reece (1995, 1996)',  loc = 'left')
plt.legend(loc='upper right')
plt.plot()






#   Figure 3: Plots of U(t), sigma(t), and 3D histograms / Edit variables to produce subplots
################################################################################
################################################################################

fig = plt.figure(figsize=(15,5))
plt.subplot(1, 2, 1)

plt.plot(mu, color = '#000000', label = u'$U(t)$, γ = 1, n = 100')
plt.xlabel('t')
plt.ylabel('w(t)')
plt.legend(loc='lower right')

plt.subplot(1, 2, 2)

plt.plot(sigma2, color = '#0000FF', label = u'$\sigma^2$, γ = 1, n = 100')
plt.xlabel('t')
plt.legend(loc='lower right')
plt.show()

fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(212, projection='3d')
nbins = 25

xvalues = numpy.zeros((0,nbins))
yvalues = numpy.zeros((0,nbins))
zvalues = numpy.zeros((0,nbins))

for row in range(simuls.shape[0]-2):
    x = numpy.ones(nbins) * row
    rsamples = simuls[row,:]
    hist, bins = numpy.histogram(rsamples, bins = nbins)
    y = (bins[:-1] + bins[1:])/2
  
    xvalues = numpy.vstack((xvalues,x))
    yvalues = numpy.vstack((yvalues,y))
    zvalues = numpy.vstack((zvalues,hist))

    c = ['r']

ax.plot_surface(xvalues, yvalues, zvalues, cmap=cm.coolwarm, linewidth=0, antialiased=False)

ax.set_xlabel('t')
ax.set_ylabel('w(t)', labelpad=10)
ax.set_zlabel(' ')
ax.set_zticklabels([])

ax.view_init(70, 30)
plt.title('Histogram of Values of $w_j$ (n = 100)', loc = 'left')

plt.show()





#   Figure 4: Histograms and Poissonness Plots / Edit variables to produce subplots
################################################################################
################################################################################
fig = plt.figure(figsize=(5,5))
plt.hist(simuls[50,:], bins = numpy.arange(21) - 0.5, label = 'n = 100')
plt.xticks(range(21))
plt.legend(loc='upper right')
plt.show()

theta = range(50)
xk = numpy.array(range(50))
poiscount = collections.Counter(simuls[50,:])
for i in theta:
    xk[i] = poiscount[i]    
theta = numpy.array(range(50))
theta = theta[:6]
xk = xk[:6]
y = numpy.log(xk) + numpy.log( scipy.misc.factorial(theta))
x = theta

fig = plt.figure(figsize=(5,5))
plt.scatter(x,y, label='n = 100')
plt.legend(loc='lower right')
plt.xlabel(u'θ')
plt.ylabel(u'ln ($w_j$) + ln θ!')
plt.legend(loc='lower right')
plt.show()





#   Figure 5: Interval Estimate of Y(t) over Time
################################################################################
################################################################################

fig = plt.figure(figsize=(13,8))
plt.fill_between(CNint[:,0],CNint[:,1]/CNarea,CNint[:,2]/CNarea, facecolor='r',alpha = 0.5, label = 'Case Nuove (n = 9)')
plt.fill_between(PIint[:,0],PIint[:,1]/PIarea,PIint[:,2]/PIarea, facecolor='blue',alpha = 0.5, label = 'Pievina (n = 34)')
plt.fill_between(MZint[:,0],MZint[:,1]/MZarea,MZint[:,2]/MZarea, facecolor='black',alpha = 0.4, label = 'Podere Marzuolo (n = 31)')
plt.fill_between(PTint[:,0],PTint[:,1]/PTarea,PTint[:,2]/PTarea, facecolor='orange',alpha = 0.5, label = 'Podere Terrato (n = 5)')
plt.xlabel('t')
plt.ylabel('coins / sq. m')
plt.title(u'Abundance of Coinage in Use (Roman Peasant Project, 2009-2014), γ = 0.1',  loc = 'left')
plt.legend(loc='upper left')
plt.plot()
    
