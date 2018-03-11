import csv
import numpy
import random
#import matplotlib.pyplot as plt
#import scipy.stats as stats    


#dates = numpy.genfromtxt(open("inputdates.csv", "rb"), delimiter=",", names=True)
#
#quant = numpy.genfromtxt(open("inputcount.csv", "rb"), delimiter=",", names=True)
##context = numpy.genfromtxt(open("context.csv", "rb"), delimiter=",", names=True)
##contextresult = numpy.genfromtxt(open("contextresult.csv", "rb"), delimiter=",", names=True)
#
#reader = csv.reader(open("inputdomain.csv", "r"), delimiter=",")
#x = list(reader)
#
#cdomain = x[1:]
#
#reader = csv.reader(open("inputcontext.csv", "r"), delimiter=",")
#y = list(reader)
#
#context2 = y[1:]

reader = csv.reader(open("outputcontext.csv", "r"), delimiter=",")
z = list(reader)

contextresult = z[1:]

samplesize = int(input('Number of Monte Carlo samples? '))

#nsample = quant.shape[0]

startdate = -200
enddate = 500
daterange = range(startdate,enddate)

reader = csv.reader(open("input.csv", "r"), delimiter=",")
rawtable = list(reader)


##################################################################################
###interval estimate

betaminusalpha = 0

assemblagesize = [36, 33, 5, 9]
sitearea = [550.49, 764.63, 636.65, 995.29]

for site in contextresult:
    
    n = assemblagesize[contextresult.index(site)]
    area = sitearea[contextresult.index(site)]
    
    resultstable = numpy.zeros((len(daterange),2))
    resultstable[:,0] = daterange
    
    finalresults = numpy.zeros((len(daterange),5))
    finalresults[:,0] = daterange #insertdate
    
    dates = numpy.zeros([0,2])
    quant = numpy.zeros([0,1])
    
    for row in range(len(rawtable)):
        if rawtable[row][0] == site[0]:
            a = int(rawtable[row][3])
            b = int(rawtable[row][4])
            dates = numpy.vstack((dates,(a,b)))
            quantint = int(rawtable[row][2])
            quant = numpy.vstack((quant,(quantint)))
    
    results2 = numpy.zeros([0,9])
    for MCcount in range(samplesize):          
        rowcount = 0
        quant2 = numpy.empty([0,1])
        dates2 = numpy.empty([0,2])
        ceramicsdomain = []
        context = []    
        
        for row in dates:
            a = row[0]
            b = row[1]
            if a < 5000:
                if b < 5000:
                    xquant = quant[rowcount][0]
                    while xquant > 0:
                            rand1 = random.randint(a,b)
                            rand2 = random.randint(a,b)
                            if rand1 < rand2:
                                alpha = rand1
                                beta = rand1 + numpy.random.geometric(0.5,1) #rand2
                                if beta > b:
                                    beta = b
                            if rand2 < rand1:
                                alpha = rand2
                                beta = rand2 + numpy.random.geometric(0.5,1)#rand1
                                if beta > b:
                                    beta = b
                            if rand1 == rand2:
                                alpha = rand1
                                beta = rand2 + 1
                            dates2 = numpy.vstack((dates2,numpy.array([alpha,beta])))
                            quant2 = numpy.vstack((quant2,numpy.array([1])))
                            #ceramicsdomain = ceramicsdomain + domain[rowcount]
                            #context = context + context2[rowcount]
                            xquant = xquant - 1
            rowcount = rowcount + 1
        rowcount = 0
        betaminusalpha = 0
        for date in dates2: 
            alpha = date[0]
            beta = date[1]
            for j in resultstable[rowcount]:
                start = numpy.where(resultstable[:,0] == alpha)[0][0]
                end = numpy.where(resultstable[:,0] == beta)[0][0]
                resultstable[start:end,1] = resultstable[start:end,1] + 1 #this adds one to each row, adding one throughout the entire interval
            betaminusalpha = betaminusalpha + abs(beta-alpha) #keeps the running sum of intervals beta/alpha
        newcol = (resultstable[:,1]/betaminusalpha) * (n / area)
        newcol2 = numpy.zeros((len(daterange),1))
        newcol2[:,0] = newcol


    mu = newcol2[:,0] #this is U(t): insert into table
    finalresults[:,1] = mu

#################################################33
#interval estimation through subsampling

    resultstable = numpy.zeros((len(daterange),2))
    resultstable[:,0] = daterange

    Rsubsampleset = numpy.zeros((len(daterange),2))
    Rsubsampleset[:,0] = daterange

###############################
###############################
    subsamplingroutine = 1000

    rmeancollected = numpy.zeros((len(daterange),0))   
    resultscollection = numpy.zeros((len(daterange),0))
    
    for r in range(subsamplingroutine): 
        resultscollection = numpy.zeros((len(daterange),0))
        samplesizeset = random.randint(1,dates.shape[0])
        randrows = numpy.random.choice(range(dates.shape[0]), samplesizeset, replace = False) #to filter out a subsample of rows
        #resultscollection = numpy.zeros((len(daterange),0)) #collect the MC estimates in cols to calculate mean and other parameters
        #print('Subsample size:')
        #print(nhatsample)
        resultstable[:,1] = 0
        betaminusalpha = 0
        for MCcount in range(samplesize):
        #generate a list of quantities - sherds (1) and corresponding random dates in dates2
        #with corresponding list of sites and domains
            rowcount = 0
            quant2 = numpy.empty([0,1])
            dates2 = numpy.empty([0,2])
            ceramicsdomain = []
            context = []
            for row in dates:
                if rowcount in randrows: #selecting from the subsample
                    a = row[0]
                    b = row[1]
                    if a < 5000:
                        if b < 5000:
                            xquant = quant[rowcount][0]
                            while xquant > 0:
                                rand1 = random.randint(a,b)
                                rand2 = random.randint(a,b)
                                if rand1 < rand2:
                                    alpha = rand1
                                    beta = rand1 + numpy.random.geometric(0.5,1) #rand2
                                    if beta > b:
                                        beta = b
                                if rand2 < rand1:
                                    alpha = rand2
                                    beta = rand2 + numpy.random.geometric(0.5,1)#rand1
                                    if beta > b:
                                        beta = b
                                if rand1 == rand2:
                                    alpha = rand1
                                    beta = rand2 + 1
                                dates2 = numpy.vstack((dates2,numpy.array([alpha,beta])))
                                quant2 = numpy.vstack((quant2,numpy.array([1])))
                                #ceramicsdomain = ceramicsdomain + domain[rowcount]
                                #context = context + context2[rowcount]
                                xquant = xquant - 1
                rowcount = rowcount + 1
            rowcount = 0
            betaminusalpha = 0
            for date in dates2: 
                alpha = date[0]
                beta = date[1]
                for j in resultstable[rowcount]:
                    start = numpy.where(resultstable[:,0] == alpha)[0][0]
                    end = numpy.where(resultstable[:,0] == beta)[0][0]
                    resultstable[start:end,1] = resultstable[start:end,1] + 1
                betaminusalpha = betaminusalpha + abs(beta-alpha)
            newcol = (resultstable[:,1]/betaminusalpha) * (n / area)
            newcol2 = numpy.zeros((len(daterange),1))
            newcol2[:,0] = newcol

        muq = newcol2 #this is U(t) according to the subsample Xq

        rmeancollected = numpy.hstack((rmeancollected,muq))

        if r/float(100) == int(r/float(100)):
            print(r) 


    for j in range(rmeancollected.shape[0]):      
        tempdata = rmeancollected[j,:]
        #tempdata = numpy.zeros((0,1)) #to remove zero values
        #for i in tempdata2:
        #    if i > 0:
        #        tempdata = numpy.vstack((tempdata,i)) #remove all zero values
        #if sum(tempdata) == 0:
        #    tempdata = numpy.zeros((1,1))
        if sum(tempdata) != 0:
            #tempdata =  tempdata[tempdata != 0]
            #tempdata = numpy.log(tempdata)
            musub = numpy.mean(tempdata)
            yset = []
            for y in tempdata:
                yset = yset + [(y - musub)**2]
            var = sum(yset) / float(tempdata.shape[0]-1)
            alpha = musub**2 / var
            beta =  var / musub 
            #var2 = alpha * (beta**2)
            #print('var1: {0}, {1}'.format(var, var2))
            
            #x1 = numpy.linspace(0,2,1000)
            #shape, loc, scale = stats.lognorm.fit(tempdata)
            #pdf = stats.lognorm.pdf(x1, shape, loc=loc, scale=scale)
            
            #fig = plt.figure()
            #plt.hist(tempdata, bins = 500, range=(0,2))
            #plt.plot(x1, pdf)#, color = '#3498db')
            #plt.show()

            finalresults[:,1] = mu        
            finalresults[j,2] = musub
            finalresults[j,3] = alpha
            finalresults[j,4] = beta

    numpy.savetxt('finalresults-{0}.csv'.format(site), finalresults, fmt='%5s', delimiter=",") #fmt="%10.50e",
    numpy.savetxt('rmeancollected-{0}.csv'.format(site), rmeancollected, fmt='%5s', delimiter=",") #fmt="%10.50e",



