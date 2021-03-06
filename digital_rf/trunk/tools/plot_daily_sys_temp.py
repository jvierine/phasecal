"""plot_daily_sys_temp.py plots multiple 24 periods of sys temp

Uses digital_rf_metadata.py

$Id: plot_daily_sys_temp.py 825 2015-11-05 18:36:39Z brideout $
"""

# standard python imports
import sys
import datetime, calendar
import time
import collections

# third party imports
import matplotlib
# set rendering
matplotlib.use('Agg')
import matplotlib.pyplot
import numpy

# Millstone imports
import digital_rf_metadata


usage = """python plot_daily_sys_temp.py <startYYYY-MM-DD> <endYYYY-MM-DD> <outputImageFile>
"""


if len(sys.argv) != 4:
    print(usage)
    sys.exit(-1)
    
t = time.time()
startDT = datetime.datetime.strptime(sys.argv[1], '%Y-%m-%d')
endDT = datetime.datetime.strptime(sys.argv[2], '%Y-%m-%d')
endDT = datetime.datetime(endDT.year, endDT.month, endDT.day, 23, 59, 59)
if startDT > endDT:
    raise ValueError, 'start %s after end %s' % (sys.argv[1], sys.argv[2])
outputImageFile = sys.argv[3]

# gather data to plot if form of ordered dict , where key = dt, and value is a tuple
# of three lists 1. seonds since UT midnight, 2. temp_misa, and 3. temp_zenith
# returns for that day
readMetaObj = digital_rf_metadata.read_digital_rf_metadata('/data0/results/rxnoise')
dataDict = collections.OrderedDict()
while startDT < endDT:
    print('working on day %s' % (str(startDT)))
    thisEndDT = datetime.datetime(startDT.year, startDT.month, startDT.day, 23, 59, 59)
    sample0 = calendar.timegm(startDT.timetuple())
    sample1 = calendar.timegm(thisEndDT.timetuple())
    thisDict = readMetaObj.read_metadata(sample0, sample1, ('temp_misa_median', 'temp_zenith_median'))
    secList = numpy.array([key % (3600*24) for key in thisDict.keys()])
    misaTempList = numpy.array([thisDict[key]['temp_misa_median'] for key in thisDict.keys()])
    zenithTempList = numpy.array([thisDict[key]['temp_zenith_median'] for key in thisDict.keys()])
    dataDict[startDT] = (secList, misaTempList, zenithTempList)
    startDT += datetime.timedelta(days=1)
print('data access took %f seconds' % (time.time() - t))
    
matplotlib.pyplot.figure(figsize=(6,12))
matplotlib.pyplot.subplot(2, 1, 1)
for i, key in enumerate(dataDict.keys()):
    label = key.strftime('Z: %Y-%m-%d')
    matplotlib.pyplot.plot(dataDict[key][0]/3600.0, dataDict[key][2]+(i*25), label=label)
matplotlib.pyplot.ylim(100, 900)
matplotlib.pyplot.xlim(0,24)
matplotlib.pyplot.legend(fontsize=8)
matplotlib.pyplot.title('Zenith Rx (K)')
matplotlib.pyplot.subplot(2, 1, 2)
for i, key in enumerate(dataDict.keys()):
    label = key.strftime('M: %Y-%m-%d')
    matplotlib.pyplot.plot(dataDict[key][0]/3600.0, dataDict[key][1]+(i*25), label=label)
matplotlib.pyplot.ylim(100, 900)
matplotlib.pyplot.xlim(0,24)
matplotlib.pyplot.legend(fontsize=8)
matplotlib.pyplot.title('Misa Rx (K)')
matplotlib.pyplot.savefig(outputImageFile)
print('plot took %f seconds' % (time.time() - t))



