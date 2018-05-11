#!/usr/bin/env python
import sys,nuri,os,numpy
import matplotlib.pyplot as plt
import matplotlib.dates as md
from datetime import datetime,timedelta

def check24hrs(date):
    """
    This operation will display the active periods for which data are
    available from every sensors.

    Parameters
    ----------

    date : str
      Year and month to display activity from. The format shoud be YYYY-MM.
    """
    # Check if date is given
    if date==None:
        print('date missing...')
        quit()
    # List all the months
    dates = numpy.empty((0,5))
    y0    = int(date.split('-')[0])
    m0    = int(date.split('-')[1])
    d0    = datetime(y0,m0,1)
    y1    = y0   if m0<12 else y0+1
    m1    = m0+1 if m0<12 else 1
    d1    = datetime(y1,m1,1)-timedelta(hours=1)
    dt    = timedelta(hours=1)
    dates = numpy.arange(d0,d1,dt)
    # Download metadata from Google Drive
    sys.stderr.write('Retrieve information from Google Drive...')
    os.system('skicka ls -r /MagneticFieldData/%s/%s/ > data'%(y0,m0))
    data = numpy.loadtxt('data',dtype=str,delimiter='\n')
    print >>sys.stderr,' done!'
    # List file path for each date and each station
    sys.stderr.write('Select active hours for each station...')
    st0,st1,st2,st3,st4 = [],[],[],[],[]
    for d in dates:
        year  = d.astype(object).year
        month = d.astype(object).month
        day   = d.astype(object).day
        hour  = d.astype(object).hour
        path  = 'MagneticFieldData/%i/%i/%i/%i/'%(year,month,day,hour)
        fname = '%i-%i-%i_%i-xx.zip'%(year,month,day,hour)
        st0.append(path+'NURI-station/'   +fname)
        st1.append(path+'NURI-station-01/'+fname)
        st2.append(path+'NURI-station-02/'+fname)
        st3.append(path+'NURI-station-03/'+fname)
        st4.append(path+'NURI-station-04/'+fname)
    st0 = numpy.array([1 if path in data else 0 for path in st0])
    st1 = numpy.array([1 if path in data else 0 for path in st1])
    st2 = numpy.array([1 if path in data else 0 for path in st2])
    st3 = numpy.array([1 if path in data else 0 for path in st3])
    st4 = numpy.array([1 if path in data else 0 for path in st4])
    print >>sys.stderr,' done!'
    # Write down information in text file
    print 'Save information in ASCII file...'
    o = open('%i-%02i.dat'%(y0,m0),'w')
    for d in dates:
        year  = d.astype(object).year
        month = d.astype(object).month
        day   = d.astype(object).day
        hour  = d.astype(object).hour
        path  = 'MagneticFieldData/%i/%i/%i/%i/'%(year,month,day,hour)
        fname = '%i-%i-%i_%i-xx.zip'%(year,month,day,hour)
        o.write('%i-%02i-%02i_%02i'%(year,month,day,hour))
        o.write('  NURI-station   ') if path+'NURI-station/'   +fname in data else o.write('  -              ')
        o.write('  NURI-station-01') if path+'NURI-station-01/'+fname in data else o.write('  -              ')
        o.write('  NURI-station-02') if path+'NURI-station-01/'+fname in data else o.write('  -              ')
        o.write('  NURI-station-03') if path+'NURI-station-01/'+fname in data else o.write('  -              ')
        o.write('  NURI-station-04') if path+'NURI-station-01/'+fname in data else o.write('  -              ')
        o.write('\n')
    o.close()
    dates = [d.astype(object) for d in dates]
    plt.rc('font', size=2, family='serif')
    plt.rc('axes', labelsize=10, linewidth=0.2)
    plt.rc('legend', fontsize=2, handlelength=10)
    plt.rc('xtick', labelsize=7)
    plt.rc('ytick', labelsize=7)
    plt.rc('lines', lw=0.2, mew=0.2)
    plt.rc('grid', linewidth=0.2)
    fig = plt.figure(figsize=(10,6))
    plt.subplots_adjust(left=0.07, right=0.95, bottom=0.1, top=0.96, hspace=0.2, wspace=0)
    print 'Plot active time for station 1...'
    ax1 = fig.add_subplot(511)
    ax1.bar(dates,st1,width=0.01,edgecolor='none',color='green')
    ax1.tick_params(direction='in')
    ax1.set_ylabel('Station 1')
    ax1.xaxis_date()
    plt.yticks([])
    ax1.grid()
    print 'Plot active time for station 2...'
    ax = fig.add_subplot(512,sharex=ax1,sharey=ax1)
    ax.bar(dates,st2,width=0.01,edgecolor='none',color='green')
    ax.tick_params(direction='in')
    ax.set_ylabel('Station 2')
    ax.xaxis_date()
    plt.yticks([])
    ax.grid()
    print 'Plot active time for station 3...'
    ax = fig.add_subplot(513,sharex=ax1,sharey=ax1)
    ax.bar(dates,st3,width=0.01,edgecolor='none',color='green')
    ax.tick_params(direction='in')
    ax.set_ylabel('Station 3')
    ax.xaxis_date()
    plt.yticks([])
    ax.grid()
    print 'Plot active time for station 4...'
    ax = fig.add_subplot(514,sharex=ax1,sharey=ax1)
    ax.bar(dates,st4,width=0.01,edgecolor='none',color='green')
    ax.tick_params(direction='in')
    ax.set_ylabel('Station 4')
    ax.xaxis_date()
    plt.yticks([])
    ax.grid()
    print 'Plot active time for station 0...'
    ax = fig.add_subplot(515,sharex=ax1,sharey=ax1)
    ax.bar(dates,st0,width=0.01,edgecolor='none',color='green')
    ax.tick_params(direction='in')
    ax.set_ylabel('Station 0')
    ax.xaxis_date()
    plt.yticks([])
    ax.grid()
    ax.set_xlabel(r'Hourly activity in %s %i (UTC)'%(d0.strftime("%B"),y0))
    ax1.xaxis.set_major_formatter(md.DateFormatter('%d'))
    ax1.xaxis.set_major_locator(md.DayLocator())
    ax1.set_xlim(d0,d1)
    ax1.set_ylim(0,1)
    plt.savefig('%i-%02i.pdf'%(y0,m0),dpi=80)
