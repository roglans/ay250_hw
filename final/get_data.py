def get_data(t0=None,t1=None,sample=False,sample_rate=3960,downfactor=1,zone='Local',seg=60,rates=[],
             station=None,verbose=False,dtype='scalar',rep='/Users/vincent/ASTRO/data/NURI/',gps='on'):
    """
    Glob all files withing user-defined period and extract data.
    
    Parameters
    ----------
    t0 : str
      Starting time in format YYYY-MM-DD-HH-MM
    t1 : str
      Ending time in format YYYY-MM-DD-HH-MM
    sample : bool
      Whether we read all the available data or only relevant files, default is False
    sample_rate : float
      Data sampling rate, default is 3960 Hz
    downfactor : float
      Down-sampling rate, default is 1
    zone : str
      Define timing zone, either Local or UTC
    rates : list
      List of successive value to decimate the data from
    station : int
      Station number
    dtype : str
      Type of output value, default is the scalar magnitude among all directions.
    rep : str
      Path in which the data are stored.
    
    Return
    ------
    ts : gwpy.timeseries.TimeSeries
      Time series data, either scalar magnitude or for individual direction.
    """
    downfactor = downfactor if len(rates)==0 else reduce(lambda x, y: x*y, rates)
    data_path = rep if sample else '%s/NURI-station-%02i'%(rep,station)
    # Check whether sample mode set or start and end timestamp defined
    if sample:
        dataset = glob.glob('%s/*time*.bin'%data_path)
    elif t0!=None and t1!=None:
        # Convert start and end dates into datetime objects
        t0 = datetime(*numpy.array(t0.split('-'),dtype=int))
        t1 = datetime(*numpy.array(t1.split('-'),dtype=int))
        # Calculate Local-UTC time difference and offset to be applied to filenames
        utc2local, offset, utc2pst = get_offset(t0,t1,zone,station,gps)
        # Initialise dataset array
        dataset = []
        # Loop through all requested hours
        for date in numpy.arange(t0+timedelta(seconds=offset),t1+timedelta(seconds=offset),timedelta(hours=1)):
            # Convert date into readable datetime format
            date = date.astype(datetime)
            # Determine timing binary filename
            path = '%s/NURI-station-%02i/'%(rep,station)+date.strftime("%Y-%-m-%-d_%-H-xx_time*.bin")
            # Search for the file and add it into storing array
            dataset += glob.glob(path)
    ts,t,x,y,z = [],[],[],[],[]
    if len(dataset)>0:
        for tfile in dataset:
            # Display file index in list
            if verbose: print(tfile.split('/')[-1])
            # Extract data from timing file
            times = get_time(tfile,downfactor) if sample else get_time(tfile,downfactor,utc2local,utc2pst,zone)
            #t = numpy.hstack((t,times))
            t.extend(times)
            # Check version of the timing
            version = 'time_v2' if 'time_v2' in tfile else 'time'
            # Extract data from X direction file
            xfile = tfile.replace(version,'rawX_uT_3960Hz')
            x.extend(read_axis(xfile,downfactor,rates,seg))
            # Extract data from Y direction file
            yfile = tfile.replace(version,'rawY_uT_3960Hz')
            y.extend(read_axis(yfile,downfactor,rates,seg))
            # Extract data from Z direction file
            zfile = tfile.replace(version,'rawZ_uT_3960Hz')
            z.extend(read_axis(zfile,downfactor,rates,seg))
        # Create individual time series if single axis requested
        if dtype=='xyz':
            x = TimeSeries(x,sample_rate=sample_rate/downfactor,epoch=t[0])
            y = TimeSeries(y,sample_rate=sample_rate/downfactor,epoch=t[0])
            z = TimeSeries(z,sample_rate=sample_rate/downfactor,epoch=t[0])
            ts = (x,y,z)
        # Create full magnitude time series is scalar type selected
        if dtype=='scalar':
            ts = numpy.vstack((x,y,z)).T
            ts = numpy.sqrt(numpy.sum(abs(ts)**2,axis=1))
            ts = TimeSeries(ts,sample_rate=sample_rate/downfactor,epoch=t[0])
    return ts
