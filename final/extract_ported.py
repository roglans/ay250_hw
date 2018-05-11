def get_offset(t0,t1,zone,station,gps):
    """
    Determine UTC to local Local offset to be applied.
    
    Parameters
    ----------
    t0 : datetime
      Starting timestamp
    t1 : datetime
      End timestamp
    zone : str
      Define timing zone, either Local or UTC
    city : str
      City where the sensor is located

    Return
    ------
    offset : datetime
      Offset time to match time in targeted filename 
    """
    # Identifying the time zone
    utc_zone = tz.gettz('UTC')
    # Format input timestamp into UTC time
    utc_epoch = t0.replace(tzinfo=utc_zone)
    # Get time in local California time
    local_epoch = utc_epoch.astimezone(tz.gettz('America/Los_Angeles'))
    # Calculate offset between UTC and PST timestamps
    utc2pst = datetime.utcoffset(local_epoch).total_seconds()
    # Consider UTC to PST offset if requested time is before fix date
    utc2pst = utc2pst if t0<datetime(2017,12,7) else 0
    # Look-up table to identify station's location over time
    locations = numpy.array([[1,datetime(2015,11,1),datetime(2017,12,3),tz.gettz('America/Los_Angeles')],
                             [1,datetime(2017,12,3),datetime.max       ,tz.gettz('America/New_York')   ],
                             [2,datetime(2015,11,1),datetime.max       ,tz.gettz('America/Los_Angeles')],
                             [3,datetime(2015,11,1),datetime(2017,10,6),tz.gettz('America/Los_Angeles')],
                             [3,datetime(2017,10,6),datetime.max       ,tz.gettz('America/New_York')   ],
                             [4,datetime(2015,11,1),datetime(2017,12,3),tz.gettz('America/Los_Angeles')],
                             [4,datetime(2017,12,3),datetime.max       ,tz.gettz('America/New_York')   ]])
    # Identify the location for requested data
    for n,start,end,loc in locations:
        if n==station and start<t0<end:
            local_zone = loc
    # Identifying the time zone
    utc_zone = tz.gettz('UTC')
    # Format input timestamp into UTC time
    utc_epoch = t0.replace(tzinfo=utc_zone)
    # Get time in local California time
    local_epoch = utc_epoch.astimezone(local_zone)
    # Calculate offset between Local and UTC timestamps
    utc2local = datetime.utcoffset(local_epoch).total_seconds()
    # Check if first version of timing data
    if t1<datetime(2016,6,10):
        # Calculate offset between provided UTC to local timestamps
        offset = -utc2local if zone=='UTC' else 0
    # Check if second version of timing data
    if t0>datetime(2016,6,10):
        # Calculate offset between provided local to UTC timestamps
        offset = -utc2local if zone=='Local' and gps=='on' else 0
    return utc2local,offset,utc2pst
 
def get_time(filename,downfactor,utc2local=0,utc2pst=0,zone='Local'):
    """
    Read timing binary file and construct full timestamp array.

    Parameters
    ----------
    filename : str
      Path to the timing binary file
    downfactor : float
      Down-sampling rate, default is 1
    offset : datetime
      Offset time to match time in targeted filename 
    zone : str
      Define timing zone, either Local or UTC
    
    Return
    ------
    tgps : numpy.array
      Reconstructed full timing array
    """
    # Read binary file and store data in array
    with open(filename,'rb') as f:
        data = f.read()
    f.close()
    # Check timing version
    if 'time_v2' in filename:
        # Define the total number records (63 bytes per record)
        size = len(data)/63
        # Unpack each record as the following succession:
        data = struct.unpack('<'+'qI?Qddcdcdd'*size,data)
        # Reshape array into 2D format
        data = numpy.array(data,dtype=object).reshape((size,11))
        # Correct for applied constant UTC/PST offset
        data[:,4] += utc2pst
        # Depending on UTC argument, convert Local timestamps into UTC format
        offset = utc2local if zone=='Local' else 0
        # Create incremented sample index array
        x = numpy.array(data[:,0]+data[:,1],dtype=float)
        # Defined timestamp array
        y = numpy.array(data[:,4]+offset)
        # Fit timestamp array
        popt,pcov = numpy.polyfit(x,y,1,cov=True)
        # Create full sample index array
        xfit = numpy.arange(x[-1])
        # Construct full timestamp array
        tgps = numpy.polyval(popt,xfit)
    else:
        # Define the total number records (28 bytes per record)
        size = len(data)/28
        # Unpack each record as the following succession:
        data = struct.unpack('<'+'qiQd'*size,data)
        # Reshape array so that each row corresponds to one record
        data = numpy.reshape(data,(size,4))
        # Convert timestamp into 
        data[:,2] = [(int(i) & 0x7fffffffffffffff)/1e7 for i in data[:,2]]
        # Correct timestamp from offset between system clock and UTC time
        data[:,2]-=(datetime(1970,1,1)-datetime(1,1,1)).total_seconds()
        # Correct for applied constant UTC/PST offset
        data[:,2] += utc2pst
        # Depending on UTC argument, convert Local timestamps into UTC format
        offset = utc2local if zone=='UTC' else 0
        # Create incremented sample index array
        x = data[:,0]+data[:,1]
        # Defined timestamp array
        y = data[:,2]+offset
        # Fit timestamp array
        popt,pcov = numpy.polyfit(x,y,1,cov=True)
        # Create full sample index array
        xfit = numpy.arange(x[-1])
        # Construct full timestamp array
        tgps = numpy.polyval(popt,xfit)
    # Downsample time and data arrays
    tgps = tgps[::downfactor] if downfactor==1 else tgps[::downfactor][:-1]
    return numpy.array(tgps)

def get_data(t0=None,t1=None,sample=False,sample_rate=3960,downfactor=1,zone='Local',seg=60,rates=[],
             station=None,verbose=False,dtype='scalar',rep='/media/shared/nuridata',gps='on'):
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
        # Create full magnitude time series if scalar type selected
        if dtype=='scalar':
            ts = numpy.vstack((x,y,z)).T
            ts = numpy.sqrt(numpy.sum(abs(ts)**2,axis=1))
            ts = TimeSeries(ts,sample_rate=sample_rate/downfactor,epoch=t[0])
    return ts

def read_axis(filename,downfactor,rates,seg):
    """
    Read the binary file of magnetic field axis.

    Parameters
    ----------
    filename : str
      Path to data file
    downfactor : float
      Resampling factor to downsample the data
    rates : list
      List of successive value to decimate the data from
    seg : int
      Number of seg to split data for faster downsampling

    Return
    ------
    data : numpy.array
      Magnetic field data.
    """
    # Read binary file and sore data in array
    with open(filename,'rb') as f:
        data = f.read()
    f.close()
    # Define the total number records (63 bytes per record)
    size = len(data)/8
    # Unpack each record
    data = numpy.array(struct.unpack('d'*size,data))
    if downfactor>1:
        # Determine number of sample in 10-minute chunk of data
        n = 3960*seg
        # Resample each 10-minute chunk of data for faster processing
        chunks = []
        for i in xrange(0, len(data), n):
            if len(data)-i>downfactor:
                chunk = data[i:i+n]
                chunk = signal.resample(chunk,len(chunk)/downfactor)
                chunks.extend(chunk)
    return data if downfactor==1 else chunks

def get_coord(filename):
    """
    Extracted longitude and latitude from version 2 timing data files.

    Parameters
    ----------
    filename : str
      Name of the timing binary file.

    Return
    ------
    lat, lon : float
      Latitute and Longitude of the magnetic sensor.
    """
    with open(filename,'rb') as f:
        data = f.read()
    f.close()
    s = len(data)/63
    data = struct.unpack('<'+'qipQddcdcdd'*s,data)
    data = numpy.reshape(data,(s,11))
    lat  = float(str(float(data[0,5])/100).split('.')[0])
    lat  = lat + (float(data[0,5])-lat*100)/60
    lat  = -lat if data[0,8]=='S' else lat
    lon  = float(str(float(data[0,7])/100).split('.')[0])
    lon  = lon + (float(data[0,7])-lon*100)/60
    lon  = -lon if data[0,6]=='W' else lon
    return lat,lon

def get_gnome(station,t0,t1,local=False,
              rep='/Users/vincent/ASTRO/data/NURI/'):
    """
    Glob all files withing user-defined period and extract data.

    Parameters
    ----------
    station : str
      Name of the station to be analysed
    t0 : str
      Starting time in format YYYY-MM-DD-HH-MM
    t1 : str
      Ending time in format YYYY-MM-DD-HH-MM
    local : bool
      Whether we extract local (GPS) or universal (UTC) time.
    rep : str
      Path in which the data are stored.

    Return
    ------
    data : numpy.array
      2D magnetic field data array.
    """
    print('Searching GNOME data files...')
    # Convert start and end dates into datetime objects
    t0 = datetime(*numpy.array(t0.split('-'),dtype=int))
    t1 = datetime(*numpy.array(t1.split('-'),dtype=int))
    dataset = []
    for date in numpy.arange(start,end,timedelta(minutes=1)):
        date = date.astype(datetime)
        fullpath = rep+'/'+station+'_'+date.strftime("%Y%m%d_%H%M*.h5")
        dataset += glob.glob(fullpath)
    if sample:
        dataset = glob.glob(rep+"/*.h5")
    gps_offset = (datetime(1980,1,6)-datetime(1970,1,1)).total_seconds()
    tdata,xdata = [],[]
    for fname in sorted(dataset):
        hfile = h5py.File(fname, "r")
        dset = hfile['MagneticFields']
        datestr = dset.attrs["Date"]
        t0str = dset.attrs["t0"]
        # Format date from extracted metadata
        instr = "%d-%d-%02dT" % tuple(map(int, datestr.split('/'))) + t0str
        # Calculate offset between UTC and GPS initial epoch time
        gps_offset = (datetime(1980,1,6)-datetime(1970,1,1)).total_seconds()
        # Calculate GPS epoch
        gps_epoch = astropy.time.Time(instr, format='isot', scale='utc').gps + gps_offset
        utc,gps     = datetime.utcfromtimestamp(gps_epoch),\
                      datetime.fromtimestamp(gps_epoch)
        delay       = round(abs((utc-gps).total_seconds())/3600.)
        start_time  = gps if local else utc
        start_time  = time.mktime(start_time.timetuple())
        sample_rate = dset.attrs["SamplingRate(Hz)"]
        end_time    = start_time + len(dset[:]) / sample_rate
        tarray      = numpy.linspace(start_time,end_time,len(dset[:]))
        t0 = time.mktime(t0.timetuple())+t0.microsecond/1e6
        t1 = time.mktime(t1.timetuple())+t1.microsecond/1e6
        idx = numpy.where(numpy.logical_and(t0-1<tarray,tarray<t1+1))[0]
        xdata       = numpy.hstack((xdata,dset[:][idx]))
        tdata       = numpy.hstack((tdata,tarray[idx]))
        hfile.close()
    data = numpy.vstack((tdata,10*xdata)).T
    return data

def read_time(t0,t1,sample=False,rep='/Users/vincent/ASTRO/data/NURI/'):
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
    rep : str
      Path in which the data are stored.
   
    Return
    ------
    tdate : numpy.array
      Extracted timing information
    """    
    t0=t1=None
    if sample:
        dataset = glob.glob(rep+"/*time*.bin")
    elif t0!=None and t1!=None:
        dataset = []
        # Convert start and end dates into datetime objects
        t0 = datetime(*numpy.array(t0.split('-'),dtype=int))
        t1 = datetime(*numpy.array(t1.split('-'),dtype=int))
        print(t0,'to',t1)
        for date in numpy.arange(t0,t1,timedelta(hours=1)):
            date = date.astype(datetime)
            path1 = rep+'/'+date.strftime("%Y-%-m-%-d_%-H-xx_time.bin")
            path2 = rep+'/'+date.strftime("%Y-%-m-%-d_%-H-xx_time_v2.bin")
            dataset += glob.glob(path1)
            dataset += glob.glob(path2)
    tdata = []
    for tfile in dataset:
        with open(tfile,'rb') as f:
            data = f.read()
        f.close()
        if 'time_v2' in tfile:
            size = len(data)/63
            data = struct.unpack('<'+'qI?Qddcdcdd'*size,data)
            data = numpy.array(data,dtype=object).reshape((size,11))
            data[:,3] = [float(i)/2533200 for i in data[:,3]]
        else:
            size = len(data)/28
            data = struct.unpack('<'+'qiQd'*size,data)
            data = numpy.reshape(data,(size,4))
            data[:,2] = [(int(i) & 0x7fffffffffffffff)/1e7 for i in data[:,2]]
        tdata = data if len(tdata)==0 else numpy.vstack((tdata,data))
    return tdata

def get_tlim(*times):
    """
    Create mask for user-defined time range.
    """
    time_list = []
    for dt in times:
        # Convert start date into timestamp
        dt = time.mktime(dt.timetuple())+dt.microsecond/1e6
        # Convert timestamp back local datetime
        now = datetime.fromtimestamp(dt)
        # Convert timestamp to UTC datetime
        utc_now = datetime.utcfromtimestamp(dt)
        # Calculate local difference with UTC time
        utc2local = (now-utc_now).total_seconds()
        # Remove UTC offset from timestamp
        dt += utc2local
        time_list.append(dt)
    return tuple(time_list)