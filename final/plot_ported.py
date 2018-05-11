import os,math,numpy,matplotlib,mlpy
import matplotlib.pyplot as plt
from datetime          import datetime
from gwpy.plotter      import SpectrogramPlot
from gwpy.timeseries   import TimeSeries
from matplotlib.ticker import LogLocator
from scipy             import fftpack,signal
from .extract          import get_tlim

def plot_spectrogram(data,tmin=None,tmax=None,fmin=None,fmax=None,vmin=None,vmax=None,
                     mode='wavelet',omega0=6,dt=1,dj=0.05,fct='morlet',
                     stride=None,nfft=None,overlap=None,scale='log',
                     funit='Hz',tunit='secs',cmap='inferno',zone='Local',fname=None):
    """
    Plot multiplot figure with time series, PSD and spectrogram.

    Parameters
    ----------
    data : TimeSeries
      Magnetic field data
    tmin, tmax : datetime
      First and last timestamps
    fmin, fmax : float
      Minimum and maximum frequencies
    vmin, vmax : float
      Minimum and maximum color values
    mode : str
      Spectrogram mode, wavelet or Fourier. Default is Fourier
    omega0 : int
      Wavelet function parameter
    dt : float
      Time step
    dj : float
      Scale resolution (smaller values of dj give finer resolution)
    fct : str
      Wavelet function (morlet,paul,dog)
    stride : float
      Length of segment
    nfft : float
      Length of FFT
    overlap : float
      Length of overlapping segment
    cmap : str
      Colormap
    scale : str
      Plotted frequency scale. Default is "log".
    funit : strg
      Frequency unit, Hz or mHz. Default is Hz.
    tunit : str
      Time unit, secs, mins or hrs. Default is mins.
    fname : str
      Output file name.
    
    Notes
    -----
    The `matplotlib.pyplot.imshow <https://matplotlib.org/api/pyplot_api.html?highlight=matplotlib%20pyplot%20imshow#matplotlib.pyplot.imshow>`_ module is
    used to plot the wavelet spectrogram. This module is usually used
    to plot raw images and assumes that the position of the cell in the
    input spectrogram array directly represents the position of the pixel
    in the raw image. That is, for an input Python array (in which rows
    are appended below previous ones), the first row in the array is
    assumed to represent the top line of pixel in the image. Therefore,
    in order to plot the spectrogram array using the imshow module, one
    needs to carefully check that the rows (which are representative of
    the frequency bands), are stored in descending order such that the
    lowest frequency is placed at the end (bottom) of the array.
    """
    if mode=='wavelet' and scale=='linear':
        print 'Warning: Wavelet mode chosen. Scale will be changed to log.'
        scale = 'log'
    # Initialise figure
    fig = plt.figure(figsize=(24,14),frameon=False)
    plt.subplots_adjust(left=0.07, right=0.95, bottom=0.1, top=0.95, hspace=0, wspace=0)
    ax1 = fig.add_axes([0.20,0.75,0.683,0.20])
    ax2 = fig.add_axes([0.20,0.10,0.683,0.64], sharex=ax1)
    ax3 = fig.add_axes([0.07,0.10,0.123,0.64])
    ax4 = fig.add_axes([0.89,0.10,0.030,0.64])
    # Prepare timing range
    tmin = data.times[0].value  if tmin==None else tmin
    tmax = data.times[-1].value if tmax==None else tmax
    mask = (data.times.value>=tmin) & (data.times.value<=tmax)
    scale_factor = 3600. if tunit=='hrs' else 60. if tunit=='mins' else 1
    times = (data[mask].times.value-tmin)/scale_factor
    # Plot time series
    ax1.plot(times,data[mask].value,alpha=0.5)
    ax1.set_ylabel('Magnetic Fields [uT]',fontsize=11)
    ax1.tick_params(bottom='off',labelbottom='off')
    ax1.set_xlim(0,(tmax-tmin)/scale_factor)
    ax1.grid(b=True, which='major', alpha=0.7, ls='--')
    if mode=='wavelet':
        # Calculate wavelet parameters
        scales = mlpy.wavelet.autoscales(N=len(data[mask].value),dt=dt,dj=dj,wf=fct,p=omega0)
        spec = mlpy.wavelet.cwt(data[mask].value,dt=dt,scales=scales,wf=fct,p=omega0)
        freq = (omega0 + numpy.sqrt(2.0 + omega0 ** 2)) / (4 * numpy.pi * scales[1:])
        freq = freq * 1000. if funit=='mHz' else freq
        spec = numpy.abs(spec)**2
        spec = spec[::-1]
        # Define minimum and maximum frequencies
        fmin_log,fmax_log = min(freq),max(freq)
        fmin_linear,fmax_linear = min(freq),max(freq)
        if fmin!=None:
            log_ratio = (numpy.log10(fmin)-numpy.log10(min(freq)))/(numpy.log10(max(freq))-numpy.log10(min(freq)))
            fmin_linear = min(freq)+log_ratio*(max(freq)-min(freq))
            fmin_log = fmin
        if fmax!=None:
            log_ratio = (numpy.log10(fmax)-numpy.log10(min(freq)))/(numpy.log10(max(freq))-numpy.log10(min(freq)))
            fmax_linear = min(freq)+log_ratio*(max(freq)-min(freq))
            fmax_log = fmax
        # Get minimum and maximum amplitude in selected frequency range
        idx = numpy.where(numpy.logical_and(fmin_log<freq[::-1],freq[::-1]<fmax_log))[0]
        vmin = vmin if vmin!=None else numpy.sort(numpy.unique(spec[idx]))[1]
        vmax = spec[idx].max() if vmax==None else vmax
        # Plot spectrogram
        img = ax2.imshow(spec,extent=[times[0],times[-1],freq[-1],freq[0]],aspect='auto',
                        interpolation='nearest',cmap=cmap,norm=matplotlib.colors.LogNorm(vmin,vmax)) 
        ax2.set_xlabel('Time [%s] from %s %s (%s)'%(tunit,datetime.utcfromtimestamp(tmin),zone,tmin),fontsize=15)
        ax2.set_xlim(0,(tmax-tmin)/scale_factor)
        ax2.set_yscale('linear')
        ax2.set_ylim(fmin_linear,fmax_linear)
        ax2.grid(False)
        # Set up axis range for spectrogram
        twin_ax = ax2.twinx()
        twin_ax.set_yscale('log')
        twin_ax.set_xlim(0,(tmax-tmin)/scale_factor)
        twin_ax.set_ylim(fmin_log,fmax_log)
        twin_ax.spines['top'].set_visible(False)
        twin_ax.spines['right'].set_visible(False)
        twin_ax.spines['bottom'].set_visible(False)
        ax2.tick_params(which='both', labelleft=False, left=False)
        twin_ax.tick_params(which='both', labelleft=False,left=False, labelright=False, right=False)
        twin_ax.grid(False)
    if mode=='fourier':
        freq, times, spec = signal.spectrogram(data[mask],fs=data.sample_rate.value,
                                           nperseg=stride,noverlap=overlap,nfft=nfft)
        # Convert time array into minute unit
        times = (numpy.linspace(data[mask].times.value[0],data[mask].times.value[-1],len(times))-tmin)/scale_factor
        # Define minimum and maximum frequencies
        freq = freq * 1000. if funit=='mHz' else freq
        fmin = freq[1]      if fmin==None    else fmin
        fmax = max(freq)    if fmax==None    else fmax
        fmin_log,fmax_log = fmin,fmax
        # Get minimum and maximum amplitude in selected frequency range
        idx = numpy.where(numpy.logical_and(fmin<=freq,freq<=fmax))[0]
        vmin = vmin if vmin!=None else numpy.sort(numpy.unique(spec[idx]))[1]
        vmax = spec[idx].max() if vmax==None else vmax
        # Plot spectrogram
        img = ax2.pcolormesh(times,freq,spec,cmap=cmap,norm=matplotlib.colors.LogNorm(vmin,vmax))
        ax2.set_xlabel('Time [%s] from %s %s (%s)'%(tunit,datetime.utcfromtimestamp(tmin),zone,tmin),fontsize=15)
        ax2.set_xlim(0,(tmax-tmin)/scale_factor)
        ax2.set_ylim(fmin,fmax)
        ax2.set_yscale(scale)
        ax2.set_ylabel('Frequency [%s]'%funit,fontsize=15,labelpad=40)
        ax2.tick_params(which='both', labelleft=False, left=False)
        ax2.grid(False)
    # Calculate Power Spectral Density
    N = len(data[mask].value)
    delta_t = 1/data.sample_rate.value
    delta_f = 1. / (N * delta_t)
    f = delta_f * numpy.arange(N / 2)
    f = f * 1000. if funit=='mHz' else f
    PSD = abs(delta_t * fftpack.fft(data[mask].value)[:N / 2]) ** 2
    psd = numpy.vstack((f,PSD)).T
    # Plot Power Spectral Density
    ticks = matplotlib.ticker.FuncFormatter(lambda v,_:("$10^{%.0f}$"%math.log(v,10)))
    ax3.loglog(psd[:,1],psd[:,0],alpha=0.5)
    ax3.invert_xaxis()
    ax3.set_ylim(fmin_log,fmax_log)
    ax3.set_ylabel('Frequency [%s]'%funit,fontsize=15)
    ax3.set_xlabel('PSD',fontsize=15)
    ax3.grid(b=True, which='major', alpha=0.7, ls='--')
    # Add color bar and save figure
    cb = fig.colorbar(img,cax=ax4)
    cb.set_ticks(LogLocator())
    cb.set_clim(vmin,vmax)
    ax4.set_ylabel('Power $|\mathrm{W}_v|^2$ $[\mu T^2/\mathrm{Hz}]$',fontsize=15)
    plt.show() if fname==None else plt.savefig(fname,frameon=False)
    plt.close(fig)

def plot_time_series(data,tmin=None,tmax=None,ymin=None,ymax=None,regions=[],
                     fname='time_series',zone='Local',tunit='secs',tbs=False):
    """
    Plot time series.

    Parameters
    ----------
    data : gwpy.timeseries.TimeSeries
      Time series data
    tmin : datetime
      First timestamp
    tmax : datetime
      Last timestamp
    fname : str
      Filename
    zone : str
      Output time zone, either UTC or Local
    tunit : str
      Scale of time axis (hrs, mins, secs)
    tbs : bool
      To Be Saved flag
    regions : list
      List of time regions to colorize
    """
    # Define timestamps if not defined by user
    tmin = data.times[0].value  if tmin==None else tmin
    tmax = data.times[-1].value if tmax==None else tmax
    # Re-determine the mask in case timestamps are defined
    mask = (data.times.value>=tmin) & (data.times.value<=tmax)
    # Estimate scale factor for time axis
    scale_factor = 3600. if tunit=='hrs' else 60. if tunit=='mins' else 1
    # Do the plotting
    fig = plt.figure(figsize=(12,5))
    plt.plot((data.times.value[mask]-tmin)/scale_factor,data.value[mask])
    if len(regions)>0:
        for text,tbeg,span in regions:
            tbeg = (get_tlim(tbeg)-tmin)/scale_factor
            tend = tbeg+span/scale_factor
            plt.axvspan(tbeg,tend,color='red',ls='dotted',alpha=0.1)
            plt.text(tbeg+(tend-tbeg)/2,max(data.value[mask]),'%s'%text,fontsize=25,color='red',ha='center')
    plt.xlabel('Time [%s] from %s %s (%i)'%(tunit,datetime.utcfromtimestamp(tmin),zone,tmin))
    plt.ylabel('Magnetic Field [$\mathrm{\mu}$T]')
    plt.xlim(0,(tmax-tmin)/scale_factor)
    if ymin!=None: plt.ylim(ymin=ymin)
    if ymax!=None: plt.ylim(ymax=ymax)
    plt.tight_layout()
    plt.savefig('%s.png'%fname,frameon=False,transparent=True) if tbs else plt.show()
    plt.close(fig)

def ts_movie(data,tmin=None,tmax=None,fname='time_series',zone='Local',tunit='secs'):
    """
    Create time series movie

    Parameters
    ----------
    data : gwpy.timeseries.TimeSeries
      Time series data
    tmin : datetime
      First timestamp
    tmax : datetime
      Last timestamp
    fname : str
      Filename
    zone : str
      Output time zone, either UTC or Local
    tunit : str
      Scale of time axis (hrs, mins, secs)
    """
    t0,t1 = data.times.value[0],data.times.value[-1]
    tmin,tmax = (t0,t1) if tmin==tmax==None else get_tlim(tmin,tmax)
    for time in numpy.arange(tmin,tmax,600):
        t0 = datetime.utcfromtimestamp(time)
        t1 = datetime.utcfromtimestamp(time+600)
        plot_time_series(data,t0,t1,fname='video_%s'%time,zone='Local',tunit='secs')
    os.system('convert -delay 60 -loop 0 -dispose previous video_*.png movie.gif')
    os.system('rm video_*.png')

def plot_wavelet(data,tmin=None,tmax=None,fmin=None,fmax=None,vmin=None,vmax=None,omega0=6,dt=1,dj=0.05,
                 fct='morlet',funit='Hz',cmap='inferno',zone='Local',tunit='secs',fname='wavelet',tbs=False):
    """
    Create wavelet spectrogram

    Parameters
    ----------
    data : gwpy.timeseries.TimeSeries
      Time series data
    tmin : datetime
      First timestamp
    tmax : datetime
      Last timestamp
    omega0 : int
      Wavelet function parameter
    dt : float
      Time step
    dj : float
      Scale resolution (smaller values of dj give finer resolution)
    fct : str
      Wavelet function (morlet,paul,dog)
    tbs : bool
      To Be Saved flag
    """
    # Define timestamps if not defined by user
    tmin = data.times[0].value  if tmin==None else tmin
    tmax = data.times[-1].value if tmax==None else tmax
    # Re-determine the mask in case timestamps are defined
    mask = (data.times.value>=tmin-1) & (data.times.value<=tmax+1)
    # Calculate wavelet parameters
    scales = mlpy.wavelet.autoscales(N=len(data[mask].value),dt=dt,dj=dj,wf=fct,p=omega0)
    spec = mlpy.wavelet.cwt(data[mask].value,dt=dt,scales=scales,wf=fct,p=omega0)
    freq = (omega0 + numpy.sqrt(2.0 + omega0 ** 2)) / (4 * numpy.pi * scales[1:])
    spec = numpy.abs(spec)**2
    spec = spec[::-1]
    # Convert time array into minute unit
    scale_factor = 3600. if tunit=='hrs' else 60. if tunit=='mins' else 1
    times = (data[mask].times.value-tmin)/scale_factor
    # Define minimum and maximum frequencies
    freq = freq * 1000. if funit=='mHz' else freq
    # Determine
    fmin_log,fmax_log = min(freq),max(freq)
    fmin_linear,fmax_linear = min(freq),max(freq)
    if fmin!=None:
        log_ratio = (numpy.log10(fmin)-numpy.log10(min(freq)))/(numpy.log10(max(freq))-numpy.log10(min(freq)))
        fmin_linear = min(freq)+log_ratio*(max(freq)-min(freq))
        fmin_log = fmin
    if fmax!=None:
        log_ratio = (numpy.log10(fmax)-numpy.log10(min(freq)))/(numpy.log10(max(freq))-numpy.log10(min(freq)))
        fmax_linear = min(freq)+log_ratio*(max(freq)-min(freq))
        fmax_log = fmax
    # Get minimum and maximum amplitude in selected frequency range
    idx = numpy.where(numpy.logical_and(fmin_log<freq[::-1],freq[::-1]<fmax_log))[0]
    vmin = spec[idx].min() if vmin==None else vmin
    vmax = spec[idx].max() if vmax==None else vmax
    # Initialise figure
    fig,ax = plt.subplots(figsize=(12,7))
    # Plot spectrogram
    img = ax.imshow(spec,extent=[times[0],times[-1],freq[-1],freq[0]],aspect='auto',
                    interpolation='nearest',cmap=cmap,norm=matplotlib.colors.LogNorm(vmin,vmax))
    ax.set_xlabel('Time [%s] from %s %s (%s)'%(tunit,datetime.utcfromtimestamp(tmin),zone,tmin),fontsize=15)
    ax.set_xlim(0,(tmax-tmin)/scale_factor)
    ax.set_yscale('linear')
    ax.set_ylim(fmin_linear,fmax_linear)
    ax.set_ylabel('Frequency [%s]'%funit,fontsize=15,labelpad=40)
    ax.grid(False)
    # Set up axis range for spectrogram
    twin_ax = ax.twinx()
    twin_ax.set_yscale('log')
    twin_ax.set_xlim(0,(tmax-tmin)/scale_factor)
    twin_ax.set_ylim(fmin_log,fmax_log)
    twin_ax.spines['top'].set_visible(False)
    twin_ax.spines['right'].set_visible(False)
    twin_ax.spines['bottom'].set_visible(False)
    ax.tick_params(which='both', labelleft=False, left=False)
    twin_ax.yaxis.tick_left()
    twin_ax.grid(False)
    # Add color bar and save figure
    cb = fig.colorbar(img,ax=ax,pad=0.01)
    cb.set_ticks(LogLocator())
    cb.set_label('Power $|\mathrm{W}_v|^2$ $[\mu T^2/\mathrm{Hz}]$',fontsize=15)
    cb.set_clim(vmin,vmax)
    plt.tight_layout()
    plt.savefig('%s.png'%fname) if tbs else plt.show()
    plt.close(fig)

def plot_psd(data,tmin=None,tmax=None,fname='psd',tbs=False):
    """
    Plot Power Spectral Density

    Parameters
    ----------
    data : gwpy.timeseries.TimeSeries
      Time series data
    tmin : datetime
      First timestamp
    tmax : datetime
      Last timestamp
    """
    # Define timestamps if not defined by user
    tmin = data.times[0].value  if tmin==None else tmin
    tmax = data.times[-1].value if tmax==None else tmax
    # Re-determine the mask in case timestamps are defined
    mask = (data.times.value>=tmin) & (data.times.value<=tmax)
    # Calculate PSD
    N = len(data[mask].value)
    dt = 1/data.sample_rate.value
    df = 1. / (N * dt)
    f = df * numpy.arange(N / 2)
    PSD = abs(dt * fftpack.fft(data[mask].value)[:N / 2]) ** 2
    psd = numpy.vstack((f,PSD)).T
    # Do the plotting
    fig = plt.figure(figsize=(12,7))
    plt.loglog(psd[:,0],psd[:,1])
    plt.xlabel('Frequency [Hz]')
    plt.ylabel('Power Spectral Density [Hz$^{-1}$]')
    plt.savefig('%s.png'%fname) if tbs else plt.show()
    plt.close(fig)

def plot_specgram(data,tmin=None,tmax=None,fmin=None,fmax=None,vmin=None,vmax=None,tbs=False,
                  funit='Hz',cmap='viridis',zone='Local',tunit='secs',fname='spectrogram',
                  stride=None,nfft=None,overlap=None):
    """
    Plot Fourier spectrogram.

    Parameters
    ----------
    data : gwpy.timeseries.TimeSeries
      Time series data
    tmin : datetime
      First timestamp
    tmax : datetime
      Last timestamp
    tbs : bool
      To Be Saved flag
    """
    # Define timestamps if not defined by user
    tmin = data.times[0].value  if tmin==None else tmin
    tmax = data.times[-1].value if tmax==None else tmax
    # Re-determine the mask in case timestamps are defined
    mask = (data.times.value>=tmin-1) & (data.times.value<=tmax+1)
    freq, times, spec = signal.spectrogram(data[mask],fs=data.sample_rate.value,
                                           nperseg=stride,noverlap=overlap,nfft=nfft)
    times = (data[mask].times.value[::int(stride-overlap)]-tmin)/scale_factor
    # Convert time array into minute unit
    scale_factor = 3600. if tunit=='hrs' else 60. if tunit=='mins' else 1
    #times = (data[mask].times.value-tmin)/scale_factor
    # Define minimum and maximum frequencies
    freq = freq * 1000. if funit=='mHz' else freq
    fmin = freq[1]      if fmin==None    else fmin
    fmax = max(freq)    if fmax==None    else fmax
    # Get minimum and maximum amplitude in selected frequency range
    idx = numpy.where(numpy.logical_and(fmin<=freq,freq<=fmax))[0]
    vmin = vmin if vmin!=None else numpy.sort(numpy.unique(spec[idx]))[1]
    vmax = spec[idx].max() if vmax==None else vmax
    # Initialise figure
    fig,ax = plt.subplots(figsize=(12,7))
    # Plot spectrogram
    img = ax.pcolormesh(times,freq,spec,cmap='inferno',norm=matplotlib.colors.LogNorm())
    ax.set_xlabel('Time [%s] from %s %s (%s)'%(tunit,datetime.utcfromtimestamp(tmin),zone,tmin),fontsize=15)
    ax.set_xlim(0,(tmax-tmin)/scale_factor)
    ax.set_ylim(fmin,fmax)
    ax.set_yscale('log')
    ax.set_ylabel('Frequency [%s]'%funit,fontsize=15,labelpad=40)
    ax.grid(False)
    # Add color bar and save figure
    cb = fig.colorbar(img,ax=ax,pad=0.01)
    cb.set_ticks(LogLocator())
    cb.set_label('Power $|\mathrm{W}_v|^2$ $[\mu T^2/\mathrm{Hz}]$',fontsize=15)
    cb.set_clim(vmin,vmax)
    plt.tight_layout()
    plt.savefig('%s.png'%fname) if tbs else plt.show()
    plt.close(fig)