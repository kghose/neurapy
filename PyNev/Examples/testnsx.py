import pylab
import nsx

lablib_fname = '/Users/kghose/RawData/dj-map-2009-03-23-03.dat'
nev_dir = '/Users/kghose/RawData/20090323/'
nev_fname = 'grfmap003.nev'
ns3_fname = 'grfmap003.ns3' 

f_nsx = open(nev_dir + ns3_fname, 'rb')
nsx_basic_header = nsx.read_basic_header(f_nsx)


lfp1 = nsx.read_electrode(f = f_nsx, 
               basic_header = nsx_basic_header,
               electrode = 14, 
               tstart_ms = 5000,
               t_dur_ms = 1000)
pylab.subplot(2,1,1)
t = pylab.arange(lfp1.size)/2000.
pylab.plot(t, lfp1)


offset = 0
t_window_ms = 500
N_wave = t_window_ms * 2
mean_ft = pylab.zeros(N_wave, dtype='float')
for n in range(100):
  lfp1 = nsx.read_electrode(f = f_nsx, 
               basic_header = nsx_basic_header,
               electrode = 14, 
               tstart_ms = t_window_ms*n + offset,
               t_dur_ms = t_window_ms)
  ft = pylab.absolute(pylab.fft(lfp1))
  mean_ft += ft

pylab.subplot(2,1,2)  
F = pylab.fftfreq(mean_ft.size, d = 1/2000.)
pylab.semilogy(F, mean_ft/50.)
pylab.xlabel('F (Hz)')
pylab.ylabel('Signal')