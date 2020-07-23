from numpy import *
from numpy.linalg import *
from scipy import signal
from scipy.signal import butter,lfilter,hilbert
from scipy.stats import ranksums
from scipy.io import savemat, loadmat, whosmat
from pylab import *
import os as os
from os import listdir
from os.path import isfile, join
from random import sample


def plot_spec2(x,fs):
 x=(x-mean(x))/std(x)
 f, Pxx_den = signal.welch(x, fs,scaling='density')
 plt.semilogy(f, Pxx_den)
 plt.xlabel('frequency [Hz]')
 plt.show()

def find_closest(A, target):
   '''
   helper function for power spectrum functions (plot and PSpec)
   '''
   #A must be sorted
   idx = A.searchsorted(target)
   idx = np.clip(idx, 1, len(A)-1)
   left = A[idx-1]
   right = A[idx]
   idx -= target - left < right - target
   return idx

def plot_spec(data,fs):

 #data=one dimensional time series

 opt=1
# #n=14 #number of bands size 5 Hz
# F=zeros((n,2))
# for i in range(n):
#  F[i]=[5*i,5*(i+1)]
 de=[1,4]# in Hz
 th=[4,8]
 al=[8,13]
 be=[13,30]
 ga=[30,60]
 hga=[60,120]

 F=[de,th,al,be,ga,hga]

 v=data

 #ro,co=shape(w18[0][0:10000])
 co=len(v)
 # Number of samplepoints
 N = co
 print co ,'N'
 # sample spacing (denominator in Hz)
 T = 1.0 / fs
 x = np.linspace(0.0, N*T, N)
 y = v#[0]
 yf = fft(y)
 xf = np.linspace(0.0, 1.0/(2.0*T), N/2)

 yff=2.0/N * np.abs(yf[0:N/2])
 print shape(yff)
 al=sum(2.0/N * np.abs(yf[0:N/2]))

# figure(1)
# plt.plot(y)

 figure(2)
 plt.semilogy(xf[0:],2.0/N * np.abs(yf[0:N/2]))
 plt.grid()
 plt.show()

 if opt==1:

  bands=zeros(len(F))
  for i in range(len(F)):
   bands[i]=sum(yff[find_closest(xf, F[i][0]):find_closest(xf, F[i][1])])
  bands=bands/sum(bands)

  figure(3)
  bar(arange(len(F)),bands,align='center',fc='none', edgecolor='k')
  print bands

def Pre2(X):
 '''
 Linear-detrend and subtract mean
 '''
 ro,co=shape(X)
 Z=zeros((ro,co))
 for i in range(ro): #maybe divide by std?
  try:
   Z[i,:]=signal.detrend((X[i,:]-mean(X[i,:]))/std(X[i,:]), axis=0)
  except ValueError:
   pass
 return Z

def Ph_rand(original_data):
 '''
 phase randomisation: multi time series in, shuffled time series out - but it has same power spectrum
 '''
 #NOTE WHEN DOING SCHREIBER_SURRO< MAYBE TRY EFFECT OF DIFFERENT PHASE DISTRIBUTIONS ON MEASURES

 surrogates = np.fft.rfft(original_data, axis=1)
 #  Get shapes
 (N, n_time) = original_data.shape
 len_phase = surrogates.shape[1]

 #  Generate random phases uniformly distributed in the
 #  interval [0, 2*Pi]
 phases = np.random.uniform(low=0, high=2 * np.pi, size=(N, len_phase))

 #  Add random phases uniformly distributed in the interval [0, 2*Pi]
 surrogates *= np.exp(1j * phases)

 #  Calculate IFFT and take the real part, the remaining imaginary part
 #  is due to numerical errors.
 return np.ascontiguousarray(np.real(np.fft.irfft(surrogates, n=n_time,axis=1)))

##############
'''
PSpec
'''
##############

#def PSpec(X):

# ro,co=shape(X)
# L=zeros((ro,5))

# fs=250.#Hz psy 250Hz, ic 250Hz
# de=[1,4]# Bands in Hz
# th=[4,8]#obs
# al=[8,15]#[8,13]
# be=[15,30]#[13,30]
# ga=[30,70]
# F=[de,th,al,be]#,ga]
# Fs=['de','th','al','be']#,'ga']
# #Sureh'[1 4; 4 8; 8 15; 15 30; 30 49;51 99]

# for i in range(ro):
#  x=X[i]
#  x=signal.detrend((x-mean(x))/std(x), axis=0)
#  Q=zeros(5) #summed up frequencies for all channels
#  f, Pxx_den = signal.welch(x, fs, nfft=2400)
#  bands=zeros(5)



#  for ii in range(5):
#   bands[ii]=sum(Pxx_den[find_closest(f, F[ii][0]):find_closest(f, F[ii][1])])

#  #print bands[i]
#  Q=bands/sum(bands) #NormC:\Temp\MATLAB\_______LUCIA_BEN\_________FINAL_Michael\Originalalization
#  L[i]=Q


# return L #mean(L,0)


def PSpec(X):
 X=Pre2(X)
 fs=250. #sampling rate in Hz

 #data=one dimensional time series

 opt=1
# #n=14 #number of bands size 5 Hz
# F=zeros((n,2))
# for i in range(n):
#  F[i]=[5*i,5*(i+1)]
 de=[1,4]# in Hz
 th=[4,8]
 al=[8,13]
 be=[13,30]
 ga=[30,60]
 hga=[60,120]

 F=[de,th,al,be]#,ga,hga]

 ro,co=shape(X)
 Q=[]

 for i in range(ro):

  v=X[i]
  co=len(v)
  # Number of samplepoints
  N = co
  # sample spacing (denominator in Hz)
  T = 1.0 / fs
  y = v
  yf = fft(y)
  xf = np.linspace(0.0, 1.0/(2.0*T), N/2)
  yff=2.0/N * np.abs(yf[0:N/2])

  bands=zeros(len(F))
  for i in range(len(F)):
   bands[i]=sum(yff[find_closest(xf, F[i][0]):find_closest(xf, F[i][1])])
  bands=bands/sum(bands)
  Q.append(bands)
 return Q



#############
'''
frequency filter
'''
#############

def butter_lowpass(cutoff, fs, order=5):
    nyq = 0.5 * fs
    normal_cutoff = cutoff / nyq
    b, a = butter(order, normal_cutoff, btype='low', analog=False)
    return b, a

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def butter_highpass(lowcut, fs, order):
    nyq = 0.5 * fs
    low = lowcut / nyq
    b, a = butter(order, low, btype='highpass')
    return b, a

def butter_highpass_filter(data, lowcut, fs, order):
    b, a = butter_highpass(lowcut, fs, order=order)
    y = lfilter(b, a, data)
    return y

def notch_iir(fs,f0,data):
    '''
    fs: Sample frequency (Hz)
    f0: Frequency to be removed from signal (Hz)
    '''

    Q = 10.# 30.0  # Quality factor
    w0 = float(f0)/(fs/2)  # Normalized Frequency
    b, a = signal.iirnotch(w0, Q)
    return lfilter(b, a, data)


##########
'''
LZc - Lempel-Ziv Complexity, column-by-column concatenation

X is continuous multidimensional time series, channels x observations
'''
##########

def cpr(string):
 '''
 Lempel-Ziv-Welch compression of binary input string, e.g. string='0010101'. It outputs the size of the dictionary of binary words.
 '''
 d={}
 w = ''
 i=1
 for c in string:
  wc = w + c
  if wc in d:
   w = wc
  else:
   d[wc]=wc
   w = c
  i+=1
 return len(d)

def str_col(X):
 '''
 Input: Continuous multidimensional time series
 Output: One string being the binarized input matrix concatenated comlumn-by-column
 '''
 ro,co=shape(X)
 TH=zeros(ro)
 M=zeros((ro,co))
 for i in range(ro):
  M[i,:]=abs(hilbert(X[i,:]))
  TH[i]=mean(M[i,:])

 s=''
 for j in xrange(co):
  for i in xrange(ro):
   if M[i,j]>TH[i]:
    s+='1'
   else:
    s+='0'

 return s


def LZc(X):
 '''
 Compute LZc and use shuffled result as normalization
 '''
 X=Pre2(X)
 SC=str_col(X)
 M=list(SC)
 shuffle(M)
 w=''
 for i in range(len(M)):
  w+=M[i]
 return cpr(SC)/float(cpr(w))


def LZs(x):

 '''
 Lempel ziv complexity of single timeseries
 '''

 #this differes from Sitt et al as they use 32 bins, not just 2.
 co=len(x)
 x=signal.detrend((x-mean(x))/std(x), axis=0)
 s=''
 r=(abs(hilbert(x)))
 th=mean(r)

 for j in xrange(co):
  if r[j]>th:
   s+='1'
  else:
   s+='0'

 M=list(s)
 shuffle(M)
 w=''
 for i in range(len(M)):
  w+=M[i]

 return cpr(s)/float(cpr(w))

############
'''
Functions for automatic renaming/copying files
'''
############


def list_files(startpath):
    '''
    get file tree printed for given path
    '''

    for root, dirs, files in os.walk(startpath):
        level = root.replace(startpath, '').count(os.sep)
        indent = ' ' * 4 * (level)
        print('{}{}/'.format(indent, os.path.basename(root)))
        subindent = ' ' * 4 * (level + 1)
        for f in files:
            print('{}{}'.format(subindent, f))


'''
useful subroutines for listing all files in a folder and renaming/moving them
'''


#def renameS(expName):
# path='/home/mic/data/kinect/'


# #expName='MH004_post_cocaine_40_60'
# oldD=path+expName
# for filename in os.listdir(oldD):
#  if filename.startswith('RID'):
#   if 'RID__' in filename:
#    os.remove(oldD+'/RID__.int16binary')

#   newD=path+'/all/%s_t%s' %(expName,filename.split('_')[1])
#   if not os.path.exists(newD):   
#    os.makedirs(newD)
#   newFN='RID.int16binary'
#   os.rename('%s/%s' %(oldD,filename),'%s/%s' %(newD,newFN))  
#   shutil.copy('%s/polyroi.csv' %oldD,'%s/polyroi.csv' %newD)
#   shutil.copy('%s/BG.mat' %oldD,'%s/BG.mat' %newD)



############
'''
Data Version saving
'''
############

def SaveDataVersions():
 '''
 Take data from Version 0 and save different versions of it
 '''
 #Version 0: minimally processed (mini)
 #Version 1: notched at multiples of 3Hz (n3Hz)
 #Version 2: notched at multiples of 10Hz (n10Hz)
 #Version 3: phase shuffled (phaseShuff)
 #Version 4: average spectrum of dark imposed (darkSpec) 

 versions=['phaseShuff']
 #versions=['phaseShuff','darkSpec']

 #path to mini
 path_dat='C:\Temp\MATLAB\_______LUCIA_BEN\_________FINAL_Michael'
 
 subjects=[2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23]
# subjects=[3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21,22,23]
 conditions=['dark','flick3','flick10']

 for subject in subjects:
  for version in versions:
   if version=='n3Hz': 
    for cond in conditions:
     d=loadmat(path_dat+'\mini\%s\%s_%s.mat' %(subject,subject,cond))['data']
     for f0 in arange(1,14)*3: #take out multiples of 3Hz up to 39Hz
      d=notch_iir(250,f0,d)
     if not os.path.exists(path_dat+version):
      os.makedirs(path_dat+version)
     savemat(path_dat+version+'\%s_%s.mat' %(subject,cond),mdict={'data':d})
   if version=='n10Hz':
    for cond in conditions:
     d=loadmat(path_dat+'mini\%s_%s.mat' %(subject,cond))['data']
     for f0 in arange(1,5)*10: #take out multiples of 3Hz up to 39Hz
      d=notch_iir(250,f0,d)
     if not os.path.exists(path_dat+version):
      os.makedirs(path_dat+version)
     savemat(path_dat+version+'\%s_%s.mat' %(subject,cond),mdict={'data':d})
   if version=='phaseShuff':
    for cond in conditions:
     d=loadmat(path_dat+'\mini\%s\%s.mat' %(subject,cond))['data']
     chans,obs=shape(d)
     dn=zeros(shape(d))
     for s in range(obs/2500):
      dn[:,s*2500:(s+1)*2500]=Ph_rand(d[:,s*2500:(s+1)*2500])
     if not os.path.exists(path_dat+version):
      os.makedirs(path_dat+version)
     savemat(path_dat+version+'\%s_%s.mat' %(subject,cond),mdict={'data':dn})
   if version=='darkSpec':
    d=loadmat(path_dat+'\mini\%s\%s.mat' %(subject,'dark'))['data']
    for cond in conditions:
     dpre=loadmat(path_dat+'\mini\%s\%s.mat' %(subject,cond))['data']
     chans,obs=shape(dpre)
     dn=zeros(shape(dpre))
     for s in range(obs/2500): #copy spectrum per segment
      R=abs(rfft(d[:,s*2500:(s+1)*2500]))
      ft=rfft(dpre[:,s*2500:(s+1)*2500])
      ft=ft/abs(ft)
      ft=nan_to_num(ft)
      dn[:,s*2500:(s+1)*2500]=irfft(R*ft)    
     if not os.path.exists(path_dat+version):
      os.makedirs(path_dat+version)
     savemat(path_dat+version+'\%s_%s.mat' %(subject,cond),mdict={'data':dn})



############
'''
Lucia EEG, sampling rate 250Hz, 64 channels times 15sec segments (40)
'''
############io.loadmat('faces_e_22.set', uint16_codec='latin1')

def LuciaV():

 func_dict ={'LZs':LZs,'LZc':LZc,'PSpec':PSpec}
 l =  2500

 subjects=[2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,21,22,23]
 #[3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,21,22,23]
 conditions=['dark','flick3','flick10']
 #['Dark','Light','Flick3','Flick10']

 #path_dat='C:\Temp\Dropbox\____My_Papers\LUCIA\EEG Data\For Michael'
 path_dat='C:\Temp\MATLAB\_______LUCIA_BEN\_________FINAL_Michael'
 path_res='C:\Temp\MATLAB\_______LUCIA_BEN\_________FINAL_Michael\Results'

 versions=['mini','n3Hz','n10Hz','phaseShuff','darkSpec']
 #versions=['mini','n3Hz','n10Hz','phaseShuff','darkSpec']
 for version in versions:
  for cond in conditions:
   for subject in subjects:

    d=loadmat(path_dat+'\%s\%s\%s.mat' %(version,subject,cond))['data']
    chs,obs=shape(d)

    print subject, shape(d)
    for measure in func_dict:
     if measure=='LZs':
      S=[]
      for t in range(obs/l):
       s0=[]
       ro,co=shape(d[:,t*l:(t+1)*l])
       for jj in range(ro):
        s0.append(func_dict[measure](d[jj,t*l:(t+1)*l]))
       print version, cond, subject, measure, t,'of',obs/l
       S.append(s0)

     if measure=='LZc':
      S=[]
      for t in range(obs/l):
       S.append(LZc(d[:,t*l:(t+1)*l]))
    
     elif measure=='PSpec':
      S=[]
      for t in range(obs/l):
       print version, cond, subject, measure, t,'of',obs/l
       S.append(PSpec(d[:,t*l:(t+1)*l]))
      S=array(S)

    
     if not os.path.exists(path_res+'%s\Sub%s' %(version,subject)):
      os.makedirs(path_res+'%s\Sub%s' %(version,subject))
     if measure=='PSpec':
      for jj in range(4):
       savetxt(path_res+'%s\Sub%s\PSpec%s_%s'%(version,subject,jj,cond),S[:,:,jj])
     else:
      print shape(S)
      savetxt(path_res+'%s\Sub%s\%s_%s' %(version,subject,measure,cond),S)



