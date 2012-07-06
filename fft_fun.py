import scipy
from scipy.fftpack import fft, ifft, fftshift, ifftshift
from scipy.signal import fftconvolve
from matplotlib import pyplot

# Useful ploting functions
show = pyplot.show
close = pyplot.close
plot = pyplot.plot
figure = pyplot.figure
xlim = pyplot.xlim
ylim = pyplot.ylim
axis = pyplot.axis
title = pyplot.title
xlabel = pyplot.xlabel
ylabel = pyplot.ylabel
legend = pyplot.legend
grid = pyplot.grid

# Useful mathematical functions
pi = scipy.pi
e = scipy.e
array = scipy.array
arange = scipy.arange
exp = scipy.exp
sqrt = scipy.sqrt
angle = scipy.angle
real = scipy.real
imag = scipy.imag
sinc = scipy.sinc
sin = scipy.sin
cos = scipy.cos
sinh = scipy.sinh
cosh = scipy.cosh
arcsin = scipy.arcsin
arccos = scipy.arccos
arcsinh = scipy.arcsinh
arccosh = scipy.arccosh


def demodulate(N,N0,Y,Inv='False'):
    """
    Removes the frequency modulation caused by having a data set where the 't=0' 
    point isn't the first value in the set.
    
    All arrays should be scipy.array objects.
    
    Inputs: Y   the descrete fourier transform
            N   the number of points in the data set
            N0  the point that represents the zero-indexed location of 't=0'
            Inv if True, this function will remodulate (undo a demodulation)
    
    Returns:    the demodulated transform
    
    Note:   The descrete fourier transform is not the same thing as the fft. It
            can be obtained by the opperation "Y = fftshift(fft(y))" where 'y'
            would be the time domain sample.
    """
    n = scipy.arange(-N0,N-N0)
    
    assert(len(Y) == len(n))
    if Inv: return Y*exp((-2.0j*scipy.pi*n/N)*(N/2))
    else: return Y*exp((-2.0j*scipy.pi*n/N)*(N/2))


def step(t,h=1.):
    """
    The step fucntion
    
        ^
      h |---------
        |
    --------------> t
        |
        
    """ 
    if t == 0.: return h*.5
    elif t > 0.: return h
    else: return 0.

def box(t,w=1.,h=1.):
    """
    The Box Function
    
           ^
      h +--|--+
        |  |  |
    --------------> t
       -w  |  w
           |
          
    """
    if abs(t) < w: return h
    elif abs(t) == w: return h*.5
    else: return 0.

def stair(t,w=1.,h=1.):
    """
    Two box functions next to the orign.
    
           ^
      h +--|
        |  |  w
    --------------> t
       -w  |  | 
           |--+ -h
          
    """
    if abs(t) < w:
        if t < 0: return -h
        elif t ==0: return 0
        else: return h
    elif abs(t) == w: return h*.5
    else: return 0

    
def gaus(t,mean=0.,sigma=1.,norm=False):
    """
    Return gaussian distribution with mean = "mean" and
    standard deviation = "sigma"
    """
    if norm: A = 1./sqrt(2.*pi*sigma**2.)
    else: A = 1.
    return A*exp(-(t - mean)**2./(2.*sigma))

#_______________________________________________________________________________
if __name__ == "__main__":
    f_s = 2.**5    # sampling frequency
    fsig = 1.      # Max signal frequency
    N = 2**8       # Number of samples (should be a power of 2)

    fo = f_s/N
    wsig = 2*pi*fsig
    n = arange(-N/2,N/2)
    t = n/f_s
    f = n*fo
    f -= scipy.median(f)
    assert(N == len(n)) # Should be true for this script

    # Calculate the time function
#    y = sin(wsig*t);
#    y = cos(wsig*t);
    y = 2*fsig*sinc(wsig*t)
#    y = array(map(lambda t:step(t,1.),t))
#    y = array(map(lambda t:gaus(t,0.,.01),t))
#    y = array(map(lambda t:box(t,1.,1.),t))
#    y = array(map(lambda t:stair(t,1.,1.),t))

    # Calculate the transfer function
#    h = sin(wsig*t);
#    h = cos(wsig*t);
#    h = 2*fsig*sinc(wsig*t)
#    h = array(map(lambda t:step(t,1.),t))
#    h = array(map(lambda t:gaus(t,0.,1,True),t))
#    h = array(map(lambda t:box(t,1.,1.),t))
#    h = array(map(lambda t:stair(t,1.,1.),t))
    
    # Calculate the transfer function starting with a filter in frequency space
#    h = array(map(lambda f:step(f,2.),f))
    h = array(map(lambda f:gaus(f,0.,1.0),f))
#    h = array(map(lambda f:box(f,1.,1.),f))
    h = ifft(ifftshift(demodulate(N,N/2,h,True))) # Convert to time space
    h = real(h)
    
    # Fourier Transform
    Y = demodulate(N,N/2,fftshift(fft(y)))
    H = demodulate(N,N/2,fftshift(fft(h)))
    
    # Convolution of y and h
    G = Y*H
    g = fftconvolve(y,h,mode='same')
#    g = ifft(ifftshift(demodulate(N,N/2,G,True)))  # Good proof of concept
    
    Yabs = abs(Y)
    Yang = angle(Y)
    Yre = real(Y)
    Yim = imag(Y)
    Ymax = max(Yabs)
    
    Habs = abs(H)
    Hang = angle(H)
    Hre = real(H)
    Him = imag(H)
    Hmax = max(Habs)
    
    Gabs = abs(G)
    Gang = angle(G)
    Gre = real(G)
    Gim = imag(G)
    Gmax = max(Gabs)
    
    
    # Views of the signals
    figure(1)
    taxis = [0. for x in t]
    a = max(abs(y))
    plot(t,y,'.-',label='y(t)',color='#3465a4')
    plot(t,real(h),'.-',label='h(t)',color='#f57900')
    
    plot(t,g,'.-',label='g',color='#73d216')
    
    plot(t,taxis,"k")
    ylim([-1.1*a,1.1*a])
    title("Time Domain")
    xlabel("t (sec)")
    legend()
    grid(True)


    # Views of the fourier transform
    figure(2)
    faxis = [0. for x in f]
    
    plot(f,Yre,'.-',label='Re(Y)',color='#3465a4')
    plot(f,Yim,'.-',label='Im(Y)',color='#cc0000')
#    plot(f,Yabs,'.-',label='Abs(Y)',color='#f57900')
#    plot(f,Yang,'.-',label='Ang(Y)',color='#75507b')

#    plot(f,Hre,'.-',label='Re(H)',color='#204a87')
#    plot(f,Him,'.-',label='Im(H)',color='#a40000')
    plot(f,Habs,'.-',label='Abs(H)',color='#ce5c00')
#    plot(f,Hang,'.-',label='Ang(H)',color='#5c3566')

#    plot(f,Gre,'.-',label='Re(G)',color='#729fcf')
#    plot(f,Gim,'.-',label='Im(G)',color='#ef2929')
#    plot(f,Gabs,'.-',label='Abs(G)',color='#fcaf3e')
#    plot(f,Gang,'.-',label='Ang(G)',color='#ad7fa8')

    plot(f,taxis,'k')
    ylim([-1.1*Ymax,1.1*Ymax])
    legend()
    title("Frequency Domain")
    xlabel("f (Hz)")
    grid(True)

    # Display plots
    show()


