import scipy
from scipy.fftpack import fft, fftshift, ifft
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


def demodulate(N,N0,Y):
    """
    Removes the frequency modulation caused by having a data set where the 't=0' 
    point isn't the first value in the set.
    
    All arrays should be scipy.array objects.
    
    Inputs: Y   the descrete fourier transform
            N   the number of points in the data set
            N0  the point that represents the zero-indexed location of 't=0'
    
    Returns:    the demodulated transform
    
    Note:   The descrete fourier transform is not the same thing as the fft. It
            can be obtained by the opperation "Y = fftshift(fft(y))" where 'y'
            would be the time domain sample.
    """
    n = scipy.arange(-N0,N-N0)
    return Y*exp((2.0j*scipy.pi*n/N)*(N/2))
    

# The box wave function
def box(t,w=1.,h=1.):
    """
    Given a value, t, will return the value of the square pulse at that time.
    """
    if abs(t) < w: return h
    else: return 0.

# Stair wave function
def stair(t,w=1.,h=1):
    if abs(t) < w:
        if t < 0: return h
        elif t ==0: return 0
        else: return -h
    else: return 0

# Trifilter function
def trifilter(t,w=1.,h=1):
    if t == w: return h
    elif t == 0: return -2.*h
    elif t == -2: return h
    else: return 0

    
# Return a Gaussian distribution
def gaus(t,mean=0.,sigma=1.):
    """
    Return gaussian distribution with mean = "mean" and
    standard deviation = "sigma"
    """
    return exp(-(t - mean)**2./(2.*sigma))/sqrt(2.*pi*sigma**2.)

#_______________________________________________________________________________
if __name__ == "__main__":
    f_s = 2.**5    # sampling frequency
    fsig = 1.      # Max signal frequency
    N = 2**10      # Number of samples (should be a power of 2)

    fo = f_s/N
    wsig = 2*pi*fsig
    n = arange(-N/2,N/2)
    t = n/f_s
    assert(N == len(n))

    # Calculate the time function
    #y = sin(wsig*t);
    #y = cos(wsig*t);
    #y = 2*fsig*sinc(wsig*t)
    #y = array(map(gaus,t))
    y = array(map(box,t))
    #y = array(map(stair,t))
    #y = array(map(trifilter,t))
    

    figure(1)
    taxis = [0. for x in t]
    a = max(abs(y))
    plot(t,y,'.-')
    plot(t,taxis,"k")
    ylim([-1.1*a,1.1*a])
    title("Time Domain")
    xlabel("t (sec)")
    ylabel("y(t)")

    # Fourier Transform
    f = n*fo
    f -= scipy.median(f)
    Y = fftshift(fft(y))
    Y = demodulate(N,N/2,Y)     # Needed to fix the "time shift"
                                # createdfrom not starting at t=0

    # Views of the transform
    Yabs = abs(Y)
    Yang = angle(Y)
    Yre = real(Y)
    Yim = imag(Y)
    A = max(Yabs)

    figure(2)
    faxis = [0. for x in f]
    plot(f,Yre,'.-',label='Re(Y)')
    plot(f,Yim,'.-',label='Im(Y)')
#    plot(f,Yabs,'.-',label='Abs(Y)')
#    plot(f,Yang,'.-',label='Ang(Y)')
    plot(f,taxis,'k')
    ylim([-1.1*A,1.1*A])
    legend()
    title("Frequency Domain")
    xlabel("f (Hz)")
    ylabel("Y(f)")

    # Display plots
    show()


