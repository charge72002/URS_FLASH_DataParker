import numpy as np
import matplotlib.pyplot as plt

def f(t):
    return np.exp(-0.5*t**2)

def g(t):
    return (t-1)*np.sin(np.pi/(t-1))

def dfdt(t,f):
    return -t * f

def dgdt(t,g):
    return g/(t-1) - np.pi/(t-1)*np.cos(np.pi/(t-1))

def EulerStep(t,y,dt,derivative):
    return y + (derivative(t, y) * dt)

def RungeKutta2(t,y,dt,derivative):
    k1 = derivative(t, y)
    k2 = derivative(t+(0.5*dt), y+ (k1*0.5*dt))
    return y + (k2 * dt)


def RungeKutta4(t,y,dt,derivative):
    k1 = derivative(t, y)
    k2 = derivative(t+(0.5*dt), y+ (k1*0.5*dt) )
    k3 = derivative(t+(0.5*dt), y+ (k2*0.5*dt) )
    k4 = derivative(t+dt, y+ (k3*dt))
    return y + ((1/6) * (k1 + (2*k2) + (2*k3) + k4) * dt)

def Integrate(step,derivative,t0,y0,dt,tf):
    tLst = np.arange(t0,tf,dt)
    n = len(tLst)
    yLst = np.zeros(n)
    yLst[0] = y0
    i = 1
    oldy = yLst[i-1]
    oldt = tLst[i-1]
    while i < n:
        yLst[i] = step(tLst[i-1],yLst[i-1],dt,derivative)
        i +=1
    return tLst,yLst

def main():
    plt.clf()
    Start = 0.0
    End = 0.99
    gInitial = 0
    fInitial = 1

    tActual = np.linspace(Start,End,1000)
    gActual = g(tActual)
    fActual = f(tActual)

    dt = 0.05
    tLstEulF, yLstEulF = Integrate(EulerStep,dfdt,Start,fInitial,dt,End)
    tLstRK2F, yLstRK2F = Integrate(RungeKutta2,dfdt,Start,fInitial,dt,End)
    tLstRK4F, yLstRK4F = Integrate(RungeKutta4,dfdt,Start,fInitial,dt,End)
    dt = 0.01
    tLstEulG, yLstEulG = Integrate(EulerStep,dgdt,Start,gInitial,dt,End)
    tLstRK2G, yLstRK2G = Integrate(RungeKutta2,dgdt,Start,gInitial,dt,End)
    tLstRK4G, yLstRK4G = Integrate(RungeKutta4,dgdt,Start,gInitial,dt,End)


    fig1 = plt.figure(1,figsize=(11,4))
    axLst1 = fig1.subplots(ncols=2,nrows=1,sharex='col')
    axLst1[0].plot(tActual,fActual,'k-',label='Analytic')
    axLst1[0].plot(tLstEulF,yLstEulF,'r.',label='Euler')
    axLst1[0].plot(tLstRK2F,yLstRK2F,'g.',label='RK2')
    axLst1[0].plot(tLstRK4F,yLstRK4F,'b.',label='RK4')

    axLst1[0].set_ylabel(r'$f(t)$')
    axLst1[0].grid()
    axLst1[0].legend()

    axLst1[1].semilogy(tLstEulF,np.abs(yLstEulF - f(tLstEulF)),'r.')
    axLst1[1].semilogy(tLstRK2F,np.abs(yLstRK2F - f(tLstRK2F)),'g.')
    axLst1[1].semilogy(tLstRK4F,np.abs(yLstRK4F - f(tLstRK4F)),'b.')

    axLst1[1].set_xlabel(r'$t$')
    axLst1[1].grid()
    axLst1[1].set_ylabel(r'$|f_\mathrm{numeric} - f_\mathrm{analytic}|$')

    fig2 = plt.figure(2,figsize=(11,4))
    axLst2 = fig2.subplots(ncols=2,nrows=1,sharex='col')
    axLst2[0].plot(tActual,gActual,'k-',label='Analytic')
    axLst2[0].plot(tLstEulG,yLstEulG,'r.',label='Euler')
    axLst2[0].plot(tLstRK2G,yLstRK2G,'g.',label='RK2')
    axLst2[0].plot(tLstRK4G,yLstRK4G,'b.',label='RK4')

    axLst2[0].set_ylabel(r'$g(t)$')
    axLst2[0].grid()
    axLst2[0].legend()

    axLst2[1].semilogy(tLstEulG,np.abs(yLstEulG - g(tLstEulG)),'r.')
    axLst2[1].semilogy(tLstRK2G,np.abs(yLstRK2G - g(tLstRK2G)),'g.')
    axLst2[1].semilogy(tLstRK4G,np.abs(yLstRK4G - g(tLstRK4G)),'b.')

    axLst2[1].set_xlabel(r'$t$')
    axLst2[1].set_ylabel(r'$|g_\mathrm{numeric} - g_\mathrm{analytic}|$')
    axLst2[1].grid()


    plt.show()



    return
main()
