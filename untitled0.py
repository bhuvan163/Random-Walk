# Calculating the MSD

import pandas
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as so
import math

# Exponential fit
def expfit(x, a, b):
    return a * np.exp( b * x )
# END def

# Fit for correlated random walk
def CorrFit(t, tau, vel):
    return ( 2 * tau * vel * vel ) * ( t - ( tau * ( 1 - np.exp( -t / tau ) ) ) )
# END of def

# Linear Fit
def LF(x, D):
    return 4*D*x
# END of def

data = pandas.read_excel('WT.xlsx', sheet_name='Tabelle1')

xcoord = np.array(data["x - coord"][1:])
ycoord = np.array(data["y - coord"][1:])
time = np.array(data["time frame"][1:])
ID = np.array(data["ID"][1:])

N_bac = int(ID[-1])  # Total number of bacteria
N_frame = time[-1]  # Total time frames

i = 0   # counter for data points
#MSD = 0 # Mean Square Displacement

# Making a different DataFrame
x_zeroarray = np.zeros((N_bac+1, N_frame+1))
y_zeroarray = np.zeros((N_bac+1, N_frame+1))

i = 0
for cellID in range(1, int(N_bac)+1):
    start = i
    while cellID == ID[i]:
        i += 1
        if i == len(xcoord):
            break
    # cell ID vs frame number
    x_zeroarray[cellID, time[start]:(time[i-1]+1) ] = xcoord[start:i]
    y_zeroarray[cellID, time[start]:(time[i-1]+1) ] = ycoord[start:i]
# END of for

# 1 px = 0.8 mu metre
#x_zeroarray = 0.8 * x_zeroarray # [mu metre]^2
#y_zeroarray = 0.8 * y_zeroarray # [mu metre]^2

# Caluculating the r(t)
r_i_t = np.zeros((N_bac+1, N_frame+1))
#r_i_t[:] = np.nan
for i in range(1,int(N_bac)+1):
    for j in range(1, int(N_frame)+1):
        tempx = np.copy(x_zeroarray[i,1])
        tempy = np.copy(y_zeroarray[i,1])
        if x_zeroarray[i,j] !=0 and y_zeroarray[i,j] !=0 :
            del_x = ( ( x_zeroarray[i,j] - tempx ) * 0.8 ) ** 2
            del_y = ( ( y_zeroarray[i,j] - tempy ) * 0.8 ) ** 2
            r_i_t[i,j] = math.sqrt( del_x + del_y )
        else:
            continue
# END of for

# Calculating the Mean-Squared displacement
rit_old = np.copy(r_i_t)
for i in range(1,int(N_bac)+1):
    for j in range(1, int(N_frame)+1):
        temp = np.copy(r_i_t[i,1])
        if r_i_t[i,j] !=0 :
            r_i_t[i,j] = ( r_i_t[i,j] - temp ) ** 2
        else:
            continue
# END of for

MSD = np.sum(r_i_t, axis = 0) / N_bac # [mu metre]^2
#MSD = 0.001 * MSD
#print(len(MSD))
"""

"""
# Plotting MSD as a function of time frame -----------------------
# dt = 100 ms = 0.1 s
plt.scatter(x = 0.1*np.array(range(len(MSD[:100]))), y = MSD[:100], s=4)
plt.xlabel('Time Frame [s]')
plt.ylabel('MSD [$\mu m^2$]')
#plt.show()
#plt.ticklabel_format(axis="x", style="sci", scilimits=(0,0))
#plt.title()

# Linear fit of the plot ------------------------------------------
popt, pcov = so.curve_fit(LF, xdata = 0.1*np.array(range(len(MSD[:100]))),
                          ydata = MSD[:100])
fitx = 0.1*np.array(range(len(MSD[:100])))
fity = LF(fitx, *popt)
plt.plot(fitx,fity, c = 'red', label = 'linear fit')
#plt.ylim(5,16)
#plt.ylim(30,40)
plt.show()
print("Linear fit parameters: ",popt)
print('covariance: ', pcov)

# Velocity correlation fit ----------------------------------------
plt.scatter(x = 0.1*np.array(range(len(MSD[:100]))), y = np.log(MSD[:100]), s=4)
plt.xlabel('Time Frame [s]')
plt.ylabel('log(MSD) [$\mu m^2$]')

popt, pcov = so.curve_fit(CorrFit, xdata = 0.1*np.array(range(len(MSD[:100]))),
                          ydata = (MSD[:100]), maxfev=18000)
fitx = 0.1*np.array(range(len(MSD[:100])))
fity = CorrFit(fitx, *popt)
plt.plot(fitx, np.log(fity), c = 'red', label = 'velocity and tau')
#plt.xlabel("Time Frame")
#plt.ylabel("$< v(t) \cdot v(t + \delta t) >$")
#plt.xlim(0, 200)
#plt.ylim(0,200)
#plt.legend()
plt.show()
print("Correlation parameters: ",popt)

# Calculating Velocity auto-correlation function -------------------
vel_t = np.zeros((N_bac+1, N_frame+1))
#vel_corr = np.zeros(N_bac+1, N_frame+1)
for i in range(1,int(N_bac)+1):
    for j in range(1, int(N_frame)+1):
        if r_i_t[i,j] !=0:
            vel_t[i,j] = ( r_i_t[i,j] - r_i_t[i,1] ) / ( (j-1) * 0.1 )
        else:
            continue
# END of for

vel_t_old = np.copy(vel_t)
for i in range(1,int(N_frame)):
    vel_t[:,i] = ( vel_t[:,i] * vel_t[:,i+1] ) # mu meter square
    #vel_t[:,i] = ( vel_t[:,i] - vel_t[:,i+1] ) ** 2 # mu meter square
# END for

vel_corr = np.sum(vel_t, axis = 0) / N_bac # [mu metre]^2

plt.scatter(0.1*np.array(range(2,N_frame+1)), vel_corr[2:], s=4)
plt.xlabel('Time [s]')
plt.ylabel("$< v(t) \cdot v(t + \delta t) >$ [$\mu m^2/s$]")
plt.show()

plt.scatter(0.1*np.array(range(2,100)), vel_corr[2:100], s=4)
plt.xlabel('Time [s]')
plt.ylabel("$< v(t) \cdot v(t + \delta t) >$ [$\mu m^2/s$]")
plt.xlim(0,10)
#plt.show()

popt, pcov = so.curve_fit(expfit, xdata = 0.1*np.array(range(2,100)),
                          ydata = vel_corr[2:100])
fitx = 0.1*np.array(range(2, len(vel_corr[2:100])))
fity = expfit(fitx, *popt)
plt.plot(fitx,(fity), c = 'red')
plt.show()
print("VAC Exponential fit paramters: ", popt)