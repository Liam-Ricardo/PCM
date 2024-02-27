
"""
Created on Thu Jan  7 17:19:49 2021


"""


import numpy as np
import matplotlib.pyplot as plt

'All constants are defined here:'
y=0.683 # $/g (Cost in Dollar/Gram)
T_melting=328  #The melting temperature of the PCM
h_sl= 1.63*10**5 #   h_l-h_s #enthalpy difference 

A= 1 #m^2
T_lucht= 303 #Air temperature at 30 degrees Celsius (India)
L= 0.01 #Thickness pcm
t_stop= 28800 # 8 hours = 28800 seconds
n= 1 
z=0 
dt= 1
k_l= 0.15 # W/(mK)
k_s= 0.24 # W/(mK)
x_melting=0 
h_lucht=8.9
e_lucht=0.95
e_pcm = 0.90
sigma = 5.7*(10**-8)
mu=-0.2 
rho= 800 # kg/m^3
I = 700 # W/M^2 # intensity irradiance
B=0.0045 # K^-1
eta_ref = 0.156
T_ref= 298 #in Kelvin
alfa = 0.9 # absorptance of PV cell, typically 0.9
T_sky= 0.0552*T_lucht**1.5
T0_eind= 0
cp = 2900
m=rho*A*L
T0 = np.zeros(t_stop+1)
T_l = np.zeros(t_stop+1)
T0[0]=T_lucht
T_l[0]=T_lucht 
T0[1] = T0[0]
T_l[1] = T_l[0]
v_melting = np.zeros(t_stop+1)
x = np.zeros(t_stop+1)
t = np.linspace(0, t_stop, t_stop+1)
t[0]=0
Tx_eind =np.zeros(t_stop+1)
eta= np.zeros(t_stop+1)
eta[0]=0
q0= np.zeros(t_stop+1)
q0[0]=0
qL= np.zeros(t_stop+1)
qL[0]=0
q=np.zeros(t_stop+1)

'All functions are defined here:'
def f_eta(T0):
    eta = eta_ref*(1-(B*(T0-T_ref)))
    if eta<0:
        return 0
    else:
        return eta
    
def flux_in(T0):
    q_conv= h_lucht*(T_lucht-T0)
    q_lw = sigma*(e_lucht*(T_sky**4)-(e_pcm*(T0**4)))
    q_sw=(1-f_eta(T0))*alfa*I
    q0 = q_conv + q_lw + q_sw 
    return q0 , q_conv , q_sw , q_lw

def flux_uit(T_l):
    q_conv=h_lucht*(T_lucht-T_l)
    q_lw = sigma*(e_lucht*(T_lucht**4)-(e_pcm*T_l**4))
    qL = q_conv + q_lw
    return qL, q_conv, q_lw


'Here is the part for the solid phase (when x = 0):'
while T0[n-1]<T_melting and z==0:
     q0[n] = flux_in(T0[n])[0]
     qL[n] = flux_uit(T0[n])[0]
     eta[n]= f_eta(T0[n-1])
     T0[n]=T0[0]+((q0[n]-qL[n])/(cp*m))*t[n]
     T_l[n]=T0[n]
     if T0[n] <= T_melting-0.0005:
         n+=1
     else:
         x[n]=0.000001
         z=1
          
'Here is the part for when the melting front starts to run(when x > 0):'
while (x[n] > 0 and x[n] < L) and z==1:#T0[n]>=T_melting: #and z==1:
    q0[n] = flux_in(T0[n-1])[0]
    qL[n] = flux_uit(T_l[n-1])[0]
    q[n]=q0[n]+qL[n]
    eta[n]= f_eta(T0[n-1])
    
    T0[n]=((q0[n]*x[n])/k_l) +T_melting
    T_l[n]=((qL[n]*(L-x[n]))/k_s) +T_melting
    
    if  x[n] > 0 and x[n] < L:
        v_melting[n] = (-k_l*(((T_melting-T0[n])/(x[n]))+k_s*((T_l[n]-T_melting)/(L-x[n])))/(rho*h_sl))
    else:
        v_melting[n]= 0
    if T_l[n] < (T_melting-0.1) and n < t_stop:
        x[n+1]=x[n]+v_melting[n]*dt
        n += 1
    else:
        T0_eind = T0[n]
        z=2
                    
dollar =0.49 # Average sale price 49 pennies per kWh
cost_zonnepaneel= 140 #Cost solar panel
'lists for the calculation'
e_1= np.zeros(t_stop+1)
opbrengst=np.zeros(t_stop+1)
average_eta= sum(eta)/len(t)
e_1 =I*average_eta
opbrengst = e_1*dollar #per day
cost_pcm= y*1000*A*L*rho
cost_total= cost_zonnepaneel+cost_pcm
        # 0= cost_total- opbrengst*s
s=cost_total/opbrengst  #time in days until breakeven point
print('tijd in dagen tot breakeven point:',s)    
print('E_PCM bespaard',e_1)

'Functions for all plots:'
plt.figure(1)   
plt.plot(t, T0, label='Stof 3' )
plt.title("Fatty Acids",fontdict={'fontname':'Comic Sans MS','fontsize':20})
plt.xlabel("time in s")
plt.ylabel('Temperature in K')
plt.legend()
plt.show()

plt.figure(2)
plt.plot(t,v_melting, label='Stof 3')
plt.title("Fatty Acids",fontdict={'fontname':'Comic Sans MS','fontsize':20})
plt.xlabel("time in s")
plt.ylabel("v_melting in m/s")
plt.legend()
plt.show()

plt.figure(3)
plt.plot(t,x, label='Stof 3')
plt.title("Fatty Acids",fontdict={'fontname':'Comic Sans MS','fontsize':20})
plt.xlabel("time in s")
plt.ylabel("x in m")
plt.legend()
plt.show()

plt.figure(4)
plt.plot(t,T_l, label='Stof 3')
plt.title("Fatty Acids",fontdict={'fontname':'Comic Sans MS','fontsize':20})
plt.xlabel("time in s")
plt.ylabel("Temperature at x=L in K")
plt.legend()
plt.show()

