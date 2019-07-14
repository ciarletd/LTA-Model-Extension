#LTA Barrier Model (Lorenzo-Trueba and Ashton, 2014)
#Recoded in Python by Daniel J. Ciarletta
print("-------------")
print("BARRIER MODEL")
print("-------------")
import numpy as np
import matplotlib.pyplot as plt

#Physical Parameters
B=0.002
Dt=15
We=800
He=2
Ae=0.02
Qow_max=100
Vd_max=0.5*He*We
K=2000
a=0.002
b=0

#Starting Conditions
A=Ae; W=We; H=He #Barrier initially in equilibrium
xt=0; xs=Dt/A; xb=xs+W; xso=xs; Z=Dt
Db=Z-B*xb #Initial Backbarrier Depth
zt=0; zs=Dt; ztop=Dt+H; zb_bot=Dt-Db

#Computational Parameters
Xmax=12000
Tmax=9000
dt=1
zfin=B*Xmax
Height=list()
Width=list()
OverwashFlux=list()
XB=list()
XT=list()
XS=list()
ZT=list()
SFslope=list()
SFFlux=list()
DepositY=list()

#Main Code
for i in range(0,Tmax,dt): #solve i for 0 to Tmax by dt 
  zdot=a
  Z=Z+zdot
  
  Vd_H=(He-H+zdot)*W
  
  if Vd_H<0: 
      Vd_H=0
  Vd_B=(We-W)*(H+Db)
  if Vd_B<0: 
      Vd_B=0
  
  Vd=Vd_H+Vd_B
  
  if Vd<Vd_max: 
      Qow_H=Qow_max*Vd_H/Vd_max
      Qow_B=Qow_max*Vd_B/Vd_max 
  else: 
      Qow_H=Qow_max*Vd_H/Vd 
      Qow_B=Qow_max*Vd_B/Vd
  
  Qow=Qow_H+Qow_B
  Qsf=K*(Ae-A)
  Hdot=Qow_H/W-zdot
  xbdot=Qow_B/(H+Db)
  xsdota=(2*Qow)/(Dt+2*H)
  xsdotb=(4*Qsf*(H+Dt))
  xsdot=xsdota-(xsdotb/((2*H+Dt)**2))
  xtdot=2*Qsf*(1/(Dt+2*H)+1/Dt)+2*zdot/A;
  if H>0 and W>0: 
      H=H+Hdot*dt
      xb=xb+xbdot
      xs=xs+xsdot
      xt=xt+xtdot
      A=Dt/(xs-xt)
      W=xb-xs
      Zb=Z
  Db=Zb-xb*B
  zt=Zb-Dt
  zs=Zb
  ztop=Zb+H
  zb_bot=Zb-Db
  ZT.append(zt)
  XT.append(xt)
  SFFlux.append(Qsf)
  OverwashFlux.append(Qow)
  Width.append(W) 
  
  #Barrier Geometry
  BarrierX=[xt,xs,xs,xb,xb,xb,xt]
  BarrierY=[zt,zs,ztop,ztop,Zb,zb_bot,zt]  
  
  #Ocean Geometry
  OceanX=[0,0,Xmax,Xmax,0]
  OceanY=[0,Z,Z,0,0]  
  
  #Water Level
  SLX=[0,Xmax,Xmax,0]
  SLY=[Z,Z,Z+0.3,Z+0.3]
  
  #Shelf Geometry
  XShelf=XT[1:i]
  XShelf.append(xt)
  XShelf.append(xb)
  XShelf.append(Xmax)
  XShelf.append(Xmax)
  YShelf=ZT[1:i]
  YShelf.append(zt)
  YShelf.append(zb_bot)
  YShelf.append(zfin)
  YShelf.append(0)
  
  #Deposit Geometry
  if zt>(xt*B):
      DepositY.append(zt)
  else:
      DepositY.append(xt*B)
      XRelict=XT[1:i]
      YRelict=DepositY[1:i]
      #XRelict.append(xt)
      XRelict.append(xb)
      XRelict.append(xt)
      #Relict.append(zt)
      YRelict.append(xb*B)
      YRelict.append(xt*B)
  
   
  #End of for loop

np.savetxt('Widthvector.txt',Width,fmt='%5.2f') #save vector in 5 digit format
xlen=len(XT) #get the length of vector
xl=0 #set limit left
xarray=np.linspace(xl,Xmax,xlen) #create an x vector

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8), sharex=True, sharey=True, constrained_layout=True) #rows by columns
plt.subplots_adjust(hspace = 1)

plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=2)
#First plot sets overall size of grid (3 x 3), start loc (0,0)-->(row,column), and then actual dimensions of the first subplot
plt.fill(OceanX,OceanY,"#98F5FF")
plt.fill(SLX,SLY,"b")
plt.fill(XShelf,YShelf,"#8B6914")
plt.fill(XRelict,YRelict,"#CDCD00")
plt.fill(BarrierX,BarrierY,"#FFFF00")
plt.plot(BarrierX,BarrierY,'k-',linewidth=0.5)
plt.plot(XT,ZT,'k-')
plt.ylabel('Elevation (m)',fontsize=14)
plt.xlabel('Time (yrs)',fontsize=14)
plt.title('Periodic Barrier Retreat, Time='+str(i+1),fontsize=16)
plt.xlim(xl, Xmax)
plt.ylim(0, 40)

plt.subplot2grid((3,1), (2,0))
plt.plot(xarray,OverwashFlux,'k-')
plt.ylabel('Overwash Flux\n($m^3/m$)',fontsize=14)
plt.xlabel('Time (yrs)',fontsize=14)
plt.title('Overwash Flux',fontsize=16)
plt.xlim(xl, Tmax)

plt.rc('xtick', labelsize=12) 
plt.rc('ytick', labelsize=12) 
plt.show()
