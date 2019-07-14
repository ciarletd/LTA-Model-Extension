#LTA Barrier Model (Lorenzo-Trueba and Ashton, 2014)
#Recoded in Python by Daniel J. Ciarletta
#This script creates an animation (GIF) of barrier evolution, with a subplot
#Below, please set the desired directory to store temporary image frames (See line 177)
print("-------------")
print("BARRIER MODEL")
print("-------------")
import matplotlib.pyplot as plt
import os
import imageio

#Physical Parameters
B=0.0015
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
Xmax=13000
Tmax=10000
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

#Ploting
atrigger=100 #Set the the interval to capture frames for animation
animate_count=0

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
  
  #Animation Frame
  animate_count=animate_count+1
  if animate_count==atrigger:  
      xl=0 #set limit left
      fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(10, 8), sharex=True, sharey=True, constrained_layout=True) #rows by columns
      plt.subplots_adjust(hspace = 1)
      #subplot 1
      plt.subplot2grid((3,3), (0,0), colspan=3, rowspan=2)
      plt.fill(OceanX,OceanY,"#98F5FF")
      plt.fill(SLX,SLY,"b")
      plt.fill(XShelf,YShelf,"#8B6914")
      plt.fill(XRelict,YRelict,"#CDCD00")
      plt.fill(BarrierX,BarrierY,"#FFFF00")
      plt.plot(BarrierX,BarrierY,'k-',linewidth=0.5)
      plt.plot(XT,ZT,'k-')
      plt.ylabel('Elevation (m)',fontsize=14)
      plt.xlabel('Distance (m)',fontsize=14)
      plt.title('Periodic Barrier Retreat, Time='+str(i+1),fontsize=16)
      plt.xlim(xl, Xmax)
      plt.ylim(0, 40)
      #subplot 2
      plt.subplot2grid((3,1), (2,0))
      plt.plot(OverwashFlux,'k-')
      plt.ylabel('Overwash Flux\n($m^3/m/yr$)',fontsize=14)
      plt.xlabel('Time (yrs)',fontsize=14)
      plt.title('Overwash Flux',fontsize=16)
      plt.xlim(xl, Tmax)
      plt.rc('xtick', labelsize=12) 
      plt.rc('ytick', labelsize=12) 
      plt.show()
      plt.savefig('zzvalue' + str(i/Tmax) + '.png')
      plt.close()
      animate_count=0
   
  #End of for loop

print("-------------------")
print("Rendering Animation")
print("-------------------")
#Compile the GIF
png_dir = 'C:\Python' #Will need to change this if using a different directory
images = []
for file_name in os.listdir(png_dir):
    if file_name.endswith('.png'):
        file_path = os.path.join(png_dir, file_name)
        images.append(imageio.imread(file_path))
imageio.mimsave('BarrierAnimation.gif', images,) #This saves the GIF

#CLEAN UP
#The next bit deletes all PNG images saved from the animation 
#Make sure you don't have any other PNG images in this directoty that you want
rmfiles = os.listdir(png_dir)
for item in rmfiles:
    if item.endswith(".png"):
        os.remove(os.path.join(png_dir, item))


