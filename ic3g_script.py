# Effect of mineralogy and morphology on thermal conductivity of sand

# 1. Input parameters:
# surface roughness of particles (assumed to be radius of sub-particle)
Ra = 0.4
# minimum and maximum CLUSTER particle size (length)
Dmin = 0.5
Dmax = 1.5
# number of different CLUSTER particle sizes
nClust = 101
# CLUSTER particle aspect ratio
ar = 0.5
# Set constituent fractions
qrtFrac=0.6
felFrac=0.2
micFrac=0.2
# initial, source & sink temperature
Tinit = 0.0
Tsrc = 10.0
Tsnk = 0.0
# Isotropic compaction stress
sigmaIso=-1e7

# 2. Simulation setup:
# packages
from yade import pack, qt, plot 
import random, math
# periodicity
O.periodic=True
O.cell.refSize=(10*Dmax,10*Dmax,10*Dmax)
# loose particle packing
sp=pack.SpherePack()
# create different sized clusters
# hexagonal packing
Ca = []
for i in range(int(nClust)):
	size=Dmin+float(i)/float(nClust-1)*(Dmax-Dmin)
	Cai = pack.regularHexa(inEllipsoid((0,0,0),(size,size,size*ar)), Ra, -0.2*Ra)
	Caj=[]
	for j in range(len(Cai)):
		Caj.append((Cai[j].state.pos,Cai[j].shape.radius))
	Ca.append(Caj)
Cb = range(int(nClust))
for i in range(int(nClust)):
	Cb[i] = pack.SpherePack(Ca[i])
sp.makeClumpCloud((0,0,0),(10*Dmax,10*Dmax,10*Dmax),Cb,periodic=True,num=5000)
sp.toSimulation()
# Compaction completion temperature start variable
compDone=False
compIter=0.0
tempStart=False
# total heat
totHeat=0.0

# 3. Materials:
# a. Quartz:
qrtMat=O.materials.append(FrictMat(density=(2.65e3),young=9.5e10,poisson=0.2,frictionAngle=0.5))
O.materials[qrtMat].specHeat=830*1e-11
O.materials[qrtMat].thrmCond=7.7
# b. Feldspar:
felMat=O.materials.append(FrictMat(density=(2.56e3),young=4.0e10,poisson=0.3,frictionAngle=0.6))
O.materials[felMat].specHeat=750*1e-11
O.materials[felMat].thrmCond=2.5
# c. Mica:
micMat=O.materials.append(FrictMat(density=(2.97e3),young=1.8e11,poisson=0.25,frictionAngle=0.3))
O.materials[micMat].specHeat=770*1e-11
O.materials[micMat].thrmCond=2.1
# d. Boundary
bndMat=O.materials.append(FrictMat(density=(2.65e3),young=9.5e10,poisson=0.2,frictionAngle=0.5))
O.materials[bndMat].specHeat=0
O.materials[bndMat].thrmCond=7.7

# 4. Assign materials:
# Identify clusters and assign materials
lClust=[]
for i in O.bodies:
	if i.isClump==True:
		nRand=random.random()
		idMat=0
		if nRand<qrtFrac: idMat=qrtMat
		elif nRand<qrtFrac+felFrac: idMat=felMat
		else: idMat=micMat
		lClust.append([i.id,idMat])
# Assign materials to individual members of clusters
for i in O.bodies: 
	if i.isClumpMember==True:
		for j in lClust:
			if i.clumpId==j[0]:
				i.material=O.materials[j[1]]

# 5. Assign initial temperature
for i in O.bodies:
	i.state.temp=Tinit
	i.state.heat=0.0

# Set time increment
mass = 4/3*math.pi*Ra**3*min(O.materials[qrtMat].density,O.materials[felMat].density,O.materials[micMat].density,O.materials[bndMat].density)
minSpec=min(O.materials[qrtMat].specHeat,O.materials[felMat].specHeat,O.materials[micMat].specHeat)
maxCond=max(O.materials[qrtMat].thrmCond,O.materials[felMat].thrmCond,O.materials[micMat].thrmCond,O.materials[bndMat].thrmCond)
Tdt=mass*minSpec/(2*Ra*maxCond)

# 6. Simulation loop
O.engines=[
   ForceResetter(label='forceReset'),
   InsertionSortCollider([Bo1_Sphere_Aabb()]),
   InteractionLoop(
      [Ig2_Sphere_Sphere_ScGeom()],
      [Ip2_FrictMat_FrictMat_FrictPhys()],
      [Law2_ScGeom_FrictPhys_CundallStrack()]
   ),
   # run initial compaction
   #PyRunner(command='checkStress()',iterPeriod=10),
   PeriTriaxController(label='triax',
      # specify target values and whether they are strains or stresses
      goal=(sigmaIso,sigmaIso,sigmaIso),stressMask=7,
      # type of servo-control
      dynCell=True,maxStrainRate=(100,100,100),
      # wait until the unbalanced force goes below this value
      maxUnbalanced=.1,relStressTol=1e-3,
      # call this function when goal is reached and the packing is stable
      doneHook='compactionFinished()'
   ),
   NewtonIntegrator(damping=.6, label='newton'),
   # run update temperatures
   PyRunner(command='updateTemp()',iterPeriod=(int(Tdt/O.dt) if int(Tdt/O.dt)>0 else 1)),
   # record data for plotting every 100 steps; addData function is defined below
   PyRunner(command='addData()',iterPeriod=100)
]

# 6. Timestep to be 1/2 of the "critical" timestep
O.dt=.5*PWaveTimeStep()

# 7. Prescribe isotropic normal deformation (constant strain rate)
# of the periodic cell
#O.cell.velGrad=Matrix3(-10,0,0, 0,-10,0, 0,0,-10)

# 8. initial compaction function
def compactionFinished():
   # set the current cell configuration to be the reference one
   O.cell.trsf=Matrix3.Identity
   # change control type: keep constant confinement in x,y, 20% compression in z
   triax.goal=(sigmaIso,sigmaIso,0)
   triax.stressMask=3
   # allow faster deformation along x,y to better maintain stresses
   triax.maxStrainRate=(1.,1.,1.)
   triax.maxUnbalanced=0.001
   global compDone
   global compIter
   #if compDone==False:
	# stress tensor as the sum of normal and shear contributions
	# Matrix3.Zero is the intial value for sum(...)
	#stress=sum(normalShearStressTensors(),Matrix3.Zero)
	# if mean stress is below (bigger in absolute value) limitMeanStress, start shearing
	#if stress.trace()/3.>sigmaIso:
	   #O.cell.velGrad=Matrix3(0,0,0, 0,0,0, 0,0,0)
	   #print ("Isostatic stress target of" + str(sigmaIso) + "reached")
   compDone=True
   newton.damping=0.9
   #O.pause()


# 9. temperature update function
def updateTemp():
   global tempStart
   global totHeat
   if (compDone==True and tempStart==False):
	# Assign initial boundary temperatures and materials
	for i in O.bodies:
	   if i.isClumpMember==True and i.state.pos[2]>O.cell.size[2]-Ra*4:
		i.material=O.materials[bndMat]
		i.state.temp=Tsrc
		i.shape.color=(1,0,0)
	   elif i.isClumpMember==True and i.state.pos[2]<Ra*4:
		i.material=O.materials[bndMat]
		i.state.temp=Tsnk
		i.shape.color=(0,0,1)
	tempStart=True
	# Set time increment
	mass = 4/3*math.pi*Ra**3*min(O.materials[qrtMat].density,O.materials[felMat].density,O.materials[micMat].density,O.materials[bndMat].density)
	minSpec=min(O.materials[qrtMat].specHeat,O.materials[felMat].specHeat,O.materials[micMat].specHeat)
	maxCond=max(O.materials[qrtMat].thrmCond,O.materials[felMat].thrmCond,O.materials[micMat].thrmCond,O.materials[bndMat].thrmCond)
	O.dt = min(mass*minSpec/(2*Ra*maxCond),.5*PWaveTimeStep())
   elif tempStart==True:
	for i in O.interactions:
	   if math.fabs(O.bodies[i.id1].state.pos[2]-O.bodies[i.id2].state.pos[2])<=0.5*O.cell.size[2]:
	   	 Reff=1/(1/i.geom.refR1+1/i.geom.refR2)
		 Rcont=(Reff*i.geom.penetrationDepth)**0.5
		 ks=2/(1/O.bodies[i.id1].material.thrmCond+1/O.bodies[i.id2].material.thrmCond)
		 O.bodies[i.id1].state.heat+=-Rcont*ks*(O.bodies[i.id1].state.temp-O.bodies[i.id2].state.temp)*O.dt
		 O.bodies[i.id2].state.heat+=Rcont*ks*(O.bodies[i.id1].state.temp-O.bodies[i.id2].state.temp)*O.dt
	totHeat=0.0
	for i in O.bodies:
	   if i.isClumpMember==True and i.state.pos[2]>O.cell.size[2]-Ra*4:
	   	 totHeat+=i.state.heat
		 i.state.heat=0.0
	   elif i.isClumpMember==True and i.material.specHeat!=0:
		 i.state.temp+=i.state.heat/(i.state.mass*i.material.specHeat)
		 i.state.heat=0.0
	totHeat/=O.dt
	
# 10. Data function
def addData():
   stress=sum(normalShearStressTensors(),Matrix3.Zero)
   sigzz=stress[2,2]
   q1T=0.0
   q1N=0
   q2T=0.0
   q2N=0
   q3T=0.0
   q3N=0
   q4T=0.0
   q4N=0
   for i in O.bodies:
      if i.isClumpMember==True and i.state.pos[2]<0.25*O.cell.size[2]:
   	q1T+=i.state.temp
   	q1N+=1
      elif i.isClumpMember==True and i.state.pos[2]<0.5*O.cell.size[2]:
   	q2T+=i.state.temp
   	q2N+=1
      elif i.isClumpMember==True and i.state.pos[2]<0.75*O.cell.size[2]:
   	q3T+=i.state.temp
   	q3N+=1
      elif i.isClumpMember==True:
   	q4T+=i.state.temp
   	q4N+=1
   plot.addData(szz=sigzz,heat=totHeat,T1=q1T/q1N,T2=q2T/q2N,T3=q3T/q3N,T4=q4T/q4N,t=O.time)

plot.plots={'t':('szz'),'t ':('heat'),'t  ':('T1','T2','T3','T4')}
plot.plot()




