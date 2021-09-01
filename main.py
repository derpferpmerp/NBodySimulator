import math
import random
import numpy as np
import json
import matplotlib.pyplot as plt

CGRAV = 6.67428 * (math.pow(10,-11))

def planetForce(m1,m2,xd,yd):
	distance = math.sqrt((xd**2)+(yd**2))
	force2reach = (CGRAV*m1*m2)/(distance**2)
	return [force2reach/m1,force2reach/m2]

def calcVirialRadius(cval=200,err=5,nparts=100):
	distances = [random.randrange(10,1000) for x in range(nparts)]
	masses = [random.randrange(10,100) for x in range(nparts)]
	good=True
	avgdensity = sum(masses)/sum(distances)
	volumes = []
	for x in np.linspace(1,1000,10000):
		if good:
			vMasses = [masses[x] for x in range(len(masses)) if distances[x] < x]
			volume = sum(vMasses)/((4/3)*(np.pi)*(x**3))
			volumes.append(round(volume,7))
			if round(volume,7) > avgdensity - 0.1 and round(volume,7) < avgdensity + 0.1:
				good=False
				outl = [[distances[g],masses[g]] for g in range(len(volumes)) if volumes[g] < x]
				return [x,outl]
				
def randFromDistance(d,F,m):
	theta = random.randrange(1,359) * (np.pi)
	theta /= 180
	plt.scatter(x=d*math.cos(theta),y=d*math.sin(theta),c='r')
	plt.plot([d*math.cos(theta),(d-50)*math.cos(theta)],[d*math.sin(theta),(d-50)*math.sin(theta)],c="black")
	return {
		"X": d * math.cos(theta),
		"Y": d * math.sin(theta),
		"FORCES": [(F/m)*math.cos(theta),(F/m)*math.sin(theta)],
		"THETA": theta
	}

def indPairs(sIndex,mx):
	outlst = []
	for x in range(mx):
		if x == sIndex: continue
		try:
			outlst.append([sIndex,x])
		except:
			continue
	return outlst

def calcEpsilon(nparts=100,cval=200,err=5,dev=False,tht=None):
	rvir, pnts = calcVirialRadius(cval=cval,err=err,nparts=nparts)
	forces = []
	epsilon = math.pow(((4*np.pi*(rvir**3))/3),(1/3)) * math.pow(nparts,(-1/3))
	for x in range(len(pnts)-1):
		ngS = []
		for g in indPairs(x,len(pnts)):
			mi, mj = [pnts[x][1],pnts[g[0]][1]]
			ri, rj = [pnts[x][0],pnts[g[1]][0]]
			nG = -CGRAV
			nG *= (mi*mj)
			nG *= (ri-rj)
			nG /= math.pow(((abs(ri-rj)**2) + (epsilon**2)),(3/2))
			ngS.append(nG)
		nG = sum(ngS)/len(ngS)
		forces.append({"MASS":mi,"D_ABS":ri,"CRDS":randFromDistance(ri,nG,mi),"ACC":nG/mi})
	if dev:
		return [(nG/mi)*math.cos(tht),(nG/mi)*math.sin(tht)]
	return forces


print(calcVirialRadius())
print(json.dumps(calcEpsilon(),indent=4))
vels = []
def lp(h=None, tvar=2, itr=0):
	h=[calcEpsilon() if h == None else h][0]
	for val in h:
		lunsquared2 = [[(val["CRDS"]["X"]-v["CRDS"]["X"])**2+(val["CRDS"]["Y"]-v["CRDS"]["Y"])**2][0] for v in h if all([v!=val])]
		accs = [v["ACC"] for v in h if v != val]
		ls = [[1 if x > 0 else -1][0] * math.sqrt(abs(x)) for x in lunsquared2]
		ls = sum(ls)/len(ls)
		accs = sum(accs)/len(accs)

		vf = ((2*accs*ls)/(abs(2*accs*ls))) * math.sqrt(abs(2*accs*ls))
		frcs = calcEpsilon(dev=True,tht=val["CRDS"]["THETA"])
		fxComp = frcs[0]/sum(frcs)
		fyComp = frcs[1]/sum(frcs)
		aX = accs*fxComp
		aY = accs*fyComp

		t=tvar

		dtX=(1/2)*(aX)*(t**2) + val["CRDS"]["X"]
		dtY=(1/2)*(aY)*(t**2) + val["CRDS"]["Y"]

		print(f"Particle #{h.index(val)+(itr*len(h))}/{(len(h)*4)} -> ({dtX},{dtY}) {vf}m/s^2")
		h[h.index(val)]["CRDS"]["X"] = dtX
		h[h.index(val)]["CRDS"]["Y"] = dtY
		plt.scatter(dtX,dtY,c=["red","green","blue","black"][t])
	return h

hh = lp(tvar=0)
hh1 = lp(h=hh, tvar=1, itr=1)
hhh = lp(h=hh1, tvar=2, itr=2)
hhhh = lp(h=hhh, tvar=3, itr=3)
plt.savefig("out.png")