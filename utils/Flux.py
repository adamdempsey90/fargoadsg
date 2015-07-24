class Vars3():
	def __init__:
		Mesh.__init__(self) #All the Mesh attributes inside Field!
        	Parameters.__init__(self) #All the Parameters attributes inside Field!
		self.dens,self.vr,self.vp,,self.sig,self.u,self.v = loadvars(self,num)
		self.omegaframe,self.mp = loadtxt(directory+'planet0.dat')[0,[5,7]]
		self.cprof = transpose(tile(self.aspectratio*self.omegaframe*pow(self.y,self.flaringindex-.5),[self.nx,1]))
		self.rsoft =  transpose(tile(self.thicknesssmoothing*self.y*self.aspectratio * pow(self.y,self.flaringindex),[self.nx,1]))


class Flux():
	def __init__:
		self.Twd = calcTwd(dat)
		self.Fp, self.drFp = calcFp(dat)
		self.Th = calcTh(dat)

def loadvars(dat,num):
	vr = fromfile('gasvy'+str(num)+'.dat').reshape(self.ny,self.nx)
	vp = fromfile('gasvx'+str(num)+'.dat').reshape(self.ny,self.nx)
	dens = fromfile('gasdens'+str(num)+'.dat').reshape(self.ny,self.nx)

	u = self.vr - transpose(tile(self.vr.mean(axis=1),[self.nx,1]))
	v = self.vp - transpose(tile(self.vp.mean(axis=1),[self.nx,1]))
	sig = self.dens - transpose(tile(self.dens.mean(axis=1),[self.nx,1]))

	return dens,vr,vp,sig,u,v

def calcTwd(dat):
	

	vr = dat.vr
	vp = dat.vp
	dens = dat.dens
	vrbar = transpose(tile(vr.mean(axis=1),[dat.nx,1]))
	vpbar = transpose(tile(vp.mean(axis=1),[dat.nx,1]))
	dbar = transpose(tile(dens.mean(axis=1),[dat.nx,1]))

	u = dat.u
	v = dat.v
	sig = dat.sig

	T,R = meshgrid(dat.x,dat.y)


	
	return 0

def calcFp(dat):

	vr = dat.vr
	vp = dat.vp
	dens = dat.dens
	vrbar = transpose(tile(vr.mean(axis=1),[dat.nx,1]))
	vpbar = transpose(tile(vp.mean(axis=1),[dat.nx,1]))
	dbar = transpose(tile(dens.mean(axis=1),[dat.nx,1]))

	u = dat.u
	v = dat.v
	sig = dat.sig

	T,R = meshgrid(dat.x,dat.y)	





	return 0


def calcTh(dat):

	T,R = meshgrid(dat.x,dat.y)
	mp = dat.mp
	rs = dat.rsoft
	dpPhi = mp*R*sin(T)*pow(R*R+1+rs*rs-2*R*cos(T),-1.5)
	
	Th = dat.sig*dpPhi

	return Th.mean(axis=1)



