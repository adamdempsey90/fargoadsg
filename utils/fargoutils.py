# Collection of useful functions to analyze fargo output
# Note numpy,matplotlib,os, & copy should be imported already 

class Vars2():
	def __init__(self,num):
		self.rho = Field('gasdens'+str(num)+'.dat')
		self.vr = Field('gasvy'+str(num)+'.dat')
		self.vp = Field('gasvx'+str(num)+'.dat')
		print 'rho,vr,vp'
		self.T = Tensor(self.rho,self.vr,self.vp)
		print 'T'
		self.dbar,self.vrbar,self.vpbar, self.sig, self.u, self.v = datsplit(self.rho,self.vr,self.vp)	
		print 'split'
		self.Tbar = Tensor(self.dbar,self.vpbar,self.vpbar)
		print 'Tbar'
		self.Tpp = Tensor(self.sig,self.u,self.v)
		print 'Tpp'
		self.Tp = pTensor(self.T,self.Tbar,self.Tpp)
		self.flux = Flux(self)
class Flux():
	def __init__(self,dat):
		self.twd = calcTwd(dat)
		self.drFp = calcFp(dat)
		self.Th = calcTh(dat)

class Vars():

	def __init__(self,num):
		self.params = Parameters()
		self.r,self.phi,self.d,self.vx,self.vy,self.dbar,self.vxbar,self.vybar,self.sig,self.u,self.v=loadvars(num,varsreturn=True)
		self.pp,self.rr = meshgrid(self.phi,self.r)
		self.x = self.rr*cos(self.pp)
		self.y = self.rr*sin(self.pp)
		self.dr= mean(diff(self.r))
		self.dphi = mean(diff(self.phi))
		self.vt =VT(self)		
		self.sig0 = self.params.sigma0*pow(self.r,-self.params.sigmaslope)
		self.vort = calcvort(self)				
class VT():
	def __init__(self,dat):
		P = dat.params
		self.Pi=zeros((P.ny,P.nx,2,2))
		
		vx = dat.vx; vy=dat.vy; d=dat.d;
		r = dat.rr; phi = dat.pp;
		om,mp = loadtxt('planet0.dat')[0,[5,7]]
		c = P.aspectratio*om*pow(r,P.flaringindex-.5)
		nu = P.nu
		
		dr = mean(diff(dat.r))
		dphi = mean(diff(dat.phi))
	
		dpvx, drvx = gradient(vx,dphi,dr)
		dpvy, drvy = gradient(vy,dphi,dr)
		dpd, drd = gradient(d,dphi,dr)


		self.Pi[:,:,0,0] = d*(-c**2 + nu*(2.0/3.0)*(2*drvx - (dpvy/r) - (vx/r) ))
		self.Pi[:,:,0,1] = nu*d*(drvy + (dpvx/r) - (vy/r))
		self.Pi[:,:,1,0] = self.Pi[:,:,0,1]
		self.Pi[:,:,1,1] = d*(-c**2+nu*(2.0/3.0)*(-drvx + 2*(dpvy/r)+2*(vx/r)))	

		self.divPi = zeros((P.ny,P.nx,2))

		self.drPi = zeros((P.ny,P.nx,2,2))
		self.dpPi = zeros((P.ny,P.nx,2,2))
		
		for i in range(2):
			for j in range(2):
				self.dpPi[:,:,i,j], self.drPi[:,:,i,j] = gradient(self.Pi[:,:,i,j],dphi,dr)

		self.divPi[:,:,0] = self.drPi[:,:,0,0] + (self.Pi[:,:,0,0] - self.Pi[:,:,1,1] + self.dpPi[:,:,0,1])/r
		self.divPi[:,:,1] = self.drPi[:,:,0,1] + (self.Pi[:,:,0,1] + self.dpPi[:,:,1,1] + self.Pi[:,:,1,0])/r 

class VTens():
	def __init__(self):
		self.rr = 0
		self.pp = 0
		self.rp = 0
		self.divr = 0
		self.divp = 0
		self.rdivr = 0
		self.rdivp = 0

def calcvort(dat):
	
	planet = loadtxt('planet0.dat')
	omega = planet[0,5]

	dr = mean(diff(dat.r))
	dphi = mean(diff(dat.phi))

	dpvx, drvx = gradient(dat.vx,dphi,dr)
	dpvy, drvy = gradient(dat.vy,dphi,dr)

	curlz = drvy + dat.vy/dat.rr - dpvx
	
	vort = curlz/dat.d
	return vort
		
def loadvars(num, plotflag=False, varsreturn=False):
	print os.getcwd()	
	P=Parameters()
	d=fromfile("gasdens"+num+".dat").reshape(P.ny,P.nx)
	vx = fromfile("gasvy"+num+".dat").reshape(P.ny,P.nx)
	vy = fromfile("gasvx"+num+".dat").reshape(P.ny,P.nx)
	
	dbar=d.mean(axis=1)
	vxbar = vx.mean(axis=1)
	vybar = vy.mean(axis=1)
	
	domain_x=loadtxt('domain_y.dat')[3:-3]
	domain_y=loadtxt('domain_x.dat')
	
	xm = 0.5*(domain_x[:-1]+domain_x[1:])
	ym = 0.5*(domain_y[:-1]+domain_y[1:])

	vyback = transpose(tile(vybar,[P.nx, 1]))
	vxback = transpose(tile(vxbar,[P.nx, 1]))
	dback = transpose(tile(dbar,[P.nx, 1]))
	
	dprime = d-dback
	vxprime = vx-vxback
	vyprime = vy - vyback


	print 'Min Mean Surface Density:  ', dbar.min() 
	
	if (plotflag):
		figure()
		plot(xm,dbar)
		figure()
		plot(xm,vxbar)
		figure()
		plot(xm,vybar)
	
	if (not(varsreturn)):
		return P,xm,ym,d,vx,vy,dbar,vxbar,vybar,dprime,vxprime,vyprime
	else:
		return xm,ym,d,vx,vy,dbar,vxbar,vybar,dprime,vxprime,vyprime
		
def calcfluxes(dat):
	P = dat.params
	dbar = swapaxes(tile(dat.d.mean(axis=1),[P.nx,1]),0,1)
	vxbar = swapaxes(tile(dat.vx.mean(axis=1),[P.nx, 1]),0,1)
	vybar = swapaxes(tile(dat.vy.mean(axis=1),[P.nx, 1]),0,1)
	Pibar = swapaxes(tile(dat.vt.Pi.mean(axis=1),[P.nx, 1,1,1]),0,1)

	rr= dat.rr	

	u = dat.u
	v = dat.v
	sig = dat.sig
	Pip = dat.vt.Pi - Pibar

	dr = mean(diff(dat.r))
	dphi = mean(diff(dat.phi))

	dpPibar=zeros(dat.vt.Pi.shape)
	drPibar=zeros(dat.vt.Pi.shape)
	dpPip=zeros(dat.vt.Pi.shape)
	drPip=zeros(dat.vt.Pi.shape)
	rdivPibar=zeros((P.ny,P.nx,2))
	rdivPip=zeros((P.ny,P.nx,2))

	print Pibar.shape, Pip.shape, rdivPibar.shape, drPibar.shape
	
	dru, dpu = gradient(u,dr,dphi)
	drv, dpv = gradient(v,dr,dphi)
	drsig,dpsig = gradient(v,dr,dphi)

	drdbar,dpdbar = gradient(dbar,dr,dphi)
	drvxbar,dpvxbar = gradient(vxbar,dr,dphi)
	drvybar,dpvybar = gradient(vybar,dr,dphi)

	for i in range(2):
		for j in range(2):
			drPibar[:,:,i,j], dpPibar[:,:,i,j] = gradient(Pibar[:,:,i,j],dr,dphi)
			drPip[:,:,i,j], dpPip[:,:,i,j] = gradient(Pip[:,:,i,j],dr,dphi)

	rdivPibar[:,:,1] = rr*drPibar[:,:,0,1] + Pibar[:,:,0,1] + Pibar[:,:,1,0]
	rdivPibar[:,:,0] = rr*drPibar[:,:,0,0] + Pibar[:,:,0,0] - Pibar[:,:,1,1]
	
	rdivPip[:,:,0] = rr*drPip[:,:,0,0] + Pip[:,:,0,0] - Pip[:,:,1,1] + dpPip[:,:,0,1]
	rdivPip[:,:,1] = rr*drPip[:,:,0,1] + Pip[:,:,0,1]+Pip[:,:,1,0] + dpPip[:,:,1,1]

 
#	nu = P.alphaviscosity * P.aspectratio * pow(rr,2*P.flaringindex - 1)

	nu = P.nu
	om,mp = loadtxt('planet0.dat')[0,[5,7]]
	c = P.aspectratio*om*pow(rr,P.flaringindex-.5)
	rs = P.thicknesssmoothing*rr*P.aspectratio * pow(rr,P.flaringindex)



	
	FFp = rr*vxbar*sig*v + rr*dbar*u*v
	Fp = FFp.mean(axis=1)
	drFp = gradient(Fp,dr)

	TTwd = rr*dbar*u*dpv - rr*sig*u*(drvybar+2*om) - (sig*sig/(dbar*dbar))*rdivPibar[:,:,1] + (sig/dbar)*rdivPip[:,:,1]   
	Twd = TTwd.mean(axis=1)

	Phi = -mp/sqrt(rr*rr + 1 + rs*rs-2*rr*cos(dat.pp))
	dpPhi = mp*rr*sin(dat.pp)*pow(rr*rr+1+rs*rs-2*rr*cos(dat.pp),-1.5)
	drPhi = mp*(rr-cos(dat.pp))*pow(rr*rr+1+rs*rs-2*rr*cos(dat.pp),-1.5)
#	figure(); imshow(Phi); colorbar();


	TTh = sig*dpPhi
	Th = TTh.mean(axis=1)

	EEp = sig*v*vxbar + sig*u*vybar
	Ep = EEp.mean(axis=1)

	

#	FFb = rr*(dbar*vxbar+sig*u)*(vybar+2*om*rr) - rr*(Pibarrphi+Pipprphi)
#	Fb = FFb.mean(axis=1)
#	drFb = gradient(Fb,dr)

#	EEb = dbar*(vxbar*vybar+u*v) - Pibarrphi + Pipprphi
#	Eb= EEb.mean(axis=1)


#	drPibar = gradient((rr*Pibarrphi).mean(axis=1),dr)

	figure()
	plot(dat.r,drFp+Th,'x',dat.r,Twd)
	xlabel('r'); legend(('drFp+Th','Twd'))

	figure()
	plot(dat.r,drFp+Th+Ep,'x',dat.r,Twd)
	xlabel('r'); legend(('drFp+Th+Ep','Twd'))
	return Fp,drFp,Th,Twd,Ep,
#,drPibar,Ep,Eb,Fb,drFb

def plotcylin(R,phi,Q, logscale=False,fname='',titlestr=''):
	rr,pp = meshgrid(R,phi)
	xx = rr*np.cos(pp)
	yy = rr*np.sin(pp)

	xx = transpose(xx)
	yy = transpose(yy)

	figure()	
	if titlestr!='':
		title(titlestr)

	if (logscale):
		pcolor(xx,yy,log10(Q)); colorbar()
	else:
		pcolor(xx,yy,Q); colorbar();

	if fname!='':
		savefig(fname,bbox_inches='tight')
		close()		


def timeseries_sig(startt, endt,numavg=50,plotflag=False):
	P=Parameters()
	
	time = [P.dt*P.ninterm*i/(2*np.pi) for i in range(startt,endt+1)]
	sig =[fromfile("gasdens"+str(i)+".dat").reshape(P.ny,P.nx).mean(axis=1).min() for i in range(startt,endt+1)]
	sigavg = sig
	sigavg[numavg:len(sig)-numavg-1] = [ mean(sig[(i-numavg):(i+numavg+1)]) for i in range(numavg,len(sig)-numavg)]

	sigavg[len(sig)-numavg:] = [mean(sig[len(sig)-numavg-1:i]) for i in range(len(sig)-numavg,len(sig))]

	print len(sigavg), len(sig)
	if (plotflag):
		figure()
		semilogy(time,sig)
		xlabel('time (orbits)')
		ylabel('$\Sigma_{min}$')
#		ylim([.0004,1])
		
		figure()
		semilogy(time,sigavg)
		xlabel('time (orbits)')
		ylabel('$<\Sigma_{min}>$')
#		ylim([.0004,1])
		
		
	return sigavg,time


def findshocks(dat,vt,plotflag=False):
	
	planet = loadtxt('planet0.dat')
	omega = planet[0,5]

	dr = mean(diff(dat.r))
	dphi = mean(diff(dat.phi))

	dpvx, drvx = gradient(dat.vx,dphi,dr)
	dpvy, drvy = gradient(dat.vy,dphi,dr)

	curlz = drvy + dat.vy/dat.rr - dpvx
	
	vort = (curlz + 2*omega)/dat.d

	dpdivP=zeros((dat.params.ny,dat.params.nx,2))
	drdivP=zeros((dat.params.ny,dat.params.nx,2))

	dpdivP[:,:,0], drdivP[:,:,0]  = gradient(vt.divPi[:,:,0],dphi,dr)

	dpdivP[:,:,1], drdivP[:,:,1]  = gradient(vt.divPi[:,:,1],dphi,dr)

	dplogd, drlogd = gradient(log(dat.d),dphi,dr)

	visc = (dat.rr*drdivP[:,:,1] + vt.divPi[:,:,1] - dpdivP[:,:,0] + vt.divPi[:,:,0]*dplogd ) / (dat.d*dat.d*dat.rr)

	res = vort - visc	
	if (plotflag):
		figure()
		pcolor(dat.rr,dat.pp,res); colorbar()
		xlabel('r'); ylabel('$\phi$')

		plotcylin(dat.r,dat.phi,res,titlestr='Residual')
		plotcylin(dat.r,dat.phi,vort,titlestr='Vortensity')
	return res,vort,visc


def makemovie(tstart,tend,stride,fbase='d'):
	
	time=[tstart+i*stride for i in range(tend-tstart+1)]

	for t in time:
		d=Vars(str(t))
		plotcylin(d.r,d.phi,d.d,logscale=True,fname=fbase+str(t),titlestr='t='+str(t))
	
			
def loadmany(dirs,times):
	
	dat=[]
	
	for i,d in enumerate(dirs):
		os.chdir(d)
		dat.append(Vars(str(times[i])))
		os.chdir('../')


	return dat



def finitediff(dx,Q,order):
	
	s = Q.shape 
	ng = order/2


	if order==2:
		coefsc=array([-.5,0, .5])
		posc = array([-1,0,1])
		coefsf=array([-1.5,2.,-.5])
		posf = array([0,1,2])
	elif order==4:
		coefsc=array([1./12, -2./3, 0,2./3, -1./12])
		posc = array([-2, -1,0, 1, 2])
		coefsf = array([-25./12,4.,-3.,4./3,-1./4])
		posf = array([0,1,2,3,4])
	elif order==6:
		coefsc=array([-1./60, 3./20, -1./4,0, 1./4, -3./20, 1./60])
		posc = array([-3, -2, -1,0, 1, 2,3])
		coefsf =array([-49./20,6.,-15./2,20./3,-15./4,6./5,-1./6])
		posfc = array([0,1,2,3,4,5,6])
	elif order==8:
		coefsc=array([1./280, -4./105, 1./5, -4./5,0, 4./5, -1./5, 4./105, -1./280])
		posc = array([-4, -3, -2, -1,0, 1, 2, 3,4])
		coefsf = array([-761./280,8.,-14.,56./3,-35./2,56./5,-14./3,8./7,-1./8])
		posf = array([0,1,2,3,4,5,6,7,8])
	elif order==10:
		coefs=array([-2., 25., -150., 600., -2100., 2100., -600., 150., -25., 2.])
		coefs /= 2520
		pos=array([-5, -4, -3, -2, -1, 1, 2, 3, 4, 5])


#	bcl = Q[:ng][::-1]
#	bcr = Q[-ng:][::-1] 

#	temp = hstack((bcl,dbar,bcr))

	out = zeros(Q.shape)

	ind= array(range(ng,len(Q)-ng))
	bcl = array(range(ng))
	bcr = array(range(len(Q)-ng,len(Q)))

	for j in range(2*ng+1):
		out[bcl] += Q[bcl+posf[j]]*coefsf[j]
		out[ind] += Q[ind+posc[j]]*coefsc[j]
		out[bcr] += Q[bcr-posf[j]]*-coefsf[j]

	out /= dx
#	temp2=zeros(temp.shape)
#	
#	ind = array(range(len(Q)))
#
#	for j in range(2*ng):
#		temp2[ind+ng] += temp[ind+ng+pos[j]]*coefs[j]
#
#	out = temp2[ind+ng]/dx
	return out


def grad(field):
	T,R = meshgrid(field.x,field.y)
	dR = R[1:,1:] - R[:-1,1:]
	dT = T[1:,1:] - T[1:,:-1]

	dRrfield = copy.deepcopy(field)
	dRrfield.nx -= 1
	dRrfield.ny -= 1
	dRrfield.x = field.x[:-1]
	dRrfield.y = field.y[:-1]
	T1,R1 = meshgrid(dRrfield.x,dRrfield.y)

	dRfield = copy.deepcopy(dRrfield)
	dTfield = copy.deepcopy(dRrfield)
	invrdTfield = copy.deepcopy(dRrfield)
	rfield = R*(field.data)
	dRrfield.data = (rfield[1:,1:]-rfield[:-1,1:])/dR
	dRfield.data = (field.data[1:,1:]-field.data[:-1,1:])/dR
	dTfield.data = (field.data[1:,1:] - field.data[1:,:-1])/dT
	invrdTfield.data = dTfield.data/R1

	
	return dRfield,dTfield,dRrfield,invrdTfield

def divV(vr,vphi):
	
	div = copy.deepcopy(vr)
	div.nx -= 1
	div.ny -= 1
	div.x = vphi.x[:-1]
	div.y = vr.y[:-1]
	print 'copy1'	
	_,_,dRrvr,_ = grad(vr)
	print 'grad1'	
	_,dTvp,_,_ = grad(vphi)
	print 'grad2'
	T,R = meshgrid(div.x,div.y)
		
	div.data = (dRrvr.data + dTvp.data)/R
	
	return div

def Tensor(rho,vr,vp):
	tens = VTens()
	print 'alloc'
	nu = rho.nu
	taurr = copy.deepcopy(vr)
	taurr.nx -= 1
	taurr.ny -= 1
	taurr.x = vp.x[:-1]
	taurr.y = vr.y[:-1]
	taupp = copy.deepcopy(taurr)
	taurp = copy.deepcopy(taurr)
	print 'copy'
	T,R = meshgrid(taurr.x,taurr.y)
	print 'mesh'
	div_v = divV(vr,vp)
	print 'div'
	div = div_v.data*(2.0/3.0)
	dRvr, dTvr, dRrvr,_ = grad(vr)
	dRvp, dTvp, dRrvp,_ = grad(vp)
	print 'grad'
	rho_corner = 0.25*(rho.data[1:,1:]+rho.data[:-1,:-1]+
                       rho.data[1:,:-1]+rho.data[:-1,1:])

	print 'here'
	taurr.data = nu*rho.data[:-1,:-1]*(2*dRvr.data - div)	
	taupp.data = nu*rho.data[:-1,:-1]*((2*dTvp.data + vr.data[1:,1:]+vr.data[:-1,1:])/R - div)
	taurp.data = nu*rho_corner*(dRvp.data + (dTvr.data - .5*(vp.data[1:,1:]+vp.data[:-1,1:]))/R)
	tens.rr = copy.deepcopy(taurr)
	tens.pp = copy.deepcopy(taupp)
	tens.rp = copy.deepcopy(taurp)
	print 'div'
	divTensor(tens)
	return tens

def divTensor(tens):
	divr = copy.deepcopy(tens.rr)
	divr.nx -= 1
	divr.ny -= 1
	divr.x = tens.rr.x[1:]
	divr.y = tens.rr.y[1:]
	divp = copy.deepcopy(divr)

	taurr = tens.rr.data
	taurp = tens.rp.data
	taupp = tens.pp.data
	Tt,Rt = meshgrid(tens.rr.x,tens.rr.y)
	dR = Rt[1:,1:] - Rt[:-1,1:]
	dT = Tt[1:,1:] - Tt[1:,:-1]
	T,R = meshgrid(divr.x,divr.y)
	divr.data = ((Rt*taurr)[1:,1:] - (Rt*taurr)[:-1,1:])/(dR*R)
	divr.data +=  (taurp[1:,1:] - taurp[1:,:-1])/(dT*R)
	divr.data -= .5*(taurr[1:,1:] + taurr[:-1,1:])/R

	divp.data = (taupp[1:,1:] - taupp[1:,:-1])/(dT*R)
	divp.data += ((Rt*Rt*taurp)[1:,1:]-(Rt*Rt*taurp)[:-1,1:])/(dR*R*R)


	tens.divr = copy.deepcopy(divr)
	tens.divp = copy.deepcopy(divp)
	print 'okay'
	divp.data *= R
	divr.data *= R
	tens.rdivr = copy.deepcopy(divr)
	tens.rdivp = copy.deepcopy(divp)
	return

def pTensor(T,Tbar,Tpp):

	Tp = copy.deepcopy(T)

	Tp.rr.data  = T.rr.data - (Tbar.rr.data +Tpp.rr.data )
	Tp.pp.data  = T.pp.data - (Tbar.pp.data +Tpp.pp.data )
	Tp.rp.data  = T.rp.data  - (Tbar.rp.data  + Tpp.rp.data )
	Tp.divr.data  = T.divr.data  - (Tbar.divr.data  + Tpp.divr.data )
	Tp.divp.data  = T.divp.data  - (Tbar.divp.data  + Tpp.divp.data )
	

	return Tp
def datsplit(rho,vr,vp):
	
	dbar = copy.deepcopy(rho)
	vrbar = copy.deepcopy(vr)
	vpbar = copy.deepcopy(vp)

	sig = copy.deepcopy(rho)
	u = copy.deepcopy(vr)
	v = copy.deepcopy(vp)
	
	db = rho.data.mean(axis=1)
	du = vr.data.mean(axis=1)
	dv = vp.data.mean(axis=1)

	dbar.data = transpose(tile(db,[rho.nx,1]))
	vrbar.data = transpose(tile(du,[vr.nx,1]))
	vpbar.data = transpose(tile(dv,[vp.nx,1]))

	sig.data = rho.data-dbar.data
	u.data = vr.data - vrbar.data
	v.data = vp.data - vpbar.data

	return dbar,vrbar,vpbar,sig,u,v

def calcTwd(dat):
	om = dat.rho.omegaframe	
	c = dat.rho.cprof 
	dbar = copy.deepcopy(dat.dbar)
	vrbar = copy.deepcopy(dat.vrbar)
	vpbar = copy.deepcopy(dat.vpbar)
	sig = copy.deepcopy(dat.sig)
	u = copy.deepcopy(dat.u)
	v = copy.deepcopy(dat.v)
	print 'copy'
	drv,_,_,_ = grad(v)
	drvpbar,_,_,_ = grad(vpbar)
	drdbar,dpdbar,_,invrdpdbar = grad(dbar)
	drsig,dpsig,_ ,invrdpsig= grad(sig)	
	print 'grad'
	sigdbar = sig.data/dbar.data
	sigdbar2 = sigdbar**2
	print 'sigdbar'
	print dat.Tbar.rdivp.data.shape
	print c.shape
	print dpdbar.data.shape
	rdivpbar = dat.Tbar.rdivp.data - c[1:-1,1:-1]**2 * dpdbar.data[:-1,:-1]
	print 'rdiv1'
	rdivpp = dat.Tp.rdivp.data - c[1:-1,1:-1]**2 * dpsig.data[:-1,:-1]
	print 'rdiv'
	dbar_corner = 0.25*(dbar.data[1:,1:]+dbar.data[:-1,:-1]+
                       dbar.data[1:,:-1]+dbar.data[:-1,1:])
	sig_corner = 0.25*(sig.data[1:,1:]+sig.data[:-1,:-1]+
                       sig.data[1:,:-1]+sig.data[:-1,1:])
	print 'corner'
	T1,R1 = meshgrid(drv.x,drv.y)
	print 'mesh'

	Twd =  R1*dbar_corner*.5*(u.data[1:,1:]+u.data[1:,:-1])*drv.data
	print '1'
	Twd += R1*sig_corner*.5*(vrbar.data[1:,1:]+vrbar.data[1:,:-1])*(drvpbar.data+2*om)
	print '2'
	Twd = Twd[:-1,:-1]
	print 'shape'
	Twd += .5*(sigdbar[1:-1,1:-1]+sigdbar[1:-1,2:])*rdivpp
	print '3'
	Twd -= .5*(sigdbar2[1:-1,1:-1]+sigdbar2[1:-1,2:])*rdivpbar
	print '4'
	return Twd.mean(axis=1)

def calcFp(dat):
	
	FFp = copy.deepcopy(dat.rho)
	T,R = meshgrid(dat.rho.x,dat.rho.y)
	FFp.data = R*dat.vrbar.data*dat.sig.data*dat.v.data
	FFp.data += R*dat.dbar.data*dat.u.data*dat.v.data

	drFp,_,_,_ = grad(FFp)
	print 'grad'
	drFp.data += (dat.sig.data*(dat.v.data*dat.vrbar.data+dat.u.data*dat.vpbar.data))[:-1,:-1]
	print 'add'

	
	return drFp.data.mean(axis=1)

def calcTh(dat):

	T,R = meshgrid(dat.rho.x,dat.rho.y)
	mp = dat.rho.mp
	rs = dat.rho.rsoft
# 	Phi = -mp/sqrt(rr*rr + 1 + rs*rs-2*rr*cos(dat.pp))
	dpPhi = mp*R*sin(T)*pow(R*R+1+rs*rs-2*R*cos(T),-1.5)
#	drPhi = mp*(rr-cos(dat.pp))*pow(rr*rr+1+rs*rs-2*rr*cos(dat.pp),-1.5)

	Th = dat.sig.data*dpPhi

	return Th.mean(axis=1)

def make_dens_mov(dir,numf,jstart,backslope,clims=(None,None)):
	if len(dir) != 0:
		if dir[-1] != '/' and dir!='':
			dir += '/'
	outdir = dir+'mov/'
	
	totchar = len(str(numf))
	dat=Field(dir+'gasdens0.dat')
	nx = dat.nx
	ny = dat.ny
	yy,xx = meshgrid(dat.x,dat.y)
	x = xx*cos(yy)
	y = xx*sin(yy)
	rhob = xx**(-backslope)

	
	
	for j in range(jstart,numf):
		rho = fromfile(dir+'gasdens'+str(j)+'.dat').reshape(ny,nx) - rhob

		fig=figure()
		ax = subplot(111)
		if clims == (None,None):
		
			ax.pcolormesh(x,y,rho)
		else:
		
			ax.pcolormesh(x,y,rho,vmin=clims[0],vmax=clims[1])
		ax.axis('tight')
	
		fig.canvas.draw()
		
		zerostr = ''
		for i in range(totchar-len(str(j))):
			zerostr += '0'
		

		
		fname = outdir+'density' + zerostr + str(j) + '.png'
	
		savefig(fname,bbox_inches='tight')
		close(fig)
		print 'Finished ' + str(j)

	return

