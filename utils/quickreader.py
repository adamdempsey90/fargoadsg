from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.animation as animation
from numpy import *
from matplotlib.pyplot import *
#from pylab import *

class fargo():
	def __init__(self,num,iso=True):
		timedat = loadtxt('planet0.dat')
		self.time = timedat[num,7]
		self.r = loadtxt('used_rad.dat')
		self.r = .5*(self.r[:-1]+self.r[1:])
		self.nrad = len(self.r)
		self.nphi = int(loadtxt('dims.dat')[-1])

		self.phi = linspace(-pi,pi,self.nphi)
		
		rr,pp = meshgrid(self.r,self.phi)
		
		self.rr = rr.transpose()
		self.pp = pp.transpose()
		self.x = self.rr*cos(self.pp)
		self.y = self.rr*sin(self.pp)
		
		self.dlr = diff(log(self.r))[0]
		self.dphi = diff(self.phi)[0]

		dens_fname = 'gasdens%d.dat' % num
		vr_fname = 'gasvrad%d.dat' % num
		vp_fname = 'gasvtheta%d.dat' % num
		temp_fname = 'gasTemperature%d.dat' % num



#		dat_rho =[fromfile(dens_fname)]
#		dat_vr = [fromfile(vr_fname)]
#		dat_vp = [fromfile(vp_fname)]
#		dat_temp = [ fromfile(temp_fname)]
#		for i in range(1,np):
#			if np<10:
#				suffix = '.0000%d' % i
#			elif np<100:
#				suffix = '.000%d' % i
#			elif np<1000:
#				suffix = '.00%d' % i
#			elif np<10000:
#				suffix = '.0%d' % i
#			else:
#				suffix = '.%d' % i		
#			dat_rho.append(fromfile(dens_fname + suffix))
#			dat_vr.append(fromfile(vr_fname + suffix))
#			dat_vp.append(fromfile(vp_fname + suffix))
#			dat_temp.append(fromfile(temp_fname + suffix))
#
#		self.rho = hstack(dat_rho).reshape(self.nrad,self.nphi)
#		self.vr = hstack(dat_vr).reshape(self.nrad,self.nphi)
#		self.vp = hstack(dat_vp).reshape(self.nrad,self.nphi)
#		self.temp = hstack(dat_temp).reshape(self.nrad,self.nphi)
		self.rho = fromfile(dens_fname).reshape(self.nrad,self.nphi)
		self.vr = fromfile(vr_fname).reshape(self.nrad,self.nphi)
		self.vp = fromfile(vp_fname).reshape(self.nrad,self.nphi)
		
		if iso:
			self.iso = True
			self.temp = ones((self.nrad,self.nphi))
		else:
			self.iso = False
			self.temp = fromfile(temp_fname).reshape(self.nrad,self.nphi)

		self.vort,self.vortens = self.vortencity()

		self.mrho = fft.rfft(self.rho)/self.nphi
		self.mvr = fft.rfft(self.vr)/self.nphi
		self.mvp = fft.rfft(self.vp)/self.nphi
		self.mtemp = fft.rfft(self.temp)/self.nphi
		self.mvort = fft.rfft(self.vort)/self.nphi
		self.mvortens = fft.rfft(self.vortens)/self.nphi
		
		self.dbar = tile((self.mrho[:,0].real).reshape(self.nrad,1),(1,self.nphi))
		self.vrbar = tile((self.mvr[:,0].real).reshape(self.nrad,1),(1,self.nphi))
		self.vpbar = tile((self.mvp[:,0].real).reshape(self.nrad,1),(1,self.nphi))
		self.Tbar = tile((self.mtemp[:,0].real).reshape(self.nrad,1),(1,self.nphi))
		self.vortbar = tile((self.mvort[:,0].real).reshape(self.nrad,1),(1,self.nphi))
		self.vortensbar = tile((self.mvortens[:,0].real).reshape(self.nrad,1),(1,self.nphi))
		
		self.power = zeros(self.mrho.shape[1])
		for i in range(self.mrho.shape[1]):
			self.power[i] = (conj(self.mrho[:,i])*self.mrho[:,i]).sum()
			

	def vortencity(self):
 		_,dpvr=gradient(self.vr,self.dlr,self.dphi)
		dlrvp,_ = gradient(self.vp,self.dlr,self.dphi)
		vort = (self.vp + dlrvp - dpvr)/self.rr
		vortens = vort/self.rho
 		return vort,vortens	
	def plot(self,q,rmax=None,xlims=None,ylims=None,logr=False,cmap='autumn',scale=1,output=False,fname=None):
		if q=='all':
			fig,((axs,axt),(axvr,axvp))=subplots(2,2,sharex='col',sharey='row',figsize=(9,7),dpi=100)
			axs.set_title('$\log_{10} \Sigma$')
			if self.iso:
				axt.set_title('$ \\log_{10}\\left|\\frac{\\nabla \\times v}{\\Sigma}\\right| $')
			else:
				axt.set_title('$ T $')
			axvr.set_title('$v_r$')
			axvp.set_title('$v_\phi$')

			ps=axs.pcolormesh(self.x,self.y,log10(self.rho))
			if self.iso:
				pt=axt.pcolormesh(self.x,self.y,log10(abs(self.vortens)))
			else:
				pt=axt.pcolormesh(self.x,self.y,self.temp)
			pr=axvr.pcolormesh(self.x,self.y,self.vr)
			pp=axvp.pcolormesh(self.x,self.y,self.vp)
			
			dividers = make_axes_locatable(axs)
			caxs = dividers.append_axes("right", size="5%", pad=0.05)
			cbars = colorbar(ps, cax=caxs)
#			axs.xaxis.set_visible(False)

			dividert = make_axes_locatable(axt)
			caxt = dividert.append_axes("right", size="5%", pad=0.05)
			cbart = colorbar(pt, cax=caxt)

			dividervr = make_axes_locatable(axvr)
			caxvr = dividervr.append_axes("right", size="5%", pad=0.05)
			cbarvr = colorbar(pr, cax=caxvr)

			dividervp = make_axes_locatable(axvp)
			caxvp = dividervp.append_axes("right", size="5%", pad=0.05)
			cbarvp = colorbar(pp, cax=caxvp)
			
			if rmax != None:
				axs.set_xlim((-rmax,rmax))
				axt.set_xlim((-rmax,rmax))
				axvr.set_xlim((-rmax,rmax))
				axvp.set_xlim((-rmax,rmax))
				axs.set_ylim((-rmax,rmax))
				axt.set_ylim((-rmax,rmax))
				axvr.set_ylim((-rmax,rmax))
				axvp.set_ylim((-rmax,rmax))

			else:
				if xlims != None:
					axs.set_xlim(xlims)
					axt.set_xlim(xlims)
					axvr.set_xlim(xlims)
					axvp.set_xlim(xlims)
				if ylims != None:
					axs.set_ylim(ylims)
					axt.set_ylim(ylims)
					axvr.set_ylim(ylims)
					axvp.set_ylim(ylims)
		if q=='rho':
			figure()
			pcolormesh(self.x,self.y,log10(self.rho))
			colorbar()
			title('$\log_{10} \Sigma$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='vr':
			figure()
			pcolormesh(self.x,self.y,self.vr)
			colorbar()
			title('$v_r$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='vp':
			figure()
			pcolormesh(self.x,self.y,self.vp)
			colorbar()
			title('$v_\phi$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='vort':
			figure()
			pcolormesh(self.x,self.y,self.vort)
			colorbar()
			title('$\\nabla \\times v$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='vortens':
			figure()
			pcolormesh(self.x,self.y,self.vortens)
			colorbar()
			title('\\frac{$\\nabla \\times v}{\\Sigma}$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='omega':
			figure()
			pcolormesh(self.x,self.y,self.omega)
			colorbar()
			title('$\Omega$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='temp':
			figure()
			pcolormesh(self.x,self.y,self.temp)
			colorbar()
			title('$T$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='mean':
			if logr==True:
				r = log10(self.r)
				xstr = '$ \log_{10} r $'
			else:
				r = self.r
				xstr = '$r$'

			fig,((axs,axt),(axvr,axvp))=subplots(2,2,sharex='col')
			axs.set_ylabel('$\\bar{\Sigma}$')
			axvr.set_ylabel('$\\bar{v}_r$')
			axvp.set_ylabel('$\\bar{v}_\phi$')
			axt.set_ylabel('$\\bar{T}$')
			axvr.set_xlabel(xstr)
			axvp.set_xlabel(xstr)
			axs.plot(r,real(self.mrho[:,0]))
			axvr.plot(r,real(self.mvr[:,0]))
			axvp.plot(r,real(self.mvp[:,0]))
			axt.plot(r,real(self.mtemp[:,0]))
		if q=='vrp':
			figure()
			pcolormesh(self.x,self.y,scale*(self.vr-self.vrbar),cmap=cmap)
			colorbar()
			title('$v_r-\\bar{v}_r$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='vpp':
			figure()
			pcolormesh(self.x,self.y,self.vp-self.vpbar)
			colorbar()
			title('$v_\\phi-\\bar{v}_\\phi$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='rhop':
			figure()
			pcolormesh(self.x,self.y,self.rho-self.dbar)
			colorbar()
			title('$\\Sigma-\\bar{\\Sigma}$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))
		if q=='Tp':
			figure()
			pcolormesh(self.x,self.y,self.temp-self.Tbar)
			colorbar()
			title('$T-\\bar{T}$')
			if rmax != None:
				xlim((-rmax,rmax))
				ylim((-rmax,rmax))	
		if q=='power':
			figure()
			semilogy(self.power[:10]/self.power[0],'-s')
			ylabel('$\\left| \\frac{\\Sigma_m}{\\bar{\\Sigma}} \\right|^2$',fontsize=20)
			xlabel('m',fontsize=20)
		if output:
			if q=='all':
				fig.savefig(fname)
			else:
				savefig(fname)
		
	def plotmode(self,m):
		fig,((axs,axt),(axvr,axvp))=subplots(2,2,sharex='col')
	
		axs.set_ylabel('$\Sigma_%d$'%m)
		axvr.set_ylabel('$v_{r,%d}$'%m)
		axvp.set_ylabel('$v_{\phi,%d}$'%m)
		if self.iso:
			axt.set_ylabel('$(\\nabla \\times v)_%d$'%m)
		else:
			axt.set_ylabel('$T_%d$'%m)
		axvr.set_xlabel('$r$')
		axvp.set_xlabel('$r$')

		axs.plot(self.r,real(self.mrho[:,m]),'-k',label=r'Real')
		axvr.plot(self.r,real(self.mvr[:,m]),'-k')
		axvp.plot(self.r,real(self.mvp[:,m]),'-k')
		if self.iso:
			axt.plot(self.r,real(self.mvort[:,m]),'-k')
		else:
			axt.plot(self.r,real(self.mtemp[:,m]),'-k')
	
		axs.plot(self.r,imag(self.mrho[:,m]),'--k',label=r'Imag')
		axvr.plot(self.r,imag(self.mvr[:,m]),'--k')
		axvp.plot(self.r,imag(self.mvp[:,m]),'--k')
		if self.iso:
			axt.plot(self.r,imag(self.mvort[:,m]),'--k')
		else:
			axt.plot(self.r,imag(self.mtemp[:,m]),'--k')
	
		axs.legend(loc='best')

		
def animate(trange,q,logy=False,rmax=None,iso=False):
	nt = len(trange)
	fld=fargo(trange[0],iso=iso)
	x = fld.x
	y = fld.y
	r = fld.r
	nr = fld.nrad
	nphi = fld.nphi
	dat = [None]*nt

	if q=='rho':
		base_name = 'gasdens'
	elif q=='vr':
		base_name = 'gasvrad'
	elif q=='vp':
		base_name = 'gasvtheta'
	elif q=='temp':
		if iso:
			print 'No temperature for Isothermal EOS'
			return
		base_name = 'gasTemperature'
	else:
		print 'Cannnot animate ' + q
		return

	
	if rmax != None:
		x = x[r <= rmax,:]
		y = y[r <= rmax,:]
	
	vmin_val = 0; vmax_val = 0;

	for i,t in enumerate(trange):
		
		fname = base_name + '%d.dat' % t
		print 'Loading %d, ' % t + fname 
		if logy:
			dat[i]=log10(fromfile(fname).reshape(nr,nphi))
		else:
			dat[i] = fromfile(fname).reshape(nr,nphi)
		if rmax != None:
			dat[i] = dat[i][r <= rmax,:]
		vmin_val = min(vmin_val,dat[i].min())
		vmax_val = max(vmax_val,dat[i].max())
#	fig=figure()
#	pc = pcolormesh(x,y,dat[0])
#	animate = lambda j: pc.set_array(dat[j].ravel())
#	animation.FuncAnimation(fig, animate, frames = range(nt), blit = False)
#	show()
#
	fig=figure()
	ax=subplot(111)	
	z = pcolormesh(x,y,dat[0])

	divider = make_axes_locatable(ax)
	caxp = divider.append_axes("right", size="5%", pad=0.05)
	cbar = colorbar(z, cax=caxp)

	show()
	for i,t in enumerate(trange[1:]):
		z.set_array(dat[i].ravel())
		z.set_clim(vmin=dat[i].min(),vmax=dat[i].max())

#		cbar=colorbar(z,cax=caxp)
		ax.set_title('%d' % t)
		draw()

	show()
	return

def modesplot(rout,maxmode,plotzero=False,showlegend=True):
	dat=loadtxt('modehistory.dat')
	pout = 2*pi*rout**(1.5)
	t = dat[:,0]/pout
	dat= dat[:,1:]
	figure();
	if plotzero:
		loglog(t,dat[:,0],label='m=0')

	for m in range(1,maxmode+1):
		loglog(t,dat[:,m],label='m=%d'%m)
	
	xlabel('t (outer period)',fontsize=20)
	ylabel('$|A_m|$',fontsize=20)
	if showlegend:
		legend(loc='best')

