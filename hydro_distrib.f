c distributed version for hydraulic diffusion part on 2d fault
c based on Aochi et al. (2014) and the following papers
c kappa_phi2 function : non-linear kappa function 
c       here depending on effective normal stress, should in m2
c injection function :  an example, should be in m3/s in main program
c
c For all the variables, please refer the equations in Aochi et al. (2014)

	implicit none
	integer imax, jmax, itmax
	integer i0, j0
        real ds, dt
c model dimension (1:20, 1:20)
	parameter (imax=20, jmax=20)
c total time step  
c	parameter (itmax=3600*24*4)
	parameter (itmax=3600)
c injection point 
	parameter (i0 = imax/2, j0 = jmax/2)
c element size (m) and time step (s)
        parameter (ds = 30., dt = 1.0)
c for above setting, physical domain (600 m, 600 m) for 4 days

	real pr(0:imax+1, 0:jmax+1), prold(0:imax+1, 0:jmax+1)
	real prmax(0:imax+1, 0:jmax+1)
        real tau_n(0:imax+1, 0:jmax+1)
	real alpha, fac, dpdtmax, dpdt, p0

	real phi(0:imax+1, 0:jmax+1), phiold(0:imax+1, 0:jmax+1)
	real phico(0:imax+1, 0:jmax+1), dphipl(0:imax+1, 0:jmax+1)
	real kappa(0:imax+1, 0:jmax+1)
	real h(0:imax+1, 0:jmax+1), hold(0:imax+1, 0:jmax+1)
	real dgamma(0:imax+1, 0:jmax+1)

	character       num*10, dir*60
	character*50    output, output2, output7

	integer i, j, k, l, m, ndir
	integer ixhypo, iyhypo, i1, j1, i2, j2, i3, j3, icos

	real    beta, dh, dh0, h0
	real	kappa0, kappamax, kappa_phi2, dphi_pl, kappamin, prevk
	real	tau0, pr0, phi0, phimax, phimin, phiinf, wc
	real	betap, betaf, rho, eta 
	real	fs, sigma_n, v0, eff, dtau, tres
	real	dkdp1, dkdp2, dp, time, ratec, vtot(itmax)
	real	diff, rheo, inject, injection
	external kappa_phi2
	external injection

c output directory to be created in advance
	dir = "hydro_distrib/"
	ndir = index(dir, ' ')-1
	output = dir(1:ndir)//"hoge1.dat"
	output2 = dir(1:ndir)//"hoge2.dat"
	output7 = dir(1:ndir)//"initial.dat"

c parameters
	betaf = 5.*10.**(-10)
	betap = betaf/10.
	beta = betap + betaf

	dh0 = 0.
	h0 = 5.0

	phi0 = 0.05
	phimin = 0.
	phimax = 0.15
	phiinf = 0.1
	dphi_pl = 0.0
	wc = 0.01

	kappamax = 10.**(-12)
	kappamin = 10.**(-16)
        kappa0 = 10.**(-8)

	rho = 1.*10.**3
	eta = 2.*10.**(-4)
c injection is given by external function
	eff = 1.0
c alpha = prmax
	alpha = 100.*10.**6
	fac = 1.0
	pr0 = 30.*10.**6

        sigma_n = 100.*10.**6
        fs = 0.65

c writing parameters
	open(11, file=output)
	write(11,*) dir
	write(11,*) ds, dt
	write(11,*) betaf, betap
	write(11,*) dh0, h0
	write(11,*) phi0, phimin, phimax, phiinf, wc
	write(11,*) kappa0, kappamax
	write(11,*) rho, eta
	write(11,*) pr0, sigma_n, fs
	close(11)

!!! initial condition hydro
	open(27, file=output7)
	do i=0, imax+1
	  do j=0, jmax+1
	    pr(i,j) = 0.0
	    prmax(i,j) = alpha - pr0

	    phi(i,j) = phi0
	    phico(i,j) = 0.0
	    dphipl(i,j) = 0.0

	    kappa(i,j) = kappamax

	    h(i,j) = h0
	    hold(i,j) = h0

            dgamma(i,j) = 0.

	    icos = 0
            tau_n(i,j) = sigma_n/10.**6 - (pr0 + pr(i,j))/10.**6
	    kappa(i,j) = kappa_phi2(icos, tau_n(i,j), 
     &		kappa0, kappamax, kappamin)

	    if (i.ge.1.and.i.le.imax) then
	      if (j.ge.1.and.j.le.jmax) then
	        write(27, '(2i5, e12.4, 4f12.4)') i, j, 
     &	        kappa(i,j), phi(i,j), pr(i,j)/10.**6, h(i,j), dgamma(i,j)
	      endif
	    endif


	  enddo
	enddo
	close(27)
c
c iteration
c
        open(22, file=output2)
        do k=1, itmax
          time = dt*k
          vtot(k) = 0.
          if (k.ge.2) vtot(k) = vtot(k-1)

c updating 1: pr and phiold
	  do i=1, imax
	    do j=1, jmax
	      prold(i,j) = pr(i,j)
	      phiold(i,j) = phi(i,j)    
	    enddo
	  enddo
c evolution with time (initially planned)
          ratec = injection(time/3600./24.)

	  do i=1, imax
	    do j=1, jmax
	      dgamma(i,j) = 0.0
	      if (i.eq.i0.and.j.eq.j0) then 
                dgamma(i,j) = ratec
    	      endif
	      vtot(k) = vtot(k) + dgamma(i,j)*dt
c updating 2
	      dh = (h(i,j) - hold(i,j))/dt
	      hold(i,j) = h(i,j)

c fluid pressure	
	      dkdp1 = (kappa(i+1,j)-kappa(i-1,j))*
     &                  (prold(i+1,j)-prold(i-1,j))
	      dkdp2 = (kappa(i,j+1)-kappa(i,j-1))*
     &                  (prold(i,j+1)-prold(i,j-1))
	      dp = prold(i+1,j) + prold(i,j+1) + prold(i-1,j) 
     &		+ prold(i,j-1) - 4*prold(i,j)

c bug fix (31/05/2011) ds**2
	      diff = 
     &		1./eta*(dkdp1 + dkdp2)/(4.*ds**2) + kappa(i,j)/eta*dp/ds**2 

     	      rheo = phiold(i,j)*dh/hold(i,j) + dphi_pl 
	      inject = dgamma(i,j)/(ds**2*hold(i,j)) 

	      dpdt = dt/(phiold(i,j)*beta)*(diff - rheo + inject)

	      dpdtmax = fac*(prmax(i,j) - prold(i,j))

	      if( dpdt.gt.dpdtmax ) then
	        dpdt = dpdtmax
      		dphipl(i,j) = (1-phi(i,j))*(diff + inject - dpdt*phi(i,j)*beta)
	      else
	        dphipl(i,j) = 0.0
	      endif
	      pr(i,j) = prold(i,j) + dpdt

	    enddo
	  enddo

c porosity and fault strength, fault width
	 
	  do i=1, imax
	    do j=1, jmax

c if already ruptured, Pmax imposed

	        phi(i,j) = phi(i,j) + phi(i,j)*betap*(pr(i,j)-prold(i,j))
     &		+ dphipl(i,j)

	      if( phi(i,j).gt.phimax) phi(i,j) = phimax
	      if( phi(i,j).lt.phimin) phi(i,j) = phimin

	      h(i,j) = h(i,j)*(1.0 + (phi(i,j)-phiold(i,j))/(1.-phi(i,j)) 
     &		- betap*(pr(i,j) - prold(i,j) ))

	      icos = 0
	      tau_n(i,j) = sigma_n/10.**6 - (pr0 + pr(i,j))/10.**6
	      kappa(i,j) = kappa_phi2(icos, tau_n(i,j), 
     &		kappa0, kappamax, kappamin)

	    enddo
	  enddo

c outer-boundary condition for spatial derivatives (kappa, pr)
	  do j=1, jmax
	    pr(0,j) = pr(1,j)
	    pr(imax+1,j) = pr(imax,j)
	  enddo
	  do i=1, imax
	    pr(i,0) = pr(i,1)
	    pr(i,jmax+1) = pr(i,jmax)
	  enddo

	  if(mod(int(time),3600).eq.0.or.k.eq.itmax) then
          write(*,*) "output", k, time,  time/3600., time/3600./24., 
     &                  ratec
          write(22,'(i10,f12.0,f12.4,f12.5,2i5,f12.4,e12.3)') 
     &          k, time, time/3600, ratec, 
     &          i0, j0, pr(i0,j0)/10.**6, kappa(i0, j0)
	  write(num, '(i4.4)') int(time)/3600
	  output = dir(1:ndir)//"snapshot"//num(1:4)//".dat"
	  open(21, file=output)
	  do i=1, imax
	    do j=1, jmax
	      write(21,'(2i5,f10.3, f10.4, e15.4, f10.4)') i, j, 
     &		pr(i,j)/10.**6, phi(i,j), kappa(i,j), h(i,j)
	    enddo
	  enddo
	  close(21)
	  endif
	
	 enddo
         close(22)

	 stop
	 END
c
	real function kappa_phi2(icos, eff, kappa0, kappamax, kappamin)
c icos :  not used
c eff : effective normal stress
c
	implicit none
	integer icos
	real eff, s0
	real kappa0, kappamax, kappamin

	s0 = 4.

	kappa_phi2 = kappa0*exp(-eff/s0)

	if (kappa_phi2.gt.kappamax) kappa_phi2 = kappamax
	if (kappa_phi2.lt.kappamin) kappa_phi2 = kappamin

	return
	END

c
c INJECTION: BALE(2016)-like 
c time below is in day, injection should be m3/s in main program
c
	REAL function injection(time)

	real time

	injection = 0.
	if (time.le.0.11) then
	  injection = 0.4*time/0.11
	else if (time.le.0.32) then
	  injection = 0.4 + (2.0 - 0.4)*(time-0.11)/(0.32-0.11)
	else if (time.le.0.74) then
	  injection = 2.0 
	else if (time.le.0.75) then
	  injection = 2.0 + (5.99-2.0)*(time-0.74)/(0.75-0.74)
	else if (time.le.0.9) then
	  injection = 5.99
	else if (time.lt.1.76) then
	  injection = 7.98
	else if (time.lt.2.15) then
	  injection = 7.98
	else if (time.lt.2.17) then
	  injection = 7.98 + (15.17-7.98)*(time-2.15)/(2.17-2.15)
	else if (time.lt.2.76) then
	  injection = 15.17
	else if (time.lt.2.78) then
	  injection = 15.17 + (19.96-15.17)*(time-2.76)/(2.78-2.76)
	else if (time.lt.2.85) then
	  injection = 19.96
	else if (time.lt.2.93) then
	  injection = 25.15
	else if (time.lt.2.94) then
	  injection = 25.15 + (30.74-25.15)*(time-2.93)/(2.94-2.93)
	else if (time.lt.3.84) then
	  injection = 30.74 + (30.34-30.74)*(time-2.94)/(3.84-2.94)
	else if (time.lt.3.86) then
	  injection = 30.34 + (53.10-30.34)*(time-3.84)/(3.86-3.84)
	else if (time.lt.3.88) then
	  injection = 53.10 + (0. - 53.10)*(time-3.86)/(3.88-3.86)
	else if (time.lt.3.93) then
      	  injection = 0. + 50.30*(time-3.88)/(3.93-3.88)
    	else if (time.lt.3.99) then
     	  injection = 50.30 + (47.51-50.30)*(time-3.93)/(3.99-3.93)
    	else if (time.lt.4.0) then
      	  injection = 47.51 + (41.52-47.51)*(time-3.99)/(4.0-3.99)
    	else if (time.lt.4.76) then
      	  injection = 41.52
    	else if (time.lt.4.81) then
      	  injection = 63.08 + (57.49-63.08)*(time-4.76)/(4.81-4.76)
    	else if (time.lt.4.88) then
      	  injection = 57.49 + (55.49-57.49)*(time-4.81)/(4.88-4.81)
    	else if (time.lt.5.12) then
      	  injection = 55.49
    	else if (time.lt.5.15) then
      	  injection = 55.49 + (52.70-55.49)*(time-5.12)/(5.15-5.12)
    	else if (time.lt.5.42) then
      	  injection = 52.7 + (53.9-52.7)*(time-5.15)/(5.42-5.15)
    	else if (time.lt.5.73) then
      	  injection = 29.94 + (30.74-29.94)*(time-5.42)/(5.73-5.42)
    	else if (time.lt.5.76) then
      	  injection = 30.74 + (0.-30.74)*(time-5.73)/(5.76-5.73)
    	else 
      	  injection = 0.
    	endif
    
    	injection = injection*0.001
    	return
    	end



