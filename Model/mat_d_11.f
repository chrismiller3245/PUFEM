      subroutine modd_11(l,nhvd,nb,u0,du0,imatd,d,T,ddu,ddn,hflag)
c
c---------------------------------------------------------TCG.28.08.2002
c
c     Discrete Finite Deformation Material Models for FEAP
c          
c cohesive zone material models (traction seperation laws)
c  imatd=1    isotropic traction Model  
c  imatd=2    transversely isotropic traction model
c  imatd=3    isotropic viscoelastic traction model
c  imatd=4    smooth traction model
c  imatd=6    Exponential normalised anisotropic Traction model            
c
c.... INPUT variables
c
c         l          no. of Gausian point
c         nc         normal vect. in cur. conf.
c         u0         Gap displacement
c         du0        Increment of Gap displacement
c         imatd      Discrete material model
c               (1) Linear isotropic Traction
c               (2) Linear transversally Traction
c         d(100)     Material parameters
c
c.... OUTPUT variables
c
c         T(3)       Traction at discontinuity 
c         ddu(3,3)      Traction moduli with respect to u0
c         ddn(3,3)     Traction moduli with respect to n
c         hn1(nhv)   History variables at Gauss points, time t_n+1
c
c--------------------------------------------------------------------71
c
c
c..... Declare variable types
      implicit none
	logical hflag
      integer imatd,l, nhvd, nn, i
c..... Declare array types
      real*8  d(*)
	real*8  t(3), ddu(3,3), ddn(3,3), u0(3), du0(3), nb(3)
      real*8  check(2), hn1(nhvd)
c
      include  'tdata.h'
      include  'hdata.h'
      include  'comblk.h'
c
c	write(*,*) hn1(1)

c.... Move history t_n-related history variables from the hr-array to hn1-array (read at the nh1 position from FEAP's hr vector)
      
      nn = nh1 + int(d(78)) + (l-1)*nhvd
      do i = 1,nhvd
        hn1(i) = hr(nn+i-1)
        check(i) = hr(nn+i-1)
      enddo

c.... Compute traction vector and traction moduli
c
      if (imatd .eq. 1) then
c isotropic traction Model  
        call Trac11_iso_pen_sept05(l,d,u0,nb,T,ddu, ddn)
	elseif (imatd.eq.2) then 
c transversely isotropic traction model
       call Trac11_transiso(l,d, u0,nb, T,ddu, ddn)
      elseif (imatd .eq. 3) then
c isotropic viscoelastic traction model
        call Trac11_viso_pen_sept05(l,d,u0,nb,T,ddu, ddn)
      elseif (imatd .eq. 4) then
c smooth traction model
	  call Trac11_iso_smooth(l,d,u0,nb,T,ddu, ddn)
      elseif (imatd .eq. 6) then 
c Linear Elastic Traction model (CM 2022) -- Crack initialisation criteria => Maximum principal stress 
        call Trac11_lin_elas(l,hn1,nhvd,d,u0,nb,T,ddu,ddn)  
      elseif (imatd .eq. 7) then 
c Exponential isotropic Traction model (CM 2022) -- Crack initialisation criteria => Maximum principal stress 
        call Trac11_iso_exp(l,hn1,nhvd,d,u0,nb,T,ddu,ddn) 
      elseif (imatd .eq. 8) then  
c Exponential Transversally isotropic Traction model (CM 2022) -- Crack initialisation criteria => Maximum normal stress in given fiber direction    
        call Trac11_norm_tran(l,hn1,nhvd,d,u0,nb,T,ddu,ddn) 
      else
	  write (*,*)
        write(*,*)'==========================================='
        write(*,*)'ERROR --> NO DISCRETE MATERIAL MODEL CALLED'
        write(*,*)'==========================================='
	  write (*,*)
      endif
      
c update history (write at the nh2 position into FEAP's hr vector)

      nn = nh2 + int(d(78)) + (l-1)*nhvd
      do i = 1,nhvd
        hr(nn+ i-1) = hn1(i) 
      enddo       
      

9999  end



       subroutine Trac11_iso_pen_sept05(l,d,u0,nc,T,ddu, ddn)
c==========================================================TCG.-28.09.2005
c  
c  Isotropic Traction model with isotropic damage
c
c    t = c u_0	; 
c   ddu = c II
c   ddn = 0
c
c   INPUT:
c         nc..........Normal vector at fictitious discontinuity
c         u0..........Gap displacement
c         uH..........history variable (max Gap)
c		d...........Material parameter vector
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n 
c                   
c
c   USED:
c         c, a,b..material parameters
c
c
c
c   History vector structure at the Gauss point level 
c	
c       entry 1 - maximal experienced upening uH = max|u0|
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'
	include  'tdata.h'

	logical dflag
	integer i,j, l, k, nhvd, lbh, n1, n2, n3, lint
      integer rot
      real*8 d(*), T(3), ddu(3,3),ddn(3,3), u0(3)
	real*8 a,b,c, te0, uh,tp(3),nc(3), u0p(3)
	real*8 un, ue, ue0, eps, te, gamma, pen,un2,twoun
	real*8 c0, cmin, u0n(3), xi, matrix(3,3), lamb(3)
      real*8 fact
      real*8 stf
 
      !stf = 1e8
      !
      !T  = K*u0
      !
      !do i= 1,3
      !  ddu(i,i) = stf 
      !enddo	  
      !
      !ddn=0.0d0      
      
      

	lbh= int(d(78))     ! lenght of bulk-related history     
      nhvd= int(d(75))    ! length of history per Gausspoint of discontinuity
     
	n1 = (l-1)*nhvd + nh1 + lbh
	n2 = (l-1)*nhvd + nh2 + lbh
	n3 = nh3 + 355   ! pointer to t0   
      
c material parameters 
c      te0 = hr(n3)
c      a   = te0/2.0d0
      te0 = d(20)
      a   = d(21)
      b   = d(22)
	ue0 = d(19)
	pen = d(30)

	c0=te0/ue0
	cmin=1.0d-5*c0
      
      
      gamma=0.0

c get ue from history or from inital compliance
        !ue =  0.0d0
        ue =  hr(n1)

C set initial stiffness
	  if (ue.le.ue0) ue=ue0  

c	  te=te0*exp(-a*(ue)**b) 

        uh = sqrt(u0(1)*u0(1)+u0(2)*u0(2)+u0(3)*u0(3))
c decide if elastic or damage
        if (uh.lt.ue) then    !elastic
	   if (uh.lt.ue0) then
	     c=c0
	   else
           c = te0*exp(-a*(ue)**b)/ue
	   endif
	   t = c*u0
c elastic tangent
          ddu = 0.0d0
          do i= 1,3
	      ddu(i,i) = c  
	    enddo	  
	  else  !damage
c traction        
c	    c = te*exp(-a*(uh**b-ue**b))/uh
          c = te0*exp(-a*(uh)**b)/uh
	    t = c*u0
          ddu=0.0d0
c elastic part of tangent
          do i= 1,3
	      ddu(i,i) = c
          enddo
          
c damage part of tangent (softening)
          If(d(24).eq.1) then
            gamma=c/uh*(a*b*uh**(b-1.0d0)+1.0d0/uh)
	      do i= 1,3
	        do j= 1,3
	         ddu(i,j) = ddu(i,j) - gamma*u0(i)*u0(j)
	        enddo
            enddo
          endif

          
c update ue
          hr(n2) = uh
      endif
            

c recompute if stiffness is too low
	  if(c.lt.cmin) then
	    c = cmin
	    t = c*u0
c elastic tangent
          ddu = 0.0d0
          do i= 1,3
	      ddu(i,i) = c  
	    enddo	  
        endif

c set Traction moduli	with respect to n to zero
         ddn=0.0d0

c check penetration
      un= u0(1)*nc(1) + u0(2)*nc(2) + u0(3)*nc(3)
	if (un.gt.0.0d0) then  ! add penalty
	  twoun = un + un
	  un2 = un*un
c traction
	  t = t + pen*un2*nc
c stiffness
        do i= 1,3
	   do j= 1,3
	     ddu(i,j) = ddu(i,j) + twoun*pen*nc(i)*nc(j)
	   enddo
	  enddo
        do i= 1,3
	    ddn(i,i) = ddn(i,i) + pen*un2
	  enddo
        do i= 1,3
	   do j= 1,3
	     ddn(i,j) = ddn(i,j) + twoun*pen*nc(i)*u0(j)
	   enddo
	  enddo
      endif

	end

       subroutine Trac11_transiso(l,d,u0,nb,T,ddu, ddn)
c==========================================================TCG.-10.04.2003
c  
c  transversely Isotropic Traction model with isotropic damage
c
c    t = c u_0	;   c = t0/delta exp(-a delta^b)
c   ddu = c (alpha II +(1-alpha)n x n)
c   ddn = c (ue.n II + n x ue)
c
c   INPUT:
c         l ..........no. of Gausian point
c         nb..........Normal vector at fictitious discontinuity
c         u0..........Gap displacement
c         uH..........history variable (max Gap)
c		d...........Material parameter vector
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n == 0
c                   
c
c   USED:
c         c, a,b..material parameters
c         lh     ..length of history vector at integration point
c         n1     ..pointer to history variable at previous time
c         n2     ..pointer to history variable at current time
c
c
c
c  d(72) Discrete Material number   (imatd)
c  d(74) History in bulk material   (nhv)
c  d(75) History in disc. material  (nhvd)
c
c
c   local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE


      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'

	logical dflag
	integer i,j, l, nhv, n1, n2, n3, lbh, nhvd
      real*8 d(*), T(3), ddu(3,3),ddn(3,3), u0(3)
	real*8 a,b,c, t0, uh
	real*8 norm, gamma, delta
	real*8 nb(3), alpha, un(3), ut(3), ua(3)
	real*8 pen, ca, coa, pencoa, cp, tcppen
	real*8 nbnb(3,3), nbu0(3,3), uh0, eps, dmd0

      data eps/1.0d-6/
      

      lbh= int(d(78))     ! lenght of bulk-related history     
      nhvd= int(d(75))    ! length of history per Gausspoint of discontinuity
     
	n1 = (l-1)*nhvd + nh1 + lbh
	n2 = (l-1)*nhvd + nh2 + lbh
	n3 = nh3 + 355   ! pointer to t0

      
c material parameters 

      t0    = d(20)

c      t0 = hr(n3)

c      write (iow,*) 't0 =', t0 
c      write (*,*) 't0 =', t0 
      
      a     = d(21)
      b     = d(22)
	alpha = d(23)
	uh0   = d(24)
      cp    = d(30)


c  update maximum gap	 if necessary and set damage evol. flag
      norm = sqrt(u0(1)*u0(1) + u0(2)*u0(2) + u0(3)*u0(3))

	 if (norm.gt.hr(n1)) then
	  dflag=.true.
	  hr(n2) = norm
	 else
	  dflag=.false.
	 endif

c  get history variable
	 uh = hr(n2)
c       dflag=.false.




c Phase I
      if (uh.lt.uh0+eps) then
c initial value
        if (uh.gt.eps) then
	    delta = uh
	  else
	    delta = eps
	  endif
      else
c Phase II
        delta = uh
	  dmd0  = uh-uh0    
	endif

c compute penetration
      pen = u0(1)*nb(1) + u0(2)*nb(2) + u0(3)*nb(3)  

c compute normal and transvers gap displacement
      un = pen*nb
	ut = u0 - un

c elastic response
c compute kinematic quantities
      do i= 1,3
	  do j= 1,3
	    nbnb(i,j) = nb(i)*nb(j)
	  enddo
	enddo
      do i= 1,3
	  do j= 1,3
	    nbu0(i,j) = nb(i)*u0(j)
	  enddo
	enddo

      if (uh.lt.uh0+eps) then
c Phase I
        c = t0/delta 
      else
c Phase II
        c = t0*exp(-a*dmd0**b)/delta 
	endif

c compute traction
        ua = (un + alpha*ut)
	  T = c*ua

c elastic tangent ddu
        ca = c*alpha
        ddu = 0.0d0
        do i= 1,3
	     ddu(i,i) = ca  
	  enddo	 
	  coa = c - ca
	  do i = 1,3
	    do j= 1,3
	      ddu(i,j) = ddu(i,j) + coa*nbnb(i,j)
	    enddo
	  enddo 
c elastic tangent ddn
        ddn=0.0d0
        pencoa = pen*coa
	  do i= 1,3
          ddn(i,i) = pencoa
	  enddo
	  do i= 1,3
	    do j= 1,3
	      ddn(i,j) = ddn(i,j) + coa*nbu0(i,j)
	    enddo
	  enddo


c add damage softening
	 if(dflag) then

c begin analytical tangent
      if (uh.lt.uh0+eps) then
c Phase I
        gamma = c/(delta*delta)
      else
c Phase II
	 gamma = c*(delta*(1.0d0+a*b*dmd0*b)-uh0)/(delta**2*dmd0)
	endif
	  
	do i= 1,3
	 do j= 1,3
	  ddu(i,j) = ddu(i,j) - gamma*ua(i)*u0(j)
	 enddo
      enddo
c end analytical tangent
	 endif
       
	
c add penalty contribution if necessary
        if (pen.gt.0.0d0) then
c add penalty to tractionc
	    t = t - cp*pen*pen*un
c add penalty to stiffness
          tcppen = 2.0d0*cp*pen
          do i= 1,3
		  do j= 1,3
	        ddu(i,j) = ddu(i,j) - tcppen*nbnb(i,j)
	        ddn(i,j) = ddn(i,j) - tcppen*nbu0(i,j)
		  enddo
           ddn(i,i) = ddn(i,i) + 0.5d0*pen
		enddo		 
	  endif


	end

       subroutine Trac11_viso_pen_sept05(l,d,u0,nc,T,ddu, ddn)
c==========================================================TCG.-28.09.2005
c  
c  Isotropic Traction model with isotropic damage
c
c    t = c u_0	; 
c   ddu = c II
c   ddn = 0
c
c   INPUT:
c         nc..........Normal vector at fictitious discontinuity
c         u0..........Gap displacement
c         uH..........history variable (max Gap)
c		d...........Material parameter vector
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n 
c                   
c
c   USED:
c         c, a,b..material parameters
c
c
c
c   local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'
	include  'tdata.h'


	logical dflag


	integer i,j, l, nhvd, n1, n2, lbh
      real*8 d(*), T(3), ddu(3,3),ddn(3,3), u0(3)
	real*8 a,b,c, te0, uh,tp(3),nc(3), u0p(3)
	real*8 un, ue, ue0, eps, te, gamma, pen,un2,twoun
	real*8 mu, u0n(3), f1, f2


      lbh= int(d(78))     ! lenght of bulk-related history     
      nhvd= int(d(75))    ! length of history per Gausspoint of discontinuity
     
	n1 = (l-1)*nhvd + nh1 + lbh
	n2 = (l-1)*nhvd + nh2 + lbh



c material parameters 

      te0 = d(20)
      a   = d(21)
      b   = d(22)
      mu  = d(23)
	ue0 = d(19)
	pen = d(30)


c get ue from history or from inital compliance
      ue =  hr(n1)
	do i= 1,3
	 u0n(i) = hr(n1+i)
	enddo

c compute displacement increment

	if (dt.ne.0.0d0) then
	  f1 = mu/dt
	  f2 = 1.0d0 + f1
	else
	  if(ttim.ne.0) then
	    write (*,100) 
	    write (iow,100)
	  endif
	endif

	if (ue.le.ue0) ue=ue0

c	te=te0*exp(-a*(ue-ue0)**b) 


      uh = sqrt(u0(1)*u0(1)+u0(2)*u0(2)+u0(3)*u0(3))
c decide if elastic or damage
      if (uh.lt.ue) then    !elastic
	  if (uh.lt.ue0) then
	    c=te0/ue0
	  else
          c = te0*exp(-a*(ue)**b)/ue
	  endif
	  t = c*(f2*u0 - f1*u0n)
c elastic tangent
         ddu = 0.0d0
         do i= 1,3
	     ddu(i,i) = f2*c  
	   enddo	  
	else  !damage
c traction        
	  c = te0*exp(-a*(uh)**b)/uh
	  t = c*(f2*u0 - f1*u0n)
        ddu=0.0d0
c elastic part of tangent
        do i= 1,3
	    ddu(i,i) = f2*c
	  enddo
c damage part of tangent (softening)
        gamma=c/uh*(a*b*(uh)**(b-1.0d0)+1.0d0/uh)
	  do i= 1,3
	    do j= 1,3
	      ddu(i,j) = ddu(i,j) - f2*gamma*u0(i)*u0(j)
     1                          + f1*gamma*u0n(i)*u0(j)
	    enddo
	  enddo
c update ue
        hr(n2) = uh
c update u0
	  do i= 1,3
	    hr(n2+i) = u0(i) 
	  enddo
	endif
            

c set Traction moduli	with respect to n	to zero
      ddn=0.0d0

c check penetration
      un= u0(1)*nc(1) + u0(2)*nc(2) + u0(3)*nc(3)


	if (un.gt.0.0d0) then  ! add penalty
	  twoun = un + un
	  un2 = un*un
c traction
	  t = t + pen*un2*nc
c stiffness
        do i= 1,3
	   do j= 1,3
	     ddu(i,j) = ddu(i,j) + twoun*pen*nc(i)*nc(j)
	   enddo
	  enddo

        do i= 1,3
	    ddn(i,i) = ddn(i,i) + pen*un2
	  enddo
        do i= 1,3
	   do j= 1,3
	     ddn(i,j) = ddn(i,j) + twoun*pen*nc(i)*u0(j)
	   enddo
	  enddo

	endif
100   format(
     1 5x, 'ERROR Trac11_viso_pen_sept05: Division by zero'/)
	end

       subroutine Trac11_iso_smooth(l,d,u0,nc,T,ddu, ddn)
c==========================================================TCG.-28.04.2006
c  
c  Smooth Isotropic Traction model with isotropic damage
c
c    t   = c u_0	;  c = e Exp[-a u_h] 
c    ddu = c II - a c/u_h (u_0 x u_0)
c    ddn = 0
c
c   INPUT:
c         nc..........Normal vector at fictitious discontinuity
c         u0..........Gap displacement
c         uH..........history variable (max Gap)
c		d...........Material parameter vector
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n 
c                   
c
c   USED:
c         c, a,b..material parameters
c
c
c
c   local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'


	logical dflag


	integer i,j, l, nhvd, n1, n2, lbh
      real*8 d(*), T(3), ddu(3,3),ddn(3,3), u0(3)
	real*8 a,b,c,e, te0, uh,tp(3),nc(3), u0p(3)
	real*8 un, ue, ue0, eps, te, gamma, pen,un2,twoun


      lbh= int(d(78))     ! lenght of bulk-related history     
      nhvd= int(d(75))    ! length of history per Gausspoint of discontinuity
     
	n1 = (l-1)*nhvd + nh1 + lbh
	n2 = (l-1)*nhvd + nh2 + lbh



c material parameters 

      e   = d(21)
      a   = d(22)
	pen   = d(30)


c get ue from history or from inital compliance
        ue =  max(hr(n1),1.0d-3)

c	ue=1.0d-3

c current state
        uh = sqrt(u0(1)*u0(1)+u0(2)*u0(2)+u0(3)*u0(3))

c decide if elastic or damage
        if (uh.lt.ue) then    !elastic
         c = e*exp(-a*(ue))
	   t = c*u0
c elastic tangent
          ddu = 0.0d0
          do i= 1,3
	      ddu(i,i) = c  
	    enddo	  
	  else  !damage


c traction        
          c = e*exp(-a*(uh))
	    T = c*u0
          ddu=0.0d0
c elastic part of tangent
          do i= 1,3
	      ddu(i,i) = c
	    enddo
c damage part of tangent (softening)
          gamma = a*c/uh
	    do i= 1,3
	      do j= 1,3
	        ddu(i,j) = ddu(i,j) - gamma*u0(i)*u0(j)
	      enddo
	    enddo
c update ue
          hr(n2) = uh
	  endif
            

c set Traction moduli	with respect to n	to zero
         ddn=0.0d0

c check penetration
      un= u0(1)*nc(1) + u0(2)*nc(2) + u0(3)*nc(3)

     

	if (un.gt.0.0d0) then  ! add penalty

c      write(iow,*) un


	  twoun = un + un
	  un2 = un*un
c traction
c	  t = t + pen*un2*nc
	  t = t + pen*un*nc
c stiffness
        do i= 1,3
	   do j= 1,3
c	     ddu(i,j) = ddu(i,j) + twoun*pen*nc(i)*nc(j)
	     ddu(i,j) = ddu(i,j) + pen*nc(i)*nc(j)
	   enddo
	  enddo

        do i= 1,3
c	    ddn(i,i) = ddn(i,i) + pen*un2
	    ddn(i,i) = ddn(i,i) + pen*un
	  enddo
        do i= 1,3
	   do j= 1,3
c	     ddn(i,j) = ddn(i,j) + twoun*pen*nc(i)*u0(j)
	     ddn(i,j) = ddn(i,j) + pen*nc(i)*u0(j)
	   enddo
	  enddo

	endif
        end
 
        
      subroutine Trac11_lin_elas(l,hn,nhvd,d,u0,nb,T,ddu,ddn)
c==========================================================CM 2022
c  
c  Linear elastic traction seperation law
c
c   INPUT:
c		d...........Material parameter vector      
c         u0..........Gap displacement
c         nb..........Normal vector at fictitious discontinuity      
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n 
c                   
c   USED:
c         K,..............Stiffness
c
c   local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'

      integer i, j, sof
      integer l, nhvd, n1, n2, lbh, n3
      
      real*8 sig0, K
      real*8 u0(3), d(*)
      real*8 T(3)
      real*8 hn(nhvd), nb(3)
      REAL*8 ddu(3,3), ddn(3,3)
      
c-----Parameters----------------------------------------------------------      
      
      K = d(21)
      
c-----Determine the Traction and stiffness--------------------------------      
      
      T  = K*u0
      
      do i= 1,3
        ddu(i,i) = K 
      enddo	  
      
      ddn=0.0d0    
      
	end          
        
              
    
      subroutine Trac11_iso_exp(l,hn,nhvd,d,u0,nb,T,ddu,ddn)
c==========================================================CM 2022
c  
c  Exponential isotropic traction seperation law
c
c   INPUT:
c		d...........Material parameter vector      
c         u0..........Gap displacement
c         nb..........Normal vector at fictitious discontinuity      
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n 
c                   
c   USED:
c         a,b,..............Shape related parameters
c         T0(1) or sig0.....Peak stress
c         sof...............Whether to include the softening in the stiffness (if nonzero)  
c
c   local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'

      integer i, j, sof
      integer l, nhvd, n1, n2, lbh, n3
      
      real*8 a, b, c, cd, sig0
      real*8 T0
      real*8 u0(3), d(*)
      real*8 u0u0(3,3)
      real*8 T(3), ddu(3,3), ddn(3,3)
      real*8 hn(nhvd), nb(3), del
      logical dflag 
      !real*8 un, ue, ue0, eps, te, gamma, pen,un2,twoun

c-----History Pointers-------------------------------------------------------------
     
      n3 = nh3 + 360   ! pointer to t0

c-----Define parameters-------------------------------------------------------------
   
      !T0 = hr(n3)
      sig0 = d(20)
      
      a = d(21)
      b = d(22)
      
c-----Determine damage evolution and whether softening part of stiffness should be used-------
      
      del = sqrt(u0(1)*u0(1)+u0(2)*u0(2)+u0(3)*u0(3))
      
      if (del.gt.hn(1)) then
        hn(1) = del
	  if (sof.eq.1) then
           dflag=.true.
        else
           dflag=.false.    
        endif 
	else
	  dflag=.false.
      endif
      
c-----Determine the normal and tangential stiffness constants-------------      
      
      c = (sig0/hn(1))*exp(-a*hn(1)**b) 
      !c = (T0/hn(1))*exp(-a*hn(1)**b) 
      
      if (dflag) then
        cd = (c/hn(1))*exp(1+a*b*hn(1)**b)*0.5d0
      endif
      
c-----compute kinematic quantities---------------------------------------------------         
      
      do i= 1,3
	  do j= 1,3
          if (dflag) then
             u0u0(i,j) = u0(i)*u0(j)
          endif
	  enddo
      enddo  
      
c-----Determine the Traction and stiffness--------------------------------
      
      T = c*u0
      
      !Stiffness wrt gap displacement (and damage if softening)
      ddu = 0.0d0
      do i= 1,3
         ddu(i,i) = c  
      enddo	       
      do i = 1,3
         do j= 1,3
            if (dflag) then
               ddu(i,j) = ddu(i,j) - cd*u0u0(i,j)
            endif           
         enddo
      enddo       
      
      !Stiffness wrt normal
      ddn = 0.0d0         
      
      end        
        
  

      subroutine Trac11_norm_tran(l,hn,nhvd,d,u0,nb,T,ddu,ddn)
c==========================================================CM 2022
c  
c  Exponential normalised transversely isotropic traction seperation law
c
c   INPUT:
c		d...........Material parameter vector      
c         u0..........Gap displacement
c         nb..........Normal vector at fictitious discontinuity      
c
c   OUTPUT:
c		T..........Traction	 
c		ddu........Traction moduli with respect to u0	 
c         ddn........Traction moduli with respect to n 
c                   
c   USED:
c         a,b,c,ds,e,f......Shape related parameters
c         fac...............Integral parameter relating to shape parameters 
c         T0(2).............Peak stress in normal and transverse directions
c         G(2)..............Fracture energy in normal and transverse directions
c         alp(2)............Normalisation parameters  
c         sof...............Whether to include the softening in the stiffness (if nonzero)  
c
c   local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'iofile.h'

      integer i, j, sof
      integer l, nhvd, n1, n2, lbh, n3
      
      real*8 a, b, c, ds, e, fac, k
      real*8 T0(2), G(2), alp(2)
      real*8 un(3), ut(3), u0(3), d(*)
      real*8 nb(3), nbnb(3,3), nbu0(3,3), tb(3), unnb(3,3), uttb(3,3)
      real*8 T(3), ddu(3,3), ddn(3,3)
      real*8 normn, normt, cn, ct, cdn, cdt
      real*8 hn(nhvd)
      real*8 check, temp, pen
      logical dflag1, dflag2  
      !real*8 un, ue, ue0, eps, te, gamma, pen,un2,twoun

c-----History Pointers-------------------------------------------------------------
     
      n3 = nh3 + 360   ! pointer to t0

c-----Define parameters-------------------------------------------------------------
   
      T0(1) = hr(n3)
      T0(2) = hr(n3+1)

      a = d(21)
      b = d(22)
      c = d(23)
      ds = d(24)
      e = d(25)
      fac = d(26)
      G(1) = d(27)
      G(2) = d(28)
      sof = int(d(29))
      
      do i=1,2
         alp(i) = G(i)/(T0(i)*fac)
      enddo

c-----normal and transverse gap displacement-----------------------------------------
      
      temp = u0(1)*nb(1)+u0(2)*nb(2)+u0(3)*nb(3)      
      un = temp*nb
      ut = u0 - un
      
c-----Determine damage evolution and whether softening part of stiffness should be used-------
      
      normn = sqrt(un(1)*un(1)+un(2)*un(2)+un(3)*un(3))
      
      if (normn.gt.hn(1)) then
        hn(1) = normn  
	  if (sof.eq.1) then
           dflag1=.true.
        else
           dflag1=.false.    
        endif 
	else
	  dflag1=.false.
      endif
      
      normt = sqrt(ut(1)*ut(1)+ut(2)*ut(2)+ut(3)*ut(3))
      
      if (normt.gt.hn(2)) then
	  hn(2) = normt
        if (sof.eq.1) then
           dflag2=.true.
        else
           dflag2=.false.    
        endif 
        tb = ut/normt
	else
	  dflag2=.false.
      endif
      
c-----Determine the normal and tangential stiffness constants-------------      
      
      cn = (T0(1)/hn(1))*((a*exp(-b*(hn(1)/alp(1))**c))
     1                  +((1-a)*exp(-ds*(hn(1)/alp(1))**e)))    
      
      ct = (T0(2)/hn(2))*((a*exp(-b*(hn(2)/alp(2))**c))
     1                  +((1-a)*exp(-ds*(hn(2)/alp(2))**e)))
      
      if (dflag1) then
        cdn = (T0(1)*hn(1)**-2d0)*(a*(b*c*(hn(1)/alp(1))**c+1)*exp(-b*c
     1  *(hn(1)/alp(1))**c) + (1-a)*(ds*e*(hn(1)/alp(1))**e+1)*exp(-ds*e
     1  *(hn(1)/alp(1))**e))
      endif
      
      if (dflag2) then
        cdt = (T0(2)*hn(2)**-2d0)*(a*(b*c*(hn(2)/alp(2))**c+1)*exp(-b*c
     1  *(hn(2)/alp(2))**c) + (1-a)*(ds*e*(hn(2)/alp(2))**e+1)*exp(-ds*e
     1  *(hn(2)/alp(2))**e))
      endif
      
c-----compute kinematic quantities---------------------------------------------------         
      
      do i= 1,3
	  do j= 1,3
	    nbnb(i,j) = nb(i)*nb(j)
          nbu0(i,j) = nb(i)*u0(j)
          if (dflag1) then
             unnb(i,j) = un(i)*nb(j)
          endif
          if (dflag2) then
             uttb(i,j) = ut(i)*tb(j)
          endif
	  enddo
      enddo  
      
c-----Determine the Traction and stiffness--------------------------------
      
      T = cn*un + ct*ut
      
      !Stiffness wrt gap displacement (and damage if softening)
      ddu = 0.0d0
      do i= 1,3
         ddu(i,i) = ct  
      enddo	       
      do i = 1,3
         do j= 1,3
            ddu(i,j) = ddu(i,j) + (cn-ct)*nbnb(i,j)
            if (dflag1) then
               ddu(i,j) = ddu(i,j) - cdn*unnb(i,j)
            endif
            if (dflag2) then
               ddu(i,j) = ddu(i,j) - cdt*uttb(i,j) 
            endif            
         enddo
      enddo       
      
      !Stiffness wrt normal
      ddn = 0.0d0
      do i = 1,3
         do j= 1,3
            ddn(i,j) = ddn(i,j) + (cn-ct)*nbu0(i,j)
         enddo
      enddo       
      
      !check = 1d0       
      
      end  


