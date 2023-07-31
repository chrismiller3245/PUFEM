      subroutine modc_11(flag,nhv,l,n,imat,d,fi,finv,detf,
     1          sig,aa,aas)
c
c--------------------------------------------------------------------71
c
c     Finite Deformation bulk Material Models for FEAP
c
c   imat=1   linear elastic material (St. Venant Kirchhoff)
c   imat=2   neo Hookean model
c   imat=3   fiber composite model with neoHookean matrix and HGO fiber model
c   imat=4   modified Delfino model
c   imat=5   fiber composite model with neoHookean matrix and GOH fiber model
c   imat=6   six parameter vessel wall model
c   imat=7   small strain Mises plasticity model (from FEAP)
c   imat=8   finite strain general pasticity model (from FEAP)
c   imat=9   Discrete nonlinear viscoelastic fiber model with damage 
c
c.... INPUT variables
c
c         flag       .true. if call from PUFEM element; .false. otherwise
c         nhv        Number of history variables at quadrature point
c         l          Quadrature point number
c         n          element number
c         imat       Mechanical material model
c         d(100)     Material parameters
c         fi(3,3)    Total deformation gradient
c         finv(3,3)  Inverse deformation gradient
c         detf       Total Jacobian determinant
c
c.... OUTPUT variables
c         sig(6)     CAUCHY stress tensor
c         aa(6,6)    CAUCHY (spatial) elastic moduli
c      
c.... History variables
c         nh1        pointer to FEAP t_n history vector 
c         nh2        pointer to FEAP t_n+1 history vector 
c         hn1(nhv)   History variables at Gauss points:
c                    at the constitutive routines' INPUT it holds variables from the last time point t_n, i.e. data that is stored in the t_n-related history of FEAP 
c                    at the constitutive routines' OUTPUT it holds variables from the present iteration, which then are stored into the t_n+1-related history of FEAP
c--------------------------------------------------------------------71
c
c..... Declare variable types
      implicit none
      logical flag
      integer nhv, ll, lint, l, nn, nn2, n
      include  'comblk.h'
      include  'hdata.h'
      
	integer i, imat,j,k
      real*8  detf, bb(6)
c..... Declare array types
      real*8  d(*),aa(6,6),aas(6,6)
	real*8  hn1(nhv), check(2400), check2(8)
      real*8  sig(9),fi(3,3),finv(3,3)
c
      include  'tdata.h'
      include 'Visco.h'

      elno = n
      gano = l
      
      call pzero(aa,36)
      call pzero(aas,36)
      
c.... Move history t_n-related history variables from the hr-array to hn1-array (read at the nh1 position from FEAP's hr vector)
      
      !hn1 = 0d0
      
      nn = nh1 + (l-1)*nhv    
      do i = 1,nhv
         hn1(i) = hr(nn+ i-1)
      enddo
      
      if (imat .eq. 1) then
c  linear elastic material (St. Venant Kirchhoff)
         call line_11(d,fi,  sig,aa)
	elseif(imat.eq.2) then
c  neo Hookean model
         call leftCauchy_11(Fi,bb) 
	   call neoh_11(d,bb,detf, sig,aa)      
	elseif(imat.eq.3) then
c  fiber composite model with neoHookean matrix and HGO fiber model
         call composite_11(d,Fi,detf,  sig,aa)
	elseif(imat.eq.4) then
c  modified Delfino model
         call delfino_11(d,Fi,detf,  sig,aa)
	elseif(imat.eq.5) then
c  fiber composite model with neoHookean matrix and GOH fiber model
         call dist_composite_11(d,Fi,detf,  sig,aa)
	elseif(imat.eq.6) then
c  six parameter vessel wall model
         call six_para_11(d,Fi,detf,  sig,aa )
      elseif(imat.eq.7) then
c  small strain Mises plasticity model (from FEAP)
         call plasticity_11(flag,d,Fi,detf, nhv, hn1,  sig,aa)   
      elseif(imat.eq.8) then
c   finite strain general pasticity model (from FEAP)
         call plasticity_11_finite(flag,d,Fi,detf, nhv, hn1,  sig,aa)  
      elseif(imat.eq.9) then
c   Discrete nonlinear viscoelastic fiber model with damage   
         call dis_nonlin_vis_dam_mod(flag,d,Fi,detf,hn1,
     1         nhv,sig,aa,aas)  
         
         !aa = 0d0
         !do i=1,6
         !    do j=1,6
         !       aa(i,j) = aas(i,j)    
         !    enddo
         !enddo
         
	else
	  write (*,*)
        write(*,*)'=================================='
        write(*,*)'ERROR --> NO MATERIAL MODEL CALLED'
        write(*,*)'=================================='
	  write (*,*)
      endif
      
c If were not using the collagen fiber damage model (imat=9) then we need an arbitrary output for aas for all other cases

      if (imat.ne.9) aas = aa      
      
c
c.... convert to Cauchy stresses and tangent moduli if necessary i.e. model formulated in kirchoff stress
c
      if((imat.eq.1).or.(imat.eq.2).or.(imat.eq.3).or.
     1   (imat.eq.4).or.(imat.eq.5).or.(imat.eq.6).or.
     1   (imat.eq.9)) then       
          do  i = 1,6
            sig(i) = sig(i)/detf
            do  j = 1,6		  
              aa(i,j) = aa(i,j)/detf
              aas(i,j) = aas(i,j)/detf
	      enddo
          enddo
      endif
      
c update history (write at the nh2 position into FEAP's hr vector)

      nn = nh2 + (l-1)*nhv  
      do i = 1,nhv
         hr(nn+ i-1) = hn1(i) 
      enddo          
      
9999  end


      subroutine plasticity_11(flag,d,Fn1,detf, nhv, Hn1,  sig,aa)
c=======================================================TCG.05.05.2021
c  
c  interface to small strain Mises plastictity model of FEAP
c
c   INPUT:
c         flag       .true. if call from PUFEM element; .false. otherwise
c		d...........Material parameter vector
c		Fn1(3,3)...... Deformation gradient at t_n+1
c	    detf.........Jacobian determinant det[Fn1]
c	    nhv.........length of History vector
c	    Hn.........History
c
c   OUTPUT:
c		sig..........Cauchy stress	 
c		aa...........Spatial Eulerian tangent	 
c
c   USED:
c		eps......strain (small)
c
c   HISTORY:
c       Hn1(1)        epp(1)  Cumulative plastic strain
c       Hn1(2)        epp(2)  
c       Hn1(3)        eps(1)  linear strain
c       Hn1(4)        eps(2)  linear strain
c       Hn1(5)        eps(3)  linear strain
c       Hn1(6)        eps(4)  linear strain
c       Hn1(7)        eps(5)  linear strain
c       Hn1(8)        eps(6)  linear strain
c
c==============================================
c..... Declare variable types

      IMPLICIT NONE

      logical flag
      integer nhv, i,j,k, istrt
      real*8 d(*), detf
      real*8  Fn1(3,3), Hn1(nhv)
      real*8 eps(6), epst(3,3)
	real*8 sig(9), aa(6,6), aar(6,6), temp1, temp2

c  compute strain
      do i=1,3
          do j= 1,3
              epst(i,j) = 0.5d0*(Fn1(i,j) + Fn1(j,i))
          enddo
      enddo
      do i= 1,3
          epst(i,i) = epst(i,i) - 1.0d0
      enddo
     
      eps(1)= epst(1,1) 
      eps(2)= epst(2,2) 
      eps(3)= epst(3,3) 
      eps(4)= 2.0*epst(1,2) 
      eps(5)= 2.0*epst(2,3) 
      eps(6)= 2.0*epst(1,3) 
     
c call FEAP plasticity model 
      istrt=0  ! Start state: 0 = elastic; 1 = last solution
c change temporarely entris of the d vector
      temp1 = d(21)
      temp2 = d(27) 
      d(21) = d(1)
      d(27) = d(2)
      call mises11(flag,d,eps,Hn1(3),Hn1(1),6,istrt, sig,aa,aar)     
      d(21) = temp1
      d(27) = temp2
      
c pass on for plotting  
      sig(7) = hn1(1)   ! cumulative plastic strain
      
1000  end      
      
 
      subroutine mises11(flag,d,eps,epsp,epp,ntm,istrt, sig,dd,dr)

c     Mises (J2) plasticity with isotropic and kinematic hardening
c     modified from FEAP

      implicit  none

      include  'counts.h'
      include  'elauto.h'
      include  'iofile.h'
      include  'pconstant.h'
      include  'sdata.h'
      include  'setups.h'
      include  'tdata.h'

      logical   conv,state, flag
      integer   i, j,istrt, ntm, count, mm,m1
      real*8    d(*),eps(6),epsp(6),epp(2), sig(6), dd(6,6),dr(6,6)
      real*8    alp(6), ep(6), en(6), xi(6)
      real*8    aa,bb,cc, press, theta, dlam, lam, xin, tolc, dot
      real*8    k, g, gbar, twog, hi,hi23,his23,hk,hk23, epss, dyld
      real*8    r0,rinf,rn,rb, beta, expb,expl, bbar, s23
      real*8    ff,phi,dphi, Tbar, vbar

      save

      data      tolc / 1.d-10 /
      
      
c     Check state for iterations

      if(niter.eq.0) then         ! First iteration in step
        if(istrt.eq.0) then       ! Elastic state requested
          state = .false.
          dyld  =  0.0d0
        else                      ! Last state requested
          state = .true.
          dyld  =  1.0d-08*d(41)
        endif
      else                        ! Not first iteration in step
        state = .true.
        dyld  =  0.0d0
        if(rank.gt.0) dyld = -1.0d-08*d(41)
      endif

c     Set parameters

      s23   = sqrt(two3)

      g     = d(27)                ! shear moduls
      k     = d(21) - 2.d0*two3*g  ! bulk modulus

      r0    = s23*d(41)  ! Initial yield stress 
      rinf  = s23*d(42)  ! Final yield stress
      beta  =     d(43)  ! Exponential hardening rate

      hi    = d(44)   ! Isotropic hardening modulus
      hk    = d(45)   ! Kinematic hardening modulus

c     Viscoplastic parameters

      if(dt.gt.0.0d0) then
        Tbar  = d(180)/dt
      else
        Tbar  = 0.0d0
      endif
      mm = nint(d(181))
      m1 = mm - 1

c     Compute constants

      hk23  = two3*hk
      hi23  = two3*hi
      his23 =  s23*hi

      twog  = 2.d0*g

      expb  = exp(-beta*epp(1))
      bbar  =  s23*beta*(rinf - r0)*expb
      gbar  = twog + hi23 + hk23

c     Compute volumetric and deviatoric strains

      theta = (eps(1) + eps(2) + eps(3))*one3
      do i = 1,3
        ep(i)   = eps(i) - theta
      end do ! i
      do i = 4,ntm
        ep(i) = eps(i)*0.5d0
      end do ! i

      theta = theta*3.0d0

c     Compute trial values

      do i = 1,ntm
        sig(i) = sig(i) + twog*(ep(i) - epsp(i))
        alp(i) = hk23*epsp(i)
        xi(i)  = sig(i) - alp(i)
      end do ! i

c     Compute trial norm of stress - back stress

      xin = sqrt(dot(xi,xi,ntm) + dot(xi(4),xi(4),ntm-3))
      rn  = rinf + his23*epp(1)
      rb  = (r0  - rinf)*expb

c     Check yield

cTCG      if(xin+dyld .gt. (rn+rb) .and. state) then
      if(xin+dyld .gt. (rn+rb) .and. state .and. not(flag)) then

c       Compute strain factor

        if(r0.gt.0.0d0) then
          epss   = r0/twog
        else
          epss   = 1.e-4
        endif

c       Compute viscoplastic or plasticity consistency

        conv   = .false.
        count  =  0

c       Plasticity

        if(d(180).eq.0.0d0)  then
          vbar = 0.0d0
          lam  = (xin - rn - rb) / (gbar + bbar)
          do while (.not.conv .and. count.le.25)
            expl = exp(-s23*beta*lam)
            dlam = (xin - rn - rb*expl - gbar*lam)/(gbar + bbar*expl)
            lam  = lam + dlam
            if(abs(r0-rinf).le.tolc*r0) then
              conv = .true.
            elseif(abs(dlam) .lt. tolc*abs(lam) .or.
     &         abs(lam)  .lt. tolc*epss)          then
              conv = .true.
            endif
            count = count + 1
          end do ! while

c       Viscoplasticity

        else
          lam = 0.0d0
c         lam = (xin - rn - rb) / (gbar + bbar)
          do while (.not.conv .and. count.le.25)
            expl = exp(-s23*beta*lam)
            ff   = xin - rn - rb*expl - gbar*lam
            if(m1.gt.0) then
              dphi  = (ff/r0)**m1 / r0
            else
              dphi  = 1.d0/r0
            endif
            phi  = dphi*ff
            dphi = dphi*mm

            dlam = (phi - Tbar*lam)/(dphi*(gbar + bbar*expl) + Tbar)
            lam  = lam + dlam
            if(abs(dlam) .lt. tolc*abs(lam) .or.
     &         abs(lam)  .lt. tolc*epss)          then
              conv = .true.
            endif
            count = count + 1
          end do ! while

c         Mask lambda to be positive

          lam  = max(0.0d0,lam)
          vbar = Tbar/dphi

        endif

c       Warning: Not converged

        if(.not.conv .and. niter.gt.0) then
          write(  *,*) '  *WARNING* No convergence in MISES'
          write(iow,*) '  *WARNING* No convergence in MISES'
          write(iow,*) '   lam ',lam,' dlam ',dlam,' count ',count
c         call plstop()
        endif

c       Set auto time stepping factor

        rmeas = max(rmeas,2.5d0*rvalu(1)*s23*lam/epss)

c       Compute normal vector

        do i = 1,ntm
          en(i) = xi(i)/xin
        end do ! i

        if(lam.gt.0.0d0) then
          bb     = twog*lam/xin
          aa     = twog*(1.d0 - bb)
          bb     = twog*(bb - twog/(gbar + bbar*expl + vbar))
          epp(2) = 1.0d0
        else
          aa     = twog
          bb     = 0.0d0
          epp(2) = 0.0d0
        endif

c       Compute plastic tangent from normal

        do i = 1,ntm
          cc = bb*en(i)
          do j = 1,ntm
            dd(i,j) = cc*en(j)
          end do ! j
        end do ! i

c     Set for elastic increment

      else

        epp(2) = 0.0d0
        do i = 1,ntm
          en(i) = 0.0d0
        end do ! i

        lam = 0.0d0
        aa  = twog

        do i = 1,ntm
          do j = 1,ntm
            dd(i,j) = 0.0d0
          end do ! j
        end do ! i

      endif

c     Compute deviator stress, plastic strain, & accumul plastic strain

      do i = 1,ntm
        sig(i)  = sig(i)  - twog*lam*en(i)
        epsp(i) = epsp(i) + lam*en(i)
      end do ! i
      epp(1) = epp(1) + s23*lam

c     Add pressure

      press = k*theta
      do i = 1,3
        sig(i) = sig(i) + press
      end do ! i

c     Compute tangent moduli

      cc = k - aa*one3
      bb = k - twog*one3

      do i = 1,3
        do j = 1,3
          dd(j,i) = dd(j,i) + cc
          dr(j,i) =           bb
        end do ! j
        dd(i  ,i  ) = dd(i  ,i  ) + aa
        dr(i  ,i  ) = dr(i  ,i  ) + twog
      end do ! i
      do i = 4,ntm
        dd(i,i) = dd(i,i) + 0.5d0*aa
        dr(i,i) = g
      end do ! i
      end

      
      subroutine plasticity_11_finite(flag,d,Fn1,detf, 
     1           nhv, Hn1,  sig,aa)
c=======================================================TCG.15.04.2021
c  
c  interface to plastictity models of FEAP
c
c   INPUT:
c         flag       .true. if call from PUFEM element; .false. otherwise
c		d...........Material parameter vector
c		Fn1(3,3)...... Deformation gradient at t_n+1
c	    detf.........Jacobian determinant det[Fn1]
c	    nhv.........length of History vector
c	    Hn.........History
c
c   OUTPUT:
c		sig..........Kirchhoff stress	 
c		aa...........Spatial Eulerian tangent	 
c
c   USED:
c		F(3,3,*)......deformation gradients and displacement gradients
c             F(3,3,1).........Deformation gradient at t_n+1
c             F(3,3,2).........Deformation gradient at t_n
c             F(3,3,3).........ref. displacement gradient at t_n+1
c             F(3,3,4).........ref. displacement gradient at t_n
c
c   HISTORY:
c       Hn1(1)        epp(1)  Cumulative plastic strain
c       Hn1(2)        epp(2)  | SIGMA_n | general plasticity
c       Hn1(3)        be(1)  Left Cauchy-Green tensor
c       Hn1(4)        be(2)  Left Cauchy-Green tensor
c       Hn1(5)        be(3)  Left Cauchy-Green tensor
c       Hn1(6)        be(4)  Left Cauchy-Green tensor
c       Hn1(7)        be(5)  Left Cauchy-Green tensor
c       Hn1(8)        be(6)  Left Cauchy-Green tensor
c       Hn1(9)        epl(1) Plastic Strain for Hardening
c       Hn1(10)       epl(2) Plastic Strain for Hardening
c       Hn1(11)       epl(3) Plastic Strain for Hardening
c       Hn1(12)       F_n(1,1) Deformation gradient
c       Hn1(13)       F_n(1,2) Deformation gradient
c       Hn1(14)       F_n(1,3) Deformation gradient
c       Hn1(15)       F_n(2,1) Deformation gradient
c       Hn1(16)       F_n(2,2) Deformation gradient
c       Hn1(17)       F_n(2,3) Deformation gradient
c       Hn1(18)       F_n(3,1) Deformation gradient
c       Hn1(19)       F_n(3,2) Deformation gradient
c       Hn1(20)       F_n(3,3) Deformation gradient
c
      
c==============================================
c..... Declare variable types

      IMPLICIT NONE

	logical state, flag
      integer nhv, i,j,k, istrt
      real*8 d(*),F(3,3,4), detf
      real*8  Fn1(3,3), Hn1(nhv)
	real*8 sig(9), aa(6,6), temp1, temp2

c  deformation gradients and displacement gradients at t_n+1 
      do  i = 1,3
          do j=1,3
              F(i,j,1) = Fn1(i,j)
              F(i,j,3) = Fn1(i,j)
          enddo
      enddo
      do i=1,3
          F(i,i,3) = F(i,i,3) - 1.0d0
      enddo

c  deformation gradient at t_n (read from history)
      k=12      
      do i=1,3
          do j=1,3
            F(i,j,2) = Hn1(k)  
            k=k+1
          enddo
      enddo

c  displacement gradient at t_n 
      do i=1,3
          do j=1,3
            F(i,j,4) =  F(i,j,2)    
          enddo
      enddo
      do i=1,3
          F(i,i,4) = F(i,i,4) - 1.0d0
      enddo
     
C call FEAP plasticity models 
      istrt=0  ! Start state: 0 = elastic; 1 = last solution
      state = .false.
c change temporarely entris of the d vector
      temp1 = d(21)
      temp2 = d(22) 
      d(21) = d(1)
      d(22) = d(2)
      call plasfd(d,detf,F, hn1(1),hn1(3),hn1(9),
     &                  6,istrt, aa,sig,3,state)
      d(21) = temp1
      d(22) = temp2
      
c update history (deformation gradient only)
      k=12      
      do i=1,3
          do j=1,3
            Hn1(k) = F(i,j,1)   
            k=k+1
          enddo
      enddo
      
c pass on for plotting  
      sig(7) = hn1(1)   ! cumulative plastic strain
      
1000  end      

      subroutine line_11(d,f,  sig,aa)
c=======================================================TCG.28.08.2002
c  
c  Linear isotropic material   St. Venant Kirchhoff
c
c   INPUT:
c		d...........Material parameter vector
c		F(3,3)......Deformation gradient
c	    jac.........Jacobian determinant jac=det[F]
c
c   OUTPUT:
c		sig..........Kirchhoff stress	 
c		aa...........Spatial Eulerian tangent	 
c
c   USED:
c         S............Second Piola Kirchhoff stress
c         CC...........Material Tangent
c         E............Young's modulus
c         b............Left Cauchy Green strain
c         b2...........Square of Left Cauchy Green strain
c         tr...........Trace of Euler Lagrange stain
c         nu...........Poisson's ratio 
c         lamb.........First Lamme constant
c         mu...........Second Lamme constant
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

	integer i,j
      real*8 d(*), F(3,3), b(6), b2(6), bb(6,6)
	real*8 sig(6), aa(6,6)
	real*8 E, nu, mu, lamb, tr

c material parameters
       E = d(1)
	 nu = d(2)

	if (nu.eq.0.5d0) nu = 0.499

c compute left Cauchy Green strain and square
	 call leftCauchy_11 (F,b)
       call quad (b,b2)
c Trace of Euler Lagrage strain tensor
       tr=0.5d0*(b(1)+b(2)+b(3)-3.0d0)
c Lamme constants
      lamb = E*nu/((1.0d0+nu)*(1.0d0-2.0d0*nu))
	 mu  = E/(2.0d0+2.0d0*nu) 
c compute Kirchhoff stress
       b2 = b2-b
	 sig = lamb*tr*b + mu*b2

c compute stiffness
	 do i=1,6
	  do j= i,6
	   bb(i,j) = b(i)*b(j) 
	  enddo
	 enddo
c part 1
       aa(1,1) = 2.0d0*bb(1,1)
       aa(1,2) = 2.0d0*bb(4,4)
       aa(1,3) = 2.0d0*bb(6,6)
       aa(1,4) = 2.0d0*bb(1,4)
       aa(1,5) = 2.0d0*bb(4,6)
       aa(1,6) = 2.0d0*bb(1,6)

       aa(2,2) = 2.0d0*bb(2,2)
       aa(2,3) = 2.0d0*bb(5,5)
       aa(2,4) = 2.0d0*bb(2,4)
       aa(2,5) = 2.0d0*bb(2,5)
       aa(2,6) = 2.0d0*bb(4,5)

       aa(3,3) = 2.0d0*bb(3,3)
       aa(3,4) = 2.0d0*bb(5,6)
       aa(3,5) = 2.0d0*bb(3,5)
       aa(3,6) = 2.0d0*bb(3,6)

       aa(4,4) = bb(1,2) + bb(4,4)
       aa(4,5) = bb(4,5) + bb(2,6)
       aa(4,6) = bb(1,5) + bb(4,6)

       aa(5,5) = bb(2,3) + bb(5,5)
       aa(5,6) = bb(3,4) + bb(5,6)

       aa(6,6) = bb(1,3) + bb(6,6)
c part 2
       aa = mu*aa + lamb*bb

c sym
       do i= 1,6
	  do j= i+1,6
	   aa(j,i) = aa(i,j)
	  enddo
	 enddo
1000  end


      subroutine neoh_11(d,bb,detf, tau,aa)
c
c--------------------------------------------------------TCG.28.08.2002
c           Compressible neo-Hookean model for FEAP
c                   (J_2/3 regularization)
c             _ _       _                    _
c           W(J,bb) = U(J) + 0.5*mu*(J^(-2/3)bb:1 - 3)
c           U(J)    = 0.5*lambda*(ln J)**2
c
c.... INPUT variables
c         d(*)       Material properties
c         bb(6)      Left Cauchy-Green tensor
c         detf       Jacobian determinant at Gauss points
c
c.... OUTPUT variables
c         tau(6)     Deviatoric Kirchhoff stress tensor
c         aa(6,6)    Deviatoric tangent moduli
c
c--------------------------------------------------------------------71
c
c..... Declare variable types
      integer i, j
      real*8  detf, two3, xj23, trbb3, taui, mub, mub2, mub3
c..... Declare array types
      real*8 d(2), bb(6), tau(6), aa(6,6)
      real*8 check
      data two3/0.6666666666666667d0/
      
      if (d(2).gt.1e3) then
         check = 1d0    
      endif      
      
c
c.... Compute J^2/3 * bb
c
      
      xj23 = detf**(- two3)
      do 100 i = 1,6
        bb(i) = xj23*bb(i)
100   continue
c
c.... Compute deviatoric Kirchhoff stress and map to tensor form
c
      trbb3  = (bb(1) + bb(2) + bb(3))/3.d0
c
      do 150 i = 1,3
        tau(i  ) = d(2)*(bb(i) - trbb3)
        tau(i+3) = d(2)*bb(i+3)
150   continue
c
c.... Rank one update: -2/3 tau x g - 2/3 g x tau
      do 210 i = 1,6
        taui = two3*tau(i)
        do 200 j = 1,3
          aa(i,j) =  aa(i,j) - taui
          aa(j,i) =  aa(j,i) - taui
200     continue
210   continue
c                       __
c.... Deviatoric term 2 mu [ I - 1/3 g x g ]
c
      mub  = d(2)*trbb3
      mub2 = mub + mub
      mub3 = mub2/3.d0
c
      do 230 i = 1,3
        aa(i  ,i  ) = aa(i  ,i  ) + mub2
        aa(i+3,i+3) = aa(i+3,i+3) + mub
        do 220 j = 1,3
          aa(i ,j ) = aa(i ,j )   - mub3
220     continue
230   continue
      end


      subroutine vol_11(d,det, sig, aa)
c
c-------------------------------------------------------TCG.28.08.2002
c
c             Compute part of the stresses and tangent
c             (for displacement FE formulation)
c
c.... variables
c
c       d      Material parameters
c       det    Jacobian determinat
c       sig    Cauchy stresses
c       aa     Eulerian moduli
c       up     partial U/ partial det
c       upp    partial^2 U/ partial det^2
c
c
c  aktiv volumetric potential
c 
c
c         Model 2. (ACTIVE)     U(J) = K*0.5*(log J)^2
c                                up  = K*(log(J) - const )/J
c                                upp = ( K/J     - up )   /J
c
c--------------------------------------------------------------------71
c
	implicit none
	real*8 d(*), det, sig(6), aa(6,6)
	real*8 up, upp,  temp, K
	integer i


	K=d(1)
	!up  = K*(log(det))/det
	!upp = ( K/det- up )/det
      up  = K*(log(det))
	upp = ( K/det- up )

c  add volumetric contribution to stress
      do i= 1,3
	 sig(i) = sig(i) + up 
	enddo
c  add volumetric contribution to tangent
      temp = up + upp*det
	aa(1,1) = aa(1,1) + temp
	aa(1,2) = aa(1,2) + temp
	aa(1,3) = aa(1,3) + temp
	aa(2,1) = aa(2,1) + temp
	aa(2,2) = aa(2,2) + temp
	aa(2,3) = aa(2,3) + temp
	aa(3,1) = aa(3,1) + temp
	aa(3,2) = aa(3,2) + temp
	aa(3,3) = aa(3,3) + temp
    
     	temp =  up 
	aa(4,4) = aa(4,4) - temp
	aa(5,5) = aa(5,5) - temp
	aa(6,6) = aa(6,6) - temp

	temp = temp + up 
	aa(1,1) = aa(1,1) - temp
	aa(2,2) = aa(2,2) - temp
	aa(3,3) = aa(3,3) - temp


      end


      subroutine composite_11(d,f,jac,  tau,dd)
c==========================================================TCG.-04,03,2002
c  
c  Large-strain composite material
c
c   INPUT:
c		d...........Material parameter vector
c			d(4)=k1
c			d(5)=k2
c		F(3,3)......Deformation gradient
c	    jac.........Jacobian determinant jac=det[F]
c
c   OUTPUT:
c		tau..........Kirchhoff stress	 tau = 2 dPsi/dg
c		dd...........Eulerian tangent	 dd = 2 dtau/dg
c
c   USED:
c		Ft(3,3)......Modified deformation gradient,  Ft=jac^(-1/3)*F
c		anl(3).......Unit fibre direction vector in lagrangian description
c		ane(3).......Unit fibre direction vector in eulerian description
c		b(6).........Left Cauchy-Green tensor
c		tau_f........Kirchhoff stress due to fibres
c		dd_f.........Eulerian tangent due to fibres
c		one3.........1/3
c		jm13.........jac**(-one3)
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

	integer i,j,ii
      real*8 d(*), F(3,3), jac, Ft(3,3)
	real*8 tau(6), dd(6,6), b(6), k1, k2
	real*8 tau_f(6), dd_f(6,6)
	real*8 tau_m(6), dd_m(6,6)
	real*8 anl(3), ane(3),one3, jm13
	data one3  /0.33333333333333333/

      include 'eldata.h'
	include  'vecdata.h'


	jm13=jac**(-one3)

 	k1= d(3)
	k2= d(4)
c initialize
 	dd   =0.0d0
 	tau  =0.0d0
	dd_m =0.0d0
	tau_m=0.0d0

c  load direction vector (unit)
      do ii= 1,2
c initialize
 	 dd_f =0.0d0
	 tau_f=0.0d0
	  do j= 1,3
         anl(j) = rvec(ii,j,n)
	  enddo

c compute modified elastic deformation gradient
       Ft=jm13*F
c compute direction vector in eulerian description
       ane = 0.0d0
       do i= 1,3
	  do j= 1,3
	  ane(i)= ane(i) + Ft(i,j)*anl(j)
	  enddo
	 enddo

c compute Kirchhoff stress and eulerian moduli of the fibre
	  call fibre_11 (k1, k2 ,ane, tau_f, dd_f)
	  tau = tau + tau_f
	  dd = dd + dd_f
	enddo
c=REM.2================================================================
c add contribution to stress and elasticity due to matrix material if 
c needed 
	call leftCauchy_11 (F,b)
      call neoh_11(d,b,jac, tau_m,dd_m)
c======================================================================
c summation
	tau = tau + tau_m
	dd  = dd + dd_m

      end


	subroutine fibre_11(k1,k2,ae, tau, dd)
c==========================================================TCG.-29.09.1999
c
c   compute the kirchhoff stress and eulerian tangent moduli 
c   of a transversally isotropic material based on a fibre-
c   reinforced composite.
c
c  Input:
c		k1.........First fiber parameter
c		k2.........Second fiber parameter
c         ae.........Eulerian fibre-vector ae_i=F_ij*al_j 
c                    (not identity norm!)
c
c  Output:
c		tau........Kirchhoff stress
c		dd.........Eulerian tangent
c
c  Used:
c         I4.........Fourth invariant I4 = ae_i*ae_i
c         w4.........First derivative of the stored energy function
c                    with respect to I4
c         w44........Second derivative of the stored energy function
c                    with respect to I4
c         axa........(ae x ae)
c         dev_axa....dev (ae x ae)
c         tempi......parts of the material tangent
c
c=======================================================================
 	IMPLICIT none
      real*8  tau(6), ae(3), dd(6,6), I4, w4, w44
	real*8 temp1(6,6), temp2(6,6), temp3(6,6), temp4(6,6),temp5(6,6)
	real*8 dev_axa(6), axa(6), one3, k1, k2
	integer i,j
	data one3  /0.33333333333333333/
c compute I4
      I4=0.0d0
      do i= 1,3
	  I4 = I4 + ae(i)*ae(i)
	enddo
c restrict computation to stretched fibres
	if (I4.lt.1.0d0)   then
	 tau=0.0d0
	 dd=0.0d0
	 goto 1000
	endif
c get w4 and w44
      call Psi_I4_11 (k1,k2,I4, w4,w44)
c derive Kirchhoff stress
      call vec_prod_11(ae,axa)
	dev_axa = axa
	call dev_euler_11(dev_axa)

	tau=2.0d0*w4*dev_axa
c derive spatial tangent
	temp1=0.0d0
      do i= 1,6
	 do j= i,6
	  temp5(i,j) = dev_axa(i)*dev_axa(j)
	 enddo
	 do j= i,3
	  temp1(i,j) = axa(i)
	 enddo
	enddo
	temp2=0.0d0
	temp3=0.0d0
      do i= 1,3
	 do j= i,6
	  temp2(i,j) = axa(j)
	 enddo
      do j= i,3
	  temp3(i,j) = 1.0d0
	enddo
	enddo
	temp4=0.0d0
	temp4(1,1)=-1.0d0
	temp4(2,2)=-1.0d0
	temp4(3,3)=-1.0d0
	temp4(4,4)=-0.5d0
	temp4(5,5)=-0.5d0
	temp4(6,6)=-0.5d0

      dd=0.0d0
	dd=(4.0d0*w4/3.0d0)*(one3*I4*temp3-temp1-temp2-I4*temp4)
	dd = dd + 4.0d0*w44*temp5 
c  symmetrize tangent
 	do i = 1,6
	 do j= i+1,6
	 dd(j,i) = dd(i,j)
	 enddo
	enddo

1000	end


	subroutine vec_prod_sym_11(a,b,axbpbxa)
c==========================================================TCG.-29.09.1999
c
c   compute the symmetric vector product of the vector a and b
c
c  Input:
c		a..........Vector
c		b..........Vector
c
c  Output:
c		axbpbxa....Vectorproduct a x b + b x a
c
c=======================================================================
 	IMPLICIT none
      real*8 a(3), b(3), axbpbxa(6), temp
	integer i

	do i= 1,3
	 temp = a(i)*b(i)
	 axbpbxa(i) = temp + temp
	enddo
	axbpbxa(4)=a(1)*b(2) + b(2)*a(1)
	axbpbxa(5)=a(2)*b(3) + b(3)*a(2)
	axbpbxa(6)=a(1)*b(3) + b(3)*a(1)

	end


	subroutine vec_prod_11 (a,axa)
c==========================================================TCG.-29.09.1999
c
c   compute the vector product of the vector a
c
c  Input:
c		a..........Vector
c
c  Output:
c		axa........Vectorproduct a x a
c
c=======================================================================
 	IMPLICIT none
      real*8 a(3), axa(6)
	integer i

	do i= 1,3
	 axa(i) = a(i)*a(i)
	enddo
	axa(4)=a(1)*a(2)
	axa(5)=a(2)*a(3)
	axa(6)=a(1)*a(3)

	end


	subroutine dev_euler_11(a)
c==========================================================TCG.-29.09.1999
c
c  compute deviator (eulerian)  from a
c
c         deva = a - 1/3 tra I
c
c  Input:
c		a........Second order tensor
c    
c  Output:
c         deva.....Deviator of a
c
c  Used:
c         tra......Trace of a
c         one3.....1/3
c
c===============================================================
 	IMPLICIT none
	real*8 a(6), tra, one3
	integer i
	data one3  /0.33333333333333333/
	tra = a(1)+a(2)+a(3)
	tra = one3*tra 
	do i= 1,3
	 a(i) = a(i) - tra
	enddo
      end	 


      subroutine Psi_I4_11(k1,k2,I4, w4,w44)
c==========================================================TCG.-29.09.1999
c
c  derivatives with respect to the fourth invariants
c
c     Wt^aniso = (k1/2*k2)*(exp(k2*(I4-1)^2)-1)
c
c  Input:
c         k1.......First material parameter
c         k2.......Second material parameter
c         I4.......Fourth invariant
c    
c  Output:
c         w4.......First derivative of the stored-energy function 
c                  with respect to the fourt invariant     
c         w44......Second derivative of the stored-energy function 
c                  with respect to the fourt invariant     
c
c  Used:
c         I4m1.....I4-1
c         temp.....exp(k2*I4m1*I4m1)
c
c===============================================================
 	IMPLICIT none
      real*8 w4, w44, I4, k1, k2 ,I4m1 , temp
      
c get material parameters

	I4m1= I4-1.0d0
	temp= exp(k2*I4m1*I4m1)

      w4  = k1*I4m1*temp    
      w44 = k1*(1.0d0+2.0d0*k2*I4m1*I4m1)*temp
      end


      subroutine delfino_11(d,f,jac,  tau,aa)
c==========================================================TCG.-26.11.2003
c  
c  Exponential isotropic material for biological tissues
c  modified Delfino potential
c
c   INPUT:
c		d...........Material parameter vector
c		F(3,3)......Deformation gradient
c	    jac.........Jacobian determinant jac=det[F]
c
c   OUTPUT:
c		tau..........Kirchhoff stress	 tau = 2 dPsi/dg
c		aa...........Eulerian tangent	 dd = 2 dtau/dg
c
c   USED:
c		bb(6)........Modified Left Cauchy-Green tensor
c         I1...........First modified invariant
c         psi1.........First derivative of PSI with respect to I1
c         psi11........Second derivative of PSI with respect to I1
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

	integer i,j
      real*8 d(*), F(3,3), jac
	real*8 tau(6), aa(6,6)
	real*8 k1, k2, mu, two3, taui
	real*8 xj23, bb(6), psi1, psi11, I1
	real*8 psi1b, psi1b2, psi1b3, dev_bb(6), trbb3
	data two3  /0.666666666666667d0/

c material parameters
       mu = d(2)
	 k1 = d(3)
	 k2 = d(4)

c compute left Cauchy Green strain
	 call leftCauchy_11 (F,bb)
 
c
c.... Compute J^2/3 * bb
c
      xj23 = jac**(- two3)
      bb = xj23*bb
c
c.... Compute first modified invariant I1 = tr_bb
c
	I1 = bb(1) + bb(2) + bb(3)
      trbb3  = I1/3.d0
c.... Compute dev_bb and map to tensor form
      do i = 1,3
        dev_bb(i  ) = bb(i) - trbb3
        dev_bb(i+3) = bb(i+3)  
	enddo
c
c.... Compute psi1 and psi11
c
      call psi_I1 (mu,k1,k2,I1,psi1,psi11)
c
c.... Compute deviatoric Kirchhoff stress and map to tensor form 
c
      psi1 = psi1 + psi1
	tau = psi1*dev_bb
c
c.... Compute elasticity moduli and map to to tensor form
c
c.... Rank one update: -2/3 tau x g - 2/3 g x tau
      do i = 1,6
        taui = two3*tau(i)
        do j = 1,3
          aa(i,j) =  aa(i,j) - taui
          aa(j,i) =  aa(j,i) - taui
        enddo
      enddo
c                       ____
c.... Deviatoric term 2 psi1 [ I - 1/3 g x g ]
c
      psi1b  = psi1*trbb3
      psi1b2 = psi1b + psi1b
      psi1b3 = psi1b2/3.d0
c
      do i = 1,3
        aa(i  ,i  ) = aa(i  ,i  ) + psi1b2
        aa(i+3,i+3) = aa(i+3,i+3) + psi1b
        do j = 1,3
          aa(i ,j ) = aa(i ,j )   - psi1b3
        enddo
      enddo

c                       
c.... Deviatoric term  4 psi11 [ dev_bb x dev_bb ]
c
	psi11 = psi11 + psi11 + psi11 + psi11
     
      do i= 1,6
	 do j= 1,6
	   aa(i,j) = aa(i,j) + psi11*dev_bb(i)*dev_bb(j)  
	 enddo
	enddo

1000  end



	subroutine psi_I1 (mu,k1,k2,I1,psi1,psi11)
c==========================================================TCG.-10.07.2001
c   Computes derivatives of strain-energy function with respect to the first
c   modified invariant
c
c  Strain-energy function        Psi(I1) = Psi_nH(I1) + Psi_exp(I1)
c
c             neo-Hookean:       Psi_nH  = 0.5 k1 (I1-3)
c             exponential:       Psi_exp = 0.5 k1/k2 (Exp[k2(I1-3)^2]-1)
c
c input variables:
c       mu........neo-Hookean parameter
c       k1........Exponential parameter 1
c       k2........Exponential parameter 2
c       I1........First modified invarianat
c output variables: 
c       psi1......First derivative
c       psi11.....Second derivative
c used:
c       c1.........Constant
c       c2.........Constant
c
c=====================================================================

      IMPLICIT NONE

      real*8 k1, k2, mu, I1, psi1, psi11, c1, c2
      
      c1 = (I1-3.0d0)
      c2 = EXP(k2*c1*c1)
	
      psi1 = 0.5d0*mu + k1*c1*c2 
	psi11= (k1+2.0d0*k1*k2*c1*c1)*c2
	
      end


      subroutine dist_composite_11(d,f,jac,  tau,dd)
c==========================================================TCG.-10.11.2004
c  
c  distributed fiberreinforced composite material
c
c   INPUT:
c		d...........Material parameter vector
c			d(4)=k1			fiberparameter
c			d(5)=k2			fiberparameter
c			d(6)=kappa		distribution parameter
c		F(3,3)......Deformation gradient
c	    jac.........Jacobian determinant jac=det[F]
c
c   OUTPUT:
c		tau..........Kirchhoff stress	 tau = 2 dPsi/dg
c		dd...........Eulerian tangent	 dd = 2 dtau/dg
c
c   USED:
c		Ft(3,3)......Modified deformation gradient,  Ft=jac^(-1/3)*F
c		anl(3).......Unit fibre direction vector in lagrangian description
c		ane(3).......Unit fibre direction vector in eulerian description
c		he(3)........Structural tensor in eulerian description
c		b(6).........Left Cauchy-Green tensor
c		tau_f........Kirchhoff stress due to fibres
c		dd_f.........Eulerian tangent due to fibres
c		one3.........1/3
c		jm13.........jac**(-one3)
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

	logical aflag

	integer i,j,ii
      real*8 d(*), F(3,3), jac, Ft(3,3)
	real*8 tau(6), dd(6,6), b(6), k1, k2
	real*8 tau_f(6), dd_f(6,6), kappa, he(6)
	real*8 tau_m(6), dd_m(6,6), bt(6)
	real*8 anl(3), ane(3),one3, jm13, I4
	data one3  /0.33333333333333333/

      include 'eldata.h'
	include  'vecdata.h'


	jm13=jac**(-one3)

 	k1		= d(3)
	k2		= d(4)
	kappa	= d(5)
c initialize
 	dd   =0.0d0
 	tau  =0.0d0
	dd_m =0.0d0
	tau_m=0.0d0

c  load direction vector (unit)
      do ii= 1,2
c initialize
 	 dd_f =0.0d0
	 tau_f=0.0d0
	  do j= 1,3
         anl(j) = rvec(ii,j,n)
	  enddo

c compute modified elastic deformation gradient
       Ft=jm13*F
c compute modified left Cauchy Green strain
	 call leftCauchy_11 (Ft,bt)
 
c compute direction vector in eulerian description
       ane = 0.0d0
       do i= 1,3
	  do j= 1,3
	  ane(i)= ane(i) + Ft(i,j)*anl(j)
	  enddo
	 enddo
c restrict anisotropy to stretched fibres
	 I4 = ane(1)*ane(1) + ane(2)*ane(2) + ane(3)*ane(3)  
	 if (I4.lt.1.0d0)   then
        aflag =.false.
	 else
	  aflag=.true.
	 endif
	  aflag=.true.
c compute structural tensor in eulerian description
       call dist_struct_11 (aflag,kappa, bt,ane, he)

c compute Kirchhoff stress and eulerian moduli of distributed fibers
	  call dist_fibre_11 (k1, k2, he, tau_f, dd_f)
	  tau = tau + tau_f
	  dd = dd + dd_f
	enddo
c=REM.2================================================================
c add contribution to stress and elasticity due to matrix material if 
c needed 
	call leftCauchy_11 (F,b)
      call neoh_11(d,b,jac, tau_m,dd_m)
c======================================================================
c summation
	tau = tau + tau_m
	dd  = dd + dd_m

      end


      subroutine six_para_11(d,f,jac,  tau,dd)
c==========================================================TCG.-31.05.2006
c  
c  6 parameter artery wall model
c
c   INPUT:
c		d...........Material parameter vector
c			d(2)=mu			neoHookean parameter
c			d(3)=k1			fiberparameter media
c			d(4)=k2			fiberparameter media
c			d(5)=kappa		distribution parameter media
c			d(6)=k1			fiberparameter adventitia
c			d(7)=k2			fiberparameter adventitia
c		F(3,3)......Deformation gradient
c	    jac.........Jacobian determinant jac=det[F]
c
c   OUTPUT:
c		tau..........Kirchhoff stress	 tau = 2 dPsi/dg
c		dd...........Eulerian tangent	 dd = 2 dtau/dg
c
c   USED:
c		Ft(3,3)......Modified deformation gradient,  Ft=jac^(-1/3)*F
c		anl(3).......Unit fibre direction vector in lagrangian description
c		ane(3).......Unit fibre direction vector in eulerian description
c		he(3)........Structural tensor in eulerian description
c		b(6).........Left Cauchy-Green tensor
c		tau_f........Kirchhoff stress due to fibres
c		dd_f.........Eulerian tangent due to fibres
c		one3.........1/3
c		jm13.........jac**(-one3)
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

	logical aflag

	integer i,j,ii
      real*8 d(*), F(3,3), jac, Ft(3,3)
	real*8 tau(6), dd(6,6), b(6), k1, k2,k3,k4
	real*8 tau_f(6), dd_f(6,6), kappam, kappaa, he(6)
	real*8 tau_m(6), dd_m(6,6), bt(6)
	real*8 anl(3), ane(3),one3, jm13, I4
	data one3  /0.33333333333333333/

      include 'eldata.h'
	include  'vecdata.h'


	jm13=jac**(-one3)

 	k1		= d(3)
	k2		= d(4)
	kappam	= d(5)
 	k3		= d(6)
	k4		= d(7)
	kappaa	= d(8)
c initialize
 	dd   =0.0d0
 	tau  =0.0d0
	dd_m =0.0d0
	tau_m=0.0d0

C======= BEGIN MEDIA =============================

c  load direction vector (unit)
      do ii= 1,2
c initialize
 	 dd_f =0.0d0
	 tau_f=0.0d0
	  do j= 1,3
         anl(j) = rvec(ii,j,n)
	  enddo

c compute modified elastic deformation gradient
       Ft=jm13*F
c compute modified left Cauchy Green strain
	 call leftCauchy_11 (Ft,bt)
 
c compute direction vector in eulerian description
       ane = 0.0d0
       do i= 1,3
	  do j= 1,3
	  ane(i)= ane(i) + Ft(i,j)*anl(j)
	  enddo
	 enddo
c restrict anisotropy to stretched fibres
	 I4 = ane(1)*ane(1) + ane(2)*ane(2) + ane(3)*ane(3)  
	 if (I4.lt.1.0d0)   then
        aflag =.false.
	 else
	  aflag=.true.
	 endif
	  aflag=.true.
c compute structural tensor in eulerian description
       call dist_struct_11 (aflag,kappam, bt,ane, he)

c compute Kirchhoff stress and eulerian moduli of distributed fibers
	  call dist_fibre_11 (k1, k2, he, tau_f, dd_f)
	  tau = tau + tau_f
	  dd = dd + dd_f
	enddo

C======= END MEDIA =============================

C======= BEGIN ADVENTITIA =============================
c  load direction vector (unit)
      do ii= 3,4
c initialize
 	 dd_f =0.0d0
	 tau_f=0.0d0
	  do j= 1,3
         anl(j) = rvec(ii,j,n)
	  enddo

c compute modified elastic deformation gradient
       Ft=jm13*F
c compute modified left Cauchy Green strain
	 call leftCauchy_11 (Ft,bt)
 
c compute direction vector in eulerian description
       ane = 0.0d0
       do i= 1,3
	  do j= 1,3
	  ane(i)= ane(i) + Ft(i,j)*anl(j)
	  enddo
	 enddo
c restrict anisotropy to stretched fibres
	 I4 = ane(1)*ane(1) + ane(2)*ane(2) + ane(3)*ane(3)  
	 if (I4.lt.1.0d0)   then
        aflag =.false.
	 else
	  aflag=.true.
	 endif
	  aflag=.true.
c compute structural tensor in eulerian description
       call dist_struct_11 (aflag,kappaa, bt,ane, he)

c compute Kirchhoff stress and eulerian moduli of distributed fibers
	  call dist_fibre_11 (k3, k4, he, tau_f, dd_f)
	  tau = tau + tau_f
	  dd = dd + dd_f
	enddo

C======= END ADVENTITIA =============================



c=REM.2================================================================
c add contribution to stress and elasticity due to matrix material if 
c needed 
	call leftCauchy_11 (F,b)
      call neoh_11(d,b,jac, tau_m,dd_m)
c======================================================================
c summation
	tau = tau + tau_m
	dd  = dd + dd_m

      end


	subroutine dist_fibre_11(k1,k2, he, tau, dd)
c==========================================================TCG.-10.11.2004
c
c   compute the kirchhoff stress and eulerian tangent moduli 
c   of a transversally isotropic material based on distributed fibre-
c   reinforced composite.
c
c  Input:
c		k1.........First fiber parameter
c		k2.........Second fiber parameter
c         he.........Eulerian structure tensor he_ij 
c
c  Output:
c		tau........Kirchhoff stress
c		dd.........Eulerian tangent
c
c  Used:
c         I4.........Fourth invariant I4 = ae_i*ae_i
c         w4.........First derivative of the stored energy function
c                    with respect to I4
c         w44........Second derivative of the stored energy function
c                    with respect to I4
c         axa........(ae x ae)
c         dev_axa....dev (ae x ae)
c         tempi......parts of the material tangent
c
c=======================================================================
 	IMPLICIT none
      real*8  tau(6), he(6), dd(6,6), I4, w4, w44
	real*8 temp1(6,6), temp2(6,6), temp3(6,6), temp4(6,6),temp5(6,6)
	real*8 dev_axa(6), axa(6), one3, k1, k2
	integer i,j
	data one3  /0.33333333333333333/
c compute I4
	I4 = he(1) + he(2) + he(3)
c get w4 and w44
      call dist_Psi_I4_11 (k1,k2,I4, w4,w44)
c derive Kirchhoff stress
      axa = he
	dev_axa = axa
	call dev_euler_11(dev_axa)

	tau=2.0d0*w4*dev_axa
c derive spatial tangent
	temp1=0.0d0
      do i= 1,6
	 do j= i,6
	  temp5(i,j) = dev_axa(i)*dev_axa(j)
	 enddo
	 do j= i,3
	  temp1(i,j) = axa(i)
	 enddo
	enddo
	temp2=0.0d0
	temp3=0.0d0
      do i= 1,3
	 do j= i,6
	  temp2(i,j) = axa(j)
	 enddo
      do j= i,3
	  temp3(i,j) = 1.0d0
	enddo
	enddo
	temp4=0.0d0
	temp4(1,1)=-1.0d0
	temp4(2,2)=-1.0d0
	temp4(3,3)=-1.0d0
	temp4(4,4)=-0.5d0
	temp4(5,5)=-0.5d0
	temp4(6,6)=-0.5d0

      dd=0.0d0
	dd=(4.0d0*w4/3.0d0)*(one3*I4*temp3-temp1-temp2-I4*temp4)
	dd = dd + 4.0d0*w44*temp5 
c  symmetrize tangent
 	do i = 1,6
	 do j= i+1,6
	 dd(j,i) = dd(i,j)
	 enddo
	enddo

1000	end

      subroutine dist_struct_11 (aflag, kappa, bt,ae, he)
c==========================================================TCG.-15.02.2005
c
c   compute structural tensor in eulerian description 
c
c==========================================================================
      implicit none

	logical aflag
      real*8 kappa, bt(6), ae(3), he(6), omk, axa(6)

      omk     = 1.0d0-3.0*kappa

      call vec_prod_11 (ae,axa)
	if (aflag) then
	 he    = kappa*bt + omk*axa
	else
	 he    = kappa*bt 
	endif

	end


      subroutine dist_Psi_I4_11(k1,k2,I4, w4,w44)
c==========================================================TCG.-29.09.1999
c
c  derivatives with respect to the fourth invariants of the 
c  distributed fiber model
c
c     Wt^aniso = (k1/2*k2)*(exp(k2*(I4-tm2k)^2)-1)
c
c  Input:
c         k1.......First material parameter
c         k2.......Second material parameter
c         I4.......Fourth invariant
c    
c  Output:
c         w4.......First derivative of the stored-energy function 
c                  with respect to the fourt invariant     
c         w44......Second derivative of the stored-energy function 
c                  with respect to the fourt invariant     
c
c  Used:
c         I4m1.....I4-tm2k
c         temp.....exp(k2*I4m1*I4m1)
c
c===============================================================
 	IMPLICIT none
      real*8 w4, w44, I4, k1, k2 ,I4m1 , temp
      
c get material parameters

	I4m1= I4-1.0d0
	temp= exp(k2*I4m1*I4m1)

      w4  = k1*I4m1*temp    
      w44 = k1*(1.0d0+2.0d0*k2*I4m1*I4m1)*temp
      end

!--------------------------------------------------------------------------------------------CHRIS'S CONSTITUTIVE MODEL------------------

      subroutine dis_nonlin_vis_dam_mod(flag,d,f,jac,hn,
     1           nhv,tau,dd,dds)
c==================================================================CJM.-21/09/2022
c      
c     COLLAGEN FIBER DAMAGE MODEL SUBROUTINE for F E A P 
c      
c     PURPOSE
c        Main subroutine for collagen fiber damage model 
c      
c     INPUT variables
c        flag.........If there is a discontinuity in the element? .true. = yes       
c	   d............Material parameter vector   
c	   f(3,3).......Deformation gradient  
c	   jac..........Jacobian determinant jac=det[F]  
c        hn...........History variables associated with the gauss point
c        nhv..........Number of history variables at each gauss point      
c      
c     OUTPUT variables
c	   tau..........Kirchhoff stress
c        dd...........Eulerian tangent  
c        dds..........Eulerian tangent including softening (determined numerically)       
c      
c     USED variables      
c        tau_f........Kirchhoff stress due to fibers 
c        dd_f.........Eulerian tangent due to fibers
c        Ft...........Modified deformation gradient,  Ft=jac**(-1/3)*F
c        l_sph........Spherical integration order, then the number of the spherical integration point (two uses)
c        lint_sph.....Total number of spherical integration points
c        sg...........Array containing unit vector components for each sph. int. point     
c        w............weighting function of spherical design 4Pi/lint
c        nhvs.........Number of history variables in a fiber (spherical integration point)
c        m0...........Unit direction vector defining referential direction of a fiber
c        m............push-forward of m0, deformed direction of fiber
c        I4...........4th invariant m.m associated with a fiber
c        lam..........stretch of a fiber lam = sqrt(I4)
c        rho..........Density associated with a fiber
c        lmin.........Minimum straightening stetch of all CFPG complexes present within the fiber
c        b............Left Cauchy-Green strain Tensor     
c        tau_m........Kirchhoff stress due to Elastin 
c        dd_m.........Eulerian tangent due to Elastin
c      
c     Local notation:	11	-	1
c					22	-	2
c					33	-	3
c					12	-	4
c					23	-	5
c					13	-	6
c      
c=================================================================================      
      
      IMPLICIT NONE  
      
      logical flag
      integer i, j, nn, nhvs, nhv
      integer l_sph, lint_sph
      real*8 d(*), f(3,3), Ft(3,3), jac, jm13
      real*8 tau(6), dd(6,6), dds(6,6)
      real*8 tau_m(6), dd_m(6,6)
      real*8 tau_f(6), dd_f(6,6), dds_f(6,6)
      real*8 hn(nhv), sg(3,240), w
      real*8 m0(3), m(3)
      real*8 I4, lam, rho, lmin, b(6)
      real*8 one3, pi4
	data one3  /0.33333333333333333/
	data pi4   /12.56637061435917/      
      
      include 'Visco.h'
      
c...  initialize stress and stiffness
      
 	tau = 0.0d0    
      dd = 0.0d0
      dds = 0.0d0
      
c=====COLLAGEN CONTRIBUTION==========================================================================================    

c...  Initialize the Collagen Kirchoff stress and Eulerian tangent

      tau_f = 0.0d0
      dd_f = 0.0d0
      dds_f = 0.0d0
      
c...  Compute modified elastic deformation gradient
      
	jm13 = jac**(-one3)
      Ft = jm13*F      
      
c...  Define integration order based on user input
    
      l_sph = d(40)
      call int_sphere(l_sph,lint_sph,sg) !sg is array containing unit vector components for each sph. int. point
      
c...  Compute integration weight (associated with each sph. int. point)
      
      w = pi4/float(lint_sph) 
      nhvs=d(43)   ! length of hist. vector at sph. int. point   
      
c...  Loop over orientations (spherical integration points)
      
      nn = 1  ! pointer to the first collagen-related history element
      
      do l_sph=1,lint_sph        
      
         spno = l_sph 
          
c...  Define direction vector to spherical integration point
      
         do i=1,3
	      m0(i) = sg(i,l_sph)
         enddo          
          
c...  Compute push forward m0 (can be sped-up when principal stretches are available)
	
	   m = 0.0d0
	   do i= 1,3
	      do j= 1,3
	         m(i)= m(i)+Ft(i,j)*m0(j)
	      enddo
         enddo          
          
         I4 = m(1)*m(1)+m(2)*m(2)+m(3)*m(3)
         lam = sqrt(I4)          
          
c...  Compute Kirchhoff stress and eulerian moduli for a fiber bundle         
         
         rho = hn(nn+0)
         lmin = d(10)      
         
         if (lam.gt.lmin) then !FIBER HAS ENTERED TENSION
            call dis_nonlin_vis_dam_fib(flag,d,I4,hn(nn),nhvs,m,
     1             tau_f,dd_f,dds_f) 
         else
            tau_f = 0d0
            dd_f = 0d0
            dds_f = 0d0
         endif  
         
         tau = tau+tau_f*rho 
	   dd = dd+dd_f*rho   
         dds = dds+dds_f*rho 
          
         nn = nn+nhvs
         
      enddo    
          
      tau = tau*w !weight stress and eulerian moduli of fibers
      dd = dd*w      
      dds = dds*w    
      ! NOTE-can do this because the weight for each sph int point contribution to stress and stiffness is the same.    
          
c=====ELASTIN CONTRIBUTION===========================================================================================          
          
	dd_m = 0.0d0
	tau_m = 0.0d0 
      
c...  Determine left Cuachy green strain tensor for elastic deformation

      call leftCauchy_11(F,b)     
      
c...  Determine contribution to stress and elasticity due to matrix material      
      
      call neoh_11(d,b,jac,tau_m,dd_m) !determine Kirchoff stress and eulerian tangent based on this deformation 
      
c...  Summation of stresses
      
	tau = tau+tau_m
	dd = dd+dd_m    
      dds = dds+dd_m    
      
      end
      
      
      
      subroutine dis_nonlin_vis_dam_fib(flag,d,I4,hn,nhvs,ane,
     1           tau,dd,dds)
c==================================================================CJM.-21/09/2022
c      
c     COLLAGEN FIBER CONTRIBUTION SUBROUTINE for F E A P 
c   
c     PURPOSE
c        Determines the 3D contribution of stress and stiffness contribution
c        resulting from a collagen fiber in a particular direction "ane"      
c          
c     INPUT variables
c        flag.........If there is a discontinuity in the element? .true. = yes             
c        d............Material parameter vector      
c        I4...........Fourth invariant, stretch in deformed fibre direction squared
c        hn...........History variables corresponding to sph int point
c        nhvs.........Number of history variables per sph int point
c        ane..........Deformed fiber direction vector (not unit)
c      
c     OUTPUT variables
c        tau..........Kirchhoff stress
c        dd...........Eulerian tangent   
c        dds..........Numerically determined Eulerian tangent, which includes softening       
c  
c     USED variables      
c        w4...........First derivative of the stored energy function with respect to I4
c        w44..........Second derivative of the stored energy function with respect to I4
c        axa..........(ae x ae)
c        dev_axa......dev (ae x ae)
c        tempi........parts of the material tangent
c      
c=================================================================================      

 	IMPLICIT NONE  
      
      logical flag
	integer i,j, nhvs
	real*8 d(*),hn(nhvs),ane(3),dev_axa(6), axa(6)
      real*8 tau(6),dd(6,6),dds(6,6),I4,w4,w44,w44_s      
      real*8 temp1(6,6),temp2(6,6),temp3(6,6),temp4(6,6),temp5(6,6)
      real*8 PB, check
      real*8 one3
	data one3  /0.33333333333333333/    
      
      include 'Visco.h'
      
c-----Get w4 and w44-----------------------------------------------------------------------------------
           
      call dis_nonlin_vis_dam_psi(flag,d,I4,hn,nhvs,w4,w44,w44_s)
      
	axa(1) = ane(1)*ane(1)
	axa(2) = ane(2)*ane(2)
	axa(3) = ane(3)*ane(3)
      axa(4) = ane(1)*ane(2)
      axa(5) = ane(2)*ane(3)
      axa(6) = ane(1)*ane(3)

c-----Derive Kirchhoff stress--------------------------------------------------------------------------

	dev_axa = axa
	call dev_euler11(dev_axa)

	tau=2.0d0*w4*dev_axa
      
      !tau=2.0d0*I4*w4*dev_axa
      
c-----Derive spatial tangent---------------------------------------------------------------------------

	temp1=0.0d0

      do i= 1,6
         
          do j= i,6
	      temp5(i,j) = dev_axa(i)*dev_axa(j)
         enddo
	 
         do j= i,3
	      temp1(i,j) = axa(i)
         enddo
         
      enddo
      
	temp2=0.0d0
	temp3=0.0d0
      
      do i= 1,3
          
         do j= i,6
	      temp2(i,j) = axa(j)
         enddo
         
         do j= i,3
	      temp3(i,j) = 1.0d0
         enddo
         
      enddo
      
	temp4=0.0d0
	temp4(1,1)=-1.0d0
	temp4(2,2)=-1.0d0
	temp4(3,3)=-1.0d0
	temp4(4,4)=-0.5d0
	temp4(5,5)=-0.5d0
	temp4(6,6)=-0.5d0

      dd=0.0d0
	dd=(4.0d0*w4/3.0d0)*(one3*I4*temp3-temp1-temp2-I4*temp4)
	dd = dd+4.0d0*w44*temp5 
      
      dds=0.0d0
	dds=(4.0d0*w4/3.0d0)*(one3*I4*temp3-temp1-temp2-I4*temp4)
	dds = dds+4.0d0*w44_s*temp5       
      
c-----Symmetrize tangent-------------------------------------------------------------------------------

 	do i = 1,6
         do j= i+1,6
	      dd(j,i) = dd(i,j)
            dds(j,i) = dds(i,j)
	   enddo
      enddo
      
      if (dds(3,3).lt.0d0) then
         check = 1d0    
      endif
      
      
      end      



	subroutine dev_euler11(a)
      
c==================================================================TCG.-29.09.1999  
c    
c              COMPUTE EULERIAN DEVIATOR SUBROUTINE for F E A P
c
c     PURPOSE 
c
c        Compute the eulerian deviator of a second order tensor. Note that the
c        tensor is inputted and then the deviator is outputted at the same variable
c      
c        deva = a - 1/3 tra I      
c      
c     INPUT variables
c      
c        a            Second order tensor
c
c     OUTPUT variables
c
c        deva         Deviator of a 
c      
c     USED variables     
c
c        tra          Trace of a
c        one3         1/3
c
c=================================================================================         
      
 	IMPLICIT none
	real*8 a(6), tra, one3
	integer i
	data one3  /0.33333333333333333/
      
	tra = a(1)+a(2)+a(3)
	tra = one3*tra 
      
	do i= 1,3
         a(i) = a(i) - tra
      enddo
      
      end


      subroutine dis_nonlin_vis_dam_psi(flag,d,I4,hn,nhvs,w4,w44,w44_s)
c==================================================================CJM.-21/09/2022
c      
c     COLLAGEN STRAIN ENERGY DERIVATIVES SUBROUTINE for F E A P
c    
c     PURPOSE
c        1D determination of the stress and stiffness of a collagen fiber
c          
c     INPUT variables
c        flag.........If there is a discontinuityin the element? .true. = yes             
c        d............Material parameter vector      
c        I4...........Fourth invariant, stretch in deformed fibre direction squared
c        hn...........History variables corresponding to sph int point
c        nhvs.........Number of history variables per sph int point
c      
c     OUTPUT variables
c        w4...........First derivative of the stored energy function with respect to I4
c        w44..........Second derivative of the stored energy function with respect to I4        
c        w44_s........Numerically determined stiffness that includes softening           
c      
c     USED variables      
c        
c=================================================================================      
      
      IMPLICIT NONE
      
      logical flag
      integer i, j, k
      integer nhvs, nf, cut, max_it
      real*8 d(*), hn(nhvs)
      real*8 kf, eta_slid, Tf_hom, eta_rec
      real*8 A, B, p
      real*8 lam1, lam2, lam3, Rmid
      real*8 I4, lam, I4_p, lam_p, dt_p
      real*8 w4, w44, w44_s, w4_p, w44_p
      real*8 lams_f, lams_f_p, PB, PB_p, lams_fp
      real*8 m, c, m_p, c_p      
      real*8 Tf, lampg_f, temp_lampg_f
      real*8 Tf_p, lampg_f_p, temp_lampg_f_p, lampg_p, temp_lampg_p
      real*8 lambar_1, lambar_2, lambar_dam, lambar_nf
      real*8 lambar_1_p, lambar_2_p, lambar_dam_p, lambar_nf_p
      real*8 check
      
      include 'tdata.h'      
      include 'Visco.h'
      
c-----Definition of Parameters-------------------------------------------------------------------------
      
      kf = d(5)   
      eta_slid = d(6)
      Tf_hom = d(7)
      eta_rec = d(8)
      
      A = d(13)
      B = d(14)
      
      if (flag) then
         B = B*1e5       
      endif
      
      lam1 = d(10)
      lam2 = d(11)
      lam3 = d(12)       

      nf = d(41)
      cut = int(d(45))
      Rmid = d(44)      
      
      !max_it = int(d(16))
      max_it = 5
      
      lam = sqrt(I4) 
      
      p = d(15)
      I4_p = I4+p
      lam_p = sqrt(I4_p)
      
      if ((I4-hn(2)).ge.1e-6) then
        dt_p = dt*((I4_p-hn(2))/(I4-hn(2)))
      else
        dt_p = dt
      endif  
 
c-----lams corresponding to previous failure (Speed up Newton method)----------------------------------     
      
      if (hn(5).lt.Rmid) then
        lams_fp = lam1+(hn(5)*(lam2-lam1)*(lam3-lam1))**0.5d0
      elseif (hn(5).ge.Rmid) then
        lams_fp = lam3-((1d0-hn(5))*(lam3-lam2)*(lam3-lam1))**0.5d0
      endif 
 
      if (hn(5).gt.0d0) then
        check = 1d0      
      endif
      
c-----Loop through all the fibrils---------------------------------------------------------------------       

      w4 = 0d0
      w44 = 0d0
      PB = 0d0
      
      w4_p = 0d0
      w44_p = 0d0
      PB_p = 0d0
      
      w44_s = 0d0
      
      do j=1,nf !========================================================================

        !=== UPDATE PG STRETCH FOR PERTURBED DEFORMATION ===  
          
        if ((kf*(lam_p/(lams(j)*hn(5+j))-1)).ge.Tf_hom) then
          lampg_p=((lams(j)*(lams(j)*Tf_hom**2d0*dt_p**2d0+2d0*lams(j)
     1    *Tf_hom*dt_p**2d0*kf-2d0*lams(j)*Tf_hom*dt_p*eta_slid
     1     *hn(5+j)+lams(j)*dt_p**2d0*kf**2d0-2d0*lams(j)*dt_p*eta_slid
     1     *kf*hn(5+j)+4d0*lam_p*dt_p*eta_slid*kf+lams(j)*eta_slid**2d0
     1     *hn(5+j)**2d0))**0.5d0-Tf_hom*dt_p*lams(j)-dt_p*kf*lams(j)
     !     +eta_slid*lams(j)*hn(5+j))/(2d0*eta_slid*lams(j))
     !!     lampg_p=hn(5+j)+(dt_p/eta_slid)*((kf*((lam_p/(lams(j)
     !!1    *hn(5+j)))-1))-Tf_hom)
        elseif ((kf*(lam_p/(lams(j)*hn(5+j))-1)).lt.Tf_hom) then
          lampg_p=(dt_p+eta_rec*hn(5+j))/(dt_p+eta_rec)
        endif

        Tf_p=kf*(lam_p/(lams(j)*lampg_p)-1)
        lampg_f_p=A*(Tf_p-Tf_hom)+B

        !=== UPDATE PG STRETCH FOR CURRENT DEFORMATION ===
        
        if((kf*(lam/(lams(j)*hn(5+j))-1)).ge.Tf_hom) then
          hn(5+j)=((lams(j)*(lams(j)*Tf_hom**2d0*dt**2d0+2d0*lams(j)
     1    *Tf_hom*dt**2d0*kf-2d0*lams(j)*Tf_hom*dt*eta_slid*hn(5+j)
     1    +lams(j)*dt**2d0*kf**2d0-2d0*lams(j)*dt*eta_slid*kf*hn(5+j)
     1    +4d0*lam*dt*eta_slid*kf+lams(j)*eta_slid**2d0*hn(5+j)**2d0))
     1    **0.5d0-Tf_hom*dt*lams(j)-dt*kf*lams(j)+eta_slid*lams(j)
     1    *hn(5+j))/(2d0*eta_slid*lams(j))
     !!     hn(5+j)=hn(5+j)+(dt/eta_slid)*((kf*((lam/(lams(j)
     !!1    *hn(5+j)))-1))-Tf_hom)          
        elseif((kf*(lam/(lams(j)*hn(5+j))-1)).lt.Tf_hom) then
          hn(5+j)=(dt+eta_rec*hn(5+j))/(dt+eta_rec)
        endif

        Tf=kf*(lam/(lams(j)*hn(5+j))-1)
        lampg_f=A*(Tf-Tf_hom)+B
        
        !=== WORK WITHIN THE SEGMENTS ===
        
        !if ((j.gt.1).and.(.not.(flag))) then ! DAMAGE CAN DEVELOP
        if (j.gt.1) then ! DAMAGE CAN DEVELOP            
              
          if ((hn(5).gt.cgam(j-1)).and.(hn(5).lt.cgam(j))) then ! Partial failure at previous time in segment
          
            ! ---=== FIRST FOR THE CURRENT DEFORMATION ===--- 
              
            lambar_1 = hn(5+j-1)*lams(j-1)
            lambar_2 = hn(5+j)*lams(j)
              
            if ((Tf.gt.Tf_hom).and.(hn(5+j-1).gt.temp_lampg_f)
     1      .and.(hn(5+j).lt.lampg_f)) then ! Partial failure at current time in segment

              !=== DETERMINE STRAIGHTENING STRETCH CORRESPONDING TO FAILURE AT CURRENT TIME ===
            
              m = (hn(5+j)-hn(5+j-1))/(lams(j)-lams(j-1))
              c = hn(5+j)-m*lams(j)
            
              lams_f = min(lams(j),lams_fp) ! Initialise for Newton iteration
              k = 0
                
              do while (k.lt.max_it)
                lams_f=lams_f-((lams_f**2d0*(c+lams_f*m)**2d0*(c-B+A*kf
     1          +lams_f*m+A*Tf_hom-(A*kf*lam)/(lams_f*(c+lams_f*m))))/(c
     1          **2d0*lams_f**2d0*m+2d0*c*lams_f**3d0*m**2d0+A*kf*lam*c
     1          +lams_f**4d0*m**3d0+2d0*A*kf*lam*lams_f*m))
                k = k+1   
              enddo            
            
              !=== DETERMINE CORRESPONDING TOTAL PERCENTAGE BROKEN AT CURRENT TIME ===
              
              if ((j-1).le.cut) then
                PB = (lams_f-lam1)**2d0/((lam2-lam1)*(lam3-lam1))
              else
                PB = 1d0-(lam3-lams_f)**2d0/((lam3-lam2)*(lam3-lam1))
              endif
              
              !=== IS CURRENT DAMAGE MORE THAN PREVIOUS DAMAGE??? ===
              
              PB = max(PB,hn(5))
              
              !=== DETERMINE THE STRESS AND STIFFNESS ===
              
              lambar_dam = (lambar_1**2d0-((PB-cgam(j-1))*(lambar_1**2d0
     1        -lambar_2**2d0))/gam(j-1))**0.5d0
              
              if (lam.gt.lambar_1) then
                if (lam.ge.lambar_2) then
                  if (lam.gt.lambar_dam) then
                    w4=w4+(kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0*lambar_2
     1              **2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)-2d0*PB
     1              *lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lambar_1**2d0-((PB-cgam(j-1))
     1              *(lambar_1**2d0-lambar_2**2d0))/gam(j-1))**0.5d0)
     1              -(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     1              *cgam(j-1)+lambar_1**2d0*gam(j-1)+lambar_2**2d0
     1              *gam(j-1)-PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0
     1              *(lambar_1**2d0-lambar_2**2d0)*lambar_2)           
                  endif
                else
                  if (lam.gt.lambar_dam) then
                    w4=w4+(kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0*lambar_2
     1              **2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)-2d0*PB
     1              *lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lambar_1**2d0-((PB-cgam(j-1))
     1              *(lambar_1**2d0-lambar_2**2d0))/gam(j-1))**(0.5d0))
     1              -(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     1              *cgam(j-1)+lam**2d0*gam(j-1)+lambar_1**2d0*gam(j-1)
     1              -PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lam**2d0)**(0.5d0))
                    w44=-(kf*(PB-cgam(j-1)+(gam(j-1)*(lam**2d0-lambar_1
     1              **2d0))/(lambar_1**2d0-lambar_2**2d0)))/(4d0*lam
     1              **3d0)                  
                  endif                
                endif
              endif
              
            else ! No failure at current time in segment
              
              !=== WE USE THE PREVIOUS DAMAGE ===
              
              PB = hn(5)
              
              !=== DETERMINE THE STRESS AND STIFFNESS ===
              
              lambar_dam = (lambar_1**2d0-((PB-cgam(j-1))*(lambar_1**2d0
     1        -lambar_2**2d0))/gam(j-1))**0.5d0
              
              if (lam.gt.lambar_1) then
                if (lam.ge.lambar_2) then
                  if (lam.gt.lambar_dam) then
                    w4=w4+(kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0*lambar_2
     1              **2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)-2d0*PB
     1              *lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lambar_1**2d0-((PB-cgam(j-1))
     1              *(lambar_1**2d0-lambar_2**2d0))/gam(j-1))**0.5d0)
     1              -(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     1              *cgam(j-1)+lambar_1**2d0*gam(j-1)+lambar_2**2d0
     1              *gam(j-1)-PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0
     1              *(lambar_1**2d0-lambar_2**2d0)*lambar_2)           
                  endif
                else
                  if (lam.gt.lambar_dam) then
                    w4=w4+(kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0*lambar_2
     1              **2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)-2d0*PB
     1              *lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lambar_1**2d0-((PB-cgam(j-1))
     1              *(lambar_1**2d0-lambar_2**2d0))/gam(j-1))**(0.5d0))
     1              -(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     1              *cgam(j-1)+lam**2d0*gam(j-1)+lambar_1**2d0*gam(j-1)
     1              -PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lam**2d0)**(0.5d0))
                    w44=-(kf*(PB-cgam(j-1)+(gam(j-1)*(lam**2d0-lambar_1
     1              **2d0))/(lambar_1**2d0-lambar_2**2d0)))/(4d0*lam
     1              **3d0)                  
                  endif                
                endif
              endif            
              
            endif 
            
            ! ---=== SECOND FOR THE PERTURBED DEFORMATION ===--- 
            
            lambar_1_p = (temp_lampg_p*lams(j-1))
            lambar_2_p = (lampg_p*lams(j)) 
       
            if ((Tf_p.gt.Tf_hom).and.(temp_lampg_p.gt.temp_lampg_f_p)
     1         .and.(lampg_p.lt.lampg_f_p)) then ! Partial failure at current time in segment 

              !=== DETERMINE STRAIGHTENING STRETCH CORRESPONDING TO FAILURE AT CURRENT TIME ===  
                
              m_p = (lampg_p-temp_lampg_p)/(lams(j)-lams(j-1))
              c_p = lampg_p - m_p*lams(j)            
           
              lams_f_p = min(lams(j),lams_fp)       
              k = 0      
              do while (k.lt.max_it)    
                lams_f_p=lams_f_p-((lams_f_p**2d0*(c_p+lams_f_p*m_p)
     1          **2d0*(c_p-B+A*kf+lams_f_p*m_p+A*Tf_hom-(A*kf*lam_p)
     1          /(lams_f_p*(c_p+lams_f_p*m_p))))/(c_p**2d0*lams_f_p**2d0
     1          *m_p+2d0*c_p*lams_f_p**3*m_p**2d0+A*kf*lam_p*c_p
     1          +lams_f_p**4*m_p**3+2d0*A*kf*lam_p*lams_f_p*m_p))
                k = k+1
              enddo              
          
              !=== DETERMINE CORRESPONDING TOTAL PERCENTAGE BROKEN AT PERTURBED DEFORMATION === 
                
              if ((j-1).le.cut) then
                PB_p=(lams_f_p-lam1)**2d0/((lam2-lam1)*(lam3-lam1))
              else
                PB_p=1-(lam3-lams_f_p)**2d0/((lam3-lam2)*(lam3-lam1))
              endif 
          
              !=== IS PERTURBED DAMAGE MORE THAN PREVIOUS DAMAGE??? ===
              
              PB_p = max(PB_p,hn(5))         
       
              !=== DETERMINE THE PERTURBED STRESS ===
              
              lambar_dam_p=(lambar_1_p**2d0-((lambar_1_p**2d0-lambar_2_p
     1        **2d0)*(PB_p-cgam(j-1)))/gam(j-1))**(0.5d0)
              
              if (lam_p.gt.lambar_1_p) then
                if (lam_p.ge.lambar_2_p) then
                  if (lam_p.gt.lambar_dam_p) then
                    w4_p=w4_p+(kf*(2d0*lambar_1_p**2d0*cgam(j-1)-2d0
     1              *lambar_2_p**2d0*cgam(j-1)+2d0*lambar_1_p**2d0
     1              *gam(j-1)-2d0*PB_p*lambar_1_p**2d0+2d0*PB_p
     1              *lambar_2_p**2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p
     1              **2d0)*(lambar_1_p**2d0-((PB_p-cgam(j-1))
     1              *(lambar_1_p**2d0-lambar_2_p**2d0))/gam(j-1))
     1              **(0.5d0))-(kf*(lambar_1_p**2d0*cgam(j-1)-lambar_2_p
     1              **2d0*cgam(j-1)+lambar_1_p**2d0*gam(j-1)+lambar_2_p
     1              **2d0*gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p
     1              **2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p**2d0)
     1              *lambar_2_p)
                  endif
                else    
                  if (lam_p.gt.lambar_dam_p) then
                    w4_p=w4_p+(kf*(2d0*lambar_1_p**2d0*cgam(j-1)-2d0
     1              *lambar_2_p**2d0*cgam(j-1)+2d0*lambar_1_p**2d0
     1              *gam(j-1)-2d0*PB_p*lambar_1_p**2d0+2d0*PB_p
     1              *lambar_2_p**2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p
     1              **2d0)*(lambar_1_p**2d0-((PB_p-cgam(j-1))
     1              *(lambar_1_p**2d0-lambar_2_p**2d0))/gam(j-1))
     1              **(0.5d0))-(kf*(lambar_1_p**2d0*cgam(j-1)-lambar_2_p
     1              **2d0*cgam(j-1)+lam_p**2d0*gam(j-1)+lambar_1_p**2d0
     1              *gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p
     1              **2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p**2d0)
     1              *(lam_p**2d0)**(0.5d0))
                  endif
                endif   
              endif          
          
            else ! No failure at current time in segment     

              !=== WE USE THE PREVIOUS DAMAGE ===
          
              PB_p = hn(5) 
          
              !=== DETERMINE THE PERTURBED STRESS ===
              
              lambar_dam_p=(lambar_1_p**2d0-((lambar_1_p**2d0-lambar_2_p
     1        **2d0)*(PB_p-cgam(j-1)))/gam(j-1))**(0.5d0)
           
              if (lam_p.gt.lambar_1_p) then
                if (lam_p.ge.lambar_2_p) then
                  if (lam_p.gt.lambar_dam_p) then
                    w4_p=w4_p+(kf*(2d0*lambar_1_p**2d0*cgam(j-1)-2d0
     1              *lambar_2_p**2d0*cgam(j-1)+2d0*lambar_1_p**2d0
     1              *gam(j-1)-2d0*PB_p*lambar_1_p**2d0+2d0*PB_p
     1              *lambar_2_p**2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p
     1              **2d0)*(lambar_1_p**2d0-((PB_p-cgam(j-1))
     1              *(lambar_1_p**2d0-lambar_2_p**2d0))/gam(j-1))
     1              **(0.5d0))-(kf*(lambar_1_p**2d0*cgam(j-1)-lambar_2_p
     1              **2d0*cgam(j-1)+lambar_1_p**2d0*gam(j-1)+lambar_2_p
     1              **2d0*gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p
     1              **2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p**2d0)
     1              *lambar_2_p)
                  endif
                else    
                  if (lam_p.gt.lambar_dam_p) then
                    w4_p=w4_p+(kf*(2d0*lambar_1_p**2d0*cgam(j-1)-2d0
     1              *lambar_2_p**2d0*cgam(j-1)+2d0*lambar_1_p**2d0
     1              *gam(j-1)-2d0*PB_p*lambar_1_p**2d0+2d0*PB_p
     1              *lambar_2_p**2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p
     1              **2d0)*(lambar_1_p**2d0-((PB_p-cgam(j-1))
     1              *(lambar_1_p**2d0-lambar_2_p**2d0))/gam(j-1))
     1              **(0.5d0))-(kf*(lambar_1_p**2d0*cgam(j-1)-lambar_2_p
     1              **2d0*cgam(j-1)+lam_p**2d0*gam(j-1)+lambar_1_p**2d0
     1              *gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p
     1              **2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p**2d0)
     1              *(lam_p**2d0)**(0.5d0))
                  endif
                endif   
              endif

            endif            

          elseif (hn(5).le.cgam(j-1)) then ! No failure at previous time in segment

            ! ---=== FIRST FOR THE CURRENT DEFORMATION ===---   
              
            lambar_1 = hn(5+j-1)*lams(j-1)
            lambar_2 = hn(5+j)*lams(j)            
            
            if ((Tf.gt.Tf_hom).and.(hn(5+j-1).gt.temp_lampg_f)
     1      .and.(hn(5+j).lt.lampg_f)) then ! Partial failure at current time in segment

              !=== DETERMINE STRAIGHTENING STRETCH CORRESPONDING TO FAILURE AT CURRENT TIME ===
            
              m = (hn(5+j)-hn(5+j-1))/(lams(j)-lams(j-1))    !!!!!!!!!!!----------------------------------------------------------------PUT A STOP MARK HERE, WHEN DAMAGE FIRST STARTS-------------------
              c = hn(5+j)-m*lams(j)
            
              lams_f = min(lams(j),lams_fp) ! Initialise for Newton iteration
              k = 0
                
              do while (k.lt.max_it)
                lams_f=lams_f-((lams_f**2d0*(c+lams_f*m)**2d0*(c-B+A*kf
     1          +lams_f*m+A*Tf_hom-(A*kf*lam)/(lams_f*(c+lams_f*m))))/(c
     1          **2d0*lams_f**2d0*m+2d0*c*lams_f**3d0*m**2d0+A*kf*lam*c
     1          +lams_f**4d0*m**3d0+2d0*A*kf*lam*lams_f*m))
                k = k+1     
              enddo            
            
              !=== DETERMINE CORRESPONDING TOTAL PERCENTAGE BROKEN AT CURRENT TIME ===
              
              if ((j-1).le.cut) then
                PB = (lams_f-lam1)**2d0/((lam2-lam1)*(lam3-lam1))
              else
                PB = 1d0-(lam3-lams_f)**2d0/((lam3-lam2)*(lam3-lam1))
              endif
              
              !=== DETERMINE THE STRESS AND STIFFNESS ===
              
              lambar_dam = (lambar_1**2d0-((PB-cgam(j-1))*(lambar_1**2d0
     1        -lambar_2**2d0))/gam(j-1))**0.5d0
              
              if (lam.gt.lambar_1) then
                if (lam.ge.lambar_2) then
                  if (lam.gt.lambar_dam) then
                    w4=w4+(kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0*lambar_2
     1              **2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)-2d0*PB
     1              *lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lambar_1**2d0-((PB-cgam(j-1))
     1              *(lambar_1**2d0-lambar_2**2d0))/gam(j-1))**0.5d0)
     1              -(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     1              *cgam(j-1)+lambar_1**2d0*gam(j-1)+lambar_2**2d0
     1              *gam(j-1)-PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0
     1              *(lambar_1**2d0-lambar_2**2d0)*lambar_2)           
                  endif
                else
                  if (lam.gt.lambar_dam) then
                    w4=w4+(kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0*lambar_2
     1              **2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)-2d0*PB
     1              *lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lambar_1**2d0-((PB-cgam(j-1))
     1              *(lambar_1**2d0-lambar_2**2d0))/gam(j-1))**(0.5d0))
     1              -(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     1              *cgam(j-1)+lam**2d0*gam(j-1)+lambar_1**2d0*gam(j-1)
     1              -PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0*(lambar_1
     1              **2d0-lambar_2**2d0)*(lam**2d0)**(0.5d0))
                    w44=-(kf*(PB-cgam(j-1)+(gam(j-1)*(lam**2d0-lambar_1
     1              **2d0))/(lambar_1**2d0-lambar_2**2d0)))/(4d0*lam
     1              **3d0)                  
                  endif                
                endif
              endif
              
            else ! No failure at current time in segment
              
              !=== DETERMINE THE STRESS AND STIFFNESS ===
              
              if (lam.gt.lambar_1) then
                if (lam.ge.lambar_2) then
                  w4=w4+-(kf*(lambar_1-lambar_2)*(cgam(j-1)*lambar_1
     1            +cgam(j-1)*lambar_2+gam(j-1)*lambar_1-PB*lambar_1-PB
     1            *lambar_2))/(2d0*lambar_1**2d0*lambar_2+2d0*lambar_2
     1            **2d0*lambar_1)
                else
                  w4=w4+(kf*(lam-lambar_1)*(lambar_1**2d0*cgam(j-1)
     1            -lambar_2**2d0*cgam(j-1)+lambar_1**2d0*gam(j-1)-PB
     1            *lambar_1**2d0+PB*lambar_2**2d0-gam(j-1)*lam
     1            *lambar_1))/(2d0*(lambar_1**2d0-lambar_2**2d0)*lam
     1            *lambar_1)
                  w44=-(kf*(PB-cgam(j-1)+(gam(j-1)*(lam**2d0-lambar_1
     1            **2d0))/(lambar_1**2d0-lambar_2**2d0)))/(4d0*lam**3d0)
                endif
              endif              
              
            endif
            
            ! ---=== SECOND FOR THE PERTURBED DEFORMATION ===---
            
            lambar_1_p = (temp_lampg_p*lams(j-1))
            lambar_2_p = (lampg_p*lams(j)) 
       
            if ((Tf_p.gt.Tf_hom).and.(temp_lampg_p.gt.temp_lampg_f_p)
     1         .and.(lampg_p.lt.lampg_f_p)) then ! Partial failure at current time in segment 

              !=== DETERMINE STRAIGHTENING STRETCH CORRESPONDING TO FAILURE AT CURRENT TIME ===  
                
              m_p = (lampg_p-temp_lampg_p)/(lams(j)-lams(j-1))
              c_p = lampg_p - m_p*lams(j)            
           
              lams_f_p = min(lams(j),lams_fp)       
              k = 0      
              do while (k.lt.max_it)    
                lams_f_p=lams_f_p-((lams_f_p**2d0*(c_p+lams_f_p*m_p)
     1          **2d0*(c_p-B+A*kf+lams_f_p*m_p+A*Tf_hom-(A*kf*lam_p)
     1          /(lams_f_p*(c_p+lams_f_p*m_p))))/(c_p**2d0*lams_f_p**2d0
     1          *m_p+2d0*c_p*lams_f_p**3*m_p**2d0+A*kf*lam_p*c_p
     1          +lams_f_p**4*m_p**3+2d0*A*kf*lam_p*lams_f_p*m_p))
                k = k+1
              enddo              
          
              !=== DETERMINE CORRESPONDING TOTAL PERCENTAGE BROKEN AT PERTURBED DEFORMATION === 
                
              if ((j-1).le.cut) then
                PB_p=(lams_f_p-lam1)**2d0/((lam2-lam1)*(lam3-lam1))
              else
                PB_p=1-(lam3-lams_f_p)**2d0/((lam3-lam2)*(lam3-lam1))
              endif        
       
              !=== DETERMINE THE PERTURBED STRESS ===
              
              lambar_dam_p=(lambar_1_p**2d0-((lambar_1_p**2d0-lambar_2_p
     1        **2d0)*(PB_p-cgam(j-1)))/gam(j-1))**(0.5d0)
              
              if (lam_p.gt.lambar_1_p) then
                if (lam_p.ge.lambar_2_p) then
                  if (lam_p.gt.lambar_dam_p) then
                    w4_p=w4_p+(kf*(2d0*lambar_1_p**2d0*cgam(j-1)-2d0
     1              *lambar_2_p**2d0*cgam(j-1)+2d0*lambar_1_p**2d0
     1              *gam(j-1)-2d0*PB_p*lambar_1_p**2d0+2d0*PB_p
     1              *lambar_2_p**2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p
     1              **2d0)*(lambar_1_p**2d0-((PB_p-cgam(j-1))
     1              *(lambar_1_p**2d0-lambar_2_p**2d0))/gam(j-1))
     1              **(0.5d0))-(kf*(lambar_1_p**2d0*cgam(j-1)-lambar_2_p
     1              **2d0*cgam(j-1)+lambar_1_p**2d0*gam(j-1)+lambar_2_p
     1              **2d0*gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p
     1              **2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p**2d0)
     1              *lambar_2_p)
                  endif
                else    
                  if (lam_p.gt.lambar_dam_p) then
                    w4_p=w4_p+(kf*(2d0*lambar_1_p**2d0*cgam(j-1)-2d0
     1              *lambar_2_p**2d0*cgam(j-1)+2d0*lambar_1_p**2d0
     1              *gam(j-1)-2d0*PB_p*lambar_1_p**2d0+2d0*PB_p
     1              *lambar_2_p**2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p
     1              **2d0)*(lambar_1_p**2d0-((PB_p-cgam(j-1))
     1              *(lambar_1_p**2d0-lambar_2_p**2d0))/gam(j-1))
     1              **(0.5d0))-(kf*(lambar_1_p**2d0*cgam(j-1)-lambar_2_p
     1              **2d0*cgam(j-1)+lam_p**2d0*gam(j-1)+lambar_1_p**2d0
     1              *gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p
     1              **2d0))/(2d0*(lambar_1_p**2d0-lambar_2_p**2d0)
     1              *(lam_p**2d0)**(0.5d0))
                  endif
                endif   
              endif          
          
            else ! No failure at current time in segment     
                         
              !=== DETERMINE THE PERTURBED STRESS ===

              if (lam_p.gt.lambar_1_p) then
                if (lam_p.ge.lambar_2_p) then
                  w4_p=w4_p+-(kf*(lambar_1_p-lambar_2_p)*(cgam(j-1)
     1            *lambar_1_p+cgam(j-1)*lambar_2_p+gam(j-1)*lambar_1_p
     1            -PB_p*lambar_1_p-PB_p*lambar_2_p))/(2d0*lambar_1_p
     1            **2d0*lambar_2_p+2d0*lambar_2_p**2d0*lambar_1_p)
                else    
                  w4_p=w4_p+(kf*(lam_p-lambar_1_p)*(lambar_1_p**2d0
     1            *cgam(j-1)-lambar_2_p**2d0*cgam(j-1)+lambar_1_p**2d0
     1            *gam(j-1)-PB_p*lambar_1_p**2d0+PB_p*lambar_2_p**2d0
     1            -gam(j-1)*lam_p*lambar_1_p))/(2d0*(lambar_1_p**2d0
     1            -lambar_2_p**2d0)*lam_p*lambar_1_p)
                endif   
              endif

            endif             
            
          endif
      
          !=== INTEGRALS RELATING TO FULL RECRUITMENT ===
        
          if (j.eq.nf) then
             
            ! ---=== FIRST FOR THE CURRENT DEFORMATION ===---  
              
            lambar_nf = hn(5+j)*lams(j)  
              
            if ((hn(5).eq.1d0)
     1      .or.((Tf.gt.Tf_hom).and.(hn(5+j).ge.lampg_f))) then ! If complete failure previously or complete failure currently
              PB=1d0    
            endif    
            
            if (lam.ge.lambar_nf) then 
              w4 = w4+(kf/2d0-(PB*kf)/2d0)/lambar_nf
     1        -(kf/2d0-(PB*kf)/2d0)/lam 
              w44 = -(kf*(PB-1))/(4d0*lam**3d0)
            endif
            
            ! ---=== SECOND FOR THE PERTURBED DEFORMATION ===---
            
            lambar_nf_p = lampg_p*lams(j)

            if ((hn(5).eq.1)
     1      .or.((Tf_p.gt.Tf_hom).and.(lampg_p.ge.lampg_f))) then
              PB_p = 1d0
            endif      
      
            if (lam_p.ge.lambar_nf_p) then  
              w4_p = w4_p + (kf/2d0-(PB_p*kf)/2d0)/lambar_nf_p
     1        -(kf/2d0-(PB_p*kf)/2d0)/lam_p
            endif    
            
          endif
          
     !!   elseif ((j.gt.1).and.(flag)) then  ! DAMAGE CAN NO LONGER DEVELOP  
     !!     
     !!     ! ---=== FIRST FOR THE CURRENT DEFORMATION ===---   
     !!       
     !!     lambar_1 = (hn(5+j-1)*lams(j-1))
     !!     lambar_2 = (hn(5+j)*lams(j))   
     !!
     !!     if ((hn(5).gt.cgam(j-1)).and.(hn(5).lt.cgam(j))) then ! Partial failure at previous update in segment     
     !!
     !!       PB = hn(5) 
     !!   
     !!       lambar_dam = (lambar_1**2d0 - ((lambar_1**2d0 - lambar_2
     !!1      **2d0)*(PB - cgam(j-1)))/gam(j-1))**(1/2d0)
     !!      
     !!       if (lam.gt.lambar_1) then
     !!         if (lam.ge.lambar_2) then
     !!           if (lam.gt.lambar_dam) then
     !!             w4 = w4 + (kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0
     !!1            *lambar_2**2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)
     !!1            -2d0*PB*lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0
     !!1            *(lambar_1**2d0-lambar_2**2d0)*(lambar_1**2d0-((PB
     !!1            -cgam(j-1))*(lambar_1**2d0-lambar_2**2d0))/gam(j-1))
     !!1            **(1/2d0))-(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     !!1            *cgam(j-1)+lambar_1**2d0*gam(j-1)+lambar_2**2d0
     !!1            *gam(j-1)-PB*lambar_1**2d0+PB*lambar_2**2d0))/(2d0
     !!1            *(lambar_1**2d0-lambar_2**2d0)*lambar_2)
     !!           endif
     !!         else    
     !!           if (lam.gt.lambar_dam) then
     !!             w4 = w4 + (kf*(2d0*lambar_1**2d0*cgam(j-1)-2d0
     !!1            *lambar_2**2d0*cgam(j-1)+2d0*lambar_1**2d0*gam(j-1)
     !!1            -2d0*PB*lambar_1**2d0+2d0*PB*lambar_2**2d0))/(2d0
     !!1            *(lambar_1**2d0-lambar_2**2d0)*(lambar_1**2d0-((PB
     !!1            -cgam(j-1))*(lambar_1**2d0-lambar_2**2d0))/gam(j-1))
     !!1            **(1/2d0))-(kf*(lambar_1**2d0*cgam(j-1)-lambar_2**2d0
     !!1            *cgam(j-1)+lam**2d0*gam(j-1)+lambar_1**2d0*gam(j-1)-PB
     !!1            *lambar_1**2d0+PB*lambar_2**2d0))/(2d0*(lambar_1**2d0
     !!1            -lambar_2**2d0)*(lam**2d0)**(1/2d0))
     !!             w44 = -(kf*(PB-cgam(j-1)+(gam(j-1)*(lam**2d0-lambar_1
     !!1            **2d0))/(lambar_1**2d0-lambar_2**2d0)))/(4d0*lam**3d0)
     !!           endif
     !!         endif   
     !!       endif       
     !!  
     !!     elseif (hn(5).le.cgam(j-1)) then ! No failure at previous update in segment  
     !!  
     !!       if (lam.gt.lambar_1) then
     !!         if (lam.ge.lambar_2) then
     !!           w4=w4+-(kf*(lambar_1-lambar_2)*(cgam(j-1)*lambar_1
     !!1          +cgam(j-1)*lambar_2+gam(j-1)*lambar_1-PB*lambar_1-PB
     !!1          *lambar_2))/(2d0*lambar_1**2d0*lambar_2+2d0*lambar_2
     !!1          **2d0*lambar_1)
     !!         else
     !!           w4=w4+(kf*(lam-lambar_1)*(lambar_1**2d0*cgam(j-1)
     !!1          -lambar_2**2d0*cgam(j-1)+lambar_1**2d0*gam(j-1)-PB
     !!1          *lambar_1**2d0+PB*lambar_2**2d0-gam(j-1)*lam
     !!1          *lambar_1))/(2d0*(lambar_1**2d0-lambar_2**2d0)*lam
     !!1          *lambar_1)
     !!           w44=-(kf*(PB-cgam(j-1)+(gam(j-1)*(lam**2d0-lambar_1
     !!1          **2d0))/(lambar_1**2d0-lambar_2**2d0)))/(4d0*lam**3d0)
     !!         endif
     !!       endif              
     !!       
     !!     endif  
     !!
     !!     if (j.eq.nf) then
     !!  
     !!       lambar_nf = hn(5+j)*lams(j)   
     !!  
     !!       if (lam.ge.lambar_nf) then 
     !!         w4 = w4+(kf/2d0-(PB*kf)/2d0)/lambar_nf
     !!1        -(kf/2d0-(PB*kf)/2d0)/lam 
     !!         w44 = -(kf*(PB-1))/(4d0*lam**3d0)
     !!       endif
     !! 
     !!     endif      
              
        endif

        !=== STORE FOR USE IN SEGMENTS ===

        temp_lampg_f = lampg_f

        temp_lampg_f_p = lampg_f_p
        temp_lampg_p = lampg_p
        
      enddo  !============================================================================ 
      
c-----Determine the stiffness accounting for softening also------------------------------------      
      
      if (w4.lt.0d0) then
        w4 = 0d0
        check = 1d0
      endif
      
      if (w4_p.gt.0d0) then
        w44_s = (w4_p-w4)/d(15)
      endif
      
      if (w44_s.lt.0) then
         check = 1d0
      endif
      
      if ((w4/lam+2*lam*w44_s).lt.0) then
         check = 1d0
      endif
            
      if ((2*w4+2*lam**2*w44_s).lt.0) then
         check = 1d0
      endif
      
      if (isnan(w4).or.isnan(w4_p)) then
         check = 1d0
      endif
      
      
c-----Update the fiber level history-----------------------------------------------------------   
      
      hn(2) = I4
      !hn(3) = w4
      hn(4) = w44
      hn(5) = PB             
      
      end    