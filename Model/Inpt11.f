
      subroutine inpt11 (d)
c-------------------------------------------------------------------
c
c  Input  Subroutine for FEAP
c
c   vector d(200) specification, see elmt11.f
c
c----------------------------------------------------------
c bulk material models
c   imat=1   linear elastic material (St. Venant Kirchhoff)
c   imat=2   neo Hookean model
c   imat=3   fiber composite model with neoHookean matrix and HGO fiber model
c   imat=4   modified Delfino model
c   imat=5   fiber composite model with neoHookean matrix and GOH fiber model
c   imat=6   six parameter vessel wall model
c   imat=7   small strain Mises plasticity model (from FEAP)
c   imat=8   finite strain general pasticity model (from FEAP)
c   imat=9   Discrete nonlinear viscoelastic fiber model with damage      
c-----------------------------------------------------------
c cohesive zone material models (traction seperation laws)
c  imatd=1    Isotropic Traction Model  
c  imatd=2    transversely isotropic Traction model
c  imatd=3    isotropic viscoelastic Traction model
c  imatd=4    smooth Traction model
c  imate=6    Exponential normalised anisotropic Traction model      
c-------------------------------------------------------------------
c Declare variable types
c
c  errck   ERROR check during FEAP input-subroutine 'dinput'
c  imat    Continuum material model
c  imatd   Discrete material model
c  FEF     FE formulation flag
c  ior     Unit number of input-file
c  iow     Unit number of output-file
c  nhv     Length of history vector per Gauss-point of bulk material
c  nhvd     Length of history vector per Gauss-point of bulk material
c  d       Vector of material parameters
c  l       order of bulk integration
c  ld      order of disc. integration
c-------------------------------------------------------------------

      implicit none
	logical errck
      common /errchk/ errck
      integer imat, imatd, nhv, sio, no
      integer i, j
      real*8  d(*),s(3,240)
      real*8  den1, den2, Rmid, R, gam_1, gam_2, R_step
      integer count, left, right   
      logical :: exist
      integer  stat  
      integer  chk(400,3), chk1, chk2, chk3
c Declare array types
      include 'iofile.h'
      include 'pvecdata.h'   
      include 'Visco.h'

c write element description
       write(iow,203) 
c--------------------------------------------------------------------
c
c    Record 1. (rho, mi)
c
c  rho     d(50) Density 
c  mi      d(51) Mass interpolator 
c                          (0) Diagonal
c                          (1) Consistent
c  st      d(52) symmetrize element tangent
c                          (0) No
c                          (1) Yes 
c--------------------------------------------------------------------
c read (rho, mi)
        call dinput(d(50),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (rho, mi, st)'
         write(iow,*) 'ERROR --> cannot read (rho, mi, st)'
	   goto 9999
	  endif
c write (rho, mi, st)
       write(iow,1001) d(50), int(d(51)), int(d(52))

c======= BEGIN CONTINUUM MATERIAL DATA ============================

c
c Read record 2 (imat,l)
c
      call dinput(d(70),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (imat,l)'
         write(iow,*) 'ERROR --> cannot read (imat,l)'
	   goto 9999
	  endif
       imat = int(d(70))
c check and write imat
	  if (imat.EQ.1) then
	   write (iow,100)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.2) then
	   write (iow,101)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.3) then
	   write (iow,102)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.4) then
	   write (iow,103)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.5) then
	   write (iow,104)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.6) then
	   write (iow,105)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.7) then
	   write (iow,106)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.8) then
	   write (iow,107)
	   write (iow,99) int(d(71))
	  elseif (imat.EQ.9) then
	   write (iow,108)
	   write (iow,99) int(d(71))         
	  else
         write(*,*)   'ERROR --> (imat,l) are not a valid numbers'
         write(iow,*) 'ERROR --> (imat,l) are not a valid numbers'
	   goto 9999
	  endif

	  if(imat.eq.1) then
c=(imat=1)============== Linear isotropic elasticity ===============71
c
c  Read record 3   (E, nu)
c--------------------------------------------------------------------
c read (E, nu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (E, nu)'
         write(iow,*) 'ERROR --> cannot read (E, nu)'
	   goto 9999
	  endif
c write (E, nu)
        write(iow,1003) d(1),d(2)
c define length of history array at each gauss point
c nhv
	  d(74) = dfloat(0)
	 elseif (imat.eq.2) then
c=(imat=2)============== neo Hookean Model=========================71
c
c Read record 3    (K, mu)
c--------------------------------------------------------------------
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1004) d(1),d(2)
c define length of history array at each gauss point
c nhv
	  d(74) = dfloat(0)
	 elseif (imat.eq.3) then
c=(imat=3)============== Fiber model ===============================71
c
c Read record 3    (K, mu)
c--------------------------------------------------------------------
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1004) d(1),d(2)
c
c Read record 4    (k1,k2)
c--------------------------------------------------------------------
c read (k1,k2)
        call dinput(d(3),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (k1, k2)'
         write(iow,*) 'ERROR --> cannot read (k1, k2)'
	   goto 9999
	  endif
c write (k1,k2)
        write(iow,1005) d(3),d(4)
c define length of history array at each gauss point
c nhv
	  d(74) = dfloat(0)
	 elseif (imat.eq.4) then
c=(imat=4)=== Modified Delfino model ===============================71
c
c Read record 3    (K, mu)
c--------------------------------------------------------------------
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1006) d(1),d(2)
c
c Read record 4    (k1,k2)
c--------------------------------------------------------------------
c read (k1,k2)
        call dinput(d(3),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (k1, k2)'
         write(iow,*) 'ERROR --> cannot read (k1, k2)'
	   goto 9999
	  endif
c write (k1,k2)
        write(iow,1007) d(3),d(4)
c define length of history array at each gauss point
c nhv
	  d(74) = dfloat(0)
	 elseif (imat.eq.5) then
c=(imat=5)============== distributed Fiber model ===============================71
c
c Read record 3    (K, mu)
c--------------------------------------------------------------------
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1008) d(1),d(2)
c
c Read record 4    (k1,k2,kappa)
c--------------------------------------------------------------------
c read (k1,k2, kappa)
        call dinput(d(3),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (k1, k2, kappa)'
         write(iow,*) 'ERROR --> cannot read (k1, k2, kappa)'
	   goto 9999
	  endif
c write (k1,k2, kappa)
        write(iow,1009) d(3),d(4), d(5)
c define length of history array at each gauss point
cc redine parameters k1 and k2
c        d(3) = d(3)/((1.0d0-(d(5)/2.0d0))**2.0)
c        d(4) = d(4)/((1.0d0-(d(5)/2.0d0))**2.0)
c nhv
	  d(74) = dfloat(0)
	 elseif (imat.eq.6) then
c=(imat=6)======== 7 parameter modelfor the arterial wall===========71
c
c Read record 3    (K, mu)
c--------------------------------------------------------------------
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1008) d(1),d(2)
c
c Read record 4    (k1,k2,kappam)
c--------------------------------------------------------------------
c read (k1,k2, kappa)
        call dinput(d(3),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (k1, k2, kappam)'
         write(iow,*) 'ERROR --> cannot read (k1, k2, kappam)'
	   goto 9999
	  endif
c write (k1,k2, kappa)
        write(iow,1009) d(3),d(4), d(5)
c Read record 5    (k3,k4,kappaa)
c--------------------------------------------------------------------
c read (k1,k2, kappa)
        call dinput(d(6),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (k3, k4,kappaa)'
         write(iow,*) 'ERROR --> cannot read (k3, k4,kappaa)'
	   goto 9999
	  endif
c write (k1,k2, kappa)
        write(iow,1010) d(6),d(7),d(8)
c define length of history array at each gauss point
c nhv
	  d(74) = dfloat(0)
      elseif (imat.eq.7) then
c=(imat=7)======== small strain von Mises elastoplasicity ===========71
c
c Read record 3    (K, mu)
c       d(1)    -  Elastic Bulk  modulus
c       d(2)    -  Elastic Shear modulus
c--------------------------------------------------------------------
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1008) d(1),d(2)
c
c Read record 4    (Y0, Yi, delta)
c--------------------------------------------------------------------
c read (Y0, Yi, delta)
c       d(41)    -  Yield stress (ep = 0)(Y0)          
c       d(42)    -  Yield stress (ep =00)(Yi)          
c       d(43)    -  Delta value  (ep = 0)(delta)
        call dinput(d(41),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Y0, Yi, delta)'
         write(iow,*) 'ERROR --> cannot read (Y0, Yi, delta)'
	   goto 9999
	  endif
c write (Y0, Yi, delta)
        write(iow,1011) d(41), d(42), d(43)

c Read record 5    (H_iso, Hk_in)
c--------------------------------------------------------------------
c read (H_iso, H_kin)
c       d(44)    -  Isotropic Hardening Modulus - H_iso 
c       d(45)    -  Kinematic Hardening Modulus - H_kin 
        call dinput(d(44),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (H_iso, H_kin)'
         write(iow,*) 'ERROR --> cannot read (H_iso, H_kin)'
	   goto 9999
	  endif
c write (Y0, Yi, delta)
        write(iow,1012) d(44), d(45)

c define length of history array at each Gauss point: nhv 
	  d(74) = dfloat(8)
      elseif (imat.eq.8) then
c=(imat=8)======== finite strain plasticity model ===========71
c
c Read record 3    (K, mu)
c       d(1)    -  Elastic Bulk  modulus
c       d(2)    -  Elastic Shear modulus
c--------------------------------------------------------------------
c
        d(46)= 1 ! von Mises plasticity (see plasfd.f)

c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1008) d(1),d(2)
c
c Read record 4    (Y0, Yi, delta)
c--------------------------------------------------------------------
c read (Y0, Yi, delta)
c       d(41)    -  Yield stress (ep = 0)(Y0)           [x sqrt(2/3)]
c       d(42)    -  Yield stress (ep =00)(Yi)           [x sqrt(2/3)]
c       d(43)    -  Delta value  (ep = 0)(delta)
        call dinput(d(41),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Y0, Yi, delta)'
         write(iow,*) 'ERROR --> cannot read (Y0, Yi, delta)'
	   goto 9999
        endif
c write (Y0, Yi, delta)
        write(iow,1013) d(41), d(42), d(43)
        d(41) = d(41)*sqrt(2.0/3.0)    ! change to fit internal FEAP definition
        d(42) = d(42)*sqrt(2.0/3.0)    ! change to fit internal FEAP definition

c Read record 5    (H_iso, Hk_in)
c--------------------------------------------------------------------
c read (H_iso, H_kin)
c       d(44)    -  Isotropic Hardening Modulus - H_iso [x sqrt(2/3)]
c       d(45)    -  Kinematic Hardening Modulus - H_kin [x sqrt(2/3)]
c      Kinematic Hardening Modulus is not availabe at finite deformations!!!!
        call dinput(d(44),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (H_iso, H_kin)'
         write(iow,*) 'ERROR --> cannot read (H_iso, H_kin)'
	   goto 9999
	  endif
c write (Y0, Yi, delta)
        write(iow,1014) d(44), d(45)
        d(44) = d(44)*sqrt(2.0/3.0)    ! change to fit internal FEAP definition
        d(45) = d(45)*sqrt(2.0/3.0)    ! change to fit internal FEAP definition

c define length of history array at each gauss point: nhv 
	  d(74) = dfloat(20)    
      elseif (imat.eq.9) then        
c=(imat=9)======== Discrete nonlinear viscoelastic fiber model with damage (CM 2022) =======================          
        
c Read record 3: Neo Hookean Matrix Material Parameters (K, mu)
c       d(1)    -  K, Bulk  modulus
c       d(2)    -  mu, Stiffness Parameter (Shear modulous)
c--------------------------------------------------------------------
c
c read (K, mu)
        call dinput(d(1),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (K, mu)'
         write(iow,*) 'ERROR --> cannot read (K, mu)'
	   goto 9999
	  endif
c write (K, mu)
        write(iow,1040) d(1),d(2)      
c
c Read record 4: Collagen Related Parameters (kf,eta_slid,T0,eta_rec)
c       d(5)    -  kf, Fibril stiffness
c       d(6)    -  eta_slid, Sliding rate
c       d(7)    -  T0, Homeostatic target stress
c       d(8)    -  eta_rec, Recovery rate         
c--------------------------------------------------------------------
c
c read (kf,eta_slid,T0,eta_rec)
        call dinput(d(5),4)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (kf,eta_slid,T0,eta_rec)'
         write(iow,*) 'ERROR --> cannot read (kf,eta_slid,T0,eta_rec)'
	   goto 9999
	  endif
c write (kf,eta_slid,T0,eta_rec)
        write(iow,1041) d(5),d(6),d(7),d(8)                   
c
c Read record 5: Collagen Failure Parameters (alpha, beta, p, no_it)
c       d(13)    -  alpha, Damage parameter 1
c       d(14)    -  beta, Damage parameter 2    
c       d(15)    -  p, fourth invariant perturberance
c       d(16)    -  max_it, Number of iterations            
c--------------------------------------------------------------------
c 
c read (alpha, beta, p, max_it)
        call dinput(d(13),4)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (alpha, beta, p, no_it)'
         write(iow,*) 'ERROR --> cannot read (alpha, beta, p, no_it)'
	   goto 9999
	  endif
c write (alpha, beta, p, max_it)
        write(iow,1042) d(13),d(14),d(15),d(16)       
c
c Read record 6: Straightening Stretch Distribution Parameters (lam1,lam2,lam3)
c       d(10)    -  lam1, Minimum straightening stretch
c       d(11)    -  lam2, Mode straightening stretch
c       d(12)    -  lam3, Maximum straightening stretch   
c--------------------------------------------------------------------
c
c read (lam1,lam2,lam3)
        call dinput(d(10),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (lam1,lam2,lam3)'
         write(iow,*) 'ERROR --> cannot read (lam1,lam2,lam3)'
	   goto 9999
	  endif
c write (lam1,lam2,lam3)
        write(iow,1043) d(10),d(11),d(12)      
c
c Read record 7: Numerical integration parameters for collagen (sio,nof)
c       d(20)    -  sio, Spherical integration order 
c       d(21)    -  nof, Number of CFPG-complexes in a fiber
c--------------------------------------------------------------------
c
c read (sio,nof)
        call dinput(d(40),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (sio,nof)'
         write(iow,*) 'ERROR --> cannot read (sio,nof)'
	   goto 9999
        endif
        
c nhv (history vector length at Gausspoint)         
        
        sio=int(d(40))
        if (ind.eq.1) then 
          d(40) = 0d0     
          sio = 0
        endif
        d(43) = 5d0 + 1d0*d(41) !density,fiber stretch,W4,W44,PB and then the proteoglycan stretch at each CFPG complex
        call int_sphere(sio,no,s)
        d(42) = no        
        
c Generate the straightening stretches that are to be used in every fiber            
        
        den1 = (d(12)-d(10))*(d(11)-d(10))
        den2 = (d(12)-d(10))*(d(12)-d(11))
        Rmid = (d(11)-d(10))**2d0/(den1)       
      
        left = nint(Rmid*(d(41)-1))
        right = (d(41)-1)-left
      
        gam_1 = Rmid/left
        gam_2 = (1-Rmid)/right
      
        R = 0d0
        count = 0
      
        do j=1,left
          count = count + 1  
          lams(count) = d(10) + sqrt(den1*R)
          cgam(count) = R
          gam(count) = gam_1
          R = R + gam_1          
        enddo    
          
        do j=1,(right+1)
          count = count + 1  
          lams(count) = d(12) - sqrt(abs(den2*(1-R)))
          cgam(count) = R
          if (count.lt.d(21)) then
            gam(count) = gam_2
          endif         
          R = R + gam_2              
        enddo           
      
        lams(count) = d(12)
      
        d(44) = Rmid
        d(45) = left        
        
c write (sio, nh, no, nof)
        write(iow,1044) d(40), d(42), d(41), d(43)
c
c Read record 8 :Collagen Density Parameter (rho_coll)
c       d(17)    -  rho_coll, Collagen density at a gauss point
c--------------------------------------------------------------------
c        
c read (rho_coll)
        call dinput(d(17),1)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (rho_coll)'
         write(iow,*) 'ERROR --> cannot read (rho_coll)'
	   goto 9999
	  endif
c write (rho_coll)
        write(iow,1045) d(17)    
c
c--------------------------------------------------------------------
c        
c Set number of history variables at Gauss point      
        
        nhv = int(d(43))*int(d(42)) ! history for collagen
	  !nhv = hist in each fiber*number of fibers
        
	  d(74) = dfloat(nhv)           
c        
c--------------------------------------------------------------------        
c                      
      endif  !imat

c        
c---------------------------------------------------------------CM (2022)    
c     Import of a boundary.txt file that provides the boundry 
c     conditions for the first three degrees of freedon (x,y,z)
c     at every node. This info is then used when using the solution
c     command "mesh,fn", meaning that instead of releasing all degrees
c     of freedom (as well as the enhanced ones) for nodes associated with
c     element that has a discontinuity, we now only release the enhanced
c     ones.
c
c     The code is also altered in umacr3.f, where the fn file is written
c
c------------------------------------------------------------------------     
c 
     !! inquire(file="boundary.txt", exist=exist)
     !!   
     !!   if (exist) then      
     !!       
     !!      write(*,*)   '***** boundary.txt successfully imported *****'
     !!      write(iow,*) '***** boundary.txt successfully imported *****'
     !!       
     !!      open(67,iostat=stat,file="boundary.txt",status="old")
     !!     
     !!      read(67,3000) bsz
     !!      do j=1,bsz
     !!         read(67,3001) boun_og(j,1),boun_og(j,2),boun_og(j,3)
     !!      enddo
     !!   
     !!      close(67)  
     !!     
     !!   else  
     !!
     !!      write(*,*)   '***** ERROR => boundary.txt file not in directo
     !!1      ry *****'
     !!      write(iow,*) '***** ERROR => boundary.txt file not in directo
     !!1      ry *****'
     !!         
     !!   endif          
c        
c--------------------------------------------------------------------        
c                           
        
c======= END CONTINUUM MATERIAL DATA =============================

        
        
        
        
        
c======= BEGIN DISCRETE MATERIAL DATA ============================
c
c Read record 4 (imatd,ld,consflag)
c
      call dinput(d(72),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (imatd,ld)'
         write(iow,*) 'ERROR --> cannot read (imatd,ld)'
	   goto 9999
	  endif
       imatd = int(d(72))

c check and write imatd
	  if (imatd.EQ.1) then
	   write (iow,150)
	   write (iow,98) int(d(73))
	  elseif (imatd.EQ.2) then
	   write (iow,151)
	   write (iow,98) int(d(73))
	  elseif (imatd.EQ.3) then
	   write (iow,152)
	   write (iow,98) int(d(73))
	  elseif (imatd.EQ.4) then
	   write (iow,153)
	   write (iow,98) int(d(73))
	  elseif (imatd.EQ.5) then
	   write (iow,154)
	   write (iow,98) int(d(73))
	  elseif (imatd.EQ.6) then
	   write (iow,157)
	   write (iow,98) int(d(73)) 
	  elseif (imatd.EQ.7) then
	   write (iow,156)
	   write (iow,98) int(d(73))   
	  elseif (imatd.EQ.8) then
	   write (iow,155)
	   write (iow,98) int(d(73))            
	  else
         write(*,*)   
	1      'ERROR --> (imatd,ld) are not  valid numbers'
         write(iow,*)   
	1      'ERROR --> (imatd,ld) are not  valid numbers'
	   goto 9999
	  endif

      if(imatd .eq. 1) then
c=(imatd=1)=========== isotropic Traction model================71
c
c Read record 5    (eps, sig_0, a, b3, uh0, cflag)
c--------------------------------------------------------------------
c read (eps, sig_0, a, b, uh0, cflag)
        call dinput(d(19),6)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (eps,sig_0,a,b,uh0,cflag)'
         write(iow,*) 'ERROR --> cannot read (eps,sig_0,a,b,uh0,cflag)'
	   goto 9999
	  endif
c write (sig_0, a, b)
        write(iow,1050) d(19),d(20),d(21),d(22),d(23)
        If (d(24).eq.1) then 
           write (iow,97) 
        else 
           write (iow,96) 
        endif
c--------------------------------------------------------------------
c Read record 6    (Pen)
c--------------------------------------------------------------------
c read (Pen)
        call dinput(d(30),1)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Pen)'
         write(iow,*) 'ERROR --> cannot read (Pen)'
	   goto 9999
	  endif
c write (Pen)
        write(iow,1052) d(30)

c define length of history array at each gauss point of disc. material
c nhvd
         d(75) = dfloat(1)

      elseif(imatd .eq. 2) then
c=(imatd=2)===== transversely isotropic Traction model=============71
c
c Read record 5    (eps, sig_0, a, b, alpha, uh0)
c--------------------------------------------------------------------
c read (sig_0, a, b, alpha)
        call dinput(d(19),6)
        if(errck) then
        write(*,*)  
	1  'ERROR --> cannot read (eps, sig_0, a, b, alpha, uh0)'
        write(iow,*)
	1  'ERROR --> cannot read (eps, sig_0, a, b, alpha, uh0)'
	   goto 9999
	  endif
c write (sig_0, a, b)
        write(iow,1051) d(19),d(20),d(21),d(22),d(23),d(24)
c Read record 6    (t0_r,t0_t,t0_z)
c--------------------------------------------------------------------
c read (t0_r,t0_t,t0_z)
        call dinput(d(25),3)
        if(errck) then
        write(*,*)  'ERROR --> cannot read (t0_r,t0_t,t0_z)'
        write(iow,*)'ERROR --> cannot read (t0_r,t0_t,t0_z)'
	   goto 9999
	  endif
c write (t0_r,t0_t,t0_z)
        write(iow,1053) d(25),d(26),d(27)
c--------------------------------------------------------------------
c Read record 7    (Pen)
c--------------------------------------------------------------------
c read (Pen)
        call dinput(d(30),1)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Pen)'
         write(iow,*) 'ERROR --> cannot read (Pen)'
	   goto 9999
	  endif
c write (Pen)
        write(iow,1052) d(30)

c define length of history array at each gauss point of disc. material
c nhvd
         d(75) = dfloat(1)

      elseif(imatd .eq. 3) then
c=(imatd=3)=========== isotropic viscoelastic Traction model================71
c
c Read record 5    (sig_0, a, b, mu)
c--------------------------------------------------------------------
c read (eps,sig_0, a, b, mu)
        call dinput(d(19),5)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (eps, sig_0, a, b, mu)'
         write(iow,*) 'ERROR --> cannot read (eps, sig_0, a, b, mu)'
	   goto 9999
	  endif
c write (sig_0, a, b, mu)
        write(iow,1060) d(19),d(20),d(21),d(22),d(23)

c--------------------------------------------------------------------
c Read record 6    (Pen)
c--------------------------------------------------------------------
c read (Pen)
        call dinput(d(30),1)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Pen)'
         write(iow,*) 'ERROR --> cannot read (Pen)'
	   goto 9999
	  endif
c write (Pen)
        write(iow,1052) d(30)
c define length of history array at each gauss point of disc. material
c nhvd
         d(75) = dfloat(4)

      elseif(imatd .eq. 4) then
c=(imatd=4)=========== smooth Traction model================71
c
c Read record 5    (c, a)
c--------------------------------------------------------------------
c read (c, a)
        call dinput(d(21),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (c, a)'
         write(iow,*) 'ERROR --> cannot read (c, a)'
	   goto 9999
	  endif
c write (c, a)
        write(iow,1070) d(21),d(22)

c write t_0, G_f
        d(20) = d(21)/(d(22)*2.71828)
        write(iow,1071) d(20),d(21)/(d(22)**2)
        

c--------------------------------------------------------------------
c Read record 6    (Pen)
c--------------------------------------------------------------------
c read (Pen)
        call dinput(d(30),1)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Pen)'
         write(iow,*) 'ERROR --> cannot read (Pen)'
	   goto 9999
	  endif
c write (Pen)
        write(iow,1052) d(30)
c define length of history array at each gauss point of disc. material
c nhvd
         d(75) = dfloat(1)

      elseif(imatd .eq. 5) then
c=(imat=5)=========== viscoelastic smooth isotropic cohesive model========71
c
c Read record 2    (c, a, mu)
c--------------------------------------------------------------------
c read (c, a)
        write(iow,1080) 
        call dinput(d(21),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (c, a, mu)'
         write(iow,*) 'ERROR --> cannot read (c, a, mu)'
	    goto 9999
	  endif
c write (c, a, mu)
        write(iow,1081) d(21),d(22),d(23)

c write t_0, G_f
        d(20) = d(21)/(d(22)*2.71828)
        write(iow,1071) d(20),d(21)/(d(22)**2)
        

c--------------------------------------------------------------------
c Read record 3    (Pen)
c--------------------------------------------------------------------
c read (Pen)
        call dinput(d(30),1)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (Pen)'
         write(iow,*) 'ERROR --> cannot read (Pen)'
	   goto 9999
	  endif
c write (Pen)
        write(iow,1052) d(30)
c define length of history array at each gauss point of disc. material
c nhv
        d(75) = dfloat(4)
 
      elseif (imatd.eq.6) then

c=(imatd=10)===== Linear Elastic model (CM 2022) =================================================     

c------------------------------------------------------------------------------
c Simple linear elastic model (t = Ku), once a maximum princiapl stress of 
c sig0 occurs. Discontinuity occurs with normal in said direction
c------------------------------------------------------------------------------
c
c Example input:
c    10,3
c    K
c
c------------------------------------------------------------------------------         
c
c Read record 10: Shape Parameters (sig0, K)
c       d(20)    -  sig0
c       d(21)    -  K, stiffness
c--------------------------------------------------------------------
c
c read (sig0, K)
        call dinput(d(20),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (sig0,K)'
         write(iow,*) 'ERROR --> cannot read (sig0,K)'
	   goto 9999
	  endif
c write (sig0, K)
        write(iow,1095) d(20),d(21)
c
c Read record 13: Additional Parameters (sig0, K)
c       d(20)    -  sig0, maximum cauchy stress stress
c       d(21)    -  K, stiffness parameter
c--------------------------------------------------------------------
c
c Number of history variables at each gauss point of the cohesive zone  
c (Dont actually need one, but specifiy one just to make life easy)     
c--------------------------------------------------------------------
c  
        d(75) = dfloat(1)         
 
      elseif (imatd.eq.7) then
        
c=(imatd=6)===== Isotropic exponential Traction model (CM 2022) ====================================       

c------------------------------------------------------------------------------
c The traction seperation law showcased in Gasser and Holzapfel 2006. The 
c discontinuity will appear once the maximum principal cauchy stress (sig0) is 
c exceeded, and the discontinuity will be normal to this direction.
c------------------------------------------------------------------------------
c
c Example input:
c    7,3
c    sig0, a, b
c    sof, del_init
c
c------------------------------------------------------------------------------   
c
c Read record 10: Shape Parameters (a,b,c,d,e,fac)
c       d(20)    -  sig0
c       d(21)    -  a, Shape parameter
c       d(22)    -  b, Shape parameter
c--------------------------------------------------------------------
c
c read (a,b)
        call dinput(d(20),3)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (sig0,a,b)'
         write(iow,*) 'ERROR --> cannot read (sig0,a,b)'
	   goto 9999
	  endif
c write (a,b)
        write(iow,1094) d(20), d(21),d(22)   
c
c Read record 13: Additional Parameters (sof, del_init)
c       d(26)    -  sof, Softening [sof=0 then no softening, softening otherwise]
c       d(27)    -  del_init, Initial damage
c--------------------------------------------------------------------
c
c read (sof,del_init)
        call dinput(d(29),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (sof,del_init)'
         write(iow,*) 'ERROR --> cannot read (sof,del_init)'
	   goto 9999
	  endif
c write (pen,sof)
        write(iow,1093) d(29),d(30)     
c        
c Number of history variables at each gauss point of the cohesive zone  
c i.e. The damage (displacement), the magnitude of the displacement vector      
c--------------------------------------------------------------------
c  
        d(75) = dfloat(1)           
        
      elseif (imatd.eq.8) then  
            
c=(imatd=8)===== Exponential normalised transversally isotropic Traction model (CM 2022)====================       

c------------------------------------------------------------------------------
c The traction seperation law showcased in Miller and Gasser 2022. The discontinuity 
c will appear once one gauss point in an element has localised i.e. the determinant
c of the acoustic tensor (detQ) in any of the referential fiber directions, specified in 
c umesh3.f becomes negative 
c------------------------------------------------------------------------------
c
c Example input:
c    8,3
c    sig0, a, b, c, d, e, fac
c    G_n, G_t
c    sof, del_init
c
c------------------------------------------------------------------------------  
c
c Read record 10: Shape Parameters (sig0,a,b,c,d,e,fac)
c       d(20)    -  sig0, peak normal stress [NOT USED BUT THERE IF WANT IN FUTURE]
c       d(21)    -  a, Shape parameter
c       d(22)    -  b, Shape parameter
c       d(23)    -  c, Shape parameter
c       d(24)    -  d, Shape parameter
c       d(25)    -  e, Shape parameter
c       d(26)    -  fac, integral shape parameter
c--------------------------------------------------------------------
c
c read (sig0,a,b,c,d,e,fac)
        call dinput(d(20),7)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (sig0,a,b,c,d,e,fac)'
         write(iow,*) 'ERROR --> cannot read (sig0,a,b,c,d,e,fac)'
	   goto 9999
	  endif
c write (sig0,a,b,c,d,e,fac)
        write(iow,1090) d(20),d(21),d(22),d(23),d(24),d(25),d(26)      
c
c Read record 12: Fracture energy (G_n,G_t)
c       d(27)    -  G_n, Fracture energy in normal direction
c       d(28)    -  G_t, Fracture energy in transversal direction
c--------------------------------------------------------------------
c
c read (G_n,G_t)
        call dinput(d(27),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (G_n,G_t)'
         write(iow,*) 'ERROR --> cannot read (G_n,G_t)'
	   goto 9999
	  endif
c write (G_n,G_t)
        write(iow,1092) d(27),d(28)      
c
c Read record 13: Additional Parameters (sof,del_init)
c       d(29)    -  sof, Softening [sof=0 then no softening, softening otherwise]
c       d(30)    -  del_init, Initial damage (normal and transverse)
c--------------------------------------------------------------------
c
c read (sof,del_init)
        call dinput(d(29),2)
        if(errck) then
         write(*,*)   'ERROR --> cannot read (sof,del_init)'
         write(iow,*) 'ERROR --> cannot read (sof,del_init)'
	   goto 9999
	  endif
c write (sof,del_init)
        write(iow,1096) d(29),d(30)      
c   
c Number of history variables at each gauss point of the cohesive zone  
c i.e. The damage (displacement), in both the normal and transverse directions        
c--------------------------------------------------------------------
c        
        d(75) = dfloat(2)        

      endif  !imatd        
        
c======= END DISCRETE MATERIAL DATA ============================

	write (iow,2000)


c Formats for output

98    format(
     1 5x, '     Discontinuity integration order:    ',i3/)
96    format(
     1 5x, '     Using inconsistent TSL tangent    '/)
97    format(
     1 5x, '     Using consistent TSL tangent    '/)
99    format(
     1 5x, '     Bulk integration order:    ',i3/)
100   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Linear isotropic Elasticity (St. Venant Model)   '/)
101   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Neo Hookean model                                '/)
102   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Fiber-reinforced composite                       '/)
103   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Modified Delfino model                           '/)
104   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Distributed fiber model                           '/)
105   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     7 parameter model for arterail wall             '/)
106   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '    FEAP small strain von Mises plasticity            '/)
107   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     FEAP finite strain plasticity                    '/)
108   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Discrete nonlinear viscoelastic fiber model with damage
     1 (CM 2022)'/)
150   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Linear isotropic Traction model                  '/)
151   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Linear transversally isotropic Traction model    '/)
152   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Linear  isotropic viscoelastic Traction model    '/)
153   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Smooth  isotropic  Traction model    '/)
154   format(
     1 5x, '------------------------------------------------------'/
     1 5x, '     Smooth  viscoelastic isotropic  Traction model   '/)
155   format(
     1 5x, '==================COHESIVE MATERIAL==================='/
     1 ''/
     1 5x, '     Exponential normalised transversally isotropic
     1 traction model (CM 2022)  '/)
156   format(
     1 5x, '==================COHESIVE MATERIAL==================='/
     1 ''/
     1 5x, '     Exponential isotropic traction Model (CM 2022)  '/)
157   format(
     1 5x, '==================COHESIVE MATERIAL==================='/
     1 ''/
     1 5x, '     Linear elastic (CM 2022)  '/)
190   format(
     1 5x, '     Eq. discontinuity in the center of the element   '/
     2 5x, '     represents a finite number of discont. with      '/
     3 10x, '      Delta L =            ', e12.5                   /)
203   format(
     1 5x, '=====================FORMULATION======================'/
     2 5x, ' '/      
     2 5x, 'Q1C1: Part. of Unity FE formulation for FEAP ver. 7.1 '/
	2 5x  '                                       (c) TCG 2002   '/
     2 5x, ' '/       
     1 5x, '====================BULK MATERIAL====================='/)
301   format(
     1 5x, 'Gauss integration order     ',i2/
     1 5x, '          (1) one point     '/
     1 5x, '          (2) 2x2x2 points  '/
     1 5x, '          (3) 3x3x3 points  '/
     1 5x, '          (4) 4x4x4 points  '/
     1 5x, '          (9) 9 points      '/)


1001  format(
     1  10x,'Density                    ',e12.5/
     2  10x,'M.-int.   (0-diag.;1-con.) ',i3/
     3  10x,'Sym.Tang. (0-No.;1-yes.)   ',i3/)



1003  format(
     1  10x,'Young Modulus              ',e12.5/
     2  10x,'Poison ratio               ',e12.5/
     5 5x, '------------------------------------------------------')

1004  format(
     1  10x,'Bulk Modulus               ',e12.5/
     2  10x,'NeoHookean parameter       ',e12.5/
     5 5x, '------------------------------------------------------')

1005  format(
     1  10x,'First fiber parameter k1   ',e12.5/
     2  10x,'Secound fiber parameter k2 ',e12.5/
     5 5x, '------------------------------------------------------')
1006  format(
     1  10x,'Bulk Modulus               ',e12.5/
     2  10x,'NeoHookean parameter       ',e12.5/
     5 5x, '------------------------------------------------------')

1007  format(
     1  10x,'First parameter k1         ',e12.5/
     2  10x,'Secound parameter k2       ',e12.5/
     5 5x, '------------------------------------------------------')

1008  format(
     1  10x,'Bulk Modulus               ',e12.5/
     2  10x,'NeoHookean parameter       ',e12.5/
     5 5x, '------------------------------------------------------')

1009  format(
     1  10x,'Media parameter k1         ',e12.5/
     2  10x,'Media parameter k2         ',e12.5/
     2  10x,'Media Distribution         ',e12.5/
     5 5x, '------------------------------------------------------')
1010  format(
     1  10x,'Adv. parameter k3          ',e12.5/
     2  10x,'Adv. parameter k4          ',e12.5/
     2  10x,'Adv. Distribution          ',e12.5/
     5 5x, '------------------------------------------------------')

1011  format(
     1  10x,'Yield stress Y0               ',e12.5/
     2  10x,'Yield stress Yi               ',e12.5/
     2  10x,'Delta parameter               ',e12.5/
     5 5x, '------------------------------------------------------')
     
1012  format(
     1  10x,'Isotropic Hardening Modulus - H_iso   ',e12.5/
     2  10x,'Kinematic Hardening Modulus - H_kin   ',e12.5/
     5 5x, '------------------------------------------------------')

1013  format(
     1  10x,'Yield stress Y0 [x sqrt(2/3)] ',e12.5/
     2  10x,'Yield stress Yi [x sqrt(2/3)] ',e12.5/
     2  10x,'Delta parameter               ',e12.5/
     5 5x, '------------------------------------------------------')
     
1014  format(
     1  10x,'Isotropic Hardening Modulus - H_iso [x sqrt(2/3)] ',e12.5/
     2  10x,'Kinematic Hardening Modulus - H_kin [x sqrt(2/3)] ',e12.5/
     5 5x, '------------------------------------------------------')


1050  format(
     1  10x,'Initial compliance         ',e12.5/
     1  10x,'Cohesive strength sig_0    ',e12.5/
     2  10x,'Damage parameter a         ',e12.5/
     3  10x,'Damage parameter b         ',e12.5/
     4  10x,'Damage ini value           ',e12.5/
     5 5x, '------------------------------------------------------')

1060  format(
     1  10x,'Initial compliance         ',e12.5/
     1  10x,'Cohesive strength sig_0    ',e12.5/
     2  10x,'Damage parameter a         ',e12.5/
     3  10x,'Damage parameter b         ',e12.5/
     4  10x,'Rayleight damping mu       ',e12.5/
     5 5x, '------------------------------------------------------')

1070  format(
     1  10x,'Initial stiffness          ',e12.5/
     1  10x,'Parameter a                ',e12.5/
     5 5x, '------------------------------------------------------')

1071  format(
     1  10x,'Cohesive strength          ',e12.5/
     1  10x,'Fracture energy            ',e12.5/
     5 5x, '------------------------------------------------------')

1080   format(
     1 10x, 'Smooth viscoelastic traction model                '/)


1081  format(
     1  10x,'Initial stiffness          ',e12.5/
     1  10x,'Parameter a                ',e12.5/
     4  10x,'Rayleight damping mu       ',e12.5/
     5  5x, '------------------------------------------------------')

1051  format(
     1  10x,'Initial compliance         ',e12.5/
     1  10x,'Cohesive strength sig_0    ',e12.5/
     2  10x,'Damage parameter a         ',e12.5/
     3  10x,'Damage parameter b         ',e12.5/
     4  10x,'Stiffness ratio  alpha     ',e12.5/
     4  10x,'Damage ini value           ',e12.5/
     5 5x, '------------------------------------------------------')
1053  format(
     1  10x,'Cohesive strength t0_r    ',e12.5/
     2  10x,'Cohesive strength t0_t    ',e12.5/
     3  10x,'Cohesive strength t0_z    ',e12.5/
     5 5x, '------------------------------------------------------')

1052  format(
     1  10x,'Penalty against penetration',e12.5/
     5 5x, '------------------------------------------------------')
      
!==================CM 2022======================================
      
c------BULK------      
      
1040  format(
     5 5x, '------------------------------------------------------'/
     1  10x,'k, Bulk  modulus                ',e12.5/
     1  10x,'mu, Shear modulous              ',e12.5/ 
     5 5x, '------------------------------------------------------')    
      
1041  format(
     1  10x,'kf, Fibril stiffness            ',e12.5/
     1  10x,'eta_slid, Sliding rate          ',e12.5/ 
     1  10x,'T0, Homeostatic target stress   ',e12.5/
     2  10x,'eta_rec, Recovery rate          ',e12.5/      
     5 5x, '------------------------------------------------------')  
      
1042  format(
     1  10x,'alpha, Damage parameter 1       ',e12.5/
     1  10x,'beta, Damage parameter 2        ',e12.5/ 
     1  10x,'p, perturberance of I4          ',e12.5/
     1  10x,'no_it, Number of NR iterations  ',e12.5/       
     5 5x, '------------------------------------------------------') 
      
1043  format(
     1  10x,'lam1, Minimum straightening stretch      ',e12.5/
     1  10x,'lam2, Mode straightening stretch         ',e12.5/
     2  10x,'lam3, Maximum straightening stretch      ',e12.5/
     5 5x, '------------------------------------------------------')  

1044  format(
     1  10x,'sio, Spherical integration order              ',e12.5/
     1  10x,'no, Number of spherical Integration points    ',e12.5/ 
     1  10x,'nof, Number of CFPG-complexes in a fiber      ',e12.5/
     2  10x,'nh, Number of history variables in a fiber    ',e12.5/ 
     5 5x, '------------------------------------------------------')   
      
1045  format(
     1  10x,'rho, Collagen density at gauss point          ',e12.5/
     1  '') 
      
c------DISCONTINUITY------            
      
1090  format(
     5 5x, '------------------------------------------------------'/  
     1  10x,'sig0, Peak stress [NOT USED]      ',e12.5/
     1  10x,'a, Shape parameter                ',e12.5/
     1  10x,'b, Shape parameter                ',e12.5/
     2  10x,'c, Shape parameter                ',e12.5/
     3  10x,'d, Shape parameter                ',e12.5/
     4  10x,'e, Shape parameter                ',e12.5/
     4  10x,'fac, integral shape parameter     ',e12.5/     
     5 5x, '------------------------------------------------------')
      
1092  format(
     1  10x,'G_n, Normal fracture energy       ',e12.5/
     1  10x,'G_t, Transverse fracture energy   ',e12.5/ 
     5 5x, '------------------------------------------------------') 
      
1093  format(
     1  10x,'sof, softening in stiffness [0=NO]        ',e12.5/
     1  10x,'del_init, Initial damage (magnitude of
     1 displacement)       ',e12.5/)   
      
1096  format(
     1  10x,'sof, softening in stiffness [0=NO]        ',e12.5/
     1  10x,'del_init, Initial damage (both in normal
     1 and transverse directions)        ',e12.5/)         
      
1094  format(
     5 5x, '------------------------------------------------------'/  
     1  10x,'sig0, Peak Stress                 ',e12.5/
     1  10x,'a, Shape parameter                ',e12.5/
     1  10x,'b, Shape parameter                ',e12.5/   
     5 5x, '------------------------------------------------------')
      
1095  format(
     5 5x, '------------------------------------------------------'/  
     1  10x,'sig0, Peak stress                 ',e12.5/
     1  10x,'k, Stiffness                      ',e12.5/)
      
2000  format(
     1 5x, '======================================================'/)

3000  format(i5)        
      
3001  format(i1,x,i1,x,i1)         
      
      
      
9999  end
c       d(5)    -  kf, Fibril stiffness
c       d(6)    -  eta_slid, Sliding rate
c       d(7)    -  T0, Homeostatic target stress
c       d(8)    -  eta_rec, Recovery rate    
c       d(13)    -  alpha, Damage parameter 1
c       d(14)    -  beta, Damage parameter 2     
c       d(10)    -  lam1, Minimum straightening stretch
c       d(11)    -  lam2, Mode straightening stretch
c       d(12)    -  lam3, Maximum straightening stretch       
      
c       d(20)    -  sio, Spherical integration order 
c       d(21)    -  nof, Number of CFPG-complexes in a fiber   
      
c       d(17)    -  rho_coll, Collagen density at a gauss point      