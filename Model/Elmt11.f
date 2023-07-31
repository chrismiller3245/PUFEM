     
      subroutine elmt11(d,ul,xl,ix,tl,s,p,ndf,ndm,nst,isw)
c
c____________________________________________________________________71
c
c     THREE DIMENSIONAL FINITE DEFORMATION ELEMENT for F E A P
c     based on the Partition of Unity Finite Element Method (PUFEM)
c
c
c  properties vector d(200) specification
c
c  d(1) .... d(49)   Material properties of bulk material and cohesive zone models.
c                    Bulk material models, see subroutine mat_c_11.f
c                    Cohesive zone models, see subroutine mat_d_11.f
c
c  d(50) Density of bulk material
c  d(51) Mass interpolator 
c                          (0) Diagonal
c                          (1) Consistent
c  d(52) Mass proportional Rayleigh damping parammeter          ! not used
c  d(53) Stiffness proportional Rayleigh damping parammeter     ! not used
c
c  d(70) Bulk material number  (imat)
c  d(71) Order of bulk material integration (l)
c  d(72) Cohesive law material number  (imatd)
c  d(73) Order of cohesive zone integration (ld)
c  d(74) Length of bulk material-related history vector at the Gausspoint (nhv; part of the time-dependent history vector)
c  d(75) Length of cohesive zone-related history vector at the Gausspoint (nhvd; part of the time-dependent history vector)
c  d(76) Length of time-dependent history vector (stored at nh1 and nh2 in FEAP)
c  d(77) Length of time-independent history vector (stored at nh3 in FEAP)
c  d(78) Length of bulk material-related history vector (part of the time-dependent history vector)
c
c  Declare variable types
c
c       d     Material parameters
c       ul    Solution array (displacements)
c       xl    (Lagrange) coordinates
c       ix    Element connections
c       tl    Temperature
c       s     Element (stiffness) matrix
c       p     Element vector
c       ndf   Degree of freedom (max) per node
c       ndm   Space dimension of mesh 
c       nst   Size of element array
c       isw   Task parameter 
c       imat  Continuum material model
c       imatd Discrete material model
c       nhv   Length of history vector per gauss point
c       fef   FE formulation flag
c       gp    Type of Gauss integration order
c       shp   Shape function and its derivatives
c       n     number of current element
c       nel   nodes per element (in common 'eldata.h')
c       l     Gausspoint id
c       lint  Number of Gauss-points lint=l*l*l
c       sg    Gauss points (1-3) and weights (4)
c       xsj   Jacobi determinante
c       dvol  xsj*sg(4,l)
c       fi    Deformation gradient
c       finv  Inverse deformation gradient
c       bei   Left Cauchy Green strain
c       xji   Jacobian determinat at Gausspoints
c       press Element pressure 
c       hel   FE formulation dependent history storage length  
c       nh1   History pointer previous time
c       nh2   History pointer current time
c       d1    Updated bulk-modulus
c       augf  Augmented Lagrangian parameter 1
c       ma    Option number (reference or current conf.) (in 'eldata.h')    
c       mct   Type of surface loading (in 'eldata.h')
c       sym   Symmetry-flag (.true.-symmetric material)
c
c   History vector structure at the element level
c
c   - time-dependent history (stored in FEAP history vector at positions nh1 and nh2, respectively)
c	
c        ---       -------                     ------
c       |       | Gauss point 1 of the bulk material |
c       |       |     entry 1                        |
c       |       |     entry 2                         > d(74)
c       |       |     entry 3                        |
c       |       |     ....                           |
c       |       |                               -----
c       |       | Gauss point 2 of the bulk material
c       | d(78)<      entry 1
c       |       |     entry 2
c       |       |     entry 3
c       |       |     ....
c       |       |                               
c       |       | Gauss point 3 of the bulk material
c       |       |     ....
c       |       |     ....
c       |       |                               
c d(76)<        -------                       ------
c       |       | Gauss point 1 of the cohesive zone |
c       |       |     entry 1                        |
c       |       |     entry 2                         > d(75)
c       |       |     entry 3                        |
c       |       |     ....                           |
c       |       |                              ------   
c       |       | Gauss point 2 of the cohesive zone
c       |       |     entry 1 
c       |       |     entry 2
c       |       |     entry 3
c       |       |     ....
c       |       |                               
c       |       | Gauss point 3 of the cohesive zone
c       |       |     ....
c       |       |     ....
c        ----   -------
c
c
c   - time-independent history (stored in FEAP history vector at positions nh3)
c        ----  
c       |      entry 1 
c d(77)<       entry 2   
c       |      entry 3       
c       |      ....       
c        ----  

c--------------------------------------------------------------------71
      IMPLICIT none
      include  'cdata.h'
      include  'eldata.h'
      include  'fdata.h'
      include  'comblk.h'
      include  'hdata.h'
      include  'augdat.h'
      include  'iofile.h'
      include  'elHist.h'
      include  'cTets.h'

	logical sym
      integer ix(nel), ndf, ndm, nst, isw, imat, nhv
	integer order, nst2, ld, nhvd
	integer i, j, k, l, hel, imatd
      integer lint, lint_sph, nn, nh, no, nof
      real*8  d(*), ul(ndf,nel,6), xl(ndm,nel), tl(nel), s(nst,nst)
      real*8 shp(10), p(nst), xu(ndm,nen), xbar(3)
	real*8 sg(5), xsj, shp_t(4,10)
      real*8 spg(3,240), rho
      integer check
      real*8  check_hist(100000)

      save

      data sg /0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.1666666666666667d0/

c recover imat,imatd,l,ld,nhv,nhvd
      imat = int(d(70))
      l    = int(d(71))
      imatd= int(d(72))
      ld   = int(d(73))
      nhv  = int(d(74))
      nhvd = int(d(75))
	sym = .true.      !symmetric stiffness of the bulk material

c Go to correct array processor depending on isw
c-0------------------------------------------------------------------
      IF(isw.EQ.0) then ! Output element description
       write(*,1000)  
c-1---Input material properties--------------------------------------
      ELSEIF (isw.eq.1) then   
        call inpt11 (d)
c recover imat,imatd,l,ld,nhv,nhvd
        l    = int(d(71))
        ld   = int(d(73))
        nhv  = int(d(74))
        nhvd = int(d(75))
         
c compute # of history elements and set pointers nh3 and nh1
	  nh3 = 370                          ! no. of time independent hist. elements
        nh1 = 2*(nhv*l*l*l) + nhvd*ld*ld   ! no. of time dependent hist. elements
                                           ! xxx twice  since every elment may split into a PU element
        d(76) = nh1                        ! store in d-vector
        d(77) = nh3
        d(78) =  2*(nhv*l*l*l)             ! lenght of bulk-related history
        d(79) = nhvd*ld*ld                 ! length of discontinuity-related history CJM-18/10/2022
c store in elHist.h variables to reallocate history when cohesive zone is inserted
c this concerns only the history of the bulk material
           iorder = int(d(71))
           nHistpGP= nhv
        
c Set 3-D Plot Sequence for 4-node tetrahedal elements
          if (nel.eq.4) then        
              call pltet4(11)
          else
c Set 3-D Plot Sequence for 10-node tetrahedal elements
             call pltet10(11)
         endif
c-14-----------------------------------------------------------------
      ELSEIF (isw.eq.14) then  !initialize history and cut element
c first zero all history
          do i=1,int(d(76)) 
              hr(nh1+ i)= 0.0d0
              hr(nh2+ i)= 0.0d0
          enddo
          do i=1,int(d(77)) 
              hr(nh3+ i)= 0.0d0
          enddo
c initialize history at Gauss points (time-dependent history)
        imat=int(d(70))
	  if (imat.eq.7) then ! small strain von Mises plasticity 
            l= int(d(71))
            l = l*l*l    ! total number of Gauss points
            j= nh1
            Do i= 1, l
              hr(j+ 0) = 0.0d0       !epp(1)  Cumulative plastic strain
              hr(j+ 1) = 0.0d0       !epp(2)  
              hr(j+ 2) = 0.0d0       !eps(1)  linear strain
              hr(j+ 3) = 0.0d0       !eps(2)  linear strain
              hr(j+ 4) = 0.0d0       !eps(3)  linear strain
              hr(j+ 5) = 0.0d0       !eps(4)  linear strain
              hr(j+ 6) = 0.0d0       !eps(5)  linear strain
              hr(j+ 7) = 0.0d0       !eps(6)  linear strain
              j= j +nhv
           enddo
	  elseif (imat.eq.8) then ! finite strain plasticity  
            l= int(d(71))
            l = l*l*l    ! total number of Gauss points
            j= nh1
            Do i= 1, l
              hr(j+ 0) = 0.0d0       !epp(1)  Cumulative plastic strain
              hr(j+ 1) = 0.0d0       !epp(2)  | SIGMA_n | general plasticity
              hr(j+ 2) = 0.0d0       !be(1)  incr. Left Cauchy-Green tensor
              hr(j+ 3) = 0.0d0       !be(2)  incr. Left Cauchy-Green tensor
              hr(j+ 4) = 0.0d0       !be(3)  incr. Left Cauchy-Green tensor
              hr(j+ 5) = 0.0d0       !be(4)  incr. Left Cauchy-Green tensor
              hr(j+ 6) = 0.0d0       !be(5)  incr. Left Cauchy-Green tensor
              hr(j+ 7) = 0.0d0       !be(6)  incr. Left Cauchy-Green tensor
              hr(j+ 8) = 0.0d0       !epl(1) Plastic Strain for Hardening
              hr(j+ 9) = 0.0d0       !epl(2) Plastic Strain for Hardening
              hr(j+10) = 0.0d0       !epl(3) Plastic Strain for Hardening
              hr(j+11) = 1.0d0       !F_n(1,1) Deformation gradient
              hr(j+12) = 0.0d0       !F_n(1,2) Deformation gradient
              hr(j+13) = 0.0d0       !F_n(1,3) Deformation gradient
              hr(j+14) = 0.0d0       !F_n(2,1) Deformation gradient
              hr(j+15) = 1.0d0       !F_n(2,2) Deformation gradient
              hr(j+16) = 0.0d0       !F_n(2,3) Deformation gradient
              hr(j+17) = 0.0d0       !F_n(3,1) Deformation gradient
              hr(j+18) = 0.0d0       !F_n(3,2) Deformation gradient
              hr(j+19) = 1.0d0       !F_n(3,3) Deformation gradient
              j= j +nhv
            enddo
        elseif (imat.eq.9) then ! Discrete nonlinear viscoelastic fiber model with damage
            
            lint = int(d(71))
            lint = 2*lint*lint*lint    ! We initialise all 16 (potentially used) Gauss points
            no=int(d(42))
            nh = int(d(43))
            nof = int(d(41))
            call int_sphere(int(d(40)),lint_sph,spg)
            
            nn = nh1
            
            !check = 0
            
            do i = 1,lint
              do j = 1,no       
                call dist07(j,spg,rho)
                !hr(nn + 0 )= d(17)/12.5663706144 ! rho per surface area
                hr(nn + 0 )= d(17)*rho ! rho per surface area
                hr(nn + 1 )= 1d0 !I4, squared fiber stretch
                hr(nn + 2 )= 27d0 !W4
                hr(nn + 3 )= 38d0 !w44
                hr(nn + 4 )= 0d0 !PB 
                !check_hist(check)=d(17)*rho
                !check_hist(check+1)=1d0
                !check_hist(check+2)=27d0
                !check_hist(check+3)=38d0 
                !check_hist(check+4)=0d0
                do k=1,nof ! Loop over all CFPG-complexes
                   hr(nn+4+((k*1)))= 1d0 !Proteoglycan stretch
                   !check_hist(check+4+((k*1)))=1d0
                enddo
                nn = nn + nh  
                !check = check + nh
              enddo  
            enddo   
            
        endif	    

        imatd=int(d(72))
        
        if (imatd.eq.7) then    

            lint = int(d(73))
            lint = lint*lint  
            nh = int(d(75))
            
            do i = 1,lint
                hr(nn + 0 )= d(30) ! Initial damage
                nn = nn + nh 
            enddo          
        
        elseif (imatd.eq.8) then
            
            lint = int(d(73))
            lint = lint*lint  
            nh = int(d(75))
            
            do i = 1,lint
                hr(nn + 0 )= d(30) ! Initial damage in the normal direction 
                hr(nn + 1 )= d(30) ! Initial damage in the transverse direction
                !check_hist(check)= d(30) ! Initial damage in the normal direction 
                !check_hist(check+1)= d(30) ! Initial damage in the transverse direction
                nn = nn + nh 
                !check = check + nh
            enddo           
            
        endif    
        
c initialize time-independent history 
          do i= 1,int(d(77))
	      hr(nh3+i-1)=0.0d0
	    enddo
	    hr(nh3+44) = 1.0d0
	    hr(nh3+48) = 1.0d0
	    hr(nh3+52) = 1.0d0
c ini res-deformation
c	    hr(nh3+360) = 1.0d0
c	    hr(nh3+364) = 1.0d0
	    !hr(nh3+368) = 978d0  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!REMEMEBER TO REMOVE!!!!!!!!!!!!!!!!!!!!!
c compute center of the element and write into history
         xbar = 0.0d0
	   do i= 1,4
	     do j= 1,3
	       xbar(j) = xbar(j) + xl(j,i)
	     enddo
	   enddo
	   xbar = 0.25d0*xbar
	   do i= 1,3
           hr(nh3+i+25) = xbar(i)
	   enddo
c compute referential element volume and write into history
        call eldata_11 (nel,ndm, nst2, order)
        call tetshp_11(sg, xl, 3, order, xsj, shp_t)
        hr(nh3+1) = xsj*sg(5)

c initialize ntet (number of cracked tetrahedrals, i.e. PUFEM elements)
         ntets = 0
c-2------------------------------------------------------------------
      ELSEIF (isw.eq.2) then    ! Check mesh for bad data in current config.
c Compute current geometry
         do i = 1,3
          do j = 1,nel
           xu(i,j) = xl(i,j) + ul(i,j,1)
          enddo
	   enddo
         call cktets ( n, ix, xl, ndm, nel, shp )
c-3------------------------------------------------------------------
	ELSEIF ((isw.eq.3).or.(isw.eq.4).or.(isw.eq.6).or.(isw.eq.8)) then  
           call kernel11 (ix,l,imat,imatd,nhv,d,ul,xl,s,p,
     1                 isw,sym)
	ELSEIF(isw.EQ.5)  then ! Compute Mass matrix (lumped/consistent)
	  write(*,*) '* *  Dummy task called !!'
	ELSEIF(isw.EQ.7)  then ! Residual and stiffness due to surface loading
        call surf_11(d, xl, ul, ma, ndf, ndm, nel, mct, nst, p,s)
      ELSEIF(isw.EQ.10) then ! Augmented Lagrangian update
       write (*,*) '* *  Dummy task called !!'
	ENDIF

1000  format(
     1 5x, '=========================================== (c) TC Gasser'/
     2 5x, '3D Partition of Unity Finite Element (PUFEM) for FEAP  '/
     3 5x, '                                     vers. 1.0  Oct. 2002'//
     3 5x, '                                     vers. 2.0  May  2021'//
     3 5x, '                                     vers. 3.0  Sept 2022'//
     8 5x, '========================================================='/)
9999  end
c
      
    
      
c
	subroutine Dist07(no,spg,rho)
      
c=================================================================================  
c    
c              ORIENTATION DISTRIBUTION SUBROUTINE for F E A P
c
c     PURPOSE 
c
c        Determine the density for a given orientation unit vector for a Von Mises
c        or Bingham Distribution      
c
c     INPUT variables
c      
c        no            spherical integration point number
c        spg           All the unit vectors
c
c     OUTPUT variables
c
c        rho          Density
c
c=================================================================================
      
	integer no, i,j,k
	real*8 spg(3,240),m_0(3), rho
      real*8 dp, m(3), A(3,3) ,Ainv(3,3)
      real*8 theta, phi
      real*8 b, den, k1, k2, c
      real*8 check
      data pi   /3.14159265359/
      
      include 'pvecdata.h'
      include 'eldata.h'
      include 'Visco.h'
      
      m_0 = spg(1:3,no)
      
      IF (ind.eq.0) then !ISOTROPY
      
         rho = 1d0/(4d0*pi) ! As we normalise relative to 4*pi
      
      ELSEIF (ind.eq.1) then !UNIAXIAL DISTRIBUTION
          
         rho = 1d0/(4d0*pi) ! As we normalise relative to 4*pi
         !rho = 1d0
         
      ELSEIF ((ind.eq.2).or.(ind.eq.3)) then !VON MISES MODELS A and B   
          
         ! Determine the angle using dot product 
          
         dp = 0d0
      
         do i=1,3
            dp = dp + m_0(i)*pvec(n,i,1)
         enddo
      
         theta = acos(dp)
         
         if (dp.lt.0d0) then
             theta = theta-pi
         endif
             
         ! Determine the density
         
         b = param(1)
         den = param(2)
         
         if (ind.eq.2) then
         
            rho = (1/pi)*sqrt(b/(2*pi))*(exp(b*(cos(2*theta)+1)))/den
            
         elseif (ind.eq.3) then
            
            rho = (1/pi)*sqrt(b/(2*pi))*(exp(b*(cos(2*theta-pi)+1)))/den
             
         endif
         
         ! Output to Fiber_Output.txt if necessary
         
         !if ((elno.eq.n_out).and.(gano.eq.g_out)) then
         !   write(77,6002) no, m_0(1), m_0(2), m_0(3), theta
         !endif                 
         
      ELSEIF (ind.eq.4) then !BINGHAM DISTRIBUTION  
          
         ! Multiply the unti vector m_0 by the rotation matrix R to change basis
          
         m(1) = bvec(n,1,1)*m_0(1) + bvec(n,1,2)*m_0(2)
     1         + bvec(n,1,3)*m_0(3)
         m(2) = bvec(n,2,1)*m_0(1) + bvec(n,2,2)*m_0(2)
     2         + bvec(n,2,3)*m_0(3)
         m(3) = bvec(n,3,1)*m_0(1) + bvec(n,3,2)*m_0(2)
     3         + bvec(n,3,3)*m_0(3) 
         
         check = sqrt(m(1)**2+m(2)**2+m(3)**2)
         
         ! Determine the azimuthal (theta) and elevation (phi) angles with new basis
         
         theta = atan(m(2)/m(1))
         phi = atan(m(3)/m(1))
         
         ! Determine the density
         
         k1 = param(1)
         k2 = param(2)  
         c = param(3)
         
         rho = 0.5d0*(1/c)*exp( (k1*((cos(phi)*cos(theta))**2d0))
     2           +(k2*((cos(phi)*sin(theta))**2d0)) )
         
         check = 1d0
         
         ! Output to Fiber_Output.txt if necessary
         
     !!     if ((elno.eq.n_out).and.(gano.eq.g_out)) then
     !!        write(77,6004) no, m_0(1), m_0(2), m_0(3), 
     !!1                      theta, phi,m(1), m(2), m(3)
     !!     endif          
         
      ENDIF      
       
       
      end