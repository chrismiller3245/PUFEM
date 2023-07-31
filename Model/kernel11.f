	subroutine kernel11(ix,l,imat,imatd,nhv,d,ul,xl,
	1                            s,p,isw,sym)
c
c--------------------------------------------------------------------71
c
c             FE routines of the PUFEM formulation 
c
c   time-independet History variables	at element level (stored in FEAP history vector at positions nh3)
c
c       hr(nh3 + 0)	   A_e
c       hr(nh3 + 1)	   Ve   
c       hr(nh3 + 2)	   fact=Ve+/Ve
c
c       hr(nh3 + 3)	   N0_1
c       hr(nh3 + 4)	   N0_1
c       hr(nh3 + 5)	   N0_1
c
c       hr(nh3 + 6)	   cini ... crack initialization flag
c								0  element intact 
c								1  element cracked 
c
c       hr(nh3 + 7)	   cn  ...  nodes of the disontinuity
c                                 3 ... triangle
c                                 4 ....4 node patch
c 
c       hr(nh3 + 8)	   conp
c       hr(nh3 + 9)	   conm
c
c       hr(nh3 + 10)	   L11
c       hr(nh3 + 11)	   L21
c       hr(nh3 + 12)	   L31
c       hr(nh3 + 13)	   L41
c
c       hr(nh3 + 14)	   L12
c       hr(nh3 + 15)	   L22
c       hr(nh3 + 16)	   L32
c       hr(nh3 + 17)	   L42
c
c       hr(nh3 + 18)	   L13
c       hr(nh3 + 19)	   L23
c       hr(nh3 + 20)	   L33
c       hr(nh3 + 21)	   L43
c
c       hr(nh3 + 22)	   L14
c       hr(nh3 + 23)	   L24
c       hr(nh3 + 24)	   L34
c       hr(nh3 + 25)	   L44
c
c       hr(nh3 + 26)	   Xb_1
c       hr(nh3 + 27)	   Xb_2
c       hr(nh3 + 28)	   Xb_3
c   
c       hr(nh3 + 30)	   Xn2d11
c       hr(nh3 + 31)	   Xn2d21
c
c       hr(nh3 + 32)	   Xn2d12
c       hr(nh3 + 33)	   Xn2d22
c
c       hr(nh3 + 34)	   Xn2d13
c       hr(nh3 + 35)	   Xn2d23
c
c       hr(nh3 + 36)	   Xn2d14
c       hr(nh3 + 37)	   Xn2d24
c
c       hr(nh3 + 38)	   sig_11
c       hr(nh3 + 39)	   sig_22
c       hr(nh3 + 40)	   sig_33
c       hr(nh3 + 41)	   sig_12
c       hr(nh3 + 42)	   sig_23
c       hr(nh3 + 43)	   sig_13
c
c       hr(nh3 + 44)	   fi_11
c       hr(nh3 + 45)	   fi_12
c       hr(nh3 + 46)	   fi_13
c       hr(nh3 + 47)	   fi_21
c       hr(nh3 + 48)	   fi_22
c       hr(nh3 + 49)	   fi_23
c       hr(nh3 + 50)	   fi_31
c       hr(nh3 + 51)	   fi_32
c       hr(nh3 + 52)	   fi_33
c
c       hr(nh3 + 53)	   nelem     ... number of neighbouring elements (max 300!!)
c       hr(nh3 + 54)	   no. elem1 
c       hr(nh3 + 55)	   no. elem2 
c        ....
c       hr(nh3 + 354)	   no. elem300
c       hr(nh3 + 355)	   t0
c
c       hr(nh3 + 356)	   Nold_1
c       hr(nh3 + 357)	   Nold_1
c       hr(nh3 + 358)	   Nold_1
c   
c       hr(nh3 + 360)	   fi_res_11
c       hr(nh3 + 361)	   fi_res_12
c       hr(nh3 + 362)	   fi_res_13
c       hr(nh3 + 363)	   fi_res_21
c       hr(nh3 + 364)	   fi_res_22
c       hr(nh3 + 365)	   fi_res_23
c       hr(nh3 + 366)	   fi_res_31
c       hr(nh3 + 367)	   fi_res_32
c       hr(nh3 + 368)	   fi_res_33
c
c
c  history associated with time (stored in FEAP history vector at positions nh1 and nh2, respectively) 
c
c     explained in subroutines mat_c_11.f and mat_d_11.f
c
c
c  Declare variable types
c
c       ix    Element connections
c       d     Material parameters
c       ul    Solution array (displacements)
c       xl    (Lagrange) coordinates
c       s     Element (stiffness) matrix
c       p     Element vector
c       ndf   Degree of freedom (max) per node
c       ndm   Space dimension of mesh 
c       nst   Size of element array
c       isw   Task parameter 
c       imat  Continuum Material model
c       imatd Discrete Material model
c       nhv   Length of history vector per gauss point
c       shp   Shape function and its derivatives
c       n     number of current element
c       nel   nodes per element
c       xu    Current geometry
c       l     Gausspoint id
c       lint  Number of Gauss-points lint=l*l*l
c       sg    Gauss points (1-3) and weights (4)
c       xsj   Jacobi determinante
c       dvol  xsj*sg(4,l)
c       fi    Deformation gradient
c       finv  Inverse deformation gradient
c       bei   Left Cauchy Green strain
c       xji   Jacobian determinat at Gausspoints
c       nn    Shiftvariable for history  
c       sigv  Cauchy stresses in Quadraturpoins
c       aa    Eulerian material moduli
c       nh1   History pointer previous time
c       nh2   History pointer current time
c       augf  Augmented Lagrangian parameter 1
c       dot   Dotproduct
c       nen   Maximum number of nodes on element (in common 'cdata')
c       rho   Mass per Gauss-point
c       cmom  Change in momentum 
c       numnp Number of nodes in problem (in 'cdata.h') 
c       nph   Pointer for nodal stress parameters (in 'prstrs.h')
c       sigl  Quantities at Gauss-points to be plotted
c       u0    Gap displacement (history)
c       Ae    Average area of disection (ref.conf.)
c       Ve    Element volume (ref.conf.)
c       av    A^e/V^e
c       N0    Normal vector of the discontinuity 
c       nb    normal vector at the fictitious discontinuity
c     cTets   test with violated crack criterion
c     ntets   number of test with violated crack criterion 
c     IXlen   Length of the IX array
c     IXpre   Precission of the IX array
c     IXflag  Flag of proper reading the IX array data
c     ng      global element nodes
c     conp    connecticity information of subelement in Omega+
c     conm    connecticity information of subelement in Omega-
c--------------------------------------------------------------------71
c
      IMPLICIT none

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
	include  'sdata.h' 
	include  'tdata.h'
	!include  'vecdata.h'
      include 'pvecdata.h'
      include  'iofile.h'

	logical sym, str_ell_flag
	integer cini, k, str_ell_case
      integer ix(nel), isw, imat, nhv
	integer i, j, l, imatd
      real*8  d(*), ul(ndf,nel,6), xl(ndm,nel),  s(nst,nst), p(nst)
	real*8  fi(3,3), finv(3,3), detF, sigv(6), sig(3,3), fpk(3,3)
      real*8  Tr(3), T0n, Tt(3), T0t, check, checkn(3), n0(3), normn
      real*8  nr(3), sig_pk, aas(6,6)
      real*8  test(3), mv(3), sig_max, temp
	
	cini  = int(hr(nh3+6))         ! crack_flag
       
	if (cini.eq.1) then ! element includes cohesive zone                                              
	  call kernel11Q1C1(nen1,ix,l,imat,imatd,nhv,d,ul,xl,
	1                            s,p,ndf,ndm,nst,isw,sym)

      else                  ! element does not include cohesive zone
        call kernel11Q1 (nen1,ix,l,imat,nhv,d,ul,xl,s,p,
     1                 ndf,ndm,nst,isw,sym, sigv,fi,aas)
        
        !Store average stress over the element in the history 
        do i= 1,6
          hr(nh3+37+i)  = sigv(i)    
        enddo   

        !Store average deformation gradient over the element in the history
        k= 0
        do i= 1,3
          do j= 1,3
            hr(nh3+44+k) = fi(i,j)    
            k = k+1
          enddo   
        enddo
        
        !If we're using the cohesive model from Miller and Gasser (2022) then we
        !need to determine the determinant of the acoustic tensor corresponding
        !to the referential principal fiber directions in the orientation density 
        !function. We use the average deviatoric stiffness of the element
        
        if (imatd.eq.8) then

          call str_ell(aas,fi)
            
        endif
        
      endif
      
c symmetrize tangent for application of symmetric solvers (loss of quadratic convergence!!!!)
      if (int(d(52)).eq.1) then
	  do i= 1, nst
	   do j= i+1,nst
	     s(i,j) = 0.5d0*(s(i,j)+s(j,i))
	     s(j,i) = s(i,j)
	   enddo
	  enddo
	endif

100   format(10x,'elem 62', e14.5,e14.5,e14.5,e14.5,e14.5,e14.5/)
1000  end


	subroutine kernel11Q1C1(nen1,ix,l,imat,imatd,nhv,d,ul,xl,
	1                            s,p,ndf,ndm,nst,isw,sym)
c
c--------------------------------------------------------------------71
c
c             Kernel of cohesive element formulation
c             based on  partition of unity
c
c
c   History variables	on element level
c
c       hr(nh3 + 0)	   A_e
c       hr(nh3 + 1)	   Ve   
c       hr(nh3 + 2)	   fact=Ve+/Ve   not used
c
c       hr(nh3 + 3)	   N0_1
c       hr(nh3 + 4)	   N0_2
c       hr(nh3 + 5)	   N0_3
c
c       hr(nh3 + 6)	   cini ... crack initialization flag
c								0  element intact 
c								1  element cracked 
c
c       hr(nh3 + 7)	   cn  ...  nodes of the disontinuity
c                                 3 ... triangle
c                                 4 ....4 node patch
c 
c       hr(nh3 + 8)	   conp
c       hr(nh3 + 9)	   conm
c
c       hr(nh3 + 10)	   L11
c       hr(nh3 + 11)	   L21
c       hr(nh3 + 12)	   L31
c       hr(nh3 + 13)	   L41
c
c       hr(nh3 + 14)	   L12
c       hr(nh3 + 15)	   L22
c       hr(nh3 + 16)	   L32
c       hr(nh3 + 17)	   L42
c
c       hr(nh3 + 18)	   L13
c       hr(nh3 + 19)	   L23
c       hr(nh3 + 20)	   L33
c       hr(nh3 + 21)	   L43
c
c       hr(nh3 + 22)	   L14
c       hr(nh3 + 23)	   L24
c       hr(nh3 + 24)	   L34
c       hr(nh3 + 25)	   L44
c
c       hr(nh3 + 26)	   Xb_1
c       hr(nh3 + 27)	   Xb_2
c       hr(nh3 + 28)	   Xb_3
c   
c       hr(nh3 + 30)	   Xn2d11
c       hr(nh3 + 31)	   Xn2d21
c
c       hr(nh3 + 32)	   Xn2d12
c       hr(nh3 + 33)	   Xn2d22
c
c       hr(nh3 + 34)	   Xn2d13
c       hr(nh3 + 35)	   Xn2d23
c
c       hr(nh3 + 36)	   Xn2d14
c       hr(nh3 + 37)	   Xn2d24
c
c       hr(nh3 + 38)	   sig_11
c       hr(nh3 + 39)	   sig_22
c       hr(nh3 + 40)	   sig_33
c       hr(nh3 + 41)	   sig_12
c       hr(nh3 + 42)	   sig_23
c       hr(nh3 + 43)	   sig_13
c
c       hr(nh3 + 44)	   fi_11
c       hr(nh3 + 45)	   fi_12
c       hr(nh3 + 46)	   fi_13
c       hr(nh3 + 47)	   fi_21
c       hr(nh3 + 48)	   fi_22
c       hr(nh3 + 49)	   fi_23
c       hr(nh3 + 50)	   fi_31
c       hr(nh3 + 51)	   fi_32
c       hr(nh3 + 52)	   fi_33
c
c       hr(nh3 + 53)	   nelem     ... number of neighbouring elements
c       hr(nh3 + 54)	   no. elem1 
c       hr(nh3 + 55)	   no. elem2 
c        ....
c       hr(nh3 + 354)	   no. elem300
c       hr(nh3 + 355)	   t0
c
c       hr(nh3 + 356)	   Nold_1
c       hr(nh3 + 357)	   Nold_1
c       hr(nh3 + 358)	   Nold_1
c   
c       hr(nh3 + 360)	   fi_res_11
c       hr(nh3 + 361)	   fi_res_12
c       hr(nh3 + 362)	   fi_res_13
c       hr(nh3 + 363)	   fi_res_21
c       hr(nh3 + 364)	   fi_res_22
c       hr(nh3 + 365)	   fi_res_23
c       hr(nh3 + 366)	   fi_res_31
c       hr(nh3 + 367)	   fi_res_32
c       hr(nh3 + 368)	   fi_res_33
c
c  history associated with time 
c
c       hr(nh1 + 0)      uH			! history variale uH = max|u0|
c
c
c  d(70) Continuum Material number  (imat)
c  d(71) order of bulk integration  (l)
c  d(72) Discrete Material number   (imatd)
c  d(73) order of disc. integration (ld)
c  d(74) History in bulk material   (nhv)
c  d(75) History in disc. material  (nhvd)
c
c     Displacement formulation  nn=0
c
c
c  Declare variable types
c
c       ix    Element connections
c       d     Material parameters
c       ul    Solution array (displacements)
c       ulc   Compatible displacements
c       ule   Enhanced displacements
c       dule  Increment of Enhanced displacements
c       xl    (Lagrange) coordinates
c       s     Element (stiffness) matrix
c       p     Element vector
c       ndf   Degree of freedom (max) per node
c       ndm   Space dimension of mesh 
c       nst   Size of element array
c       isw   Task parameter 
c       imat  Continuum Material model
c       imatd Discrete Material model
c       nhv   Length of history vector per gauss point
c       shp   Shape function and its derivatives
c       shpt  shape function and derivatives of triangular discontinuity
c       shpq  shape function and derivatives of 4 node discontinuity
c       n     number of current element
c       nel   nodes per element
c       xu    Current geometry
c       l     Order of Gauss integration
c       lint  Number of Gauss-points for volume integrals
c      lint2d Number of Gauss-points for surface integrals
c       sg    Gauss points  and weights 
c       xsj    determinant of the Jacobin 
c       dvol  current volume element
c       da    current area element
c       fi    Deformation gradient (based on compatible displacements)
c       finv  Inverse deformation gradient (based on compatible displacements)
c       xji   determinat the deformation gradient
c       nn    Shiftvariable for history  
c       sigv  Cauchy stresses in Quadraturpoins
c       aa    Eulerian material moduli
c       nh1   History pointer previous time
c       nh2   History pointer current time
c       dot   Dotproduct
c       nen   Maximum number of nodes on element (in common 'cdata')
c       rho   Mass per Gauss-point
c       numnp Number of nodes in problem (in 'cdata.h') 
c       nph   Pointer for nodal stress parameters (in 'prstrs.h')
c       sigl  Quantities at Gauss-points to be plotted
c       Ae    Area of the discontinuity (ref.conf.)
c       Ve    Element volume (ref.conf.)
c       N0    Normal vector of the discontinuity (ref.conf)
c       nb    normal vector at the fictitious discontinuity (curr. conf)
c       xid   Volume coordinates at the discontinuity corners
c       sgt   Gauss points and weights of triangular discontinuity
c       sgq   Gauss points and weights of 4 node discontinuity
c       fib   deformation gradient at the discontinuity
c       u0    enhanced displacement 
c       fact  fact= Ve+/Ve
c     conp    connecticity information of subelement in Omega+
c     conm    connecticity information of subelement in Omega-
c             1.     character .... number of nodes of subbody
c             2.-7.  character .... local node numbers of subbody
c
c   numbering convention:
c
c                 1   -
c				2    |  local tetrahedral nodes
c				3	 |
c				4	_|
c				5    |
c				6    |  nodes on discontinuity
c				7    |
c				8  (not defined for triangular disc)
c
c       nst2  no. of enhanced degrees of freedom
c       order element order     
c       xl2d  Lagrage coordinates of disc element   
c       sg2d  Gauss points and weights for 2d integration  
c
c--------------------------------------------------------------------71
c
      IMPLICIT none

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
      include  'cdata.h'
      include  'prstrs.h'
      include  'rdata.h'
      include  'iofile.h'
      include  'tdata.h'

	logical sym, hflag
      integer ix(nel), ndf, ndm, nst, isw, imat, nhv,  lint, ll, ltot
	integer i, j, l, k, nn, imatd, lint2d, cn, nst2, nen1, pom, nhvd
	integer idg, conm(7), conp(7), order
      real*8  d(*), ul(ndf,nel,6), xl(ndm,nel)
	real*8  shp(4,10), sg(5,256), sgrel(5,256), xsj
	real*8  dvol, da, ulc(ndm,nel), ule(ndm,nel), dule(ndm,nel)
	real*8  alc(ndm,nel), ale(ndm,nel)
	real*8  fi(ndm,ndm),  invfi(ndm,ndm), xji
	real*8  rho, sigl(20),mflag
	real*8  du0(3), N0(3), t(3), ddu(3,3)
	real*8  Ae, ddn(3,3), nc(3), ttemp(3)
	real*8  fact
      real*8  s(nst,nst), p(nst)
	real*8  Kb, nb(3), sgd(5,36)
	real*8  xid(4,4), shpd(4,10)
	real*8  ut(ndm,nel),u0(3)
	real*8  xt(ndm,nel),xc(ndm,nel)
      real*8  aa(6,6), sigv(9), jsub, aas(6,6)
	real*8  xlsubp(4,6), xlsubm(4,6), nodes(4,8)
	real*8  xl2d(2,4), sg2d(4,36), xsj2d
	real*8  ua(ndm,nel), xa(ndm,nel)
	real*8  fid(ndm,ndm), invfid(ndm,ndm), xjid, temp
      real*8  sigvtmp(3)

c	write(*,*) 'compute PU-Element'
	rho   = d(50)
	mflag = d(51)

c compute element data
      call eldata_11 (nel,ndm, nst2, order)

c get solution data
	call solution_11 (ul(1,1,1), ul(1,1,5), ul(1,1,2), ndf, ndm,  
	1                  nel, ulc, ule, alc, ale, dule)

c compute total displacement

      ut = ulc + ule              ! total displacement
	ua = ulc + 0.5d0*ule        ! average displacement
c ******************************diesen teil vielleicht auslagern
c Retrieve element- and directional data
      do i= 1,3
	  N0(i) = hr(nh3+i+2)		!ref. normal vector of disc.
	enddo
	Ae   = hr(nh3)				!ref. area of the disc.
	fact = hr(nh3+2)			!volume ratio Ve^+/Ve
      cn   = hr(nh3+7)			!number of nodes of disc.

c ... retrieve connectivity data
      call real_to_array_11(hr(nh3+8),conp) ! connectivity in omega+
      call real_to_array_11(hr(nh3+9),conm) ! connectivity in omega-
c ... retrieve volume coordinates of the cornes of the disc.
      k= 1
      do i= 1,cn
	  do j= 1,4
	    xid(j,i) = Hr(nh3+9+k)
	    k = k+1
	  enddo
	enddo
c ... retrieve 2d coordinates of the cornes of the disc.
      k= 1
      do i= 1,cn
	  do j= 1,2
	    xl2d(j,i) = Hr(nh3+29+k)
	    k = k+1
	  enddo
	enddo
c ... assemble volume subelements
c ... store volume coordinates of all nodes in nodes
c part a - element nodes
      nodes=0.0d0
      do i= 1,4
	  nodes(i,i) = 1.0d0
	enddo
c part b - disc nodes
      do i=1,4
	  do j= 5,8
	    nodes(i,j) = xid(i,j-4)
	  enddo
	enddo

c ...Omega+
	do i=1,4
        do j=1,conp(1)
	    xlsubp(i,j) = nodes(i,conp(j+1))
	  enddo
	enddo
c ...Omega-
	do i=1,4
        do j=1,conm(1)
	    xlsubm(i,j) = nodes(i,conm(j+1))
	  enddo
	enddo
c ******************************diesen teil vielleicht auslagern
c Compute current geometries
      do i = 1,3
        do j = 1,nel
         xt(i,j) = xl(i,j) + ut(i,j)   ! total eulerian coordinate
         xa(i,j) = xl(i,j) + ua(i,j)   ! average eulerian coordinate
         xc(i,j) = xl(i,j) + ulc(i,j)  ! compatible eulerian coordinate
        enddo
      enddo
      
      do ll= 1,2   ! loop over subelements
c retrieve integration points and weights of bulk 
	 lint=0
	 l=int(d(71))       ! order of bulk integration
       ltot = (l*l*l)   ! Number of gauss points, in history, belonging to a subelement...8 (CM 21/09/2022)
       if (ll.eq.1) then ! positive subelement
         call pu_sub_tint3d_11(l,conp(1),xlsubp, lint,sg,sgrel)
	   idg = 1
	 else
         call pu_sub_tint3d_11(l,conm(1),xlsubm, lint,sg,sgrel)
	   idg = -1
	 endif
      
       do l= 1,lint   ! loop over integration points of bulk 
c integration over Omega-
        if (idg.eq.-1) then !  Omega-
c compute determonant of subtransformation
	    call subJ_11 (conm(1),xlsubm,sg(1,l), sgrel(1,l), Jsub)
c Get shape functions and derivatives at compatible config 
          call tetshp_11( sg(1,l), xc, ndm, order, xsj, shp)
c compute compatible kinematics
          call kine_11(ndm, nel, shp, ulc, fi,invfi,xji)
          dvol = Jsub*abs(xsj)*sg(5,l)  !dvc
        elseif (idg.eq.1) then ! Omega+
c compute determonant of subtransformation
	    call subJ_11 (conp(1),xlsubp,sg(1,l), sgrel(1,l), Jsub)
c Get shape functions and derivatives at current config 
          call tetshp_11( sg(1,l), xt, ndm, order, xsj, shp)
c compute total kinematics
          call kine_11(ndm, nel, shp, ut, fi,invfi,xji)
          dvol = Jsub*abs(xsj)*sg(5,l)  !dv
        endif

c ======= BEGINN BULK RESPONSE ====================================70
c compute stress and tangent 
        sigv = 0.0d0
	  aa   = 0.0d0
        !call modc_11(.true.,nhv,l,n,imat,d,fi,invfi,xji,sigv,aa)
        pom = (ll-1)*ltot+l ! CM 21/09/2022
        ! pom is the pointer to the gauss point of interest 
        !call modc_11(.true.,nhv,pom,n,imat,d,fi,invfi,xji,sigv,aa,aas) !Damage turned off
        call modc_11(.false.,nhv,pom,n,imat,d,fi,invfi,xji,sigv,aa,aas)!Damage turned on
        if ((imat.eq.2).or.(imat.eq.3).or.(imat.eq.4).or.
     1      (imat.eq.5).or.(imat.eq.6).or.(imat.eq.9)) then ! add volumetric stress if needed
            call vol_11(d, xji, sigv, aa)
        endif

	  if ((isw.eq.3).or.(isw.eq.6)) then
c.... Compute the R.H.S and element stiffness
c
c Multiply tangent moduli and stress by volume element
c and integration factors
        do  i = 1,6
          sigv(i) = sigv(i)*dvol
          do  j = 1,6
            aa(i,j) = aa(i,j)*dvol
          enddo
        enddo

c compute current density multiplied by the reference volume

c set a limit for rho!! to limit the time step
c        rho = d(50)*dvol/xji
        rho = d(50)*max(dvol,	hr(nh3+1)/float(10*lint))/xji

c compute residual of bulk material
c.... Part 1: due to stress field
        call resi_stre3d_11(nst2,shp,sigv,idg, p)   
	    
c.... Part 2: due to inertia (needed for transient solution)
        call resi_mass3d_11(idg,nst2,alc,ale,mflag,rho,shp, p) 

c compute the tangent stiffness matrix if isw=3
        if(isw.eq.3) then
c
c.... Part 1: geometric stiffness
          call stiff_geom3d_11(nst2,shp,sigv,idg, s)     

c.... Part 2: tangent modulus
          call stiff_mate3d_11(nst2,sym,shp,aa,idg, s)
     
c.... Part 3: mass matrix (needed for transient solution)
          call stiff_mass3d_11(idg,nst2,mflag,rho,shp, s)     
        endif	 ! if((isw.eq.3)
	endif

c plot stresses
       If(isw.eq.8) then
c Load stress for plotting
         if((imat.eq.1).or.(imat.eq.2).or.(imat.eq.3).or.
     1   (imat.eq.4).or.(imat.eq.5).or.(imat.eq.6).or.(imat.eq.9)) then
c  principal stresses           
           call prinb10(sigv,sigv(7)) 
           do i=1,9
	       sigl(i) = sigv(i)
           enddo 
         endif

         if((imat.eq.7).or.(imat.eq.8)) then
           do i=1,6
	       sigl(i) = sigv(i)
           enddo
c max principal stress           
           call prinb10(sigv,sigl(7)) 
c von Mises stress
           sigl(8)= 0.7071067811*sqrt((sigv(1)-sigv(2))**2 +
     1                                (sigv(1)-sigv(3))**2 + 
     2                                (sigv(2)-sigv(3))**2 +
     3            6.0*(sigv(4)**2+ sigv(5)**2+ sigv(6)**2))
c cumulative plastic strain
           sigl(9) = sigv(7) 
         endif

        call stcn10tetra(ix,xt,dvol,shp,sigl,Hr(nph), Hr(nph+numnp)
     1                  ,ndm,nel,numnp)
	 endif 
	enddo ! loop over integration points of bulk    
	enddo ! loop over subelements  
c ======= END BULK RESPONSE =======================================70



c ======= BEGIN TRACTION RESPONSE =================================70

	if ((isw.eq.3).or.(isw.eq.6)) then
      l=int(d(73))    ! order of disc. integration
c ... get Gauss points at discontinuity
        if (cn.eq.3) then
c three noded discontinuity
c compute volume coordinates at disc. Gausspoints 
          call inttria_11(l,xid, lint2d, sgd, sg2d)
	  elseif (cn.eq.4) then
c   four noded discontinuity
c compute volume coordinates at disc. Gausspoints 
          call intquad_11(l,xid, lint2d, sgd, sg2d)
	  elseif(cn.eq.0) then
          write (*,*) '* * zero node discontinuity'
	  else
          write (*,*) '* * Undefined discontinuity'
	  endif
c loop over Gauss points on disc.
	  do l= 1,lint2d          
c compute shape funtions and derivative at Gauss points on the disc
c average deformation 
          call tetshp_11( sgd(1,l), xa, ndm, order, xsj, shpd)
c compute xsj2d 
          if (cn.eq.3) then
            call trishp_11(2,xl2d, xsj2d)  ! braechte nur einmal gemacht werden
	    elseif (cn.eq.4) then
            call quadshp_11(2,sg2d(1,l),xl2d, xsj2d)
	    endif
  
c compute average deformation gradient and inverse	
          call kine_11(ndm, nel, shpd, ua, fid,invfid,xjid)
c compute kinematical quantities at discontinuity	
          call kine_d11 (nel, fid,invfid,xjid,ule,dule,shpd, N0, 
	1                   u0,du0,Kb,nb)
c  call discrete mat model
          nhvd = int(d(75))
          call modd_11(l,nhvd,nb,u0,du0,imatd,d,t,ddu,ddn,hflag)
c define current area element
c          da = Kb*xsj2d*sgd(5,l) 
          da = xsj2d*sgd(5,l) 
	    do i= 1,ndm
	      t(i) = t(i)*da
		  do j= 1,ndm
	        ddu(i,j) = ddu(i,j)*da
	        ddn(i,j) = ddn(i,j)*da
		  enddo
		enddo 

c compute contribution of the traction to the residuum 
            call resi_trac3d_11(nst2,t,shpd, p)
c cohesive contribution
            call stiff_trac3d_11(nst2,shpd,ddu,ddn,t,nb, s)
        enddo
c ======= END TRACTION RESPONSE ===================================70
c      call plot_tangent11(s)
 	endif   ! if((isw.eq.3).or.(isw.eq.6) then

	if(isw.eq.4) then  !print quantities 
c	   call pr_quant (1,n,l,imat,0.0,Hr(nh2),hel,
c     1   fi,xji,sigv,Hr(nn+nh1-nhv),nhv,Hr(nn+nh2-nhv))
	endif


2000  format(i4,e12.5,e12.5,e12.5,e12.5,e12.5,e12.5)
1000  end



	subroutine kernel11Q1 (nen1,ix,l,imat,nhv,d,ul,xl,s,p,
     1                     ndf,ndm,nst,isw,sym, sig,fi,aas)
c
c--------------------------------------------------------------------71
c
c             Kernel of displacement formulation for
c             3-d tetrahedral element
c
c
c  d(70) Continuum Material number  (imat)
c  d(71) order of bulk integration  (l)
c  d(72) Discrete Material number   (imatd)
c  d(73) order of disc. integration (ld)
c  d(74) History in bulk material   (nhv)
c  d(75) History in disc. material  (nhvd)
c
c
c  Declare variable types
c
c       ix    Element connections
c       d     Material parameters
c       ul    Solution array (displacements)
c       xl    (Lagrange) coordinates
c       s     Element (stiffness) matrix
c       p     Element vector
c       ndf   Degree of freedom (max) per node
c       ndm   Space dimension of mesh 
c       nst   Size of element array
c       isw   Task parameter 
c       imat  Material model
c       nhv   Length of history vector per gauss point
c       shp   Shape function and its derivatives
c       n     number of current element
c       nel   nodes per element
c       xu    Current geometry
c       l     Gausspoint id
c       lint  Number of Gauss-points lint=l*l*l
c       sg    Gauss points (1-3) and weights (4)
c       xsj   Jacobi determinante
c       dvol  xsj*sg(4,l)
c       fi    Deformation gradient
c       finv  Inverse deformation gradient
c       bei   Left Cauchy Green strain
c       xji   Jacobian determinat at Gausspoints
c       nn    Shiftvariable for history  
c       sigv  Cauchy stresses in Quadraturpoins
c       aa    Eulerian material moduli
c       nh1   History pointer previous time
c       nh2   History pointer current time
c       augf  Augmented Lagrangian parameter 1
c       dot   Dotproduct
c       nen   Maximum number of nodes on element (in common 'cdata')
c       rho   Mass per Gauss-point
c       cmom  Change in momentum 
c       numnp Number of nodes in problem (in 'cdata.h') 
c       nph   Pointer for nodal stress parameters (in 'prstrs.h')
c       sigl  Quantities at Gauss-points to be plotted
c       nst2  no. of enhanced degrees of freedom
c       order element order        
c       ng    global element nodes
c       idd   node id
c       fiavg  average deformation gradient over the whole element
c       sigavg  average stress over the whole element
c       evol total element volume
c--------------------------------------------------------------------71
c
      IMPLICIT none

      include  'eldata.h'
      include  'comblk.h'
      include  'hdata.h'
      include  'cdata.h'
      include  'prstrs.h'
      include  'iofile.h'
      include  'nodes_id.h'
	include  'vecdata.h'
      include  'tdata.h'

	logical sym
      integer ix(nel), ndf, ndm, nst, isw, imat, nhv,  lint, nst2
	integer i, j, l, nn, order, nen1, k, nn1, nn2
      real*8  d(*), ul(ndf,nel,6), xl(ndm,nel),  s(nst,nst), p(nst)

	real*8  shp(4,10), sig_max, temp
	
	real*8 xu(ndm,nel), sg(5,87), xsj, dvol
	real*8  fi(ndm,ndm),  finv(ndm,ndm), xji
	real*8  aa(6,6), sigv(9), rho, sigl(20),mflag, aas(6,6)
	real*8  sig(6),ulc(ndm,nel),alc(ndm,nel)
	real*8  n0(3), sig_res(6)
	real*8  nc(3), ttemp(3), sig_max_old, sig_tmp(9)
      
      real*8 fiavg(ndm,ndm), sigavg(6), evol, aasavg(6,6)
      real*8 check,detaa, checkh(100)

	integer idd(nel), ng
	integer IXpoint, IXlen, IXpre, str_ell_case
      logical IXflag
      logical str_ell_flag



	rho=d(50)
	mflag=d(51)

	sig=0.0d0



cc read residual deformation from history
c        k= 0
c        do i= 1,3
c          do j= 1,3
c            fi_res(i,j)=hr(nh3+360+k)    
c  	      k = k+1
c    	    enddo   
c	  enddo   

c compute element data
      call eldata_11 (nel,ndm, nst2, order)


c read node of element
c  get pointer and length of the IX     
	   call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)

c retrieve element global nodes
	   if(ixflag) then
          do i= 1,nel
	      ng =  mr(ixpoint-1+i+(n-1)*nen1)
c get node id
            idd(i) = nid(ng)
	    enddo
	   endif

     
c get displacements
      do i= 1,nel
	 do j= 1,ndm
	   ulc(j,i) = ul(j,i,1)   !compatible displacements
	   if (idd(i).eq.1) then
	     ulc(j,i) = ulc(j,i)  +  ul(j+3,i,1) 
	   endif
	 enddo
	enddo

c get accelerations
      do i= 1,nel
	 do j= 1,ndm
	   alc(j,i) = ul(j,i,5)   !compatible acceleration
	   if (idd(i).eq.1) then
	     alc(j,i) = alc(j,i)  +  ul(j+3,i,5) 
	   endif
	 enddo
	enddo

c Compute current geometry
      do i = 1,nel
        do j = 1,ndm
         xu(j,i) = xl(j,i) + ulc(j,i)  
        enddo
	enddo
c initialize stress for crack criterion
      sig =0.0d0
      sig_max =0.0d0
      sig_max_old =0.0d0

	lint=0
	l=int(d(71))
      call tint3d_11(l,lint,sg)

      fiavg =0.0
       sigavg=0.0
       aasavg=0.0
       evol=0.0
       do  l = 1,lint  
c Get shape functions and derivatives at current time 
        call tetshp_11( sg(1,l), xu, ndm, order, xsj, shp)
        dvol = abs(xsj)*sg(5,l)
        evol = evol + dvol

c compute deformation measures
        call kine_11(ndm, nel, shp, ulc, fi,finv,xji)
        fiavg = fiavg + fi*dvol

        if (xji.lt.0d0) then
           check = 1d0
        endif
        
c Loop over Gauss points to compute (remaining) stiffness and R.H.S
c compute deviatoric response

        sigv = 0.0d0
	  aa   = 0.0d0
        !call modc_11(.true.,nhv,l,n,imat,d,fi,finv,xji, sigv,aa,aas)  !Damage turned off
        call modc_11(.false.,nhv,l,n,imat,d,fi,finv,xji, sigv,aa,aas)  !Damed turned on
        if ((imat.eq.2).or.(imat.eq.3).or.(imat.eq.4).or.
     1      (imat.eq.5).or.(imat.eq.6).or.(imat.eq.9)) then ! add volumetric stress if needed
     !!       if ((imat.eq.9).and.(.not.str_ell_flag)) then
     !!       ! Check to see if the loss of strong ellipticity condition has been
     !!       ! satisfied, and if so, the flag (logical) and case (integer) variables
     !!       ! then define the discontinuity. When the flag becomes true following
     !!       ! execution, we no longer check as there is no longer any further
     !!       ! need i.e. A gauss point has localised in the element.
     !!!!       call str_ell(sigv,aas,fi,finv,xji,str_ell_flag,str_ell_case
     !!!!1      ,detH)
     !!!!         sig_tmp = sigv    
     !!!!         call vol_11(d, xji, sig_tmp, aas)
     !!!!         call str_ell(sig_tmp,aas,fi,finv,xji,
     !!!!1                         str_ell_flag,str_ell_case)
     !!       endif 

            call vol_11(d, xji, sigv, aa)
            !call vol_11(d, xji, sig_tmp, aas)
           
        endif
      

	 if ((isw.eq.3).or.(isw.eq.6)) then
c.... Compute the R.H.S and element stiffness
c
c Multiply tangent moduli and stress by volume element
c and integration factors
        do  i = 1,6
          sigv(i) = sigv(i)*dvol
          sigavg(i) = sigavg(i) + sigv(i)
          do  j = 1,6
            aa(i,j) = aa(i,j)*dvol
            aas(i,j) = aas(i,j)*dvol
            aasavg(i,j) = aasavg(i,j) + aas(i,j)
          enddo
        enddo

c compute current density
       rho = d(50)*dvol/xji

c compute residua of bulk material
c.... Part 1: due to stress field
	 call resi_stre3dQ1_11(idd,nst2,shp,sigv, p) 
	    
c.... Part 2: due to inertia (needed for transient solution)
	 call resi_mass3dQ1_11(idd,nst2,alc,mflag,rho,shp, p) 

c compute the tangent stiffness matrix if isw=3
c
        if(isw.eq.3) then
c
c.... Part 1: geometric stiffness
        call stiff_geom3dQ1_11(idd,nst2,shp,sigv, s)     

c.... Part 2: tangent modulus
        call stiff_mate3dQ1_11(idd,nst2,sym,shp,aa, s)     

c.... Part 3: mass matrix (needed for transient solution)
        call stiff_mass3dQ1_11(idd,nst2,mflag,rho,shp, s)     
     
        endif	 ! if((isw.eq.3)
 	 endif   ! if((isw.eq.3).or.(isw.eq.6) then

c plot quantities
       If(isw.eq.8) then
c Load stress for plotting
         if((imat.eq.1).or.(imat.eq.2).or.(imat.eq.3).or.
     1   (imat.eq.4).or.(imat.eq.5).or.(imat.eq.6).or.(imat.eq.9)) then
c  principal stresses           
           call prinb10(sigv,sigv(7)) 
           do i=1,9
	       sigl(i) = sigv(i)
           enddo 
         endif

         if((imat.eq.7).or.(imat.eq.8)) then
           do i=1,6
	       sigl(i) = sigv(i)
           enddo
c max principal stress           
           call prinb10(sigv,sigl(7)) 
c von Mises stress
           sigl(8)= 0.7071067811*sqrt((sigv(1)-sigv(2))**2 +
     1                                (sigv(1)-sigv(3))**2 + 
     2                                (sigv(2)-sigv(3))**2 +
     3            6.0*(sigv(4)**2+ sigv(5)**2+ sigv(6)**2))
c cumulative plastic strain
           sigl(9) = sigv(7) 
         endif
         call stcn10tetra(ix,xu,dvol,shp,sigl,Hr(nph), Hr(nph+numnp)
     1                  ,ndm,nel,numnp)
	 endif
      enddo    ! Gauss loop

	 if(isw.eq.4) then  !print quantities 
c	  call pr_quant (1,n,l,imat,0.0,Hr(nh2),hel,
c     1   fi,xji,sigv,Hr(nn+nh1-nhv),nhv,Hr(nn+nh2-nhv))
      endif

c computer average element stress and deformaton gardient
      If (evol.ne.0.0) then  
          sig = sigavg/evol
          !fi = fi/evol
          fi = fiavg/evol
          aas = aasavg/evol
          check = 1d0
      else
	  write (*,*)
        write(*,*)'=================================='
        write(*,*)'ERROR --> zero element volume'
        write(*,*)'          failed to average'
        write(*,*)'=================================='
	  write (*,*)
      endif

c Move the history for the discontinuity (Only needed if its defined as non zero at the beginning) CJM - 18/10/2022
c Note that there is history for the bulk gauss points not used that is lost becuase it isn't moved.
      
      nn1 = nh1 + int(d(78))
      nn2 = nh2 + int(d(78))
      do i = 1,int(d(79))
        !checkh(i) = hr(nn1+i-1)  
        hr(nn2+i-1) = hr(nn1+i-1)
        hr(nn2+i-1) = hr(nn1+i-1)
      enddo

	end


