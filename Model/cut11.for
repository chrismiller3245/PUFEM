	subroutine cut_tetra_11(elmtnr,xl,N0,xbs,Ae,Ve,fact,Lcn,
	1					node_f,xn_f,id,k_f,rconp,rconm)
c---------------------------------------------------------------------71
c.... Compute Ae, Ve, fact, Lcn for tetraeder cutted by plane defined
c     due to N0 and xb
c
c  Declare variable types
c
c       Ae          A^e
c       Ve          V^e
c       fact        Ve+/Ve
c       N0          normal vector N to plane
c       xb          Coordinates of point on plane
c       xl          Lagrange Coordinates at nodes
c       fact        fact= Ve+/Ve
c       elmtnr      Element number
c       sg          natural coordinates at center (xi=0)
c       shp         shape functions and derivatives
c       xsj         det[J] at center (xi=0)
c       node        node at discontinuity
c       n           element number
c       k_f         No. of corners of disc.
c       conp		  connectivity of subelement in omega+
c       conm		  connectivity of subelement in omega-
c---------------------------------------------------------------------71
	implicit none
	logical flag, fflag, mflag
	integer i,j,  k, k_f, elmtnr
	integer ctemp(7),conp(7),conm(7)
	integer ii, id(4), sgn, it


	real*8 anl(2,3), N0(3), Ae, s(5), xl(3,4), shp(4,4)
	real*8 Ve, xsj, max, temp, xb(3), xi(3)
	real*8 dxx(3), dxb(3), AA(3,3), det
	real*8 node(2,6), invAA(3,3)
	real*8 node_f(2,4),eps, norm
	real*8 x3d_f(3,4),xi3d_f(3,4)
	real*8 tr_anl(3,3), xc(3)
	real*8 mind, tempv(3)

      real*8 xn(3,6), xn_f(3,4), xlsub(3,6)
	real*8 fact, Lcn(4,4), xbs(3)
!tcg 2018-02-18	real*8 jp(8), jm(8),rconp(7),rconm(7)
	real*8 jp(8), jm(8),rconp,rconm

      include  'comblk.h'
      include  'hdata.h'
      include  'eldata.h'
      include  'iofile.h'

      data s /0.25d0, 0.25d0, 0.25d0, 0.25d0, 0.1666666666666667d0/

c write message
         write(*,*)   'ACTIVATE PU Element No.',elmtnr
         write(iow,*) 'ACTIVATE PU Element No.',elmtnr

      xb = xbs
      mflag=.false.
c normalize N0
	norm = sqrt(N0(1)*N0(1) + N0(2)*N0(2) + N0(3)*N0(3))
      If(norm.GT.0.0) then
	  N0=N0/norm
      else
        write(*,*) 'Divison by zero in cut11'
        write(iow,*) 'Divison by zero in cut11'
      endif


c decide local coordinate
!c trial directions
!      tr_anl(1,1)= 0.0d0
!      tr_anl(2,1)= N0(3)
!      tr_anl(3,1)=-N0(2)
!
!
!      tr_anl(1,2)=-N0(3)
!      tr_anl(2,2)= 0.0d0
!      tr_anl(3,2)= N0(1)
!
!c      tr_anl(1,3)= N0(2)
!c      tr_anl(2,3)=-N0(1)
!c      tr_anl(3,3)= 0.0d0
!
!      
!      tr_anl(1,3)= 0
!      tr_anl(2,3)=-N0(3)
!      tr_anl(3,3)= N0(2)
!
!c filter out largest vector
!      max=0.0d0
!      do i= 1,3
!       temp = sqrt(tr_anl(1,i)*tr_anl(1,i)+tr_anl(2,i)*tr_anl(2,i)+
!	1          tr_anl(3,i)*tr_anl(3,i))
!	 if (max.lt.temp) then 
!	   max = temp
!	   k = i
!	 endif
!	enddo
!      do i= 1,3
!	  anl(1,i) = tr_anl(i,k)/max
!	enddo
      
c decide upon local coordinate
      If((abs(N0(1)).ge.abs(N0(2))).and.
     1   (abs(N0(1)).ge.abs(N0(3)))) then
        anl(1,1) = -N0(2)
        anl(1,2) =  N0(1)
        anl(1,3) =  0.0

      elseif((abs(N0(2)).ge.abs(N0(1))).and.
     1       (abs(N0(2)).ge.abs(N0(3)))) then
        anl(1,1) =  N0(2)
        anl(1,2) =  -N0(1)
        anl(1,3) =  0.0

      elseif((abs(N0(3)).ge.abs(N0(1))).and.
     1       (abs(N0(3)).ge.abs(N0(2)))) then
        anl(1,1) =  N0(3)
        anl(1,2) =  0.0
        anl(1,3) =-N0(1)
      endif
c normalize
      temp = sqrt(anl(1,1)**2+ anl(1,2)**2+ anl(1,3)**2)    
      if (temp.gt.0.0) then
          Do i=1,3
             anl(1,i) = anl(1,i)/temp
          enddo
      else
          write(*,*) 'Divison by zero in cut11'
          write(iow,*) 'Divison by zero in cut11'         
      endif

c compute second base vector N0 x anl
      anl(2,1) =  anl(1,3)*N0(2) - anl(1,2)*N0(3)
	anl(2,2) = -anl(1,3)*N0(1) + anl(1,1)*N0(3)
	anl(2,3) =  anl(1,2)*N0(1) - anl(1,1)*N0(2)
      
      
c compute Element volume in reference configuration
!tcg 20180-02-19      call tetshp(s, xl, 3, 1, xsj, shp)
      call tetshp_11(s, xl, 3, 1, xsj, shp)
      Ve = xsj*s(5)
c define eps
      eps = 1.0d-12*Ve**0.333333333d0
c area_11 of dissection Ae
10    flag=.false.
	k=0
	do i=1,6  ! edges
	 if(i.eq.1) then
	  do j=1,3
	    dxx(j)  = xl(j,2) - xl(j,1)
	    dxb(j) = xb(j)    - xl(j,1)
	  enddo 
	 elseif(i.eq.2) then
	  do j=1,3
	    dxx(j)  = xl(j,3) - xl(j,1)
	    dxb(j) = xb(j)    - xl(j,1)
	  enddo 
	 elseif(i.eq.3) then
	  do j=1,3
	    dxx(j)  = xl(j,4) - xl(j,1)
	    dxb(j) = xb(j)    - xl(j,1)
	  enddo 
	 elseif(i.eq.4) then
	  do j=1,3
	    dxx(j)  = xl(j,3) - xl(j,2)
	    dxb(j) = xb(j)    - xl(j,2)
	  enddo 
	 elseif(i.eq.5) then
	  do j=1,3
	    dxx(j)  = xl(j,4) - xl(j,2)
	    dxb(j) = xb(j)    - xl(j,2)
	  enddo 
	 elseif(i.eq.6) then
	  do j=1,3
	    dxx(j)  = xl(j,4) - xl(j,3)
	    dxb(j) = xb(j)    - xl(j,3)
	  enddo 
	 endif
c set up linear system of eq.
	 do j= 1,3
        AA(j,1) =  dxx(j)
	  AA(j,2) = -anl(1,j)
	  AA(j,3) = -anl(2,j)
	 enddo
c solve system of eq.
c compute det AA
       call det_AA_11(AA, det)

c	 if (abs(det).gt.1.0d-8) then
c        call invert(AA,3,3)
!tcg 2018-02-18       CALL DLINRG (3, AA, 3, invaa, 3)
         invAA=AA
c         call invert(invAA,3,3)
         call invert3(invAA, det)
      !tcg 2018-02-18

	   do j= 1,3
	     xi(j) = invAA(j,1)*dxb(1) + invAA(j,2)*dxb(2) + 
     1           invAA(j,3)*dxb(3)
	   enddo
c check point
         if((xi(1).GE.0.0d0).AND.(xi(1).LE.1.0d0)) then
	     flag=.true.
           k = k+1
c 2d-setting
	     node(1,k) = xi(2) 
	     node(2,k) = xi(3)
c 3d-setting
           do ii= 1,3
             xn(ii,k) = xb(ii) + xi(2)*anl(1,ii) + xi(3)*anl(2,ii)
           enddo
         endif
c       endif
       If(abs(det).lt.eps) then
	   write(*,*) 
     1		'CRITICAL WARNING ** cannot cut edge and plane'
	   write(iow,*) 
     1		'CRITICAL WARNING ** cannot cut edge and plane'
       endif
       enddo
c plausibility check
       if ((k.gt.4).or.(k.lt.3)) then
          write(*,*) 'ERROR ** unexpected solution in cut11'
          write(iow,*) 'ERROR ** unexpected solution in cut11'
      endif
       
       
c filter out identical nodes
      node_f=0.0d0
	xn_f = 0.0d0
	k_f=0
	do i=1,k
        fflag =.false.
	  do j= i+1,k
	   if ((abs(node(1,i)-node(1,j)).LT.eps).and.
	1       (abs(node(2,i)-node(2,j)).LT.eps)) fflag =.true.  
	  enddo
	  if(.not.fflag) then
	    k_f = k_f+1
	    if (k_f.gt.4) then
	      write(*,*) 'CRITICAL WARNING ** k_f too large in cut11'
	      write(iow,*) 'CRITICAL WARNING ** k_f too large in cut11'
	    endif
c 2d-setting
         do ii= 1,2
	    node_f(ii,k_f) = node(ii,i)
	   enddo
c 3d-setting
         do ii= 1,3
	    xn_f(ii,k_f) = xn(ii,i)
	   enddo          
	  endif
      enddo

 
c compute Ae and sort nodes (counterclockwise)
      if (k_f.eq.0) then
	 flag=.false. 
	else
	 flag=.true.
	endif

      if (not(flag)) then
        
c find closest node
         k=0
	   mind = 1.0d20
         do i= 1,4
	     do j= 1,3
	       tempv(j) = xl(j,i)-xb(j)
	     enddo
	     temp = tempv(1)*N0(1) + tempv(2)*N0(2) + tempv(3)*N0(3) 
	     if(abs(temp).lt.mind) then
	       k = i
	       mind = abs(temp)
	     endif
	   enddo
c compute center
         xc =0.0d0
         do i= 1,4
	     do j= 1,3
	       xc(j) = xc(j) + 0.25d0*xl(j,i)
	     enddo
	   enddo
c change xb
         do i= 1,3
	     xb(i)   = xl(i,k)
	    tempv(i) = xc(i) - xl(i,k)
	   enddo

c move xb into element
	   xb = xb + 1.0d-4*tempv
	   mflag = .true.

         write(*,*)   'WARNING --> discontinuity shift required'
         write(*,*)   'Old and new front point '
         write(iow,*)   'WARNING --> discontinuity shift required'
         write(iow,1001) (xbs(i),i=1,3), (xb(i),i=1,3)
	else
	  if (k_f.gt.2) then 
	    call area_11 (k_f,node_f,x3d_f, xi3d_f, Ae)
	  else
	    Ae   = 0.0d0
	    xn_f = 0.0d0
	  endif
c 3d-setting
         do i = 1, k_f
           do ii= 1,3
             xn_f(ii,i) = xb(ii) + node_f(1,i)*anl(1,ii) + 
	1                             node_f(2,i)*anl(2,ii)
	     enddo
	   enddo
         mflag = .false.

      endif

      if (mflag) goto 10

c compute volume coordinates of corner nodes
      Lcn = 0.0d0
	call volume_coord_11(Ve,xl,k_f,xn_f,Lcn)
c compute Ve+ and Ve-
c compute sign of nodes
      do i = 1,4
        temp = N0(1)*(xl(1,i)-xb(1)) +
	1         N0(2)*(xl(2,i)-xb(2)) + N0(3)*(xl(3,i)-xb(3))
	  if(temp.gt.0.0d0) then
	    id(i) = -1
	  else
	    id(i) = +1
	  endif
	enddo

	jp=0.0d0
	jm=0.0d0

      if (k_f.eq.3) then   ! triangular discontinuity
c assemble subelement = tetra and store connectivity in history
        xlsub =0.0d0
        call asseble_tetra_11(id, xn_f,xl,xlsub,ctemp,sgn)
c store connectivity
        if (sgn.eq.1) then 
	    conp = ctemp
	  else
	    conm = ctemp
	  endif
c assemble remaining subelement = prisma
        xlsub =0.0d0
	  it =-1*sgn
        call asseble_prisma_b_11(it,id, xn_f,xl,xlsub,ctemp)
c store connectivity
        if (it.eq.1) then 
	    conp = ctemp
	  else
	    conm = ctemp
	  endif

	elseif (k_f.eq.4) then   ! patch discontinuity
c assemble positiv subelement = prisma
        xlsub =0.0d0
         call asseble_prisma_11(+1,id, xn_f,xl,xlsub,conp)

c assemble negative subelement = prisma
        xlsub =0.0d0
        call asseble_prisma_11(-1,id, xn_f,xl,xlsub,conm)
	endif






c transform conp and conm into real type
      call array_to_real_11(conp,rconp)
      call array_to_real_11(conm,rconm)



1000  format(
     1  10x,'Tet data   '/
     2  10x,'node 1      ',e12.5,e12.5,e12.5/
     2  10x,'node 2      ',e12.5,e12.5,e12.5/
     2  10x,'node 3      ',e12.5,e12.5,e12.5/
     2  10x,'node 4      ',e12.5,e12.5,e12.5/)
1001  format(
     2  10x,'Plane data              '/
     2  10x,'node        ',e12.5,e12.5,e12.5/
     2  10x,'normal      ',e12.5,e12.5,e12.5/)
	end





      subroutine area_11 (nn,x,x3d,xi3d, Ae)
c--------------------------------------------------------------------71
c.... compute area_11 included by a convex polygon of n nodes
c
c  Declare variable types
c
c       nn         Number of nodes  
c       x          Coordinates of node in 2d
c       Ae         Included area_11  
c       xm         midpoint node
c       xi3d_f     nodes in a 3d natural setting
c       x3d_f      nodes in a 3d cartesian setting
c--------------------------------------------------------------------71
	implicit none
	logical flag
	integer nn,i,j
      real*8 x(2,4), Ae, xm(2), temp(2)
	real*8 r(nn), phi(nn)
	real*8 x3d(3,4), xi3d(3,4), temp3d(3)

c compute mitpoint node
      xm=0.0d0
      do i=1,nn
	 xm(1) = xm(1) + x(1,i)
	 xm(2) = xm(2) + x(2,i)
	enddo
	xm = xm/dfloat(nn)

c transver problem into polar coordinates
        do i = 1,nn
	    r(i) = sqrt((x(1,i)-xm(1))*(x(1,i)-xm(1))+
	1                (x(2,i)-xm(2))*(x(2,i)-xm(2)))
	    phi(i)=dacos((x(1,i)-xm(1))/r(i))
	    if((x(2,i)-xm(2)).LE.0.0d0) phi(i) = -phi(i)
	  enddo 

c sort nodes with repect to phi   
c      flag=.true.
c      do while (flag)
c	  flag=.false.
c 	  do i= 1,nn-1
c	   j = i+1
c	   if (phi(j).ge.phi(i)) then
c
c	    flag=.true.
c	    temp(1) = x(1,i) 
c	    temp(2) = x(2,i)
c		 
c  	    x(1,i)  = x(1,j)
c	    x(2,i)  = x(2,j)
c
c	    x(1,j)= temp(1) 
c	    x(2,j)= temp(2)
cc 3d cartesian setting
c
c		temp3d(1) = x3d(1,i) 
c		temp3d(2) = x3d(2,i) 
c		temp3d(3) = x3d(3,i) 
c
c	    x3d(1,i) = x3d(1,j)
c	    x3d(2,i) = x3d(2,j)
c	    x3d(3,i) = x3d(3,j)
c
c	    x3d(1,j) = temp3d(1)
c	    x3d(2,j) = temp3d(2)
c	    x3d(3,j) = temp3d(3)
c
cc 3d natural setting
c
c		temp3d(1) = xi3d(1,i) 
c		temp3d(2) = xi3d(2,i) 
c		temp3d(3) = xi3d(3,i) 
c
c	    xi3d(1,i) = xi3d(1,j)
c	    xi3d(2,i) = xi3d(2,j)
c	    xi3d(3,i) = xi3d(3,j)
c
c	    xi3d(1,j) = temp3d(1)
c	    xi3d(2,j) = temp3d(2)
c	    xi3d(3,j) = temp3d(3)
cc  polar
c
c	    temp(1) = phi(i) 
c	    temp(2) = r(i) 
c  
c  	    phi(i)  = phi(j)
c	    r(i)    = r(j)
c
c	    phi(j)= temp(1) 
c	    r(j)= temp(2)
c	   endif
c	  enddo
c	enddo




 	  do j= 1,nn
         do i =j+1,nn
	   if (phi(j).ge.phi(i)) then
	    temp(1) = x(1,i) 
	    temp(2) = x(2,i)
		 
  	    x(1,i)  = x(1,j)
	    x(2,i)  = x(2,j)

	    x(1,j)= temp(1) 
	    x(2,j)= temp(2)
c 3d cartesian setting

		temp3d(1) = x3d(1,i) 
		temp3d(2) = x3d(2,i) 
		temp3d(3) = x3d(3,i) 

	    x3d(1,i) = x3d(1,j)
	    x3d(2,i) = x3d(2,j)
	    x3d(3,i) = x3d(3,j)

	    x3d(1,j) = temp3d(1)
	    x3d(2,j) = temp3d(2)
	    x3d(3,j) = temp3d(3)

c 3d natural setting

		temp3d(1) = xi3d(1,i) 
		temp3d(2) = xi3d(2,i) 
		temp3d(3) = xi3d(3,i) 

	    xi3d(1,i) = xi3d(1,j)
	    xi3d(2,i) = xi3d(2,j)
	    xi3d(3,i) = xi3d(3,j)

	    xi3d(1,j) = temp3d(1)
	    xi3d(2,j) = temp3d(2)
	    xi3d(3,j) = temp3d(3)
c  polar

	    temp(1) = phi(i) 
	    temp(2) = r(i) 
  
  	    phi(i)  = phi(j)
	    r(i)    = r(j)

	    phi(j)= temp(1) 
	    r(j)= temp(2)
	   endif
	   enddo
	  enddo

	  Ae=0.0d0
        do i=1,nn-2
	   call area_11_tr (x(1,1),x(1,i+1),x(1,i+2),ae) 
	  enddo
	end

	subroutine area_11_tr (x,y,z,a) 
c--------------------------------------------------------------------71
c.... compute area_11 of a triangle x,y,z (2d)
c
c  Declare variable types
c
c       x,y,z Cordinates of nodes  
c       A     Included area_11  
c
c--------------------------------------------------------------------71
	implicit none
	real*8 x(2),y(2),z(2), a, ymx(2), zmx(2)

	ymx = y - x
	zmx = z - x
	 
	a = a + 0.5d0*Abs(ymx(1)*zmx(2) - ymx(2)*zmx(1)) 
      end







	subroutine volume_coord_11(Ve,x,k_f,xn,Lcn)
c--------------------------------------------------------------------71
c.... compute Volume coordinates of cornernodes
c
c  Declare variable types
c       Ve          volume of the element
c       x           Lagrange Coordinates at element nodes
c       k_f         Number of corner nodes
c       xn          Lagrange Coordinates at corners of the disc.
c       Lcn         Volume Coordinates at corners of the disc.
c
c   
c      |Lcn(1)|     | a(1)  b(1) c(1) d(1)|  |  1  |  
c      |Lcn(2)|     | a(2)  b(2) c(2) d(2)|  |xn(1)|     
c      |Lcn(3)|  =  | a(3)  b(3) c(3) d(3)|  |xn(2)|     
c      |Lcn(4)|     | a(4)  b(4) c(4) d(4)|  |xn(3)|     
c--------------------------------------------------------------------71
      implicit none

      integer i, j, k_f
	real*8 x(3,4), xn(3,4), Lcn(4,4)
	real*8 a(4), b(4), c(4), d(4), Ve, Ve6

c compute constant a_1, a_2, a_3, a_4
      a(1) = -x(1, 4)*x(2, 3)*x(3, 2) + x(1, 3)*x(2, 4)*x(3, 2) + 
	1		x(1, 4)*x(2, 2)*x(3, 3) - x(1, 2)*x(2, 4)*x(3, 3) - 
     2		x(1, 3)*x(2, 2)*x(3, 4) + x(1, 2)*x(2, 3)*x(3, 4)

	a(2) =  x(1, 4)*x(2, 3)*x(3, 1) - x(1, 3)*x(2, 4)*x(3, 1) - 
	1		x(1, 4)*x(2, 1)*x(3, 3) + x(1, 1)*x(2, 4)*x(3, 3) + 
     2		x(1, 3)*x(2, 1)*x(3, 4) - x(1, 1)*x(2, 3)*x(3, 4)

	a(3) = -x(1, 4)*x(2, 2)*x(3, 1) + x(1, 2)*x(2, 4)*x(3, 1) + 
	1		x(1, 4)*x(2, 1)*x(3, 2) - x(1, 1)*x(2, 4)*x(3, 2) - 
     2		x(1, 2)*x(2, 1)*x(3, 4) + x(1, 1)*x(2, 2)*x(3, 4)

	a(4) =  x(1, 3)*x(2, 2)*x(3, 1) - x(1, 2)*x(2, 3)*x(3, 1) - 
	1		x(1, 3)*x(2, 1)*x(3, 2) + x(1, 1)*x(2, 3)*x(3, 2) + 
     2		x(1, 2)*x(2, 1)*x(3, 3) - x(1, 1)*x(2, 2)*x(3, 3)

c compute constant b_1, b_2, b_3, b_4
      b(1) =  x(2, 4)*(-x(3, 2) + x(3, 3)) + 
	1		x(2, 3)*( x(3, 2) - x(3, 4)) + 
     2		x(2, 2)*(-x(3, 3) + x(3, 4))

	b(2) =  x(2, 4)*( x(3, 1) - x(3, 3)) + 
	1		x(2, 1)*( x(3, 3) - x(3, 4)) + 
     2		x(2, 3)*(-x(3, 1) + x(3, 4))

	b(3) =  x(2, 4)*(-x(3, 1) + x(3, 2)) + 
	1		x(2, 2)*( x(3, 1) - x(3, 4)) + 
     2		x(2, 1)*(-x(3, 2) + x(3, 4))

	b(4) =  x(2, 3)*( x(3, 1) - x(3, 2)) + 
	1		x(2, 1)*( x(3, 2) - x(3, 3)) + 
     2		x(2, 2)*(-x(3, 1) + x(3, 3))

c compute constant c_1, c_2, c_3, c_4
      c(1) =  x(1, 4)*( x(3, 2) - x(3, 3)) + 
	1		x(1, 2)*( x(3, 3) - x(3, 4)) + 
     2		x(1, 3)*(-x(3, 2) + x(3, 4))

	c(2) =  x(1, 4)*(-x(3, 1) + x(3, 3)) + 
	1		x(1, 3)*( x(3, 1) - x(3, 4)) + 
     2		x(1, 1)*(-x(3, 3) + x(3, 4))

	c(3) =  x(1, 4)*( x(3, 1) - x(3, 2)) + 
	1		x(1, 1)*( x(3, 2) - x(3, 4)) + 
     2		x(1, 2)*(-x(3, 1) + x(3, 4))

	c(4) =  x(1, 3)*(-x(3, 1) + x(3, 2)) + 
	1		x(1, 2)*( x(3, 1) - x(3, 3)) + 
     2		x(1, 1)*(-x(3, 2) + x(3, 3))

c compute constant d_1, d_2, d_3, d_4
      d(1) =  x(1, 4)*(-x(2, 2) + x(2, 3)) + 
	1		x(1, 3)*( x(2, 2) - x(2, 4)) + 
     2		x(1, 2)*(-x(2, 3) + x(2, 4))

      d(2) =  x(1, 4)*( x(2, 1) - x(2, 3)) + 
     1        x(1, 1)*( x(2, 3) - x(2, 4)) + 
     2		x(1, 3)*(-x(2, 1) + x(2, 4)) 

	d(3) =  x(1, 4)*(-x(2, 1) + x(2, 2)) + 
	1		x(1, 2)*( x(2, 1) - x(2, 4)) + 
     2		x(1, 1)*(-x(2, 2) + x(2, 4))

	d(4) =  x(1, 3)*( x(2, 1) - x(2, 2)) + 
	1		x(1, 1)*( x(2, 2) - x(2, 3)) + 
     2		x(1, 2)*(-x(2, 1) + x(2, 3))
c compute volume coordinates
      Ve6 = 6.0d0*Ve

      do i= 1,k_f
	  do j= 1,4
	    Lcn(j,i) = (a(j) + b(j)*xn(1,i) + c(j)*xn(2,i) + 
	1					   d(j)*xn(3,i))/Ve6
	  enddo
	enddo

      end