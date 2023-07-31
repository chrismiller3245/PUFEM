c$Id: eig3.f,v 1.1 1999/08/25 21:54:03 rlt Exp $
      subroutine eigen_11(v,d,rot)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-1999: Robert L. Taylor

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute eigenvalues/vectors for 3x3 symmetric matrix

c      Inputs:
c         v(3,3) - matrix with initial values (only upper half used)

c      Outputs:
c         v(3,3) - matrix of eigenvectors (by column)
c         d(3)   - eigenvalues associated with columns of v
c         rot    - number of rotations to diagonalize
c-----[--.----+----.----+----.-----------------------------------------]
c     Storage done as follows:

c       | v(1,1) v(1,2) v(1,3) |     |  d(1)  a(1)  a(3)  |
c       | v(2,1) v(2,2) v(2,3) |  =  |  a(1)  d(2)  a(2)  |
c       | v(3,1) v(3,2) v(3,3) |     |  a(3)  a(2)  d(3)  |

c        Transformations performed on d(i) and a(i) and v(i,j) become
c        eigenvectors.  Thus, original array is destroyed.

c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   rot, its, i,j,k
      real*8    g,h, aij, sm,thresh, t, c,s,tau
      real*8    v(3,3), d(3), a(3), b(3), z(3)

      save

c     Move array into one-d arrays

      a(1) = v(1,2)
      a(2) = v(2,3)
      a(3) = v(1,3)

      do i = 1,3
        d(i) = v(i,i)
        b(i) = d(i)
        z(i) = 0.0d0
        do j = 1,3
          v(i,j) = 0.0d0
        end do ! j
        v(i,i) = 1.d0
      end do ! i

      rot = 0

      do its = 1,50

c       Set convergence test and threshold

        sm = abs(a(1)) + abs(a(2)) + abs(a(3))
        if (sm.eq.0.0d0) return

        if(its.lt.4) then
          thresh = 0.011d0*sm
        else
          thresh = 0.0d0
        end if

c       Perform sweeps for rotations

        do i = 1,3
          j = mod(i,3) + 1
          k = mod(j,3) + 1

          aij  = a(i)
          g    = 100.d0*abs(aij)
          if(abs(d(i))+g.ne.abs(d(i)) .or.
     &       abs(d(j))+g.ne.abs(d(j))) then

            if(abs(aij).gt.thresh) then
              a(i) = 0.0d0
              h    = d(j) - d(i)
              if(abs(h)+g.eq.abs(h)) then
                t = aij/h
              else
                t = sign(2.d0,h/aij)/(abs(h/aij)+sqrt(4.d0+(h/aij)**2))
              endif

c             Set rotation parameters

              c    = 1.d0/sqrt(1.d0+t*t)
              s    = t*c
              tau  = s/(1.d0+c)

c             Rotate diagonal terms

              h    = t*aij
              z(i) = z(i) - h
              z(j) = z(j) + h
              d(i) = d(i) - h
              d(j) = d(j) + h

c             Rotate off-diagonal terms

              h    = a(j)
              g    = a(k)
              a(j) = h + s*(g - h*tau)
              a(k) = g - s*(h + g*tau)

c             Rotate eigenvectors

              do k = 1,3
                g      = v(k,i)
                h      = v(k,j)
                v(k,i) = g - s*(h + g*tau)
                v(k,j) = h + s*(g - h*tau)
              end do ! k

              rot = rot + 1

            end if
          else
            a(i) = 0.0d0
          end if
        end do ! i

c       Update diagonal terms

        do i = 1,3
          b(i) = b(i) + z(i)
          d(i) = b(i)
          z(i) = 0.0d0
        end do ! i

      end do ! its

      end


      subroutine max_lamb_11 (lamb,k)
c--------------------------------------------------TCG.06.08.2004----71
c
c       sort eigenvalues
c
c--------------------------------------------------------------------71
      implicit none

	integer k
	real*8 lamb(3)

	k=1   ! default value if lamb(1)=lamb(2)=lamb(3)

	if ((lamb(1).gt.lamb(2)).and.(lamb(1).gt.lamb(3))) k=1
	if ((lamb(2).gt.lamb(1)).and.(lamb(2).gt.lamb(3))) k=2
	if ((lamb(3).gt.lamb(1)).and.(lamb(3).gt.lamb(2))) k=3

	end



      subroutine tetshp_11( xi, xl, ndm, order, xsj, shp)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2002: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Compute 3-d tetrahedral element shape
c               functions and their derivatives w/r x,y,z

c      Inputs:
c         xi(4)     - Natural volume coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c         order     - Interpolation order: 1=linear; 2=quadratic

c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx
c                     shp(2,i) = dN_i/dy
c                     shp(3,i) = dN_i/dz
c                     shp(4,i) =  N_i
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

      include 'iofile.h'
      include 'eldata.h'

      integer  ndm, order, i, j, k, l
      real*8   xsj, detr, xi(4), xl(ndm,nel), shp(4,nel), a(4,3), xip(4)
      real*8   qsh(4,10), xxi(3,4)

c     Linear

      if(order.eq.1) then

c       Compute determinants for transformation minors

        do i = 1,4
          j      = mod(i,4) + 1
          k      = mod(j,4) + 1
          l      = mod(k,4) + 1
          a(i,1) = xl(2,j)*(xl(3,k) - xl(3,l))
     &           + xl(2,k)*(xl(3,l) - xl(3,j))
     &           + xl(2,l)*(xl(3,j) - xl(3,k))

          a(i,2) = xl(3,j)*(xl(1,k) - xl(1,l))
     &           + xl(3,k)*(xl(1,l) - xl(1,j))
     &           + xl(3,l)*(xl(1,j) - xl(1,k))

          a(i,3) = xl(1,j)*(xl(2,k) - xl(2,l))
     &           + xl(1,k)*(xl(2,l) - xl(2,j))
     &           + xl(1,l)*(xl(2,j) - xl(2,k))
          xip(i) = 256.d0*xi(j)*xi(k)*xi(l)
        end do ! i



c       Correct signs on determinants

        do i = 1,3
          a(1,i) = -a(1,i)
          a(3,i) = -a(3,i)
        end do ! i

c       Determinant for element volume

        xsj  = (xl(1,1)*a(1,1) + xl(1,2)*a(2,1)
     &        + xl(1,3)*a(3,1) + xl(1,4)*a(4,1))

        if(xsj.ne.0.0d0) then
          detr = 1.d0/xsj
        else
          write(iow,*) ' TETSHP: Determinant =',xsj
          detr = 1.d0
        endif

c       Linear and bubble mode shape functions

        shp(1,5) = 0.0d0
        shp(2,5) = 0.0d0
        shp(3,5) = 0.0d0
        shp(4,5) = 256.d0*xi(1)*xi(2)*xi(3)*xi(4)
        do i = 1,4
          shp(1,i) = a(i,1)*detr
          shp(2,i) = a(i,2)*detr
          shp(3,i) = a(i,3)*detr
          shp(4,i) = xi(i)

          shp(1,5) = shp(1,5) + shp(1,i)*xip(i)
          shp(2,5) = shp(2,5) + shp(2,i)*xip(i)
          shp(3,5) = shp(3,5) + shp(3,i)*xip(i)
        end do ! i

c     Quadratic

      elseif(order.eq.2) then

c       Shape functions and natural derivatives

        do j = 1,10
          do i = 1,4
            qsh(i,j) = 0.0d0
          end do ! i
        end do ! j

c       Vertex values

        do i = 1,4
          shp(4,i) = 2.d0*xi(i)*xi(i) - xi(i)
          qsh(i,i) = 4.d0*xi(i) - 1.d0
        end do ! i

c       Edge values

        do i = 1,3
          j          = mod(i,3) + 1
          shp(4,i+4) = 4.d0*xi(i)*xi(j)
          qsh(i,i+4) = 4.d0*xi(j)
          qsh(j,i+4) = 4.d0*xi(i)

          shp(4,i+7) = 4.d0*xi(i)*xi(4)
          qsh(i,i+7) = 4.d0*xi(4)
          qsh(4,i+7) = 4.d0*xi(i)
        end do ! i

c       Compute natural derivatives of coordinates

        do i = 1,3
          xxi(i,1) = xl(i,1)*qsh(1,1) + xl(i, 5)*qsh(1, 5)
     &             + xl(i,7)*qsh(1,7) + xl(i, 8)*qsh(1, 8)
          xxi(i,2) = xl(i,2)*qsh(2,2) + xl(i, 5)*qsh(2, 5)
     &             + xl(i,6)*qsh(2,6) + xl(i, 9)*qsh(2, 9)
          xxi(i,3) = xl(i,3)*qsh(3,3) + xl(i, 6)*qsh(3, 6)
     &             + xl(i,7)*qsh(3,7) + xl(i,10)*qsh(3,10)
          xxi(i,4) = xl(i,4)*qsh(4,4) + xl(i, 8)*qsh(4, 8)
     &             + xl(i,9)*qsh(4,9) + xl(i,10)*qsh(4,10)
        end do ! i




c       Compute determinants for transformation minors

        do i = 1,4
          j      = mod(i,4) + 1
          k      = mod(j,4) + 1
          l      = mod(k,4) + 1
          a(i,1) = xxi(2,j)*(xxi(3,k) - xxi(3,l))
     &           + xxi(2,k)*(xxi(3,l) - xxi(3,j))
     &           + xxi(2,l)*(xxi(3,j) - xxi(3,k))

          a(i,2) = xxi(3,j)*(xxi(1,k) - xxi(1,l))
     &           + xxi(3,k)*(xxi(1,l) - xxi(1,j))
     &           + xxi(3,l)*(xxi(1,j) - xxi(1,k))

          a(i,3) = xxi(1,j)*(xxi(2,k) - xxi(2,l))
     &           + xxi(1,k)*(xxi(2,l) - xxi(2,j))
     &           + xxi(1,l)*(xxi(2,j) - xxi(2,k))
          xip(i) = 256.d0*xi(j)*xi(k)*xi(l)
        end do ! i



c       Correct signs on determinants

        do i = 1,3
          a(1,i) = -a(1,i)
          a(3,i) = -a(3,i)
        end do ! i


c       Determinant for element volume


        xsj  = (xxi(1,1)*a(1,1) + xxi(1,2)*a(2,1)
     &        + xxi(1,3)*a(3,1) + xxi(1,4)*a(4,1))

        if(xsj.ne.0.0d0) then
          detr = 1.d0/xsj
        else
          write(iow,*) ' TETSHP: Determinant =',xsj
          detr = 1.d0
        endif

        do j = 1,3
          do i = 1,4
            a(i,j) = a(i,j)*detr
          end do ! i
        end do ! j



c       Compute global derivatives

        do j = 1,4
          do i = 1,3
            shp(i,j) = qsh(j,j)*a(j,i)
          end do ! i
        end do ! j

        do j = 1,3
          do i = 1,3
            k = mod(i,3) + 1
            shp(i,j+4) = qsh(i,i+4)*a(i,i) + qsh(j,j+4)*a(j,i)
            shp(i,j+7) = qsh(i,i+7)*a(i,i) + qsh(4,j+7)*a(4,i)
          end do ! i
        end do ! j

      do i= 1,10
	  do j= 1,3
	      shp(j,i) = qsh(1,i)*a(1,j) + qsh(2,i)*a(2,j) + 
	1                 qsh(3,i)*a(3,j) + qsh(4,i)*a(4,j) 
	  enddo
	enddo



c      xsj = -xsj

c     Error - Not coded

      else

        write(iow,2000) order
        write(  *,2000) order
        call plstop()

      endif

c     Format

2000  format(/' *ERROR* TRISHP not coded for order =',i4)

      end





      subroutine tint3d_11(ll,lint,s)

c--------------------------------------------------TCG.26.02.2003----71
c      Purpose: Gauss poionts and weights for 3-d tetrahedral element
c
c
c      Inputs:
c         ll       - Type of quadrature

c      Outputs:
c         lint     - Number of quadrature points
c         s(5,*)   - Values of volume coordinates and weights
c
c    ll=1      1 pt. Gauss quadrature O(1)
c    ll=2      4 pt. Gauss quadrature O(2)
c    ll=3      5 pt. Gauss quadrature O(3)
c    ll=4     11 pt. ? quadrature O(?)
c    ll=5     10 pt. Newton Cotes quadrature O(2)
c    ll=6      8 pt. Newton Cotes quadrature O(3)
c    ll=7     35 pt. Newton Cotes quadrature O(4)
c    ll=8     16 pt. Gauss quadrature O(5)
c    ll=11    87 pt. Gauss quadrature O(11)
c
c
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer ll, lint
      real*8  s(5,87)
	real*8  a, b, c, d, e, f, g, h, m
	real*8 	w0, a0, w11, a11, b11, w12, a12, b12, w13, a13
	real*8  b13, w14, a14, b14, w15, a15, b15, w21, a21, b21, w31 
	real*8  a31, b31, c31, w32, w35, a35, b35, c35
	real*8  a32, b32, c32, w33, a33, b33, c33, w34, a34, b34, c34

      save


      if(ll.eq.1) then ! 1 pt. quadrature O(h^1)
c set no. of integration points
        lint = 1

c define constants
	  a = 0.25d0
	  b = 0.1666666666666666d0

c load integration points
        s(1,1) = a
        s(2,1) = a
        s(3,1) = a
        s(4,1) = a
        s(5,1) = b


      elseif(ll.eq.2) then  ! 4 pt. quadrature O(h^2)
c set no. of integration points
        lint = 4

c define constants
	  a = 0.5854101966249658d0
	  b = 0.1381966011250105d0
	  c = 0.0416666666666666d0

c load integration points
        s(1,1:4) = (/ a, b, b, b /)
        s(2,1:4) = (/ b, a, b, b /)
        s(3,1:4) = (/ b, b, a, b /)
        s(4,1:4) = (/ b, b, b, a /)
        s(5,1:4) = (/ c, c, c, c /)


      elseif(ll.eq.3) then   ! 5 pt. quadrature O(h^3)
c set no. of integration points
        lint = 5

c define constants
        a = 0.25d0
	  b = -0.13333333333333d0
	  c = 0.5d0
	  d = 0.16666666666666d0
	  e = 0.075d0
c load integration points
        s(1,1:5) = (/ a, c, d, d, d /)
        s(2,1:5) = (/ a, d, c, d, d /)
        s(3,1:5) = (/ a, d, d, c, d /)
        s(4,1:5) = (/ a, d, d, d, c /)
        s(5,1:5) = (/ b, e, e, e, e /)


      elseif(ll.eq.4) then  ! 11 pt. quadrature O(?)
c set no. of integration points
        lint = 11
      
c define constants
        a = 0.25d0
	  b = -0.013155555555555d0
	  c = 0.071428571428571d0
	  d = 0.785714285714286d0
	  e = 0.007622222222222d0
	  f = 0.399403576166799d0
	  g = 0.100596423833201d0
	  h = 0.024888888888888d0

c load integration points
	  s(1,1:11) = (/a,	c,	c,	c,	d,	f,	f,	g,	g,	g,	f /)
	  s(2,1:11) = (/a,	c,	c,	d,	c,	f,	g,	g,	f,	f,	g /)
	  s(3,1:11) = (/a,	c,	d,	c,	c,	g,	f,	f,	f,	g,	g /)
	  s(4,1:11) = (/a,	d,	c,	c,	c,	g,	g,	f,	g,	f,	f /)
	  s(5,1:11) = (/b,	e,	e,	e,	e,	h,	h,	h,	h,	h,	h /)

      elseif(ll.eq.5) then   ! 10 point NEWTON COTES O(2)

c set no. of integration points
        lint = 10
      
c define constants
        a =  1.0d0
        b =  0.5d0
        c = -0.0083333333333333d0
        d =  0.0333333333333333d0
        e =  0.0d0
c load integration points
        s(1,1:10) =   (/ a, e, e, e, b, e, e, b, b, e /)
        s(2,1:10) =   (/ e, a, e, e, e, b, e, b, e, b /)
        s(3,1:10) =   (/ e, e, a, e, e, e, b, e, b, b /)
        s(4,1:10) =   (/ e, e, e, a, b, b, b, e, e, e /)
        s(5,1:10) =   (/ c, c, c, c, d, d, d, d, d, d /)


      elseif(ll.eq.6) then   ! 8 point NEWTON COTES O(3)

c set no. of integration points
        lint = 8
      
c define constants
        a = 1.0D+00
        b = 0.0041666666666667d0
        c = 0.3333333333333333d0
        d = 0.0375d0
        e = 0.0D+00

        s(1,1:8) =   (/ e, a, e, e, c, c, e, c /)
        s(2,1:8) =   (/ e, e, a, e, c, e, c, c /)
        s(3,1:8) =   (/ e, e, e, a, e, c, c, c /)
        s(3,1:8) =   (/ a, e, e, e, c, c, c, e /)
        s(5,1:8) =   (/ b, b, b, b, d, d, d, d /)


      elseif(ll.eq.7) then   ! 35 point NEWTON COTES O(4)

c set no. of integration points
        lint = 35

c define constants
        a =   0.25D+00
        b =   0.50D+00
        c =   0.75D+00
        d =   1.00D+00
        e =  -5.0D+00 / (6.0d0*420.0D+00)
        f = -12.0D+00 / (6.0d0*420.0D+00)
        g =  16.0D+00 / (6.0d0*420.0D+00)
        h = 128.0D+00 / (6.0d0*420.0D+00)
        m =   0.0D+00


        s(1,1:35) =   (/ m, d, m, m, a, m, m, c, c, c, m, a, m, m, a, 
	1				   m, b, m, m, b, b, m, a, b, a, a, b, m, b, m, 
     2				   a, a, m, a, a /)
        s(2,1:35) =   (/ m, m, d, m, m, a, m, m, a, m, c, c, c, m, m, 
	1                   a, m, b, m, b, m, b, a, a, b, m, m, a, a, b, 
     2                   b, m, a, a, a /)
        s(3,1:35) =   (/ m, m, m, d, m, m, a, m, m, a, m, m, a, c, c, 
	1                   c, m, m, b, m, b, b, m, m, m, a, a, a, a, a, 
     2                   a, b, b, b, a /)
        s(4,1:35) =   (/ d, m, m, m, c, c, c, a, m, m, a, m, m, a, m, 
	1                   m, b, b, b, m, m, m, b, a, a, b, a, b, m, a, 
     2                   m, a, a, m, a /)

        s(5,1:35) =   (/ e, e, e, e, g, g, g, g, g, g, g, g, g, g, g,  
	1                   g, f, f, f, f, f, f, g, g, g, g, g, g, g, g,  
     2                   g, g, g, g, h /)
  
      elseif(ll.eq.8) then   ! 16 point Gauss O(5)
        lint = 16
        a = 0.7611903264425430d-01
	  b = 0.7716429020672371d+00
	  c = 0.8395632516687135d-02

	  d = 0.1197005277978019d+00
	  e = 0.7183164526766925d-01
	  f = 0.4042339134672644d+00
	  g = 0.1109034477221540d-01 


        s(1,1:16)=(/ b, a, a, a, d, f, f, e, d, f, e, f, d, e, f, f/)  
        s(2,1:16)=(/ a, b, a, a, e, d, f, f, f, d, f, e, f, d, e, f/)  
        s(3,1:16)=(/ a, a, b, a, f, e, d, f, e, f, d, f, f, f, d, e/)  
        s(4,1:16)=(/ a, a, a, b, f, f, e, d, f, e, f, d, e, f, f, d/)  
        s(5,1:16)=(/ c, c, c, c, g, g, g, g, g, g, g, g, g, g, g, g/)  


      elseif(ll.eq.11) then   ! 87 point Gauss O(11)
        lint = 87
	  w0  = -0.3229059250896649271190741d0/6.0d0
	  w11 =-0.3831136086645949490922834d0/6.0d0
	  w12 = 0.1259876832639002206343925d0/6.0d0
	  w13 = 0.7772656110490364391225834d-2/6.0d0
	  w14 = 0.4475842042017354585386039d-5/6.0d0
	  w15 = 0.3076630972851224623855561d-1/6.0d0
	  w21 = 0.2230322290225118827491926d-1/6.0d0
	  w31 = 0.5167456484634155355210072d-3/6.0d0
	  w32 = 0.1484538986489890388750796d0/6.0d0
	  w33 = 0.9330967352789100478689519d-3/6.0d0
	  w34 = 0.9319130804165715418184301d-2/6.0d0
	  w35 = 0.1272850504266610340365066d-1/6.0d0

	  a0  = 0.25d0

	  a11 = 0.3197881306061907190476732d0
	  b11 = 0.4063561442097275838489758d-1

	  a12 = 0.2745875432484354948800959d0
	  b12 = 0.1762373724529911216791656d0

	  a13 = 0.4902463231623282344872756d-1
	  b13 = 0.8529260850827039660680600d0

	  a14 =-0.5889205032331655032726359d-1
	  b14 = 0.1176676123352849073840473d+1

	  a15 = 0.1436980650803076374165338d0
	  b15 = 0.5689057952549437662582103d0

	  a21 = 0.4334059320676971741712688d0
	  b21 = 0.6659406793230282582873115d-1

	  a31 = 0.5031834294032451116078467d0
	  b31 =-0.6735432578129541713779964d-1
	  c31 = 0.6098746697480519392210609d-1

	  a32 = 0.2944561694949265009906537d0
	  b32 = 0.3926655492603751889759375d-1
	  c32 = 0.3718211060841094791210986d0

	  a33 = 0.0d0
	  b33 = 0.1383898530902673690519883d0
	  c33 = 0.8616101469097326309480116d0

	  a34 = 0.1455031635850380766054879d0
	  b34 = 0.6926735250835180200459205d0
	  c34 = 0.1632014774640582674310343d-1

	  a35 = 0.4385453179269500725262452d-1
	  b35 = 0.2775959971470881592271749d0
	  c35 = 0.6346949392675218262675760d0


	  s(1,1:87)=(/ a0,b11,a11,a11,a11,b12,a12,a12,a12,b13,a13,a13,
	1  a13,b14,a14,a14,a14,b15,a15,a15,a15,a21,a21,b21,b21,b21,a21,
     2  a31,a31,a31,a31,a31,a31,b31,b31,b31,c31,c31,c31,a32,a32,a32,
     3  a32,a32,a32,b32,b32,b32,c32,c32,c32,a33,a33,a33,a33,a33,a33,
     4  b33,b33,b33,c33,c33,c33,a34,a34,a34,a34,a34,a34,b34,b34,b34,
     5  c34,c34,c34,a35,a35,a35,a35,a35,a35,b35,b35,b35,c35,c35,c35/)

	  s(2,1:87)=(/ a0,a11,b11,a11,a11,a12,b12,a12,a12,a13,b13,a13,
	1  a13,a14,b14,a14,a14,a15,b15,a15,a15,a21,b21,b21,a21,a21,b21,
     2  b31,b31,c31,c31,a31,a31,a31,a31,c31,a31,a31,b31,b32,b32,c32,
     3  c32,a32,a32,a32,a32,c32,a32,a32,b32,b33,b33,c33,c33,a33,a33,
     4  a33,a33,c33,a33,a33,b33,b34,b34,c34,c34,a34,a34,a34,a34,c34,
     5  a34,a34,b34,b35,b35,c35,c35,a35,a35,a35,a35,c35,a35,a35,b35/)

	  s(3,1:87)=(/ a0,a11,a11,b11,a11,a12,a12,b12,a12,a13,a13,b13,
	1  a13,a14,a14,b14,a14,a15,a15,b15,a15,b21,a21,a21,a21,b21,b21,
     2  c31,a31,b31,a31,b31,c31,c31,a31,a31,b31,a31,a31,c32,a32,b32,
     3  a32,b32,c32,c32,a32,a32,b32,a32,a32,c33,a33,b33,a33,b33,c33,
     4  c33,a33,a33,b33,a33,a33,c34,a34,b34,a34,b34,c34,c34,a34,a34,
     5  b34,a34,a34,c35,a35,b35,a35,b35,c35,c35,a35,a35,b35,a35,a35/)

	  s(4,1:87)=(/ a0,a11,a11,a11,b11,a12,a12,a12,b12,a13,a13,a13,
     1  b13,a14,a14,a14,b14,a15,a15,a15,b15,b21,b21,a21,b21,a21,a21,
     2  a31,c31,a31,b31,c31,b31,a31,c31,a31,a31,b31,a31,a32,c32,a32,
     3  b32,c32,b32,a32,c32,a32,a32,b32,a32,a33,c33,a33,b33,c33,b33,
     4  a33,c33,a33,a33,b33,a33,a34,c34,a34,b34,c34,b34,a34,c34,a34,
     5  a34,b34,a34,a35,c35,a35,b35,c35,b35,a35,c35,a35,a35,b35,a35/)

	  s(5,1:87)=(/ w0,w11,w11,w11,w11,w12,w12,w12,w12,w13,w13,w13,
	1  w13,w14,w14,w14,w14,w15,w15,w15,w15,w21,w21,w21,w21,w21,w21,
     2  w31,w31,w31,w31,w31,w31,w31,w31,w31,w31,w31,w31,w32,w32,w32,
     3  w32,w32,w32,w32,w32,w32,w32,w32,w32,w33,w33,w33,w33,w33,w33,
     4  w33,w33,w33,w33,w33,w33,w34,w34,w34,w34,w34,w34,w34,w34,w34,
     5  w34,w34,w34,w35,w35,w35,w35,w35,w35,w35,w35,w35,w35,w35,w35/)

	else
      endif


      end





      subroutine tint3d_11_sav(ll,lint,s)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-1999: Robert L. Taylor

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Gauss quadrature for 3-d tetrahedral element

c      Inputs:
c         ll       - Type of quadrature

c      Outputs:
c         lint     - Number of quadrature points
c         s(5,*)   - Values of volume coordinates and weights
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer   i, j, ll, lint
      real*8    s(5,*)

      save

c     1 pt. quadrature O(h^2)

      if(ll.eq.1) then
        lint = 1
        do i = 1,4
          s(i,1) = 0.25d0
        end do ! i
        s(5,1) = 1.0d0/6.d0

c     4 pt. quadrature O(h^3)

      elseif(ll.eq.2) then
        lint = 4
        s(5,4) = 0.25d0/6.d0
        do i = 1,4
          do j = 1,4
            s(i,j) = 0.1381966011250105d+00
          end do ! j
          s(i,i) = 0.5854101966249658d+00
          s(5,i) = s(5,4)
        end do ! i

c     11 pt. quadrature O(h^4)

      elseif(ll.eq.3) then
        lint = 11
        do i = 1,3
          do j = 1,10
            s(i,j) = 0.0d0
          end do ! j
          s(i, i  ) = 1.00d0
          s(i, i+4) = 0.50d0
          s(i, i+7) = 0.50d0
          s(i, 11 ) = 0.25d0
        end do ! i
        s(2, 5) = 0.50d0
        s(3, 6) = 0.50d0
        s(1,10) = 0.50d0
        do j = 1,4
          s(5,j) = 1.d0/360.d0
        end do ! j
        do j = 5,10
          s(5,j) = 1.d0/90.d0
        end do ! j
        s(5,11) = 4.d0/45.d0
      else

c     16 pt. quadrature O(h^5)

        lint = 16
        s(5,4) = 0.8395632516687135d-02
        do i = 1,3
          do j = 1,4
            s(i,j) = 0.7611903264425430d-01
          end do ! j
          s(i,i) = 0.7716429020672371d+00
          s(5,i) = s(5,4)
        end do ! i
        do i = 5,16
          s(5,i) = 0.1109034477221540d-01
        end do ! i

        s(1, 5) = 0.1197005277978019d+00
        s(2, 5) = 0.7183164526766925d-01
        s(3, 5) = 0.4042339134672644d+00
        s(1, 6) = 0.4042339134672644d+00
        s(2, 6) = 0.1197005277978019d+00
        s(3, 6) = 0.7183164526766925d-01
        s(1, 7) = 0.4042339134672644d+00
        s(2, 7) = 0.4042339134672644d+00
        s(3, 7) = 0.1197005277978019d+00
        s(1, 8) = 0.7183164526766925d-01
        s(2, 8) = 0.4042339134672644d+00
        s(3, 8) = 0.4042339134672644d+00

        s(1, 9) = 0.1197005277978019d+00
        s(2, 9) = 0.4042339134672644d+00
        s(3, 9) = 0.7183164526766925d-01
        s(1,10) = 0.4042339134672644d+00
        s(2,10) = 0.1197005277978019d+00
        s(3,10) = 0.4042339134672644d+00
        s(1,11) = 0.7183164526766925d-01
        s(2,11) = 0.4042339134672644d+00
        s(3,11) = 0.1197005277978019d+00
        s(1,12) = 0.4042339134672644d+00
        s(2,12) = 0.7183164526766925d-01
        s(3,12) = 0.4042339134672644d+00

        s(1,13) = 0.1197005277978019d+00
        s(2,13) = 0.4042339134672644d+00
        s(3,13) = 0.4042339134672644d+00
        s(1,14) = 0.7183164526766925d-01
        s(2,14) = 0.1197005277978019d+00
        s(3,14) = 0.4042339134672644d+00
        s(1,15) = 0.4042339134672644d+00
        s(2,15) = 0.7183164526766925d-01
        s(3,15) = 0.1197005277978019d+00
        s(1,16) = 0.4042339134672644d+00
        s(2,16) = 0.4042339134672644d+00
        s(3,16) = 0.7183164526766925d-01

      endif

c     Compute fourth points

      do j = 1,lint
        s(4,j) = 1.d0 - (s(1,j) + s(2,j) + s(3,j))
      end do ! j

      end




      subroutine pu_sub_tint3d_11(ll,nsub,xlsub, lint,s,srel)


c-----[--.----+----.----+----.----------------------TCG 05.03.2003----]
c      generates integration points for subelement 
c      Inputs:
c         ll        Order of Gauss quadrature
c         nsub      Nodes of SubElementflag
c                   4  ... tetrahedral subelement
c                   6  ... prismatic subelement
c         xlsub     volumetric coordinates of subelement
c         

c      Outputs:
c         lint      Number of quadrature points
c         s         volume coordinates and weights of Gauss points
c         srel      relative volume coordinates and weights of Gauss points
c-----[--.----+----.----+----.-----------------------------------------]

      implicit  none

      integer ll, nsub, i, j, l, k, lint

	real*8 xlsub(4,6), s(5,256), srel(5,256), sgtet(5,256)
	real*8 sgpris(4,256), shppris(6,256)
   
      save


      if (nsub.eq.4) then    ! tetrahedral subelement

c get rel Gausspoints of tet
        lint=0
        call tint3d_11(ll,lint,sgtet)
c transform sgtet in volume coordinates of global tet
        s = 0.0d0
        do i= 1,lint			! loop over gausspoints
	    do j= 1,4			! loop over volume coordinates
	      do k= 1,nsub		! loop over subelement nodes
	        s(j,i) = s(j,i) + xlsub(j,k)*sgtet(k,i)
	       enddo
	    enddo
c weights
          s(5,i) = sgtet(5,i)
 	  enddo
c store rel Gausspoints
        do l= 1,lint
	    do j= 1,5
	      srel(j,l) = sgtet(j,l)
	    enddo
	  enddo
	elseif  (nsub.eq.6) then    ! prismatic subelement

c get Gausspoints of prism
        lint=0
        call gauss_prism_11(ll,lint,sgpris)
	  do l= 1,lint
          call shp_prism_11m (sgpris(1,l), shppris(1,l))
	  enddo
      
c transform sgpris in volume coordinates of global tet
        s=0.0d0
        do i= 1,lint        ! loop over gausspoints
	    do j= 1,4         ! loop over volume coordinates
	      do k= 1,nsub    ! loop over subelement nodes
	        s(j,i) = s(j,i) + xlsub(j,k)*shppris(k,i)
	       enddo
	    enddo
c weights
          s(5,i)=sgpris(4,i)
   	  enddo
c store rel Gausspoints
        do l= 1,lint
	    do j= 1,4
	      srel(j,l) = sgpris(j,l)
	    enddo
	  enddo
	else
	  write (1000,*)
	endif



1000  format('  *ERROR* Illegal subelement')
      end









      subroutine kine_11(ndm,nel, shp,ul, fi,finv,det)
c
c--------------------------------------------------------------------71
c
c             Compute 3D total deformation gradient
c
c.... INPUT variables
c         lint    Number of integration points
c         ndm     No. of dimension
c         nel     Number of nodes per element
c         shp     Isoparametric shape functions (4-Node only)
c         ul      Array of total and relative displacements
c                 ul(ndf,1)     Total
c                 ul(ndf,1+nen) Relative
c
c.... OUTPUT variables
c         fi      Total deformation gradient
c         finv    Inverse deformation gradient
c         det     Total Jacobian determinant
c
c--------------------------------------------------------------------71
c
c..... Declare variable types
      integer nel, k, ndm
c..... Declare array types
      real*8 shp(4,nel),ul(ndm,nel),fi(9),finv(9),det
c..... Intrinsic
c
c.... Compute [minus] the current displacement gradient 
c
      call pzero(finv,9)
      call pzero(fi  ,9)
      do 110 k = 1,nel
           finv(1) = finv(1) - ul(1,k)*shp(1,k)
           finv(2) = finv(2) - ul(2,k)*shp(1,k)
           finv(3) = finv(3) - ul(3,k)*shp(1,k)
           finv(4) = finv(4) - ul(1,k)*shp(2,k)
           finv(5) = finv(5) - ul(2,k)*shp(2,k)
           finv(6) = finv(6) - ul(3,k)*shp(2,k)
           finv(7) = finv(7) - ul(1,k)*shp(3,k)
           finv(8) = finv(8) - ul(2,k)*shp(3,k)
           finv(9) = finv(9) - ul(3,k)*shp(3,k)
 110  continue
c
c.... Compute deformation gradient [F**-1 = I - grad u]
c
        finv(1) = 1.d0 + finv(1)
        finv(5) = 1.d0 + finv(5)
        finv(9) = 1.d0 + finv(9) 
c
c.... Compute Jacobian determinant and deformation gradient
c     [F = ( F**-1)**-1]
c
        det = 1.d0/(finv(1)*finv(5)*finv(9)
     1  +finv(4)*finv(8)*finv(3)+finv(7)*finv(2)*finv(6)
     2  -finv(3)*finv(5)*finv(7)-finv(6)*finv(8)*finv(1)
     3  -finv(9)*finv(2)*finv(4))
c
        fi(1) = (finv(5)*finv(9)-finv(6)*finv(8))*det
        fi(4) =-(finv(4)*finv(9)-finv(6)*finv(7))*det
        fi(7) = (finv(4)*finv(8)-finv(5)*finv(7))*det
        fi(2) =-(finv(2)*finv(9)-finv(3)*finv(8))*det
        fi(5) = (finv(1)*finv(9)-finv(3)*finv(7))*det
        fi(8) =-(finv(1)*finv(8)-finv(2)*finv(7))*det
        fi(3) = (finv(2)*finv(6)-finv(3)*finv(5))*det
        fi(6) =-(finv(1)*finv(6)-finv(3)*finv(4))*det
        fi(9) = (finv(1)*finv(5)-finv(2)*finv(4))*det
c
      end




      subroutine resi_stre3d_11(nst2,shp,sigv,idg,p)     
c--------------------------------------------------------------------71
c.... Compute the continuum residuum due to stress field  for
c     partition of unity finite element
c
c
c  Declare variable types
c
c       nst2    enhanced degrees of freedom
c       idg     id of integration point 
c                +1 point is in omega+
c                -1 point is in omega-
c       sigv    Cauchy stress	times dv
c       shp     shape functions and its deerivatives 
c               (spatial formulation)
c       p       residuum consists of compatible and enhanced loads
c
c--------------------------------------------------------------------71

	Implicit none

      include  'cdata.h'
      include  'eldata.h'
      include  'sdata.h'

	integer i,j, j1, j0, idg, nst2
	real*8 shp(4,nel), sigv(6)
	real*8 p(nst)

      real*8 a(nst2)

c compatible part

c
c.... Compute the residual due to stress field
c
      j1 = 1
      do j = 1,nel
       a(j1)   = - shp(1,j)*sigv(1)-shp(2,j)*sigv(4)
     1              - shp(3,j)*sigv(6)
       a(j1+1) = - shp(1,j)*sigv(4)-shp(2,j)*sigv(2)
     1              - shp(3,j)*sigv(5)
       a(j1+2) = - shp(1,j)*sigv(6)-shp(2,j)*sigv(5)
     1              - shp(3,j)*sigv(3)

       j1 = j1 + 3
      enddo


c store in right format

      j1=0
	j0=0
      do i= 1,nel
	 do j= 1,ndm
	  if (idg.eq.-1) then     ! Omega-
	    p(j1+j)  = p(j1+j)   +  a(j0+j) 
	  elseif (idg.eq.1) then  ! Omega+
	    p(j1+j)  = p(j1+j)   + a(j0+j) 
	    p(j1+j+3)= p(j1+j+3) + a(j0+j) 
	  endif
	 enddo
	 j0 = j0 + ndm
	 j1 = j1 + ndf
	enddo

	end




      subroutine resi_mass3d_11(idg,nst2,ulc,ule,mflag,rho,shp, p)     
c--------------------------------------------------------------------71
c.... Compute the residuum due to inertia parts for
c     displacement finite element and store in special format
c
c
c  Declare variable types
c
c       idg     id of integration point 
c                +1 point is in omega+
c                -1 point is in omega-
c       nst2    enhanced degrees of freedom
c       shp     shape functions and its deerivatives 
c               (spatial formulation)
c       rho     current density times current volume =
c               ref. density times ref. volume = mass
c       mflag   Massinterpolation flag
c                mflag = 0 ...... lumped mass interpolation
c                mflag = 1 ...... variationally consitent mass interpolation
c       ulc     compaible accelerations
c       ule     enhanced accelerations
c       p       residuum consists of compatible and enhanced loads
c
c--------------------------------------------------------------------71

	Implicit none

      include  'cdata.h'
      include  'eldata.h'
      include  'sdata.h'

	integer i,j, j1, j0, nst2,idg
	real*8 rho, shp(4,nel)
	real*8 cmomc(3),cmome(3), p(nst)
	real*8 mflag, rhol, rhom
	real*8 ulc(3,nel), ule(3,nel)

      real*8 bc(nst2), be(nst2)

     

	rhol = rho
	rhom = rhol*mflag
	rhol = rhol- rhom
	do i = 1,3
	 cmomc(i) = 0.0d0
	 cmome(i) = 0.0d0
	 do j = 1,nel
	   cmomc(i) = cmomc(i) + shp(4,j)*ulc(i,j)
	   cmome(i) = cmome(i) + shp(4,j)*ule(i,j)
       enddo
	 cmomc(i) = rhom*cmomc(i)
	 cmome(i) = rhom*cmome(i)
      enddo
c
c.... Compute the residual due to inertia parts
c
      j1 = 1
      do j = 1,nel

	 bc(j1)   = - shp(4,j)*(cmomc(1) + rhol*ulc(1,j))
	 bc(j1+1) = - shp(4,j)*(cmomc(2) + rhol*ulc(2,j))
	 bc(j1+2) = - shp(4,j)*(cmomc(3) + rhol*ulc(3,j))

	 be(j1)   = - shp(4,j)*(cmome(1) + rhol*ule(1,j))
	 be(j1+1) = - shp(4,j)*(cmome(2) + rhol*ule(2,j))
	 be(j1+2) = - shp(4,j)*(cmome(3) + rhol*ule(3,j))

       j1 = j1 + 3
      enddo


c store in right format


      j1=0
	j0=0
      do i= 1,nel
	 do j= 1,ndm
	  if (idg.eq.-1) then     ! Omega-
	    p(j1+j)  = p(j1+j)   +  bc(j0+j) 
	  elseif (idg.eq.1) then  ! Omega+
	    p(j1+j)  = p(j1+j)   + bc(j0+j) 
	    p(j1+j+3)= p(j1+j+3) + be(j0+j) 
	  endif
	 enddo
	 j0 = j0 + ndm
	 j1 = j1 + ndf
	enddo

	end








      subroutine stiff_mate3d_11(nst2,sym,shp,aa,idg, s)     
c--------------------------------------------------------------------71
c.... Compute the continuum material stiffness matrix and add to s
c     for partition of unity finite element
c  Declare variable types
c
c       nst2    enhanced degrees of freedom
c       aa      Elasticity matrix	times dv
c       idg     id of integration point 
c                +1 point is in omega+
c                -1 point is in omega-
c       shp     shape functions and its derivatives 
c               (spatial formulation)
c       sym     symmetry flag  
c                 sym=.true.   aa .eq. aa^T
c                 sym=.false.  aa .ne. aa^T
c       s       stiffness matrix
c
c--------------------------------------------------------------------71
	implicit none
      include  'sdata.h'
      include  'eldata.h'

	logical sym
	integer jj, i1, i, ii, j1, k, j, j0, i0, idg, nst2
	real*8 shp(4,nel), aa(6,6), bbd(3,6), s(nst,nst)
	real*8 ss(nst2,nst2)


      ss = 0.0d0
          i1 = 0
          do  i  = 1,nel
c
c.... Compute bmat-t * aa * dvol
            do jj = 1,6
              bbd(1,jj) = shp(1,i)*aa(1,jj)+shp(2,i)*aa(4,jj)
     1                  + shp(3,i)*aa(6,jj)
              bbd(2,jj) = shp(1,i)*aa(4,jj)+shp(2,i)*aa(2,jj)
     1                  + shp(3,i)*aa(5,jj)
              bbd(3,jj) = shp(1,i)*aa(6,jj)+shp(2,i)*aa(5,jj)
     1                  + shp(3,i)*aa(3,jj)
             enddo
c
            j1 = 0
c ... Sym check 
	      If(sym) then
	       k = i
	      else
	        k = nel
	      endif

            do  j  = 1,k
c
c.... Compute mechanics part of tangent stiffness
              do  jj = 1,3
                ss(i1+jj,j1+1) =  bbd(jj,1)*shp(1,j)
     1                 + bbd(jj,4)*shp(2,j) + bbd(jj,6)*shp(3,j)
                ss(i1+jj,j1+2) =  bbd(jj,4)*shp(1,j)
     1                 + bbd(jj,2)*shp(2,j) + bbd(jj,5)*shp(3,j)
                ss(i1+jj,j1+3) =  bbd(jj,6)*shp(1,j)
     1                 + bbd(jj,5)*shp(2,j) + bbd(jj,3)*shp(3,j)
              enddo 
			  
              j1 = j1 + 3
	      enddo
            i1 = i1 + 3
	    enddo
c  symmetrisieren
      If(sym) then
         do j = 1,nst2
           do i = 1,j
             ss(i,j) = ss(j,i)
           end do
         end do	  
	endif



      i1=0
      i0=0
      do i = 1,nel
	 j1=0
	 j0=0
       do j= 1,nel
	  do ii = 1,ndm
	    do jj= 1,ndm
	     if (idg.eq.-1) then    ! Omega-
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)     + ss(i0+ii,j0+jj)
	     elseif (idg.eq.1) then ! Omega+
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)     + ss(i0+ii,j0+jj)
	       s(i1+ii,j1+jj+3) = s(i1+ii,j1+jj+3)   + ss(i0+ii,j0+jj)
	       s(i1+ii+3,j1+jj) = s(i1+ii+3,j1+jj)   + ss(i0+ii,j0+jj)
	       s(i1+ii+3,j1+jj+3)= s(i1+ii+3,j1+jj+3)+ ss(i0+ii,j0+jj)
	     endif
	    enddo
	   enddo
	   j0 = j0 + ndm
	   j1 = j1 + ndf
	  enddo
	 i0 = i0 + ndm
	 i1 = i1 + ndf
	enddo



      end


       subroutine stiff_geom3d_11(nst2,shp,sigv,idg, s)     
c--------------------------------------------------------------------71
c.... Compute the continuum geometric stiffness matrix and add to s
c     for partition of unity finite element
c
c  Declare variable types
c
c       nst2    enhanced degrees of freedom
c       sigv    Cauchy stress	times dv
c       idg     id of integration point 
c                +1 point is in omega+
c                -1 point is in omega-
c       shp     shape functions and deerivatives 
c               (spatial formulation)
c       s       stiffness matrix
c
c--------------------------------------------------------------------71

	Implicit none

      include  'sdata.h'
      include  'cdata.h'
      include  'eltran.h'
      include  'eldata.h'

	integer jj,ii,j0,i0, i1, i, j1, j, idg, nst2
	real*8 sigv(6), shp(4,nel), tb(4)
	real*8 bdb, ada, s(nst,nst), dot
	real*8  bb(nst2,nst2)


	bb=0.0d0
     
	i1 = 0
      do i = 1,nel
         j1 = 0

         tb(1) = shp(1,i)*sigv(1) + shp(2,i)*sigv(4)
     1           +shp(3,i)*sigv(6)
         tb(2) = shp(1,i)*sigv(4) + shp(2,i)*sigv(2)
     1           +shp(3,i)*sigv(5)
         tb(3) = shp(1,i)*sigv(6) + shp(2,i)*sigv(5)
     1           +shp(3,i)*sigv(3)

		  do j = 1,i
             bdb     = dot(tb(1),shp(1,j),3)
             ada     = tb(4)*shp(4,j)
             do jj = 1,3
c....  Geometric contribution
               bb(i1+jj,j1+jj) =  bdb
             enddo
            j1 = j1 + 3
           enddo
            i1 = i1 + 3
      enddo
c  symmetrisieren
         do j = 1,nst2
           do i = 1,j
             bb(i,j) = bb(j,i)
           end do
         end do	  


      i1=0
      i0=0
      do i = 1,nel
	 j1=0
	 j0=0
       do j= 1,nel
	  do ii = 1,ndm
	    do jj= 1,ndm
	     if (idg.eq.-1) then    ! Omega-
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)   +  bb(i0+ii,j0+jj)
	     elseif(idg.eq.1) then  ! Omega+
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)   +  bb(i0+ii,j0+jj)
	       s(i1+ii,j1+jj+3) = s(i1+ii,j1+jj+3) +  bb(i0+ii,j0+jj)
	       s(i1+ii+3,j1+jj) = s(i1+ii+3,j1+jj) +  bb(i0+ii,j0+jj)
	       s(i1+ii+3,j1+jj+3)= s(i1+ii+3,j1+jj+3)+bb(i0+ii,j0+jj)
	     endif
	    enddo
	   enddo
	   j0 = j0 + ndm
	   j1 = j1 + ndf
	  enddo
	  i0 = i0 + ndm
	  i1 = i1 + ndf
	enddo

	end






       subroutine stiff_mass3d_11(idg,nst2,mflag,rho,shp, s)     
c--------------------------------------------------------------------71
c.... Compute the weighted mass matrix and add to s
c     for displacement finite element and store in special format
c
c  Declare variable types
c
c       idg     id of integration point 
c                +1 point is in omega+
c                -1 point is in omega-
c       nst2    enhanced degrees of freedom
c       shp     shape functions and its deerivatives 
c               (spatial formulation)
c       rho     density
c       mflag   Massinterpolation flag
c                mflag = 0 ...... lumped mass interpolation
c                mflag = 1 ...... variationally consitent mass interpolation
c       ctan(3) weighting factor dependent on the solution method
c       s       stiffness matrix
c
c--------------------------------------------------------------------71

	Implicit none

      include  'sdata.h'
      include  'cdata.h'
      include  'eltran.h'
      include  'eldata.h'

	integer jj,ii,j0,i0, i1, i, j1, j, nst2, idg
	real*8 shp(4,nel), mass, rho,am 
	real*8 ada, mflag, s(nst,nst), rhom
	real*8 aa(nst2,nst2)

      aa=0.0d0

	rho = rho*ctan(3)
     
	i1 = 0
      do i = 1,nel
         j1 = 0
         am = shp(4,i)*rho
c....  Lumped mass factor
	   rhom  = am*(1.d0 - mflag)
	   do jj = 1,3
c....   Add diagonal mass effects
            aa(i1+jj,i1+jj) =  rhom
         enddo

c....  Consistent mass factor
	    mass = am*mflag

		  do j = 1,i
             ada     = mass*shp(4,j)
             do jj = 1,3
c....  Consistent mass matrix
	         aa(i1+jj,j1+jj) =  aa(i1+jj,j1+jj) + ada
             enddo
            j1 = j1 + 3
           enddo
            i1 = i1 + 3
      enddo
c  symmetrisieren
      do j = 1,nst2
        do i = 1,j
          aa(i,j) = aa(j,i)
        end do
      end do	  

      i1=0
      i0=0
      do i = 1,nel
	 j1=0
	 j0=0
       do j= 1,nel
	  do ii = 1,ndm
	    do jj= 1,ndm
	     if (idg.eq.-1) then    ! Omega-
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)   +  aa(i0+ii,j0+jj)
	     elseif(idg.eq.1) then  ! Omega+
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)   +  aa(i0+ii,j0+jj)
c	       s(i1+ii,j1+jj+3) = s(i1+ii,j1+jj+3) +  aa(i0+ii,j0+jj)
c	       s(i1+ii+3,j1+jj) = s(i1+ii+3,j1+jj) +  aa(i0+ii,j0+jj)
	       s(i1+ii+3,j1+jj+3)= s(i1+ii+3,j1+jj+3)+aa(i0+ii,j0+jj)
	     endif
	    enddo
	  enddo
	  j0 = j0 + ndm
	  j1 = j1 + ndf
	 enddo
	 i0 = i0 + ndm
	 i1 = i1 + ndf
	enddo

	end







      subroutine inttria_11(l,Lc, lint, Lg, sg2d)
c--------------------------------------------------------------------71
c.... Compute volume coordinates of Gausspoints at 
c     triangular discontinuities
c
c  Declare variable types
c
c       lint    number of Gausspoints at disc.
c       l       order of integration
c       Lc      Volume coordinates at disc. corners 
c       Lg      Volume coordinates at Gausspoints and weights
c       shp2d   shapefunctions (without derivatives)
c       
c
c--------------------------------------------------------------------71
      implicit none

      integer lint,l,i,j,k

	real*8 Lc(4,3), Lg(5,36), sg2d(4,36)
	real*8 shp2d(3,36)
      

	lint=0
	call tint2d_11(l,lint,sg2d)
c compute shape functions      
      do i=1,lint
	  shp2d(1,i) = sg2d(1,i)
	  shp2d(2,i) = sg2d(2,i)
	  shp2d(3,i) = sg2d(3,i)
	enddo

c compute volume coordinates at gausspoints
      Lg=0.0d0
      do i= 1,lint
	  do j= 1,4
	    do k= 1,3
	      Lg(j,i) = Lg(j,i) + Lc(j,k)*shp2d(k,i)
	    enddo
	  enddo
	enddo
c store weigths in Lg
      do i= 1,lint
	   Lg(5,i) = sg2d(4,i)
	enddo
       
      end





      subroutine intquad_11(l,Lc, lint, Lg, sg2d)
c--------------------------------------------------------------------71
c.... Compute volume coordinates of Gausspoints for
c     4 noded discontinuities
c
c  Declare variable types
c
c       lint    number of Gausspoints at disc.
c       l       order of integration
c       Lc      Volume coordinates at disc. corners 
c       Lg      Volume coordinates at Gausspoints and weights
c       shp2d   shapefunctions (without derivatives)
c       
c
c--------------------------------------------------------------------71
      implicit none

      integer lint,l,i,j,k

	real*8 Lc(4,4), Lg(5,36), sg2d(3,36)
	real*8 shp2d(4,36)
      

	lint=0
	call int2d_11(l,lint,sg2d)
c compute shape functions      
      do i=1,lint        ! kann man besser machen!!
	  shp2d(1,i) = 0.25d0*(1.0d0 - sg2d(1,i))*(1.0d0 - sg2d(2,i))
	  shp2d(2,i) = 0.25d0*(1.0d0 + sg2d(1,i))*(1.0d0 - sg2d(2,i))
	  shp2d(3,i) = 0.25d0*(1.0d0 + sg2d(1,i))*(1.0d0 + sg2d(2,i))
	  shp2d(4,i) = 0.25d0*(1.0d0 - sg2d(1,i))*(1.0d0 + sg2d(2,i))
	enddo



c compute volume coordinates at gausspoints
      Lg=0.0d0
      do i= 1,lint
	  do j= 1,4
	    do k= 1,4
	      Lg(j,i) = Lg(j,i) + Lc(j,k)*shp2d(k,i)
	    enddo
	  enddo
	enddo
c store weigths in Lg
      do i= 1,lint
	   Lg(5,i) = sg2d(3,i)
	enddo
       
      end


      subroutine kine_d11 (nel, fb, invfb, Jb,ue, due, shp, N0, 
     1                   u0,du0, Kb,nb)
c--------------------------------------------------------------------71
c compute kinematical quantities at the fictitious discontinuity 
c
c  Declare variable types
c
c       nel    node per element
c       fb     average deformation gradient
c       invfb  inverse average deformation gradient
c       Jb     volume ratio at the discontinuity
c       ue     enhanced degrees of freedom
c       due    increment of enhanced degrees of freedom
c       shp    shape function
c       N0     normal vector in the ref. conf.
c       Kb     area ratio at the discontinuity
c       nb     normal vector in the curr. conf.
c       u0     enhanced displacement at gauss point
c       du0    increment of enhanced displacement at gauss point
c
c--------------------------------------------------------------------71
    
      implicit none

	integer i,k, nel

	real*8 N0(3), Jb, Kb, nb(3), ue(3,nel), due(3,nel), shp(4,nel)
	real*8 Fb(3,3), invfb(3,3), norm, u0(3), du0(3)

c compute u0 and du0
      u0  =0.0d0
      du0 =0.0d0
	do i= 1,3
        do k= 1,nel
	    u0(i)  = u0(i)  + shp(4,k)*ue(i,k)
	    du0(i) = du0(i) + shp(4,k)*due(i,k)
	  enddo
      enddo

c compute nb
      nb=0.0d0
	do i= 1,3
	  do k= 1,3
	    nb(i) = nb(i) + N0(k)*invFb(k,i)
c	    nb(i) = nb(i) + N0(k)*invFb(i,k)
	  enddo
	enddo
	norm = sqrt(nb(1)*nb(1)+nb(2)*nb(2)+nb(3)*nb(3))
	nb = nb/norm
c compute area ratio Kb
      Kb = Jb*norm



	end


      subroutine resi_trac3d_11(nst2,t,shp,p)     
c--------------------------------------------------------------------71
c.... Compute the discrete residuum for
c     partition of unity finite element
c
c
c  Declare variable types
c
c       nst2    enhanced degrees of freedom
c       t       traction times da
c       shp     shape functions and its deerivatives (not needed)
c               (spatial formulation)
c       p       residuum consists of compatible and enhanced loads
c
c--------------------------------------------------------------------71

	Implicit none

      include  'cdata.h'
      include  'eldata.h'
      include  'sdata.h'

	integer i,j, j1, j0, jj, nst2
	real*8 shp(4,nel), t(3)
	real*8 p(nst)

      real*8 a(nst2)

c compatible part


c.... Compute the residum N^T t
c
      j1 = 0
      do j = 1,nel
	  do jj= 1,ndm
          a(j1+jj)   =  -shp(4,j)*t(jj)
        enddo
       j1 = j1 + ndm
      enddo

c store in right format
      j1=0
	j0=0
      do i= 1,nel
	 do j= 1,ndm
	  p(j1+j+3)= p(j1+j+3) + a(j0+j)
	 enddo
	 j0 = j0 + ndm
	 j1 = j1 + ndf
	enddo

c      call plot_residuum_small11(p)


	end






      subroutine stiff_trac3d_11(nst2,shp,cu,cn,t,nb, s)
c---------------------------------------------------TCG 09.04.2003---71
c.... Computes the cohesive contribution to the element stiffness.
c     K_ue_ue = rrrr + 0.5*ssss 
c     K_ue_uc = ssss 
c     
c  Declare variable types
c
c       nst2    enhanced degrees of freedom
c       shp     shape functions and its derivatives 
c       cu      cohesive stiffness with respect to ue times da
c       cn      cohesive stiffness with respect to n times da
c       t       cohesive traction times da
c       nb      normal onto fictitious discontinuity
c       s       stiffness matrix
c       shpshp  N^i N^j
c       cntk    stiffness
c       kkb     weighted interpolation of linearization of area ratio
c       nnb     interpolation of linearization of area ratio normal
c
c--------------------------------------------------------------------71
	implicit none
      include  'sdata.h'
      include  'eldata.h'

	integer i, j, ii, jj, k, l, nst2
	integer i0, j0, i1, j1
	real*8  shp(4,nel), s(nst,nst)
	real*8  shpshp(nel,nel)
	real*8  cu(3,3), cn(3,3), t(3), kkb(3,nel), nb(3), alpha
	real*8  nnb(3,3,nel), cntk(3,3,nel)
	real*8  rrrr(nst2,nst2), ssss(nst2,nst2)

c define constants 
      do i= 1,nel
	  do j= i, nel
	    shpshp(i,j) = shp(4,i)*shp(4,j)
          shpshp(j,i) = shpshp(i,j)  ! sym
	  enddo
	enddo
c compute interpolations
c weighted interpolation of area ratio
c for ndm=3!
      do k= 1,nel
	  alpha = nb(1)*shp(1,k) + nb(2)*shp(2,k) +nb(3)*shp(3,k) 
	  do i= 1,3
	    kkb(i,k) = shp(i,k) - alpha*nb(i)
	  enddo
	enddo
c  interpolation of current normal vector
      do k= 1,nel
        do i= 1,ndm
	    do j= 1,ndm
	      nnb(i,j,k) = -kkb(i,k)*nb(j)
	    enddo
	  enddo
	enddo

c sum up stiffness

      cntk =0.0d0
      do i= 1,nel
	  do ii= 1,ndm
	    do jj= 1,ndm
	      do l= 1,ndm
	        cntk(ii,jj,i) = cntk(ii,jj,i) + cn(ii,l)*nnb(l,jj,i)
		  enddo  

C contribution due to change in area (for Cauchy traction only)
c            cntk(ii,jj,i) = cntk(ii,jj,i) + t(ii)*kkb(jj,i)
	    enddo
	  enddo
	enddo


c  ssss
      ii = 0
	do i= 1,nel
	  jj = 0
	  do j= 1,nel
         do k= 1,ndm
	     do l= 1,ndm
             ssss(ii+k,jj+l) = cntk(k,l,j)*shp(4,i)
	     enddo
	   enddo
         jj = jj + ndm

        enddo
	  ii = ii + ndm
	enddo

c  rrrr (symmetric)
      ii = 0
	do i= 1,nel
	  jj = 0
c	  jj = ii
	  do j= 1,nel
c	  do j= i,nel
         do k= 1,ndm
	     do l= 1,ndm
c	     do l= k,ndm
			rrrr(ii+k,jj+l) = cu(k,l)*shpshp(i,j) 
c minor symmetry
c			rrrr(ii+l,jj+k) = rrrr(ii+k,jj+l)  
	     enddo
	   enddo
         jj = jj + ndm

        enddo
	  ii = ii + ndm
	enddo
c major symmetry
c      do i= 1,nst2
c        do j= i+1,nst2
c	    rrrr(j,i) = rrrr(i,j)
c	  enddo
c	enddo




c compute cohesive stiffness 
      
      i0 = 0
      i1 = 0
      do i= 1,nel
	  j0 = 0
	  j1 = 0
	  do j= 1,nel
	    do ii= 1,ndm
	      do jj=1,ndm
	        s(i1+ii+3,j1+jj) = s(i1+ii+3,j1+jj) + ssss(i0+ii,j0+jj)
	        s(i1+ii+3,j1+jj+3) = s(i1+ii+3,j1+jj+3) + 
	1			rrrr(i0+ii,j0+jj) + 0.5d0*ssss(i0+ii,j0+jj)
	      enddo
	    enddo
	    j0 = j0 + ndm
	    j1 = j1 + ndf
        enddo
	  i0 = i0 + ndm
	  i1 = i1 + ndf
	enddo

      end





      subroutine prisshp_11(xi, xl,  xsj, shp )


c----------------------------------------------------TCG.17.01.2003--]
c               Compute 3-d prismatic element shape
c               functions and their derivatives w/r x,y,z
c
c      Inputs:
c         xi(4)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh
c
c      Outputs:
c         xsj       - Jacobian determinant at point
c         shp(4,*)  - Shape functions and derivatives at point
c                     shp(1,i) = dN_i/dx
c                     shp(2,i) = dN_i/dy
c                     shp(3,i) = dN_i/dz
c                     shp(4,i) =  N_i
c
c   declared variables
c        xs       Transpose Jacobian matrix
c        xs       Inverse Jacobian matrix
c-----[--.----+----.----+----.-----------------------------------------]
      implicit none

	integer i, j ,k
	
	real*8 xi(3), xl(3,6), xsj, shp(4,6)
	real*8 r,s,t, one8, one4
	real*8 r8, s8, t8, rs8, rt8, st8, rst8
	real*8 r4, s4, t4, rs4, rt4, st4, rst4
	real*8 xs(3,3), rxsj, ad(3,3), c1, c2, c3

      data one8/0.125d0/ 
      data one4/0.25d0/ 

	r = xi(1)
	s = xi(2)
	t = xi(3)

	r8 = one8*r
	s8 = one8*s
	t8 = one8*t

	rs8 = r8*s
	rt8 = r8*t
	st8 = s8*t

	rst8 = rs8*t
	 
      r4 = r8 + r8
      s4 = s8 + s8
      t4 = t8 + t8

	rs4 = rs8 + rs8
	rt4 = rt8 + rt8
	st4 = st8 + st8

	rst4 = rst8 + rst8

c compute shapefunctions
      shp(4,1) =  one8 - r8 - s8 + rs8 - t8 + rt8 + st8 - rst8
      shp(4,2) =  one8 + r8 - s8 - rs8 - t8 - rt8 + st8 + rst8
      shp(4,3) =  one8 + r8 + s8 + rs8 - t8 - rt8 - st8 - rst8
      shp(4,4) =  one8 - r8 + s8 - rs8 - t8 + rt8 - st8 + rst8
      shp(4,5) =  one4 - r4 + t4 - rt4
      shp(4,6) =  one4 + r4 + t4 + rt4

c compute derivatives with respect to r
      shp(1,1) =  -one8 + s8 + t8 - st8
      shp(1,2) =  -shp(1,1)
	shp(1,3) =   one8 + s8 - t8 - st8
	shp(1,4) =  -shp(1,3)
	shp(1,5) =  -one4 - t4
	shp(1,6) =  -shp(1,5) 

c compute derivatives with respect to s
      shp(2,1) =  -one8 + r8 + t8 - rt8
	shp(2,2) =  -one8 - r8 + t8 + rt8 
      shp(2,3) =  -shp(2,2)
	shp(2,4) =  -shp(2,1)
	shp(2,5) =   0.0d0
	shp(2,6) =   0.0d0

c compute derivatives with respect to t
      shp(3,1) =  -one8 + r8 + s8 - rs8
	shp(3,2) =  -one8 - r8 + s8 + rs8
	shp(3,3) =  -one8 - r8 - s8 - rs8
	shp(3,4) =  -one8 + r8 - s8 + rs8
	shp(3,5) =   one4 - r4 
	shp(3,6) =   one4 + r4

c compute determinant of jacobian matrix 

c     Compute jacobian transformation
      
      do j = 1,3
        xs(j,1) = (xl(j,2) - xl(j,1))*shp(1,2)
     &          + (xl(j,4) - xl(j,3))*shp(1,4)
     &          + (xl(j,6) - xl(j,5))*shp(1,6)

        xs(j,2) = (xl(j,4) - xl(j,1))*shp(2,4)
     &          + (xl(j,3) - xl(j,2))*shp(2,3)

        xs(j,3) = xl(j,1)*shp(3,1) 
     &          + xl(j,2)*shp(3,2)
     &          + xl(j,3)*shp(3,3)
     &          + xl(j,4)*shp(3,4)
     &          + xl(j,5)*shp(3,5)
     &          + xl(j,6)*shp(3,6)

      enddo


c     Compute adjoint to jacobian

      ad(1,1) = xs(2,2)*xs(3,3) - xs(2,3)*xs(3,2)
      ad(1,2) = xs(3,2)*xs(1,3) - xs(3,3)*xs(1,2)
      ad(1,3) = xs(1,2)*xs(2,3) - xs(1,3)*xs(2,2)

      ad(2,1) = xs(2,3)*xs(3,1) - xs(2,1)*xs(3,3)
      ad(2,2) = xs(3,3)*xs(1,1) - xs(3,1)*xs(1,3)
      ad(2,3) = xs(1,3)*xs(2,1) - xs(1,1)*xs(2,3)

      ad(3,1) = xs(2,1)*xs(3,2) - xs(2,2)*xs(3,1)
      ad(3,2) = xs(3,1)*xs(1,2) - xs(3,2)*xs(1,1)
      ad(3,3) = xs(1,1)*xs(2,2) - xs(1,2)*xs(2,1)

c     Compute determinant of jacobian

      xsj  = xs(1,1)*ad(1,1) + xs(1,2)*ad(2,1) + xs(1,3)*ad(3,1)
      rxsj = 1.d0/xsj

c     Compute jacobian inverse

      do j = 1,3
        do i = 1,3
          xs(i,j) = ad(i,j)*rxsj
        end do
      end do

c     Compute derivatives with repect to global coords.

      do k = 1,6

        c1 = shp(1,k)*xs(1,1) + shp(2,k)*xs(2,1) + shp(3,k)*xs(3,1)
        c2 = shp(1,k)*xs(1,2) + shp(2,k)*xs(2,2) + shp(3,k)*xs(3,2)
        c3 = shp(1,k)*xs(1,3) + shp(2,k)*xs(2,3) + shp(3,k)*xs(3,3)

        shp(1,k) = c1
        shp(2,k) = c2
        shp(3,k) = c3

      end do
      end




      subroutine resi_stre3dQ1_11(idd,nst2,shp,sigv, p)     
c--------------------------------------------------------------------71
c.... Compute the continuum residuum due to stress field  for
c     displacement finite element and store in special format
c
c
c  Declare variable types
c
c       idd     node id
c       nst2    enhanced degrees of freedom
c       sigv    Cauchy stress	times dv
c       shp     shape functions and its deerivatives 
c               (spatial formulation)
c       p       residuum consists of compatible and enhanced loads
c
c--------------------------------------------------------------------71

	Implicit none

      include  'cdata.h'
      include  'eldata.h'
      include  'sdata.h'

	integer i,j, j1, j0, nst2, idd(nel)
	real*8 shp(4,nel), sigv(6)
	real*8 p(nst)

      real*8 a(nst2)
c compatible part

c
c.... Compute the residual due to stress field
c
      j1 = 1
      do j = 1,nel
       a(j1)   = - shp(1,j)*sigv(1)-shp(2,j)*sigv(4)
     1              - shp(3,j)*sigv(6)
       a(j1+1) = - shp(1,j)*sigv(4)-shp(2,j)*sigv(2)
     1              - shp(3,j)*sigv(5)
       a(j1+2) = - shp(1,j)*sigv(6)-shp(2,j)*sigv(5)
     1              - shp(3,j)*sigv(3)

       j1 = j1 + 3
      enddo
c store in right format

      j1=0
	j0=0
      do i= 1,nel
	 do j= 1,ndm
	  p(j1+j)  = p(j1+j)  + a(j0+j)
	  if (idd(i).eq.1) then
	    p(j1+j+3)  = p(j1+j+3)  + a(j0+j)
	  endif 
	 enddo
	 j0 = j0 + ndm
	 j1 = j1 + ndf
	enddo

	end



      subroutine resi_mass3dQ1_11(idd,nst2,ul,mflag,rho,shp, p)     
c--------------------------------------------------------------------71
c.... Compute the residuum due to inertia parts for
c     displacement finite element and store in special format
c
c
c  Declare variable types
c
c       idd     node id
c       nst2    enhanced degrees of freedom
c       shp     shape functions and its deerivatives 
c               (spatial formulation)
c       rho     current density times current volume =
c               ref. density times ref. volume = mass
c       mflag   Massinterpolation flag
c                mflag = 0 ...... lumped mass interpolation
c                mflag = 1 ...... variationally consitent mass interpolation
c       ul      solution array (accelerations are needed)
c       p       residuum consists of compatible and enhanced loads
c
c--------------------------------------------------------------------71

	Implicit none

      include  'cdata.h'
      include  'eldata.h'
      include  'sdata.h'

	integer i,j, j1, j0, nst2,idd(nel)
	real*8 rho, shp(4,nel)
	real*8 ul(3,nel),cmom(3), p(nst)
	real*8 mflag, rhol, rhom

      real*8 b(nst2)


     


	rhol = rho
	rhom = rhol*mflag
	rhol = rhol- rhom
	do i = 1,3
	 cmom(i) = 0.0d0
	 do j = 1,nel
	   cmom(i) = cmom(i) + shp(4,j)*ul(i,j)
       enddo
	 cmom(i) = rhom*cmom(i)
      enddo
c
c.... Compute the residual due to inertia parts
c
      j1 = 1
      do j = 1,nel

	 b(j1)   = - shp(4,j)*(cmom(1) + rhol*ul(1,j))
	 b(j1+1) = - shp(4,j)*(cmom(2) + rhol*ul(2,j))
	 b(j1+2) = - shp(4,j)*(cmom(3) + rhol*ul(3,j))

       j1 = j1 + 3
      enddo

c store in right format


      j1=0
	j0=0
      do i= 1,nel
	 do j= 1,ndm
	  p(j1+j)  = p(j1+j)  + b(j0+j) 
	  if (idd(i).eq.1) then
	    p(j1+j+3)  = p(j1+j+3)  + b(j0+j)
	  endif 
	 enddo
	 j0 = j0 + ndm
	 j1 = j1 + ndf
	enddo

	end






      subroutine stiff_mate3dQ1_11(idd,nst2,sym,shp,aa, s)     
c--------------------------------------------------------------------71
c.... Compute the continuum material stiffness matrix and add to s
c     for displacement finite element and store in special format
c  Declare variable types
c
c       idd     node id
c       nst2    enhanced degrees of freedom
c       aa      Elasticity matrix	times dv
c       shp     shape functions and its derivatives 
c               (spatial formulation)
c       sym     symmetry flag  
c                 sym=.true.   aa .eq. aa^T
c                 sym=.false.  aa .ne. aa^T
c       s       stiffness matrix
c
c--------------------------------------------------------------------71
	implicit none
      include  'sdata.h'
      include  'eldata.h'

	logical sym
	integer jj, i1, i, ii, j1, k, j, j0, i0, nst2, idd(nel)
	real*8 shp(4,nel), aa(6,6), bbd(3,6), s(nst,nst)
	real*8 ss(nst2,nst2)


      ss = 0.0d0
          i1 = 0
          do  i  = 1,nel
c
c.... Compute bmat-t * aa * dvol
            do jj = 1,6
              bbd(1,jj) = shp(1,i)*aa(1,jj)+shp(2,i)*aa(4,jj)
     1                  + shp(3,i)*aa(6,jj)
              bbd(2,jj) = shp(1,i)*aa(4,jj)+shp(2,i)*aa(2,jj)
     1                  + shp(3,i)*aa(5,jj)
              bbd(3,jj) = shp(1,i)*aa(6,jj)+shp(2,i)*aa(5,jj)
     1                  + shp(3,i)*aa(3,jj)
             enddo
c
            j1 = 0
c ... Sym check 
	      If(sym) then
	       k = i
	      else
	        k = nel
	      endif

            do  j  = 1,k
c
c.... Compute mechanics part of tangent stiffness
              do  jj = 1,3
                ss(i1+jj,j1+1) =  bbd(jj,1)*shp(1,j)
     1                 + bbd(jj,4)*shp(2,j) + bbd(jj,6)*shp(3,j)
                ss(i1+jj,j1+2) =  bbd(jj,4)*shp(1,j)
     1                 + bbd(jj,2)*shp(2,j) + bbd(jj,5)*shp(3,j)
                ss(i1+jj,j1+3) =  bbd(jj,6)*shp(1,j)
     1                 + bbd(jj,5)*shp(2,j) + bbd(jj,3)*shp(3,j)
              enddo 
			  
              j1 = j1 + 3
	      enddo
            i1 = i1 + 3
	    enddo
c  symmetrisieren
      If(sym) then
         do j = 1,nst2
           do i = 1,j
             ss(i,j) = ss(j,i)
           end do
         end do	  
	endif


c store in right format

      i1=0
      i0=0
      do i = 1,nel
	  j1=0
	  j0=0
        do j= 1,nel
	    do ii = 1,ndm
	      do jj= 1,ndm
	       s(i1+ii,j1+jj) = s(i1+ii,j1+jj)   + ss(i0+ii,j0+jj)
	       if (idd(i).eq.1) then
	         s(i1+ii+3,j1+jj)=s(i1+ii+3,j1+jj)+ss(i0+ii,j0+jj)
	       endif
	       if (idd(j).eq.1) then
	         s(i1+ii,j1+jj+3)=s(i1+ii,j1+jj+3)+ss(i0+ii,j0+jj)
	       endif
	       if ((idd(i).eq.1).and.(idd(j).eq.1)) then
	         s(i1+ii+3,j1+jj+3)=s(i1+ii+3,j1+jj+3)+ss(i0+ii,j0+jj)
	       endif 
	      enddo
	     enddo
	     j0 = j0 + ndm
	     j1 = j1 + ndf
	   enddo
	   i0 = i0 + ndm
	   i1 = i1 + ndf
	enddo



      end


       subroutine stiff_geom3dQ1_11(idd,nst2,shp,sigv, s)     
c--------------------------------------------------------------------71
c.... Compute the continuum geometric stiffness matrix and add to s
c     for displacement finite element and store in special format
c
c  Declare variable types
c
c       idd     node id
c       nst2    enhanced degrees of freedom
c       sigv    Cauchy stress	times dv
c       shp     shape functions and derivatives 
c               (spatial formulation)
c       s       stiffness matrix
c
c--------------------------------------------------------------------71

	Implicit none

      include  'sdata.h'
      include  'cdata.h'
      include  'eltran.h'
      include  'eldata.h'

	integer jj,ii,j0,i0, i1, i, j1, j, nst2, idd(nel)
	real*8 sigv(6),  shp(4,nel), tb(4)
	real*8 bdb, ada, s(nst,nst), dot
	real*8 bb(nst2,nst2)


	bb=0.0d0
     
	i1 = 0
      do i = 1,nel
         j1 = 0

         tb(1) = shp(1,i)*sigv(1) + shp(2,i)*sigv(4)
     1           +shp(3,i)*sigv(6)
         tb(2) = shp(1,i)*sigv(4) + shp(2,i)*sigv(2)
     1           +shp(3,i)*sigv(5)
         tb(3) = shp(1,i)*sigv(6) + shp(2,i)*sigv(5)
     1           +shp(3,i)*sigv(3)

		  do j = 1,i
             bdb     = dot(tb(1),shp(1,j),3)
             ada     = tb(4)*shp(4,j)
             do jj = 1,3
c....  Geometric contribution
               bb(i1+jj,j1+jj) =  bdb
             enddo
            j1 = j1 + 3
           enddo
            i1 = i1 + 3
      enddo
c  symmetrisieren
         do j = 1,nst2
           do i = 1,j
             bb(i,j) = bb(j,i)
           end do
         end do	  

c store in right format

      i1=0
      i0=0
      do i = 1,nel
	  j1=0
	  j0=0
        do j= 1,nel
	    do ii = 1,ndm
	      do jj= 1,ndm
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)   + bb(i0+ii,j0+jj)
	       if (idd(i).eq.1) then
	         s(i1+ii+3,j1+jj)=s(i1+ii+3,j1+jj)+bb(i0+ii,j0+jj)
	       endif
	       if (idd(j).eq.1) then
	         s(i1+ii,j1+jj+3)=s(i1+ii,j1+jj+3)+bb(i0+ii,j0+jj)
	       endif
	       if ((idd(i).eq.1).and.(idd(j).eq.1)) then
	         s(i1+ii+3,j1+jj+3)=s(i1+ii+3,j1+jj+3)+bb(i0+ii,j0+jj)
	       endif 
	      enddo
	     enddo
	     j0 = j0 + ndm
	     j1 = j1 + ndf
	   enddo
	   i0 = i0 + ndm
	   i1 = i1 + ndf
	enddo

	end






       subroutine stiff_mass3dQ1_11(idd,nst2,mflag,rho,shp, s)     
c--------------------------------------------------------------------71
c.... Compute the weighted mass matrix and add to s
c     for displacement finite element and store in special format
c
c  Declare variable types
c
c       idd     node id
c       nst2    enhanced degrees of freedom
c       shp     shape functions and its deerivatives 
c               (spatial formulation)
c       rho     density
c       mflag   Massinterpolation flag
c                mflag = 0 ...... lumped mass interpolation
c                mflag = 1 ...... variationally consitent mass interpolation
c       ctan(3) weighting factor dependent on the solution method
c       s       stiffness matrix
c
c--------------------------------------------------------------------71

	Implicit none

      include  'sdata.h'
      include  'cdata.h'
      include  'eltran.h'
      include  'eldata.h'

	integer jj,ii,j0,i0, i1, i, j1, j, nst2, idd(nel)
	real*8 shp(4,nel), mass, rho,am 
	real*8 ada, mflag, s(nst,nst), rhom
	real*8 aa(nst2,nst2)



      aa=0.0d0

	rho = rho*ctan(3)
     
	i1 = 0
      do i = 1,nel
         j1 = 0
         am = shp(4,i)*rho
c....  Lumped mass factor
	   rhom  = am*(1.d0 - mflag)
	   do jj = 1,3
c....   Add diagonal mass effects
            aa(i1+jj,i1+jj) =  rhom
         enddo

c....  Consistent mass factor
	   mass = am*mflag

	   do j = 1,i
           ada     = mass*shp(4,j)
           do jj = 1,3
c....  Consistent mass matrix
	       aa(i1+jj,j1+jj) =  aa(i1+jj,j1+jj) + ada
           enddo
           j1 = j1 + 3
         enddo
         i1 = i1 + 3
      enddo
c  symmetrisieren
         do j = 1,nst2
           do i = 1,j
             aa(i,j) = aa(j,i)
           end do
         end do	  

c store in right format

      i1=0
      i0=0
      do i = 1,nel
	  j1=0
	  j0=0
        do j= 1,nel
	    do ii = 1,ndm
	      do jj= 1,ndm
	       s(i1+ii,j1+jj)   = s(i1+ii,j1+jj)   + aa(i0+ii,j0+jj)
	       if (idd(i).eq.1) then
	         s(i1+ii+3,j1+jj)=s(i1+ii+3,j1+jj)+aa(i0+ii,j0+jj)
	       endif
	       if (idd(j).eq.1) then
	         s(i1+ii,j1+jj+3)=s(i1+ii,j1+jj+3)+aa(i0+ii,j0+jj)
	       endif
	       if ((idd(i).eq.1).and.(idd(j).eq.1)) then
	         s(i1+ii+3,j1+jj+3)=s(i1+ii+3,j1+jj+3)+aa(i0+ii,j0+jj)
	       endif 
	      enddo
	     enddo
	     j0 = j0 + ndm
	     j1 = j1 + ndf
	   enddo
	   i0 = i0 + ndm
	   i1 = i1 + ndf
	enddo

	end





      subroutine asseble_tetra_11(id, xn,xl,xlsub,con,sgn)
c---------------------------------------------------TCG 19.02.2003--71
c.... assembles tetrahedral subelement (the discontinuity is 
c     a triangle) 
c
c  Declare variable types
c
c       id      node id (+1,-1)               
c       xn      nodes on cutted face (triangle)               
c       xl      nodes of tet
c       xlsub   nodes of the prismatic element in right format
c       con     connectivity of subtet
c            _________________________________________
c           |    |         |        |        |        |
c           |  4 |  node 1 | node 2 | node 3 | node 4 |
c           |____|_________|________|________|________|
c
c   
c  element nodes:         1   -
c						2    |  local tetrahedral nodes
c						3	 |
c						4	_|
c						5    |
c						6    |  nodes on discontinuity
c						7   _|
c
c       flag     node flag
c                 .false.   tetraeder node is allready node of subtet
c
c--------------------------------------------------------------------71
      implicit none

	real*8 xn(3,3), xl(3,4), xlsub(3,4)
	real*8 x(3), y(3), z(3), temp, crit

	integer id(4), i,j, sgn, con(7), itemp



      con(1) = 4
	do i= 1,3
	  con(i+1) = i+4
	  do j= 1,3
	    xlsub(i,j) = xn(i,j)
	  enddo
	enddo
c add fourth node to subtet
      if(abs(id(2)+id(3)+id(4)).eq.3) then 
	  con(5) = 1
	  do i= 1,3
	    xlsub(i,4) = xl(i,1)
	    sgn = id(1)
	  enddo
      elseif(abs(id(1)+id(3)+id(4)).eq.3) then	    
	  con(5) = 2
        do i= 1,3
	    xlsub(i,4) = xl(i,2)
	    sgn = id(2)
	  enddo
      elseif(abs(id(1)+id(2)+id(4)).eq.3) then	    
	  con(5) = 3
	  do i= 1,3
	    xlsub(i,4) = xl(i,3)
	    sgn = id(3)
	  enddo
      elseif(abs(id(1)+id(2)+id(3)).eq.3) then	    
	  con(5) = 4
	  do i= 1,3
	    xlsub(i,4) = xl(i,4)
	    sgn = id(4)
	  enddo
	endif
 
c check right node numbering
      
	do i= 1,3
	  x(i) = xlsub(i,2)-xlsub(i,1)
	  y(i) = xlsub(i,3)-xlsub(i,1)
	  z(i) = xlsub(i,4)-xlsub(i,1)
	enddo 
c compute spatproduct 
      crit = (x(2)*y(3)-y(2)*x(3))*z(1) +
	1       (y(1)*x(3)-x(1)*y(3))*z(2) + 
     2	   (x(1)*y(2)-y(1)*x(2))*z(3)		   
 
 
 
      if (crit.lt.0.0d0) then  ! interchange node 2 and node 3
        itemp  = con(3)
	  con(3) = con(4)
	  con(4) = itemp
        do i= 1,3
	    temp       = xlsub(i,2) 
	    xlsub(i,2) = xlsub(i,3)
	    xlsub(i,3) = temp
	  enddo	  
      endif
      end




      subroutine asseble_prisma_11(oflag,id, xn,xl,xlsub,con)
c---------------------------------------------------TCG 19.02.2003--71
c.... assembles prismatic subelement, where the discontinuity is 
c     a 4-noded patch 
c
c  Declare variable types
c
c       id      node id (+1,-1)               
c       xn      nodes on cutted face (triangle)               
c       xl      nodes of tet
c       xlsub   nodes of the prismatic element in right format
c       con     connectivity of prisma
c            ___________________________________________________________
c           |    |         |        |        |        |        |        |
c           |  6 |  node 1 | node 2 | node 3 | node 4 | node 5 | node 6 |
c           |____|_________|________|________|________|________|________|
c
c   
c  element nodes:         1   -
c						2    |  local tetrahedral nodes
c						3	 |
c						4	_|
c						5    |
c						6    |  nodes on discontinuity
c						7    |
c						8   _|
c
c       flag     node flag
c                 .false.   tetraeder node is allready node of prisma
c
c--------------------------------------------------------------------71
      implicit none

	real*8 xn(3,4), xl(3,4), xlsub(3,6)
	real*8 y(3), x(3), eps, crit
	real*8 z(3), temp, t1, t2, t3
	real*8 a(3,2)

	integer id(4), i,j, k, kk,oflag, con(7)

	logical flag(4)

      data eps/1.0d-10/

      
	con(1) = 6
c define node 5 and 6 and trial chooce a1 and a2
      k=5
	kk=1
	do j= 1,4
        if (id(j).eq.oflag) then
	    con(k+1) = j
	    do i= 1,3
	      xlsub(i,k) = xl(i,j)
	    enddo
	    k=k+1
	  else
	    do i= 1,3
	      a(i,kk) = xl(i,j)
	    enddo
	    kk = kk+1
	  endif
      enddo

c define a1 and a2 via right coordiante system
      do i= 1,3
	  x(i) = xlsub(i,6) - a(i,1)
	  y(i) = xlsub(i,5) - a(i,1)
	  z(i) = a(i,2) - a(i,1)
	enddo
c  compute criteria
      crit = (x(2)*y(3)-y(2)*x(3))*z(1) +
	1       (y(1)*x(3)-x(1)*y(3))*z(2) + 
     2	   (x(1)*y(2)-y(1)*x(2))*z(3)		   

      if (crit.lt.0.0d0) then  ! interchange a1 and a2
        do i= 1,3
	    temp   = a(i,1) 
	    a(i,1) = a(i,2)
	    a(i,2) = temp
	  enddo	  
      endif

c define nodes on interface plane
      flag=.true.
c node 4
      do j= 1,4
        if(flag(j)) then
          do i= 1,3
	      x(i) = xlsub(i,5) - xn(i,j)
	      y(i) =     a(i,1) - xn(i,j)
	    enddo
c  compute criteria
          t1 = x(2)*y(3)-y(2)*x(3)
	    t2 = y(1)*x(3)-x(1)*y(3)
          t3 = x(1)*y(2)-y(1)*x(2)		   
	    crit = t1*t1 + t2*t2 + t3*t3
	    if(abs(crit).lt.eps) then ! node j = node 4
	      con(5) = j + 4
            do i= 1,3
	        xlsub(i,4) = xn(i,j)
	      enddo	
	      flag(j) =.false.
		endif  
	  endif
      enddo


c node 1
      do j= 1,4
        if(flag(j)) then
          do i= 1,3
	      x(i) = xlsub(i,5) - xn(i,j)
	      y(i) =     a(i,2) - xn(i,j)
	    enddo
c  compute criteria
          t1 = x(2)*y(3)-y(2)*x(3)
	    t2 = y(1)*x(3)-x(1)*y(3)
          t3 = x(1)*y(2)-y(1)*x(2)		   
	    crit = t1*t1 + t2*t2 + t3*t3
	    if(abs(crit).lt.eps) then ! node j = node 1
	      con(2) = j + 4
            do i= 1,3
	        xlsub(i,1) = xn(i,j)
	      enddo	
	      flag(j) =.false.
		endif  
	  endif
      enddo



c node 3
      do j= 1,4
        if(flag(j)) then
          do i= 1,3
	      x(i) = xlsub(i,6) - xn(i,j)
	      y(i) =     a(i,1) - xn(i,j)
	    enddo
c  compute criteria
          t1 = x(2)*y(3)-y(2)*x(3)
	    t2 = y(1)*x(3)-x(1)*y(3)
          t3 = x(1)*y(2)-y(1)*x(2)		   
	    crit = t1*t1 + t2*t2 + t3*t3
	    if(abs(crit).lt.eps) then ! node j = node 3
	      con(4) = j + 4
            do i= 1,3
	        xlsub(i,3) = xn(i,j)
	      enddo	
	      flag(j) =.false.
		endif  
	  endif
      enddo




c node 2 is the remaining node
      do j= 1,4
        if(flag(j)) then
	      con(3) = j + 4
            do i= 1,3
	        xlsub(i,2) = xn(i,j)
	      enddo	
	      flag(j) =.false.
	  endif
      enddo

	end







      subroutine asseble_prisma_b_11(oflag,id, xn,xl,xlsub,con)
c---------------------------------------------------TCG 19.02.2003--71
c.... assembles prismatic subelement, where the discontinuity is 
c     a trinagle 
c
c  Declare variable types
c
c       id      node id (+1,-1)               
c       xn      nodes on cutted face (triangle)               
c       xl      nodes of tet
c       xlsub   nodes of the prismatic element in right format
c       con     connectivity of prisma
c            ___________________________________________________________
c           |    |         |        |        |        |        |        |
c           |  6 |  node 1 | node 2 | node 3 | node 4 | node 5 | node 6 |
c           |____|_________|________|________|________|________|________|
c
c   
c  element nodes:         1   -
c						2    |  local tetrahedral nodes
c						3	 |
c						4	_|
c						5    |
c						6    |  nodes on discontinuity
c						7    |
c						8   _|
c
c       flag     node flag
c                 .false.   tetraeder node is allready node of prisma
c
c--------------------------------------------------------------------71
      implicit none

	real*8 xn(3,4), xl(3,4), xlsub(3,6)
	real*8 y(3), x(3), eps, crit
	real*8 z(3), temp, t1, t2, t3, itemp
	real*8 a(3)

	integer id(4), i,j, oflag, con(7)

	logical flag(4)

      data eps/1.0d-10/

      
	con(1) = 6
	flag=.true.
c define node a
      do j= 1,4
	if(id(j).ne.oflag) then   ! node j is node a
	  do i= 1,3
	    a(i) = xl(i,j)
	    flag(j)= .false.
	  enddo
	endif
	enddo

c define nodes on disc  node 1, 2 and 5
c trial choice

      do i= 1,3
	  xlsub(i,1) = xn(i,1)
	  con(2) = 5
	  xlsub(i,4) = xn(i,2)
	  con(5) = 6
	  xlsub(i,5) = xn(i,3)
	  con(6) = 7
	enddo  
	    
c define 1,2 and 5 via right coordiante system
      do i= 1,3
	  x(i) = xlsub(i,5) - xlsub(i,1)
	  y(i) = xlsub(i,4) - xlsub(i,1)
	  z(i) =       a(i) - xlsub(i,1)
	enddo
c  compute criteria
      crit = (x(2)*y(3)-y(2)*x(3))*z(1) +
	1       (y(1)*x(3)-x(1)*y(3))*z(2) + 
     2	   (x(1)*y(2)-y(1)*x(2))*z(3)		   

      if (crit.lt.0.0d0) then  ! interchange 4 and 5
        do i= 1,3
	    temp       = xlsub(i,4) 
	    xlsub(i,4) = xlsub(i,5)
	    xlsub(i,5) = temp
	  enddo	  
	  itemp  = con(5)
	  con(5) = con(6)
	  con(6) = itemp
      endif



c node 2      
      do j= 1,4
        if(flag(j)) then
          do i= 1,3
	      x(i) = xl(i,j) - xlsub(i,1)
	      y(i) =    a(i) - xlsub(i,1)
	    enddo
c  compute criteria
          t1 = x(2)*y(3)-y(2)*x(3)
	    t2 = y(1)*x(3)-x(1)*y(3)
          t3 = x(1)*y(2)-y(1)*x(2)		   
	    crit = t1*t1 + t2*t2 + t3*t3
	    if(abs(crit).lt.eps) then ! node j = node 2
	      con(3) = j
            do i= 1,3
	        xlsub(i,2) = xl(i,j)
	      enddo	
	      flag(j) =.false.
		endif  
	  endif
      enddo



c node 3      
      do j= 1,4
        if(flag(j)) then
          do i= 1,3
	      x(i) = xl(i,j) - xlsub(i,4)
	      y(i) =    a(i) - xlsub(i,4)
	    enddo
c  compute criteria
          t1 = x(2)*y(3)-y(2)*x(3)
	    t2 = y(1)*x(3)-x(1)*y(3)
          t3 = x(1)*y(2)-y(1)*x(2)		   
	    crit = t1*t1 + t2*t2 + t3*t3
	    if(abs(crit).lt.eps) then ! node j = node 4
	      con(4) = j
            do i= 1,3
	        xlsub(i,3) = xl(i,j)
	      enddo	
	      flag(j) =.false.
		endif  
	  endif
      enddo


c node 6      remaining node
      do j= 1,4
        if(flag(j)) then
	    con(7) = j
          do i= 1,3
	      xlsub(i,6) = xl(i,j)
	    enddo	
	    flag(j) =.false.
	  endif
      enddo

	end







      subroutine gauss_prism_11(ll,lint,s)
c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: Gauss quadrature for 3-d prismatic element

c      Inputs:
c         ll     - Number of points/direction

c      Outputs:
c         lint   - Total number of quadrature points
c         s(4,*) - Gauss points (1-3) and weights (4)
c
c
c
c note: numbers are generated by Mathematica 
c
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   ll,lint
      real*8    s(4,*)


c     1 pt. quadrature

      if(ll.eq.1) then
        lint = 1
c  locations
        s(1,1) = 0.0d0
        s(2,1) = 0.0d0
        s(3,1) = 0.0d0
c  weigths
        s(4,1) = 8.0d0

c     2 x 2 x 2 pt. quadrature

      elseif(ll.eq.2) then
        lint = 8
  
c  location 1 point
        s(1,1) = -0.5773502691896258
        s(2,1) = -0.4553418012614796
        s(3,1) = -0.5773502691896258
c  weigth 1 point
        s(4,1) = 1.0d0

c  location 2 point
        s(1,2) =  0.5773502691896258
        s(2,2) = -0.4553418012614796
        s(3,2) = -0.5773502691896258
c  weigth 2 point
        s(4,2) = 1.0d0

c  location 3 point
        s(1,3) =  0.5773502691896258
        s(2,3) =  0.4553418012614796
        s(3,3) = -0.5773502691896258
c  weigth 3 point
        s(4,3) = 1.0d0

c  location 4 point
        s(1,4) = -0.5773502691896258
        s(2,4) =  0.4553418012614796
        s(3,4) = -0.5773502691896258
c  weigth 4 point
        s(4,4) = 1.0d0

c  location 5 point
        s(1,5) = -0.5773502691896258
        s(2,5) = -0.12200846792814621
        s(3,5) =  0.5773502691896258
c  weigth 5 point
        s(4,5) = 1.0d0

c  location 6 point
        s(1,6) =  0.5773502691896258
        s(2,6) = -0.12200846792814621
        s(3,6) =  0.5773502691896258
c  weigth 6 point
        s(4,6) = 1.0d0

c  location 7 point
        s(1,7) =  0.5773502691896258
        s(2,7) =  0.12200846792814621
        s(3,7) =  0.5773502691896258
c  weigth 7 point
        s(4,7) = 1.0d0

c  location 8 point
        s(1,8) = -0.5773502691896258
        s(2,8) =  0.12200846792814621
        s(3,8) =  0.5773502691896258
c  weigth 8 point
        s(4,8) = 1.0d0

c     3 x 3 x 3 pt. quadrature
      elseif(ll.eq.3) then
        lint = 27
	s(   1   ,   1    )=   -0.7745966692414834
	s(   2   ,   1    )=   -0.6872983346207417
	s(   3   ,   1    )=   -0.7745966692414834
	s(   4   ,   1    )=   0.1714677640603567
	s(   1   ,   2    )=   0.
	s(   2   ,   2    )=   -0.6872983346207417
	s(   3   ,   2    )=   -0.7745966692414834
	s(   4   ,   2    )=   0.2743484224965706
	s(   1   ,   3    )=   0.7745966692414834
	s(   2   ,   3    )=   -0.6872983346207417
	s(   3   ,   3    )=   -0.7745966692414834
	s(   4   ,   3    )=   0.1714677640603567
	s(   1   ,   4    )=   -0.7745966692414834
	s(   2   ,   4    )=   0.
	s(   3   ,   4    )=   -0.7745966692414834
	s(   4   ,   4    )=   0.2743484224965706
	s(   1   ,   5    )=   0.
	s(   2   ,   5    )=   0.
	s(   3   ,   5    )=   -0.7745966692414834
	s(   4   ,   5    )=   0.4389574759945128
	s(   1   ,   6    )=   0.7745966692414834
	s(   2   ,   6    )=   0.
	s(   3   ,   6    )=   -0.7745966692414834
	s(   4   ,   6    )=   0.2743484224965706
	s(   1   ,   7    )=   -0.7745966692414834
	s(   2   ,   7    )=   0.6872983346207417
	s(   3   ,   7    )=   -0.7745966692414834
	s(   4   ,   7    )=   0.1714677640603567
	s(   1   ,   8    )=   0.
	s(   2   ,   8    )=   0.6872983346207417
	s(   3   ,   8    )=   -0.7745966692414834
	s(   4   ,   8    )=   0.2743484224965706
	s(   1   ,   9    )=   0.7745966692414834
	s(   2   ,   9    )=   0.6872983346207417
	s(   3   ,   9    )=   -0.7745966692414834
	s(   4   ,   9    )=   0.1714677640603567
	s(   1   ,   10   )=   -0.7745966692414834
	s(   2   ,   10   )=   -0.3872983346207417
	s(   3   ,   10   )=   0.
	s(   4   ,   10   )=   0.2743484224965706
	s(   1   ,   11   )=   0.
	s(   2   ,   11   )=   -0.3872983346207417
	s(   3   ,   11   )=   0.
	s(   4   ,   11   )=   0.4389574759945128
	s(   1   ,   12   )=   0.7745966692414834
	s(   2   ,   12   )=   -0.3872983346207417
	s(   3   ,   12   )=   0.
	s(   4   ,   12   )=   0.2743484224965706
	s(   1   ,   13   )=   -0.7745966692414834
	s(   2   ,   13   )=   0.
	s(   3   ,   13   )=   0.
	s(   4   ,   13   )=   0.4389574759945128
	s(   1   ,   14   )=   0.
	s(   2   ,   14   )=   0.
	s(   3   ,   14   )=   0.
	s(   4   ,   14   )=   0.7023319615912202
	s(   1   ,   15   )=   0.7745966692414834
	s(   2   ,   15   )=   0.
	s(   3   ,   15   )=   0.
	s(   4   ,   15   )=   0.4389574759945128
	s(   1   ,   16   )=   -0.7745966692414834
	s(   2   ,   16   )=   0.3872983346207417
	s(   3   ,   16   )=   0.
	s(   4   ,   16   )=   0.2743484224965706
	s(   1   ,   17   )=   0.
	s(   2   ,   17   )=   0.3872983346207417
	s(   3   ,   17   )=   0.
	s(   4   ,   17   )=   0.4389574759945128
	s(   1   ,   18   )=   0.7745966692414834
	s(   2   ,   18   )=   0.3872983346207417
	s(   3   ,   18   )=   0.
	s(   4   ,   18   )=   0.2743484224965706
	s(   1   ,   19   )=   -0.7745966692414834
	s(   2   ,   19   )=   -0.08729833462074169
	s(   3   ,   19   )=   0.7745966692414834
	s(   4   ,   19   )=   0.1714677640603567
	s(   1   ,   20   )=   0.
	s(   2   ,   20   )=   -0.08729833462074169
	s(   3   ,   20   )=   0.7745966692414834
	s(   4   ,   20   )=   0.2743484224965706
	s(   1   ,   21   )=   0.7745966692414834
	s(   2   ,   21   )=   -0.08729833462074169
	s(   3   ,   21   )=   0.7745966692414834
	s(   4   ,   21   )=   0.1714677640603567
	s(   1   ,   22   )=   -0.7745966692414834
	s(   2   ,   22   )=   0.
	s(   3   ,   22   )=   0.7745966692414834
	s(   4   ,   22   )=   0.2743484224965706
	s(   1   ,   23   )=   0.
	s(   2   ,   23   )=   0.
	s(   3   ,   23   )=   0.7745966692414834
	s(   4   ,   23   )=   0.4389574759945128
	s(   1   ,   24   )=   0.7745966692414834
	s(   2   ,   24   )=   0.
	s(   3   ,   24   )=   0.7745966692414834
	s(   4   ,   24   )=   0.2743484224965706
	s(   1   ,   25   )=   -0.7745966692414834
	s(   2   ,   25   )=   0.08729833462074169
	s(   3   ,   25   )=   0.7745966692414834
	s(   4   ,   25   )=   0.1714677640603567
	s(   1   ,   26   )=   0.
	s(   2   ,   26   )=   0.08729833462074169
	s(   3   ,   26   )=   0.7745966692414834
	s(   4   ,   26   )=   0.2743484224965706
	s(   1   ,   27   )=   0.7745966692414834
	s(   2   ,   27   )=   0.08729833462074169
	s(   3   ,   27   )=   0.7745966692414834
	s(   4   ,   27   )=   0.1714677640603567
c     4 x 4 x 4 pt. quadrature
      elseif(ll.eq.4) then
        lint = 64
       s(   1   ,   1    )=   -0.8611363115940526
       s(   2   ,   1    )=   -0.8013460293699308
       s(   3   ,   1    )=   -0.8611363115940526
       s(   4   ,   1    )=   0.042091477490531534
       s(   1   ,   2    )=   -0.3399810435848563
       s(   2   ,   2    )=   -0.8013460293699308
       s(   3   ,   2    )=   -0.8611363115940526
       s(   4   ,   2    )=   0.07891151579507064
       s(   1   ,   3    )=   0.3399810435848563
       s(   2   ,   3    )=   -0.8013460293699308
       s(   3   ,   3    )=   -0.8611363115940526
       s(   4   ,   3    )=   0.07891151579507064
       s(   1   ,   4    )=   0.8611363115940526
       s(   2   ,   4    )=   -0.8013460293699308
       s(   3   ,   4    )=   -0.8611363115940526
       s(   4   ,   4    )=   0.042091477490531534
       s(   1   ,   5    )=   -0.8611363115940526
       s(   2   ,   5    )=   -0.31637553273470814
       s(   3   ,   5    )=   -0.8611363115940526
       s(   4   ,   5    )=   0.07891151579507064
       s(   1   ,   6    )=   -0.3399810435848563
       s(   2   ,   6    )=   -0.31637553273470814
       s(   3   ,   6    )=   -0.8611363115940526
       s(   4   ,   6    )=   0.14794033605678134
       s(   1   ,   7    )=   0.3399810435848563
       s(   2   ,   7    )=   -0.31637553273470814
       s(   3   ,   7    )=   -0.8611363115940526
       s(   4   ,   7    )=   0.14794033605678134
       s(   1   ,   8    )=   0.8611363115940526
       s(   2   ,   8    )=   -0.31637553273470814
       s(   3   ,   8    )=   -0.8611363115940526
       s(   4   ,   8    )=   0.07891151579507064
       s(   1   ,   9    )=   -0.8611363115940526
       s(   2   ,   9    )=   0.31637553273470814
       s(   3   ,   9    )=   -0.8611363115940526
       s(   4   ,   9    )=   0.07891151579507064
       s(   1   ,   10   )=   -0.3399810435848563
       s(   2   ,   10   )=   0.31637553273470814
       s(   3   ,   10   )=   -0.8611363115940526
       s(   4   ,   10   )=   0.14794033605678134
       s(   1   ,   11   )=   0.3399810435848563
       s(   2   ,   11   )=   0.31637553273470814
       s(   3   ,   11   )=   -0.8611363115940526
       s(   4   ,   11   )=   0.14794033605678134
       s(   1   ,   12   )=   0.8611363115940526
       s(   2   ,   12   )=   0.31637553273470814
       s(   3   ,   12   )=   -0.8611363115940526
       s(   4   ,   12   )=   0.07891151579507064
       s(   1   ,   13   )=   -0.8611363115940526
       s(   2   ,   13   )=   0.8013460293699308
       s(   3   ,   13   )=   -0.8611363115940526
       s(   4   ,   13   )=   0.042091477490531534
       s(   1   ,   14   )=   -0.3399810435848563
       s(   2   ,   14   )=   0.8013460293699308
       s(   3   ,   14   )=   -0.8611363115940526
       s(   4   ,   14   )=   0.07891151579507064
       s(   1   ,   15   )=   0.3399810435848563
       s(   2   ,   15   )=   0.8013460293699308
       s(   3   ,   15   )=   -0.8611363115940526
       s(   4   ,   15   )=   0.07891151579507064
       s(   1   ,   16   )=   0.8611363115940526
       s(   2   ,   16   )=   0.8013460293699308
       s(   3   ,   16   )=   -0.8611363115940526
       s(   4   ,   16   )=   0.042091477490531534
       s(   1   ,   17   )=   -0.8611363115940526
       s(   2   ,   17   )=   -0.5769531667393063
       s(   3   ,   17   )=   -0.3399810435848563
       s(   4   ,   17   )=   0.07891151579507064
       s(   1   ,   18   )=   -0.3399810435848563
       s(   2   ,   18   )=   -0.5769531667393063
       s(   3   ,   18   )=   -0.3399810435848563
       s(   4   ,   18   )=   0.1479403360567813
       s(   1   ,   19   )=   0.3399810435848563
       s(   2   ,   19   )=   -0.5769531667393063
       s(   3   ,   19   )=   -0.3399810435848563
       s(   4   ,   19   )=   0.1479403360567813
       s(   1   ,   20   )=   0.8611363115940526
       s(   2   ,   20   )=   -0.5769531667393063
       s(   3   ,   20   )=   -0.3399810435848563
       s(   4   ,   20   )=   0.07891151579507064
       s(   1   ,   21   )=   -0.8611363115940526
       s(   2   ,   21   )=   -0.22778407679095214
       s(   3   ,   21   )=   -0.3399810435848563
       s(   4   ,   21   )=   0.14794033605678134
       s(   1   ,   22   )=   -0.3399810435848563
       s(   2   ,   22   )=   -0.22778407679095214
       s(   3   ,   22   )=   -0.3399810435848563
       s(   4   ,   22   )=   0.27735296695391276
       s(   1   ,   23   )=   0.3399810435848563
       s(   2   ,   23   )=   -0.22778407679095214
       s(   3   ,   23   )=   -0.3399810435848563
       s(   4   ,   23   )=   0.27735296695391276
       s(   1   ,   24   )=   0.8611363115940526
       s(   2   ,   24   )=   -0.22778407679095214
       s(   3   ,   24   )=   -0.3399810435848563
       s(   4   ,   24   )=   0.14794033605678134
       s(   1   ,   25   )=   -0.8611363115940526
       s(   2   ,   25   )=   0.22778407679095214
       s(   3   ,   25   )=   -0.3399810435848563
       s(   4   ,   25   )=   0.14794033605678134
       s(   1   ,   26   )=   -0.3399810435848563
       s(   2   ,   26   )=   0.22778407679095214
       s(   3   ,   26   )=   -0.3399810435848563
       s(   4   ,   26   )=   0.27735296695391276
       s(   1   ,   27   )=   0.3399810435848563
       s(   2   ,   27   )=   0.22778407679095214
       s(   3   ,   27   )=   -0.3399810435848563
       s(   4   ,   27   )=   0.27735296695391276
       s(   1   ,   28   )=   0.8611363115940526
       s(   2   ,   28   )=   0.22778407679095214
       s(   3   ,   28   )=   -0.3399810435848563
       s(   4   ,   28   )=   0.14794033605678134
       s(   1   ,   29   )=   -0.8611363115940526
       s(   2   ,   29   )=   0.5769531667393063
       s(   3   ,   29   )=   -0.3399810435848563
       s(   4   ,   29   )=   0.07891151579507064
       s(   1   ,   30   )=   -0.3399810435848563
       s(   2   ,   30   )=   0.5769531667393063
       s(   3   ,   30   )=   -0.3399810435848563
       s(   4   ,   30   )=   0.1479403360567813
       s(   1   ,   31   )=   0.3399810435848563
       s(   2   ,   31   )=   0.5769531667393063
       s(   3   ,   31   )=   -0.3399810435848563
       s(   4   ,   31   )=   0.1479403360567813
       s(   1   ,   32   )=   0.8611363115940526
       s(   2   ,   32   )=   0.5769531667393063
       s(   3   ,   32   )=   -0.3399810435848563
       s(   4   ,   32   )=   0.07891151579507064
       s(   1   ,   33   )=   -0.8611363115940526
       s(   2   ,   33   )=   -0.28418314485474633
       s(   3   ,   33   )=   0.3399810435848563
       s(   4   ,   33   )=   0.07891151579507064
       s(   1   ,   34   )=   -0.3399810435848563
       s(   2   ,   34   )=   -0.28418314485474633
       s(   3   ,   34   )=   0.3399810435848563
       s(   4   ,   34   )=   0.1479403360567813
       s(   1   ,   35   )=   0.3399810435848563
       s(   2   ,   35   )=   -0.28418314485474633
       s(   3   ,   35   )=   0.3399810435848563
       s(   4   ,   35   )=   0.1479403360567813
       s(   1   ,   36   )=   0.8611363115940526
       s(   2   ,   36   )=   -0.28418314485474633
       s(   3   ,   36   )=   0.3399810435848563
       s(   4   ,   36   )=   0.07891151579507064
       s(   1   ,   37   )=   -0.8611363115940526
       s(   2   ,   37   )=   -0.11219696679390419
       s(   3   ,   37   )=   0.3399810435848563
       s(   4   ,   37   )=   0.14794033605678134
       s(   1   ,   38   )=   -0.3399810435848563
       s(   2   ,   38   )=   -0.11219696679390419
       s(   3   ,   38   )=   0.3399810435848563
       s(   4   ,   38   )=   0.27735296695391276
       s(   1   ,   39   )=   0.3399810435848563
       s(   2   ,   39   )=   -0.11219696679390419
       s(   3   ,   39   )=   0.3399810435848563
       s(   4   ,   39   )=   0.27735296695391276
       s(   1   ,   40   )=   0.8611363115940526
       s(   2   ,   40   )=   -0.11219696679390419
       s(   3   ,   40   )=   0.3399810435848563
       s(   4   ,   40   )=   0.14794033605678134
       s(   1   ,   41   )=   -0.8611363115940526
       s(   2   ,   41   )=   0.11219696679390419
       s(   3   ,   41   )=   0.3399810435848563
       s(   4   ,   41   )=   0.14794033605678134
       s(   1   ,   42   )=   -0.3399810435848563
       s(   2   ,   42   )=   0.11219696679390419
       s(   3   ,   42   )=   0.3399810435848563
       s(   4   ,   42   )=   0.27735296695391276
       s(   1   ,   43   )=   0.3399810435848563
       s(   2   ,   43   )=   0.11219696679390419
       s(   3   ,   43   )=   0.3399810435848563
       s(   4   ,   43   )=   0.27735296695391276
       s(   1   ,   44   )=   0.8611363115940526
       s(   2   ,   44   )=   0.11219696679390419
       s(   3   ,   44   )=   0.3399810435848563
       s(   4   ,   44   )=   0.14794033605678134
       s(   1   ,   45   )=   -0.8611363115940526
       s(   2   ,   45   )=   0.28418314485474633
       s(   3   ,   45   )=   0.3399810435848563
       s(   4   ,   45   )=   0.07891151579507064
       s(   1   ,   46   )=   -0.3399810435848563
       s(   2   ,   46   )=   0.28418314485474633
       s(   3   ,   46   )=   0.3399810435848563
       s(   4   ,   46   )=   0.1479403360567813
       s(   1   ,   47   )=   0.3399810435848563
       s(   2   ,   47   )=   0.28418314485474633
       s(   3   ,   47   )=   0.3399810435848563
       s(   4   ,   47   )=   0.1479403360567813
       s(   1   ,   48   )=   0.8611363115940526
       s(   2   ,   48   )=   0.28418314485474633
       s(   3   ,   48   )=   0.3399810435848563
       s(   4   ,   48   )=   0.07891151579507064
       s(   1   ,   49   )=   -0.8611363115940526
       s(   2   ,   49   )=   -0.05979028222412169
       s(   3   ,   49   )=   0.8611363115940526
       s(   4   ,   49   )=   0.042091477490531534
       s(   1   ,   50   )=   -0.3399810435848563
       s(   2   ,   50   )=   -0.05979028222412169
       s(   3   ,   50   )=   0.8611363115940526
       s(   4   ,   50   )=   0.07891151579507064
       s(   1   ,   51   )=   0.3399810435848563
       s(   2   ,   51   )=   -0.05979028222412169
       s(   3   ,   51   )=   0.8611363115940526
       s(   4   ,   51   )=   0.07891151579507064
       s(   1   ,   52   )=   0.8611363115940526
       s(   2   ,   52   )=   -0.05979028222412169
       s(   3   ,   52   )=   0.8611363115940526
       s(   4   ,   52   )=   0.042091477490531534
       s(   1   ,   53   )=   -0.8611363115940526
       s(   2   ,   53   )=   -0.02360551085014816
       s(   3   ,   53   )=   0.8611363115940526
       s(   4   ,   53   )=   0.07891151579507064
       s(   1   ,   54   )=   -0.3399810435848563
       s(   2   ,   54   )=   -0.02360551085014816
       s(   3   ,   54   )=   0.8611363115940526
       s(   4   ,   54   )=   0.14794033605678134
       s(   1   ,   55   )=   0.3399810435848563
       s(   2   ,   55   )=   -0.02360551085014816
       s(   3   ,   55   )=   0.8611363115940526
       s(   4   ,   55   )=   0.14794033605678134
       s(   1   ,   56   )=   0.8611363115940526
       s(   2   ,   56   )=   -0.02360551085014816
       s(   3   ,   56   )=   0.8611363115940526
       s(   4   ,   56   )=   0.07891151579507064
       s(   1   ,   57   )=   -0.8611363115940526
       s(   2   ,   57   )=   0.02360551085014816
       s(   3   ,   57   )=   0.8611363115940526
       s(   4   ,   57   )=   0.07891151579507064
       s(   1   ,   58   )=   -0.3399810435848563
       s(   2   ,   58   )=   0.02360551085014816
       s(   3   ,   58   )=   0.8611363115940526
       s(   4   ,   58   )=   0.14794033605678134
       s(   1   ,   59   )=   0.3399810435848563
       s(   2   ,   59   )=   0.02360551085014816
       s(   3   ,   59   )=   0.8611363115940526
       s(   4   ,   59   )=   0.14794033605678134
       s(   1   ,   60   )=   0.8611363115940526
       s(   2   ,   60   )=   0.02360551085014816
       s(   3   ,   60   )=   0.8611363115940526
       s(   4   ,   60   )=   0.07891151579507064
       s(   1   ,   61   )=   -0.8611363115940526
       s(   2   ,   61   )=   0.05979028222412169
       s(   3   ,   61   )=   0.8611363115940526
       s(   4   ,   61   )=   0.042091477490531534
       s(   1   ,   62   )=   -0.3399810435848563
       s(   2   ,   62   )=   0.05979028222412169
       s(   3   ,   62   )=   0.8611363115940526
       s(   4   ,   62   )=   0.07891151579507064
       s(   1   ,   63   )=   0.3399810435848563
       s(   2   ,   63   )=   0.05979028222412169
       s(   3   ,   63   )=   0.8611363115940526
       s(   4   ,   63   )=   0.07891151579507064
       s(   1   ,   64   )=   0.8611363115940526
       s(   2   ,   64   )=   0.05979028222412169
       s(   3   ,   64   )=   0.8611363115940526
       s(   4   ,   64   )=   0.042091477490531534
	else
        write(iow,2000) ll
        if(ior.lt.0) then
          write(*,2000) ll
        endif
        call plstop()
      endif

c     Format

2000  format('  *ERROR* Illegal quadrature order =',i16)

      end





      subroutine jac08 (F, jac)
c--------------------------------------------------------------------71
c
c    compute determinat of 3x3 tensor
c
c--------------------------------------------------------------------71
      real*8 f(3,3), jac

       jac  = f(1,1)*f(2,2)*f(3,3) + f(1,2)*f(2,3)*f(3,1)
     1      + f(1,3)*f(2,1)*f(3,2) - f(3,1)*f(2,2)*f(1,3)
     2      - f(3,2)*f(2,3)*f(1,1) - f(3,3)*f(2,1)*f(1,2)
	end



      subroutine pushv10 (fi,vect)

c-------------------------------------------------------------------
c     Push-forward a 1 cova vector.  fi = invF
c     pull-back a 1 cova vector.     fi = F
c-------------------------------------------------------------------
 
      integer i,j
      real*8  fi(3,3),vect(3),temp(3)


      do i =1,3
        temp(i) = 0.d0
        do j = 1,3
          temp(i) = temp(i) + fi(j,i)*vect(j)
        end do
      end do

      do i = 1,3
        vect(i) = temp(i)
      end do

      end

      subroutine pushv10_contra (fi,vect)

c-------------------------------------------------------------------
c     Push-forward a 1 cova vector.  fi = invF
c     pull-back a 1 cova vector.     fi = F
c-------------------------------------------------------------------
 
      integer i,j
      real*8  fi(3,3),vect(3),temp(3)


      do i =1,3
        temp(i) = 0.d0
        do j = 1,3
          temp(i) = temp(i) + fi(i,j)*vect(j)
        end do
      end do

      do i = 1,3
        vect(i) = temp(i)
      end do

      end



      subroutine prinb10(bb,bpr)
c
c--------------------------------------------------------------------71
c
c    Principal Values of a 3D Symmetric Positive Definite Tensor FEAP
c
c.... INPUT variables
c        bb(6)  Symmetric 3x3 tensor (compressed storage)
c
c.... OUTPUT variables
c        bpr(3) Principal values 
c
c--------------------------------------------------------------------71
c
c..... Declare variable types
      integer i
      real*8 pi23, tol, al, b1, b2, b3, c1, c2,c3,temp
c..... Declare array types
      real*8  bd(6),bb(6),bpr(3)
c..... Intrinsics
      intrinsic atan2, cos, sqrt
c     data pi4,pi23/0.7853981633974483d0,2.0943951023931955d0/
      data pi23/2.0943951023931955d0/ , tol /1.d-10/
c
c.... Compute mean and deviatoric (upper trianglular part) tensors
c
      b1  = (bb(1) + bb(2) + bb(3))/3.
      do 100 i = 1,6
        bd(i) = bb(i)
100   continue
      do 110 i = 1,3
        bd(i) = bd(i) - b1
110   continue
c
c.... Compute 2nd and 3rd invariants of deviator
c
      c1 = bd(4)*bd(4)
      c2 = bd(5)*bd(5)
      c3 = bd(6)*bd(6)
      b2 = 0.5d0*(bd(1)*bd(1)+bd(2)*bd(2)+bd(3)*bd(3))
     1   + c1 + c2 + c3
      if(b2.le.tol*b1*b1) then
        bpr(1) = b1
        bpr(2) = b1
        bpr(3) = b1
      else
        b3 = bd(1)*bd(2)*bd(3)+(bd(4)+bd(4))*bd(5)*bd(6)
     1     + bd(1)*(c1-c2) + bd(2)*(c1-c3)
c
c.... Set constants
c
        c1 = 2.d0*sqrt(b2/3.d0)
        c2 = 4.d0*b3
        c3 = c1*c1*c1
        temp = sqrt(abs(c3*c3 - c2*c2))
        al = atan2(temp,c2)/3.d0
c
c.... Set principal values
c
        bpr(1) = b1 + c1*cos(al)
        bpr(2) = b1 + c1*cos(al-pi23)
        bpr(3) = b1 + c1*cos(al+pi23)
      endif
c
      end



      subroutine stcn10tetra(ix,x,xsj,shp,sigl,dt,st,ndm,nel,nnp)

c--------------------------------------------------------------------71
c.... load (stress) quantities for plotting for tetraeder element
c     (suitable for 1 Gausspoint)
c
c  Declare variable types
c
c       ix    Element connections
c       x     spatial coordinates of element nodes
c       sigl  Vector of quantities to be plotted
c       shp   shape function    
c       xsj   dv    
c       dt    
c       st    Vector of nodal quantities to be plotted
c       ndm   Number of dimension
c       nel   Nodes per element
c       nnp   Number of nodal points
c--------------------------------------------------------------------71
      implicit none 
	 
	integer ndm,nel,nnp, i,k,ll

      integer ix(nel)
      real*8  dt(nnp),st(nnp,10),x(ndm,nel)
      real*8  shp(4,10),sigl(20), xj
      real*8  xsj


      save


      do i= 1,nel
        ll = iabs(ix(i))
	  xj =  xsj*shp(4,i)
	  dt(ll) = dt(ll) + xj
	  do k=1,10 
         st(ll,k) = st(ll,k) + sigl(k)*xj
	  enddo
	enddo


      end

	subroutine leftCauchy_11(F,b)
	integer i, j, k 
	real*8 f(3,3), bb(3,3), b(6)
c left cauchy green tensor
      bb=0.0d0
	do i= 1,3
	do j= 1,3
	 do k= 1,3
	bb(i,j) = bb(i,j) + f(i,k)*f(j,k)
	 enddo
	enddo
	enddo
	b(1) = bb(1,1)
	b(2) = bb(2,2)
	b(3) = bb(3,3)
	b(4) = bb(1,2)
	b(5) = bb(2,3)
	b(6) = bb(1,3)
	end

	subroutine leftCauchy_mat_11(F,bb)

	integer i,j,k
	real*8 f(3,3), bb(3,3)
c left cauchy green tensor
      bb=0.0d0
	do i= 1,3
	do j= 1,3
	 do k= 1,3
	  bb(i,j) = bb(i,j) + f(i,k)*f(j,k)
	 enddo
	enddo
	enddo
	end


	subroutine quad (a,aa)
c
c-------------------------------------------------------------------
c
c  Coomputes quadrat of a(6)
c   
c  Declare variable types
c    a(6)        arbitrary sym. tensor
c    aa(6)       a*a
c
c-------------------------------------------------------------------

	real*8 a(6), aa(6)

	aa(1) = a(1)*a(1) + a(4)*a(4) + a(6)*a(6)
	aa(4) = a(1)*a(4) + a(4)*a(2) + a(6)*a(5)
	aa(6) = a(1)*a(6) + a(4)*a(5) + a(6)*a(3)
	aa(2) = a(4)*a(4) + a(2)*a(2) + a(5)*a(5)
	aa(5) = a(4)*a(6) + a(2)*a(5) + a(5)*a(3)
	aa(3) = a(6)*a(6) + a(5)*a(5) + a(3)*a(3)

	end



      subroutine det_AA_11(A, det)
c--------------------------------------------------------------------71
c.... compute determinant of 3x3 matrix a
c
c  Declare variable types
c
c       det   det a  
c
c--------------------------------------------------------------------71
	implicit none
      real*8 det, A(3,3)

	det = a(1,1)*(a(2,2)*a(3,3) - a(3,2)*a(2,3)) -
	1      a(1,2)*(a(2,1)*a(3,3) - a(3,1)*a(2,3)) +
	2      a(1,3)*(a(2,1)*a(3,2) - a(3,1)*a(2,2))
	end





      subroutine flag_convert10(temp,lnode_flag)
c--------------------------------------------------TCG.28.08.2002----71
c
c    converts integernode flag into logical node flag
c
c--------------------------------------------------------------------71
      implicit none

	logical  lnode_flag(4)
	integer node_flag, i, tt, iitemp, temp

	node_flag = temp
	lnode_flag=.false.

      do i=4,1,-1
	 iitemp = 10**(i-1) 
	 tt=node_flag/iitemp
	 if ((tt.eq.1))  then
	   lnode_flag(i)=.true.
	   node_flag = node_flag - iitemp 
	 endif  
	enddo

	end


      subroutine kine_11m (ndm, shp, ul, f,fi,detfi)
c-----------------------------------------------------TCG.03.02.2003-71
c
c             Compute 3D compatible deformation gradient
c
c
c         shp     Shape functions and derivatives of compatible modes
c         ul      Solution array (displacements)
c         f       Deformation gradient
c         fi      Inverse deformation gradient
c         detfi   Determinat of f
c         deti    1/detfi 
c         
c
c--------------------------------------------------------------------71
c
      include  'eldata.h'

      integer i,j,k, ndm

      real*8  deti,detfi
      real*8  shp(4,nel),ul(ndm,nel)
      real*8  f(3,3),fi(3,3)

c.... compute compatible deformation gradient 
c........F = I + GRAD u
 
      do i = 1,ndm
        do j = 1,ndm
          f(i,j) = 0.0d0
          do k = 1,nel
            f(i,j) = f(i,j) + ul(i,k)*shp(j,k)
          end do
        end do
        f(i,i) = f(i,i) + 1.0d0
      end do


c...  Invert F

      detfi  = f(1,1)*f(2,2)*f(3,3) + f(1,2)*f(2,3)*f(3,1)
     1       + f(1,3)*f(2,1)*f(3,2) - f(3,1)*f(2,2)*f(1,3)
     2       - f(3,2)*f(2,3)*f(1,1) - f(3,3)*f(2,1)*f(1,2)
      deti   = 1.d0/detfi

      fi(1,1) = (f(2,2)*f(3,3) - f(3,2)*f(2,3))*deti
      fi(1,2) =-(f(1,2)*f(3,3) - f(3,2)*f(1,3))*deti
      fi(1,3) = (f(1,2)*f(2,3) - f(2,2)*f(1,3))*deti
      fi(2,1) =-(f(2,1)*f(3,3) - f(3,1)*f(2,3))*deti
      fi(2,2) = (f(1,1)*f(3,3) - f(3,1)*f(1,3))*deti
      fi(2,3) =-(f(1,1)*f(2,3) - f(2,1)*f(1,3))*deti
      fi(3,1) = (f(2,1)*f(3,2) - f(3,1)*f(2,2))*deti
      fi(3,2) =-(f(1,1)*f(3,2) - f(3,1)*f(1,2))*deti
      fi(3,3) = (f(1,1)*f(2,2) - f(2,1)*f(1,2))*deti

      end



      subroutine pushv_11 (fi,vect)

c-------------------------------------------------------------------
c     Push-forward a 1 cova vector.
c-------------------------------------------------------------------
 
      integer i,j
      real*8  fi(3,3),vect(3),temp(3)


      do i =1,3
        temp(i) = 0.d0
        do j = 1,3
          temp(i) = temp(i) + fi(j,i)*vect(j)
        end do
      end do

      do i = 1,3
        vect(i) = temp(i)
      end do

      end



      subroutine vc_int_subtet_11 (flag,l,xid,id,fact,s,srel,idg,lint)
c-----------------------------------------------------TCG.11.02.2003-71
c
c             Compute volume coordinates of Guasspoints in subtet
c             with respect to global tetrahedral element
c
c
c         falg     +1 for Omega+ , -1 for Omega-
c         l        integration order
c         xid      Volume cvoordiantes of corner nodes of diec
c         id       id of nodes of global tet
c         fact     V+/V-
c         s        Gausspoints and weights
c         srel     relative Gausspoints and weights
c         idg      id of Gausspoints
c         Ltet     Volume coordinates of subtetcorners
c         sgtet    Gausspoints and weights in subtet
c         
c
c--------------------------------------------------------------------71
c
      implicit none
	 
	integer flag, i, j, k, id(4), l, lint, idg(4)

!	real*8  fact, s(5,4), srel(5,4), sgtet(5,4), ltet(4,4)
	real*8  fact, s(5,4), srel(5,4), sgtet(5,87), ltet(4,4)
	real*8  xid(4,3)

c get Gausspoints of tet
      lint=0
      call tint3d_11(l,lint,sgtet)
c assemble volume coordinates of subtet
      Ltet=0.0d0
      do i= 1,3       ! loop over points on disc
        do j= 1,4        ! loop over volumecoordiantes
	    Ltet(j,i) = xid(j,i)
	  enddo
	enddo
c add corner node
      do i= 1,4        ! loop over nodes
	  if(id(i).eq.flag) then
	    Ltet(i,4) = 1.0d0
  	   endif
	enddo
c transform sgtet in volume coordinates of global tet
      do i= 1,lint  ! loop over gausspoints
	  do j= 1,4      ! loop over volume coordinates
	    do k= 1,4    ! loop over subelement nodes
	      s(j,i) = s(j,i) + Ltet(j,k)*sgtet(k,i)
	     enddo
	  enddo
c weights
        s(5,i) = sgtet(5,i)
	  idg(i) = flag
 	enddo

      end




      subroutine vc_int_subpris_11(flag,l,xid,id,fact,cn, s,idg,lint)
c-----------------------------------------------------TCG.11.02.2003-71
c
c             Compute volume coordinates of Guasspoints in subpris
c             with respect to global tetrahedral element
c
c
c         falg     +1 for Omega+ , -1 for Omega-
c         l        integration order
c         cn       number of corners of disc
c         xid      Volume cvoordiantes of corner nodes of diec
c         id       id of nodes of global tet
c         fact     V+/V-
c         s        Gausspoints and weights
c         idg      id of Gausspoints
c         Lpris    Volume coordinates of subpriscorners
c         sgpris   Gausspoints and weights in subpris
c         
c
c--------------------------------------------------------------------71
c
      implicit none
	 
	integer flag, i, j, k, id(4), l, lint, idg(8), cn

	real*8  fact, s(5,8), sgpris(4,8), shppris(6,8), Lpris(4,6)
	real*8  xid(4,cn)

c get Gausspoints of prism
      lint=0
      call gauss_prism_11(l,lint,sgpris)
	do l= 1,lint
        call shp_prism_11m (sgpris(1,l), shppris(1,l))
	enddo
      
c assemble volume coordinates of subpris
      Lpris=0.0d0
      do i= 1,cn       ! loop over points on disc
        do j= 1,4        ! loop over volumecoordiantes
	    Lpris(j,i) = xid(j,i)
	  enddo
	enddo
c add corner nodes
      k = cn+1
      do i= 1,4        ! loop over nodes
	  if(id(i).eq.flag) then
	    Lpris(i,k) = 1.0d0
	    k = k+1
  	   endif
	enddo
c transform sgtet in volume coordinates of global tet
      do i= 1,lint  ! loop over gausspoints
	  do j= 1,4      ! loop over volume coordinates
	    do k= 1,6    ! loop over subelement nodes
	      s(j,i) = s(j,i) + Lpris(j,k)*shppris(k,i)
	     enddo
	  enddo
c weights
         s(5,i)=sgpris(4,i)
	  idg(i) = flag
 	enddo

      end

      subroutine shp_prism_11m (xi, shp)
c-----------------------------------------------------TCG.11.02.2003-71
c
c             Compute shapefunction of prismatic element
c             no derivatives are computed!!
c
c
c         xi    local coordinates
c         shp   shape functions (without derivatives)
c         
c
c--------------------------------------------------------------------71

      implicit none
	real*8  r,s,t, one8, one4, xi(3), shp(6)
	real*8  r8, s8, t8, rs8, rt8, st8, rst8
	real*8  r4, s4, t4, rs4, rt4, st4, rst4

      data one8/0.125d0/ 
      data one4/0.25d0/ 

c compute shapefunctions of prism
	  r = xi(1)
	  s = xi(2)
	  t = xi(3)

	  r8 = one8*r
	  s8 = one8*s
	  t8 = one8*t

	  rs8 = r8*s
	  rt8 = r8*t
	  st8 = s8*t

	  rst8 = rs8*t
	 
        r4 = r8 + r8
        s4 = s8 + s8
        t4 = t8 + t8

	  rs4 = rs8 + rs8
	  rt4 = rt8 + rt8
	  st4 = st8 + st8

	  rst4 = rst8 + rst8

c compute shapefunctions
        shp(1) =  one8 - r8 - s8 + rs8 - t8 + rt8 + st8 - rst8
        shp(2) =  one8 + r8 - s8 - rs8 - t8 - rt8 + st8 + rst8
        shp(3) =  one8 + r8 + s8 + rs8 - t8 - rt8 - st8 - rst8
        shp(4) =  one8 - r8 + s8 - rs8 - t8 + rt8 - st8 + rst8
        shp(5) =  one4 - r4 + t4 - rt4
        shp(6) =  one4 + r4 + t4 + rt4


      end


      subroutine array_to_real_11(con,rcon)
c-----------------------------------------------------TCG.13.02.2003-71
c
c             Compute connectivity reals
c
c
c      3    1 2 3 4 5 6
c    nodes
c
c
c--------------------------------------------------------------------71
      implicit none

	integer con(7), icon,  i
	real*8  rcon



c compute integer connectivity
	icon = 0
	do i= 1,7
	  icon = icon + con(8-i)*10**(i-1)
	enddo
      
c change data type
	rcon = Dfloat(icon)  
	end



      subroutine real_to_array_11(rcon, con)
c-----------------------------------------------------TCG.13.02.2003-71
c
c             Compute connectivity array
c
c
c      3    1 2 3 4 5 6
c    nodes
c
c
c--------------------------------------------------------------------71
      implicit none

	integer con(7), icon,  i
	real*8  rcon

c change data type
      icon = int(rcon)

c compute integer connectivity
	do i= 1,7
        con(i) = icon/10**(7-i)
	  icon = icon - con(i)*10**(7-i)   
      enddo
	end
      
	
	
	
	subroutine subJ_11 (nsub,xlsub,sg,sgrel, Jsub)
c-----------------------------------------------------TCG.17.02.2003-71
c
c             Compute determinant of subtransformation
c
c
c         nsub    nodes of subelement
c                    4    ... tetrahedral subelement
c                    6    ... pismatic subelement
c         xlsub    volume coordinates of subelement
c         sg       volume coordiante of Gausspoint
c         sgrel    relative coordiantes of Gausspoint
c         Jsub     determinant of subtransformation         
c
c--------------------------------------------------------------------71
c
      implicit none

	integer nsub

	real*8  xlsub(4,6), sg(4), sgrel(4)
	real*8  Jsub

c decide type of subelement
      if (nsub.eq.4) then       ! tetrahedral subelement
	   call subJ_tet_11(xlsub, sg, Jsub)
	elseif (nsub.eq.6) then   ! pismatic subelement
	   call subJ_pris_11(xlsub, sg, sgrel, Jsub)
	else
	  write (1000,*)
	endif

1000  format('  *ERROR* Illegal subelement')

      end



	subroutine subJ_tet_11(xi, sg, detJ)
c-----------------------------------------------------TCG.17.02.2003-71
c
c             Compute determinant of subtransformation for 
c             tetrahedral subelement
c
c
c         xlsub    volume coordinates of subelement
c         sg       volume coordiante of Gausspoint
c         Jsub     determinant of subtransformation         
c
c--------------------------------------------------------------------71
c
      implicit none


	real*8  xi(4,6), sg(4)
	real*8  detJ, J(3,3)

c assembly J
c first row
      J(1,1) = -xi(2,1) + xi(2,2)
      J(1,2) = -xi(3,1) + xi(3,2)
      J(1,3) = -xi(4,1) + xi(4,2)

c secound row
      J(2,1) = -xi(2,1) + xi(2,3)
      J(2,2) = -xi(3,1) + xi(3,3)
      J(2,3) = -xi(4,1) + xi(4,3)

c third row
      J(3,1) = -xi(2,1) + xi(2,4)
      J(3,2) = -xi(3,1) + xi(3,4)
      J(3,3) = -xi(4,1) + xi(4,4)

c compute determinat
      call det_AA_11(J, detJ)

      end




	subroutine subJ_pris_11(xi, sg, sgrel, detJ)
c-----------------------------------------------------TCG.17.02.2003-71
c
c             Compute determinant of subtransformation for 
c             prismatic subelement
c
c
c         xlsub    volume coordinates of subelement
c         sg       volume coordiante of Gausspoint
c         Jsub     determinant of subtransformation         
c
c--------------------------------------------------------------------71
c
      implicit none


	real*8  xi(4,6), sg(4), sgrel(3), r, s, t, one8, one4
	real*8  detJ, J(3,3)
      real*8 omsomt,opsomt,opt,omromt,opromt,omroms,oproms
	real*8 oprops,omrops,opr,omr


      data one8 /0.125d0/
      data one4 /0.25d0/


      r = sgrel(1)
      s = sgrel(2)
      t = sgrel(3)

c define constants
      omsomt = one8*(1.0d0-s)*(1.0d0-t)
      opsomt = one8*(1.0d0+s)*(1.0d0-t)
	opt    = one4*(1.0d0+t)

	omromt = one8*(1.0d0-r)*(1.0d0-t)
	opromt = one8*(1.0d0+r)*(1.0d0-t)

	omroms = one8*(1.0d0-r)*(1.0d0-s)
	oproms = one8*(1.0d0+r)*(1.0d0-s)
	oprops = one8*(1.0d0+r)*(1.0d0+s)
	omrops = one8*(1.0d0-r)*(1.0d0+s)
	opr    = one4*(1.0d0+r)
	omr    = one4*(1.0d0-r)



c assembly J
c first row
      J(1,1) = omsomt*(-xi(2,1)+xi(2,2)) + opsomt*(xi(2,3)-xi(2,4)) +
	1            opt*(-xi(2,5)+xi(2,6))

      J(1,2) = omsomt*(-xi(3,1)+xi(3,2)) + opsomt*(xi(3,3)-xi(3,4)) +
	1            opt*(-xi(3,5)+xi(3,6))

      J(1,3) = omsomt*(-xi(4,1)+xi(4,2)) + opsomt*(xi(4,3)-xi(4,4)) +
	1            opt*(-xi(4,5)+xi(4,6))

c secound row
      J(2,1) = omromt*(-xi(2,1)+xi(2,4)) + opromt*(-xi(2,2)+xi(2,3))
      J(2,2) = omromt*(-xi(3,1)+xi(3,4)) + opromt*(-xi(3,2)+xi(3,3))
      J(2,3) = omromt*(-xi(4,1)+xi(4,4)) + opromt*(-xi(4,2)+xi(4,3))

c third row
      J(3,1) = -omroms*xi(2,1) - oproms*xi(2,2) - oprops*xi(2,3) -
	1          omrops*xi(2,4) +    omr*xi(2,5) +    opr*xi(2,6)
      J(3,2) = -omroms*xi(3,1) - oproms*xi(3,2) - oprops*xi(3,3) -
	1          omrops*xi(3,4) +    omr*xi(3,5) +    opr*xi(3,6)
      J(3,3) = -omroms*xi(4,1) - oproms*xi(4,2) - oprops*xi(4,3) -
	1          omrops*xi(4,4) +    omr*xi(4,5) +    opr*xi(4,6)

c compute determinat
      call det_AA_11(J, detJ)

      end



      subroutine eldata_11 (nel,ndm, nst2, order)
c-----------------------------------------------------TCG.20.02.2003-71
c
c             Compute element data
c
c
c         nel     no. of nodes per element
c         ndm     no. of dimensions
c         nst2    no. of enhanced degrees of freedom
c         order   element order        
c
c--------------------------------------------------------------------71
      implicit none

      include  'iofile.h'

      integer nel, ndm, nst2, order

c retrieve element order
      if (nel.eq.4) then
	  order = 1   ! linear 
	elseif(nel.eq.10) then
        order = 2    ! quadratic
	else
        write(iow,1000) 
        if(ior.lt.0) then
          write(*,1000)
        endif
        call plstop()
     	endif


c degrees of enhanced freedoms
      nst2 = nel*ndm



c     Format

1000  format('  *ERROR* Element order not coded')


	end





      subroutine trishp_11(ndm,xl, xsj)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: compute jacobian determinant of linear Triangular element


c      Inputs:
c         xl(ndm,*) - Nodal coordinates for element
c         ndm       - Spatial dimension of mesh

c      Outputs:
c         xsj       - Jacobian determinant at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm
      real*8    x1,x2,x3, y1,y2,y3
      real*8    xl(ndm,*), xsj


c     Form Jacobian terms

      x1 = xl(1,1)
      x2 = xl(1,2)
      x3 = xl(1,3)

      y1 = xl(2,1)
      y2 = xl(2,2)
      y3 = xl(2,3)


      xsj  = x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)
      xsj  = 0.5d0*xsj


      end




      subroutine quadshp_11(ndm,xi,xl,xsj)

c-----[--.----+----.----+----.-----------------------------------------]
c     Purpose: compute jacobian determinant of bi-linear quadrilaterals


c      Inputs:
c         ndm       - Spatial dimension of mesh
c         xi(2)     - Natural coordinates of point
c         xl(ndm,*) - Nodal coordinates for element

c      Outputs:
c         xsj       - Jacobian determinant at point
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      integer   ndm
	real*8    xi(2)
      real*8    xo,xs,xt, yo,ys,yt
      real*8    s,t, xsj,xl(ndm,4)


      s = xi(1)
	t = xi(2)


c     Set up natural coordinate functions (times 4)

      xo =  xl(1,1)-xl(1,2)+xl(1,3)-xl(1,4)
      xs = -xl(1,1)+xl(1,2)+xl(1,3)-xl(1,4) + xo*t
      xt = -xl(1,1)-xl(1,2)+xl(1,3)+xl(1,4) + xo*s
      yo =  xl(2,1)-xl(2,2)+xl(2,3)-xl(2,4)
      ys = -xl(2,1)+xl(2,2)+xl(2,3)-xl(2,4) + yo*t
      yt = -xl(2,1)-xl(2,2)+xl(2,3)+xl(2,4) + yo*s

c     Compute jacobian (times 16)

      xsj = xs*yt - xt*ys

c     Divide jacobian by 16 (multiply by .0625)

      xsj = 0.0625d0*xsj

      end




      subroutine tint2d_11(l,lint,ss)

c-----------------------------------------------------TCG.30.04.2003-71
c      Purpose: Set gauss points and weights for triangular elements

c      Inputs:
c         l       - Number of gauss points indicator
c
c      Outputs:
c         lint    - Total number of points
c         s(4,*) - Area coordinate points and weights for quadrature
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   l, lint
      real*8    ss(4,*)
	real*8    a,b,c,d,e,f,g,h,p,q,r,s,t,u,v,w,x,y
	real*8    w1,w2,w3,w4,w5,w6,w7,w8


c  1 point, precision 1.

      if (l.eq.1) then

c define constants
        a = 0.333333333333333d0
        w = 1.0d0

c define total number of integration points
        lint = 1
c load integration points and weights
        ss(1,1) = a
        ss(2,1) = a
        ss(3,1) = a
        ss(4,1) = w


c  4 points, precision 3, Strang and Fix formula #3.

      elseif (l.eq.2) then

c define constants
        a =   0.3333333333333333d0
        b =   0.6d0
        c =   0.2d0
        d =  -0.5625d0
        e =   0.5208333333333333d0

c define total number of integration points
        lint = 4
c load integration points and weights
        ss(1,1:4) = (/ a,	b,	c,	c/) 
        ss(2,1:4) = (/ a,	c,	b,	c/) 
        ss(3,1:4) = (/ a,	c,	c,	b/) 
        ss(4,1:4) = (/ d,	e,	e,	e/) 



c  7 points, precision 5, Strang and Fix formula #7, Stroud T2:5-1

      elseif (l.eq.3) then

c define constants
        a = 0.3333333333333333d0
        b = 0.7974269853530873d0
        c = 0.1012865073234563d0
        d = 0.0597158717897698d0
        e = 0.4701420641051150d0
        u = 0.225d0
        v = 0.1259391805448271d0
        w = 0.1323941527885061d0

c define total number of integration points
        lint = 7
c load integration points and weights
        ss(1,1:7) = (/ a, b, c, c, d, e, e /)
        ss(2,1:7) = (/ a, c, b, c, e, d, e /)
        ss(3,1:7) = (/ a, c, c, b, e, e, d /)
        ss(4,1:7) = (/ u, v, v, v, w, w, w /)



c  13 points, precision 7, Strang and Fix, formula #10.

      elseif ( l.eq.4 ) then

c define constants
        a =  0.479308067841923d0
        b =  0.260345966079038d0
        c =  0.869739794195568d0
        d =  0.065130102902216d0
        e =  0.638444188569809d0
        f =  0.312865496004875d0
        g =  0.048690315425316d0
        h =  0.333333333333333d0
        t =  0.175615257433204d0
        u =  0.053347235608839d0
        v =  0.077113760890257d0
        w = -0.149570044467670d0

c define total number of integration points
        lint = 13

c load integration points and weights
        ss(1,1:13) = (/ a, b, b, c, d, d, e, e, f, f, g, g, h /)
        ss(2,1:13) = (/ b, a, b, d, c, d, f, g, e, g, e, f, h /)
        ss(3,1:13) = (/ b, b, a, d, d, c, g, f, g, e, f, e, h /)
        ss(4,1:13) = (/ t, t, t, u, u, u, v, v, v, v, v, v, w /)




c  19 points, precision 9.

      elseif ( l.eq.5 ) then

c define constants
        a = 0.333333333333333d0
        b = 0.02063496160252593d0
        c = 0.4896825191987370d0
        d = 0.1258208170141290d0
        e = 0.4370895914929355d0
        f = 0.6235929287619356d0
        g = 0.1882035356190322d0
        r = 0.9105409732110941d0
        s = 0.04472951339445297d0
        t = 0.7411985987844980d0
        u = 0.03683841205473626d0
        v = 0.22196288916076574d0

        w1 = 0.09713579628279610d0
        w2 = 0.03133470022713983d0
        w3 = 0.07782754100477543d0
        w4 = 0.07964773892720910d0
        w5 = 0.02557767565869810d0
        w6 = 0.04328353937728940d0

c define total number of integration points
        lint = 19

c load integration points and weights
       ss(1,1:19) = (/  a,  b,  c,  c,  d,  e,  e,  f,  g,  g,  r,
	1                  s,  s, t, t, u, u, v, v /)
       ss(2,1:19) = (/  a,  c,  b,  c,  e,  d,  e,  g,  f,  g,  s, 
	1                  r,  s, u, v, t, v, t, u /)
       ss(3,1:19) = (/  a,  c,  c,  b,  e,  e,  d,  g,  g,  f,  s, 
	1                  s,  r, v, u, v, t, u, t /)
       ss(4,1:19) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, 
	1                w5, w5, w6, w6, w6, w6, w6, w6 /)



c  28 points, precision 11.

      elseif ( l.eq.6 ) then

c define constants
        a = 0.333333333333333d0
        b = 0.9480217181434233d0
        c = 0.02598914092828833d0
        d = 0.8114249947041546d0
        e = 0.09428750264792270d0
        f = 0.01072644996557060d0
        g = 0.4946367750172147d0
        p = 0.5853132347709715d0
        q = 0.2073433826145142d0
        r = 0.1221843885990187d0
        s = 0.4389078057004907d0
        t = 0.6779376548825902d0
        u = 0.04484167758913055d0
        v = 0.27722066752827925d0
        w = 0.8588702812826364d0
        x = 0.0d0
        y = 0.1411297187173636d0

        w1 = 0.08797730116222190d0
        w2 = 0.008744311553736190d0
        w3 = 0.03808157199393533d0
        w4 = 0.01885544805613125d0
        w5 = 0.07215969754474100d0
        w6 = 0.06932913870553720d0
        w7 = 0.04105631542928860d0
        w8 = 0.007362383783300573d0

c define total number of integration points
        lint = 28

c load integration points and weights
        ss(1,1:28) =   (/  a, b, c, c, d, e, e, f, g, g, p, q, q, r,
	1                  s, s, t, t, u, u, v, v, w, w, x, x, y, y /)
        ss(2,1:28) =   (/  a, c, b, c, e, d, e, g, f, g, q, p, q, s,
	1                  r, s, u, v, t, v, t, u, x, y, w, y, w, x /)
        ss(3,1:28) =   (/  a, c, c, b, e, e, d, g, g, f, q, q, p, s, 
	1                  s, r, v, u, v, t, u, t, y, x, y, w, x, w /)

        ss(4,1:28) = (/ w1, w2, w2, w2, w3, w3, w3, w4, w4, w4, w5, 
	1                w5, w5, w6, w6, w6, w7, w7, w7, w7, w7, w7, 
     2               w8, w8, w8, w8, w8, w8 /)




c     Unspecified quadrature specified

      else
        write(  *,2000) l
        write(iow,2000) l
        write(ilg,2000) l
        lint    = -1
      endif

c     Format

2000  format(' *ERROR* TINT2D: Wrong quadrature, l =',i3)

      end







      subroutine int2d_11(l,lint,sg)

c-----------------------------------------------------TCG.02.05.2003-71
c      Purpose: Form Gauss points and weights for two dimensions

c      Inputs:
c         l       - Number of gauss points indicator
c
c      Outputs:
c         lint    - Total number of points
c         sg(3,*) - Array of points and weights
c
c
c note: numbers are generated by Mathematica file "Gausspoints_2d.nb"
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'iofile.h'

      integer   l,lint
      real*8    sg(3,*)

c  1 point, precision 1.
      if(l.eq.1) then

c     Set number of total points
        lint = 1

c load integration points and weights
        sg(1,1) = 0.0d0
        sg(2,1) = 0.0d0
        sg(3,1) = 4.0d0

c  4 point, precision 3.
      elseif(l.eq.2) then

c     Set number of total points
        lint = 4

c load integration points and weights
       sg(   1   ,   1   )=   -0.5773502691896257
       sg(   2   ,   1   )=   -0.5773502691896257
       sg(   3   ,   1   )=   0.9999999999999998
       sg(   1   ,   2   )=   0.5773502691896257
       sg(   2   ,   2   )=   -0.5773502691896257
       sg(   3   ,   2   )=   0.9999999999999998
       sg(   1   ,   3   )=   -0.5773502691896257
       sg(   2   ,   3   )=   0.5773502691896257
       sg(   3   ,   3   )=   0.9999999999999998
       sg(   1   ,   4   )=   0.5773502691896257
       sg(   2   ,   4   )=   0.5773502691896257
       sg(   3   ,   4   )=   0.9999999999999998

c  9 point, precision 5.
      elseif(l.eq.3) then

c     Set number of total points
        lint = 9

c load integration points and weights
       sg(   1   ,   1   )=   -0.7745966692414834
       sg(   2   ,   1   )=   -0.7745966692414834
       sg(   3   ,   1   )=   0.308641975308642
       sg(   1   ,   2   )=   0.0d0
       sg(   2   ,   2   )=   -0.7745966692414834
       sg(   3   ,   2   )=   0.493827160493827
       sg(   1   ,   3   )=   0.7745966692414834
       sg(   2   ,   3   )=   -0.7745966692414834
       sg(   3   ,   3   )=   0.308641975308642
       sg(   1   ,   4   )=   -0.7745966692414834
       sg(   2   ,   4   )=   0.0d0
       sg(   3   ,   4   )=   0.493827160493827
       sg(   1   ,   5   )=   0.0d0
       sg(   2   ,   5   )=   0.0d0
       sg(   3   ,   5   )=   0.790123456790123
       sg(   1   ,   6   )=   0.7745966692414834
       sg(   2   ,   6   )=   0.0d0
       sg(   3   ,   6   )=   0.493827160493827
       sg(   1   ,   7   )=   -0.7745966692414834
       sg(   2   ,   7   )=   0.7745966692414834
       sg(   3   ,   7   )=   0.308641975308642
       sg(   1   ,   8   )=   0.0d0
       sg(   2   ,   8   )=   0.7745966692414834
       sg(   3   ,   8   )=   0.493827160493827
       sg(   1   ,   9   )=   0.7745966692414834
       sg(   2   ,   9   )=   0.7745966692414834
       sg(   3   ,   9   )=   0.308641975308642

c  16 point, precision 7.
      elseif(l.eq.4) then

c     Set number of total points
        lint = 16

c load integration points and weights
       sg(   1   ,   1    )=   -0.8611363115940526
       sg(   2   ,   1    )=   -0.8611363115940526
       sg(   3   ,   1    )=   0.12100299328560216
       sg(   1   ,   2    )=   -0.3399810435848563
       sg(   2   ,   2    )=   -0.8611363115940526
       sg(   3   ,   2    )=   0.22685185185185194
       sg(   1   ,   3    )=   0.3399810435848563
       sg(   2   ,   3    )=   -0.8611363115940526
       sg(   3   ,   3    )=   0.22685185185185194
       sg(   1   ,   4    )=   0.8611363115940526
       sg(   2   ,   4    )=   -0.8611363115940526
       sg(   3   ,   4    )=   0.12100299328560216
       sg(   1   ,   5    )=   -0.8611363115940526
       sg(   2   ,   5    )=   -0.3399810435848563
       sg(   3   ,   5    )=   0.22685185185185194
       sg(   1   ,   6    )=   -0.3399810435848563
       sg(   2   ,   6    )=   -0.3399810435848563
       sg(   3   ,   6    )=   0.42529330301069407
       sg(   1   ,   7    )=   0.3399810435848563
       sg(   2   ,   7    )=   -0.3399810435848563
       sg(   3   ,   7    )=   0.42529330301069407
       sg(   1   ,   8    )=   0.8611363115940526
       sg(   2   ,   8    )=   -0.3399810435848563
       sg(   3   ,   8    )=   0.22685185185185194
       sg(   1   ,   9    )=   -0.8611363115940526
       sg(   2   ,   9    )=   0.3399810435848563
       sg(   3   ,   9    )=   0.22685185185185194
       sg(   1   ,   10   )=   -0.3399810435848563
       sg(   2   ,   10   )=   0.3399810435848563
       sg(   3   ,   10   )=   0.42529330301069407
       sg(   1   ,   11   )=   0.3399810435848563
       sg(   2   ,   11   )=   0.3399810435848563
       sg(   3   ,   11   )=   0.42529330301069407
       sg(   1   ,   12   )=   0.8611363115940526
       sg(   2   ,   12   )=   0.3399810435848563
       sg(   3   ,   12   )=   0.22685185185185194
       sg(   1   ,   13   )=   -0.8611363115940526
       sg(   2   ,   13   )=   0.8611363115940526
       sg(   3   ,   13   )=   0.12100299328560216
       sg(   1   ,   14   )=   -0.3399810435848563
       sg(   2   ,   14   )=   0.8611363115940526
       sg(   3   ,   14   )=   0.22685185185185194
       sg(   1   ,   15   )=   0.3399810435848563
       sg(   2   ,   15   )=   0.8611363115940526
       sg(   3   ,   15   )=   0.22685185185185194
       sg(   1   ,   16   )=   0.8611363115940526
       sg(   2   ,   16   )=   0.8611363115940526
       sg(   3   ,   16   )=   0.12100299328560216

c  25 point, precision 9.
      elseif(l.eq.5) then

c     Set number of total points
        lint = 25

c load integration points and weights
       sg(   1   ,   1    )=   -0.906179845938664
       sg(   2   ,   1    )=   -0.906179845938664
       sg(   3   ,   1    )=   0.05613434886242878
       sg(   1   ,   2    )=   -0.5384693101056831
       sg(   2   ,   2    )=   -0.906179845938664
       sg(   3   ,   2    )=   0.1134000000000001
       sg(   1   ,   3    )=   0.0d0
       sg(   2   ,   3    )=   -0.906179845938664
       sg(   3   ,   3    )=   0.13478507238752108
       sg(   1   ,   4    )=   0.5384693101056831
       sg(   2   ,   4    )=   -0.906179845938664
       sg(   3   ,   4    )=   0.1134000000000001
       sg(   1   ,   5    )=   0.906179845938664
       sg(   2   ,   5    )=   -0.906179845938664
       sg(   3   ,   5    )=   0.05613434886242878
       sg(   1   ,   6    )=   -0.906179845938664
       sg(   2   ,   6    )=   -0.5384693101056831
       sg(   3   ,   6    )=   0.1134000000000001
       sg(   1   ,   7    )=   -0.5384693101056831
       sg(   2   ,   7    )=   -0.5384693101056831
       sg(   3   ,   7    )=   0.2290854042239909
       sg(   1   ,   8    )=   0.0d0
       sg(   2   ,   8    )=   -0.5384693101056831
       sg(   3   ,   8    )=   0.2722865325507506
       sg(   1   ,   9    )=   0.5384693101056831
       sg(   2   ,   9    )=   -0.5384693101056831
       sg(   3   ,   9    )=   0.2290854042239909
       sg(   1   ,   10   )=   0.906179845938664
       sg(   2   ,   10   )=   -0.5384693101056831
       sg(   3   ,   10   )=   0.1134000000000001
       sg(   1   ,   11   )=   -0.906179845938664
       sg(   2   ,   11   )=   0.0d0
       sg(   3   ,   11   )=   0.13478507238752108
       sg(   1   ,   12   )=   -0.5384693101056831
       sg(   2   ,   12   )=   0.0d0
       sg(   3   ,   12   )=   0.2722865325507506
       sg(   1   ,   13   )=   0.0d0
       sg(   2   ,   13   )=   0.0d0
       sg(   3   ,   13   )=   0.32363456790123457
       sg(   1   ,   14   )=   0.5384693101056831
       sg(   2   ,   14   )=   0.0d0
       sg(   3   ,   14   )=   0.2722865325507506
       sg(   1   ,   15   )=   0.906179845938664
       sg(   2   ,   15   )=   0.0d0
       sg(   3   ,   15   )=   0.13478507238752108
       sg(   1   ,   16   )=   -0.906179845938664
       sg(   2   ,   16   )=   0.5384693101056831
       sg(   3   ,   16   )=   0.1134000000000001
       sg(   1   ,   17   )=   -0.5384693101056831
       sg(   2   ,   17   )=   0.5384693101056831
       sg(   3   ,   17   )=   0.2290854042239909
       sg(   1   ,   18   )=   0.0d0
       sg(   2   ,   18   )=   0.5384693101056831
       sg(   3   ,   18   )=   0.2722865325507506
       sg(   1   ,   19   )=   0.5384693101056831
       sg(   2   ,   19   )=   0.5384693101056831
       sg(   3   ,   19   )=   0.2290854042239909
       sg(   1   ,   20   )=   0.906179845938664
       sg(   2   ,   20   )=   0.5384693101056831
       sg(   3   ,   20   )=   0.1134000000000001
       sg(   1   ,   21   )=   -0.906179845938664
       sg(   2   ,   21   )=   0.906179845938664
       sg(   3   ,   21   )=   0.05613434886242878
       sg(   1   ,   22   )=   -0.5384693101056831
       sg(   2   ,   22   )=   0.906179845938664
       sg(   3   ,   22   )=   0.1134000000000001
       sg(   1   ,   23   )=   0.0d0
       sg(   2   ,   23   )=   0.906179845938664
       sg(   3   ,   23   )=   0.13478507238752108
       sg(   1   ,   24   )=   0.5384693101056831
       sg(   2   ,   24   )=   0.906179845938664
       sg(   3   ,   24   )=   0.1134000000000001
       sg(   1   ,   25   )=   0.906179845938664
       sg(   2   ,   25   )=   0.906179845938664
       sg(   3   ,   25   )=   0.05613434886242878


c  36 point, precision 11.
      elseif(l.eq.6) then

c     Set number of total points
        lint = 36

c load integration points and weights
       sg(   1   ,   1    )=   -0.932469514203152
       sg(   2   ,   1    )=   -0.932469514203152
       sg(   3   ,   1    )=   0.029352081688980326
       sg(   1   ,   2    )=   -0.6612093864662645
       sg(   2   ,   2    )=   -0.932469514203152
       sg(   3   ,   2    )=   0.06180729337238326
       sg(   1   ,   3    )=   -0.2386191860831969
       sg(   2   ,   3    )=   -0.932469514203152
       sg(   3   ,   3    )=   0.08016511731780655
       sg(   1   ,   4    )=   0.2386191860831969
       sg(   2   ,   4    )=   -0.932469514203152
       sg(   3   ,   4    )=   0.08016511731780655
       sg(   1   ,   5    )=   0.6612093864662645
       sg(   2   ,   5    )=   -0.932469514203152
       sg(   3   ,   5    )=   0.06180729337238326
       sg(   1   ,   6    )=   0.932469514203152
       sg(   2   ,   6    )=   -0.932469514203152
       sg(   3   ,   6    )=   0.029352081688980326
       sg(   1   ,   7    )=   -0.932469514203152
       sg(   2   ,   7    )=   -0.6612093864662645
       sg(   3   ,   7    )=   0.06180729337238326
       sg(   1   ,   8    )=   -0.6612093864662645
       sg(   2   ,   8    )=   -0.6612093864662645
       sg(   3   ,   8    )=   0.1301489125881675
       sg(   1   ,   9    )=   -0.2386191860831969
       sg(   2   ,   9    )=   -0.6612093864662645
       sg(   3   ,   9    )=   0.16880536708758792
       sg(   1   ,   10   )=   0.2386191860831969
       sg(   2   ,   10   )=   -0.6612093864662645
       sg(   3   ,   10   )=   0.16880536708758792
       sg(   1   ,   11   )=   0.6612093864662645
       sg(   2   ,   11   )=   -0.6612093864662645
       sg(   3   ,   11   )=   0.1301489125881675
       sg(   1   ,   12   )=   0.932469514203152
       sg(   2   ,   12   )=   -0.6612093864662645
       sg(   3   ,   12   )=   0.06180729337238326
       sg(   1   ,   13   )=   -0.932469514203152
       sg(   2   ,   13   )=   -0.2386191860831969
       sg(   3   ,   13   )=   0.08016511731780655
       sg(   1   ,   14   )=   -0.6612093864662645
       sg(   2   ,   14   )=   -0.2386191860831969
       sg(   3   ,   14   )=   0.16880536708758792
       sg(   1   ,   15   )=   -0.2386191860831969
       sg(   2   ,   15   )=   -0.2386191860831969
       sg(   3   ,   15   )=   0.2189434501672968
       sg(   1   ,   16   )=   0.2386191860831969
       sg(   2   ,   16   )=   -0.2386191860831969
       sg(   3   ,   16   )=   0.2189434501672968
       sg(   1   ,   17   )=   0.6612093864662645
       sg(   2   ,   17   )=   -0.2386191860831969
       sg(   3   ,   17   )=   0.16880536708758792
       sg(   1   ,   18   )=   0.932469514203152
       sg(   2   ,   18   )=   -0.2386191860831969
       sg(   3   ,   18   )=   0.08016511731780655
       sg(   1   ,   19   )=   -0.932469514203152
       sg(   2   ,   19   )=   0.2386191860831969
       sg(   3   ,   19   )=   0.08016511731780655
       sg(   1   ,   20   )=   -0.6612093864662645
       sg(   2   ,   20   )=   0.2386191860831969
       sg(   3   ,   20   )=   0.16880536708758792
       sg(   1   ,   21   )=   -0.2386191860831969
       sg(   2   ,   21   )=   0.2386191860831969
       sg(   3   ,   21   )=   0.2189434501672968
       sg(   1   ,   22   )=   0.2386191860831969
       sg(   2   ,   22   )=   0.2386191860831969
       sg(   3   ,   22   )=   0.2189434501672968
       sg(   1   ,   23   )=   0.6612093864662645
       sg(   2   ,   23   )=   0.2386191860831969
       sg(   3   ,   23   )=   0.16880536708758792
       sg(   1   ,   24   )=   0.932469514203152
       sg(   2   ,   24   )=   0.2386191860831969
       sg(   3   ,   24   )=   0.08016511731780655
       sg(   1   ,   25   )=   -0.932469514203152
       sg(   2   ,   25   )=   0.6612093864662645
       sg(   3   ,   25   )=   0.06180729337238326
       sg(   1   ,   26   )=   -0.6612093864662645
       sg(   2   ,   26   )=   0.6612093864662645
       sg(   3   ,   26   )=   0.1301489125881675
       sg(   1   ,   27   )=   -0.2386191860831969
       sg(   2   ,   27   )=   0.6612093864662645
       sg(   3   ,   27   )=   0.16880536708758792
       sg(   1   ,   28   )=   0.2386191860831969
       sg(   2   ,   28   )=   0.6612093864662645
       sg(   3   ,   28   )=   0.16880536708758792
       sg(   1   ,   29   )=   0.6612093864662645
       sg(   2   ,   29   )=   0.6612093864662645
       sg(   3   ,   29   )=   0.1301489125881675
       sg(   1   ,   30   )=   0.932469514203152
       sg(   2   ,   30   )=   0.6612093864662645
       sg(   3   ,   30   )=   0.06180729337238326
       sg(   1   ,   31   )=   -0.932469514203152
       sg(   2   ,   31   )=   0.932469514203152
       sg(   3   ,   31   )=   0.029352081688980326
       sg(   1   ,   32   )=   -0.6612093864662645
       sg(   2   ,   32   )=   0.932469514203152
       sg(   3   ,   32   )=   0.06180729337238326
       sg(   1   ,   33   )=   -0.2386191860831969
       sg(   2   ,   33   )=   0.932469514203152
       sg(   3   ,   33   )=   0.08016511731780655
       sg(   1   ,   34   )=   0.2386191860831969
       sg(   2   ,   34   )=   0.932469514203152
       sg(   3   ,   34   )=   0.08016511731780655
       sg(   1   ,   35   )=   0.6612093864662645
       sg(   2   ,   35   )=   0.932469514203152
       sg(   3   ,   35   )=   0.06180729337238326
       sg(   1   ,   36   )=   0.932469514203152
       sg(   2   ,   36   )=   0.932469514203152
       sg(   3   ,   36   )=   0.029352081688980326


c     Error

      else

        write(ilg,2000) l
        write(iow,2000) l
        if(ior.lt.0) then
          write(*,2000) l
        endif
        call plstop()

      endif

c     Format

2000  format(' *ERROR* INT2D: Illegal quadrature order =',i16)

      end

	subroutine rcurv (icurv, nr, nreac, ndispl, nd, rflag)
c__________________________________________________TCG. 13.03.2002___71
c
c        read data needed for curve plotting from CURVE.dat
c
c  Declare variable types
c
c     nr          Nodes of needed reactions
c     nd          Nodes of needed average displ
c     rflag       Propper reading flag
c     icurv       reading UNIT number
c     nreac       Number of nodes	of needed reactions
c____________________________________________________________________71

      implicit none

	integer nreac, ndispl, i, io, icurv
	integer nr(1000), nd(5)

	logical rflag
	
     


      rflag=.false.
      OPEN(UNIT=icurv,FILE='CURVE.dat',IOSTAT=io)
	if(io.EQ.0) rflag=.true. 


	if (rflag) then
	 read(icurv,*,IOSTAT=io) nreac
	 read(icurv,*,IOSTAT=io)
	 do i= 1,nreac
	 	read(icurv,*,IOSTAT=io) nr(i) 
	 enddo

	 read(icurv,*,IOSTAT=io)

	 read(icurv,*,IOSTAT=io) ndispl 
	 if (ndispl.gt.5) ndispl=5
	 read(icurv,*,IOSTAT=io)
	 do i= 1,ndispl
	 	read(icurv,*,IOSTAT=io) nd(i) 
	 enddo

	endif
	 
	CLOSE(icurv)
	end


	subroutine printdata (icurv, numnp,ndf, Ulen,displ,DRlen,
	1                 reac,nreac,nr,ndispl, nd)

c__________________________________________________TCG. 13.03.2002___71
c
c        write data for curve plotting into CURVE.txt
c
c  Declare variable types
c
c     ndf         ndf
c     nr          Nodes of needed reactions
c     nd          Node of needed displ
c     rflag       Propper reading flag
c     icurv       reading UNIT number
c     nreac       Number of nodes	of needed reactions
c____________________________________________________________________71
	implicit none

	integer nreac, ndispl, i, j, icurv
	integer nd(ndispl), numnp, ndf, nr(nreac)

	integer  DRlen, Ulen

	real*8 reac(DRlen), displ(Ulen)
	real*8 u(ndispl,ndf), force(ndf)
      
      include 'tdata.h'  

c recover needed displacement
      do j= 1,ndispl
        do i= 1,ndf
	    u(j,i)=displ((nd(j)-1)*ndf+i)
	  enddo
	enddo

c sum over reactions
      force=0.0d0
	do j=1,nreac
      do i= 1,ndf
	  force(i)=force(i) - reac((nr(j)-1)*ndf+i)
	enddo
      enddo 

      if (ndispl.eq.1) then
      
c append line to file
	write(icurv,10) ttim, dt, 
     1                u(1,1),u(1,2),u(1,3),
     1                u(1,4),u(1,5),u(1,6),
     1                force(1),force(2),force(3),
     1                force(4),force(5),force(6) 
      
      elseif (ndispl.eq.2) then 
          
c append line to file
	write(icurv,20) ttim, dt, 
     1                u(1,1),u(1,2),u(1,3),
     1                u(1,4),u(1,5),u(1,6),
     1                u(2,1),u(2,2),u(2,3),
     1                u(2,4),u(2,5),u(2,6),      
     1                force(1),force(2),force(3),
     1                force(4),force(5),force(6) 
      
      endif

10    format(f18.8,f18.8,
     1       f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,
     1       f18.8,f18.8,f18.8,f18.8,f18.8,f18.8)    
      
20    format(f18.8,f18.8,
     1       f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,
     1       f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,      
     1       f18.8,f18.8,f18.8,f18.8,f18.8,f18.8)          
      
!c append line to file
!	write(icurv,10) u(1,1),u(2,1),u(3,1),u(4,1),u(5,1),force(1),
!	1                u(1,2),u(2,2),u(3,2),u(4,2),u(5,2),force(2), 
!	1                u(1,3),u(2,3),u(3,3),u(4,3),u(5,3),force(3), 
!	1                u(1,4),u(2,4),u(3,4),u(4,4),u(5,4),force(4), 
!	1                u(1,5),u(2,5),u(3,5),u(4,5),u(5,5),force(5), 
!	1                u(1,6),u(2,6),u(3,6),u(4,6),u(5,6),force(6) 
!
!10	format(1x,f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,'  , ',
!     1          f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,'  , ',
!     1          f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,'  , ',
!     1          f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,'  , ',
!     1          f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,'  , ',
!     1          f18.8,f18.8,f18.8,f18.8,f18.8,f18.8,'  , ')

	end


      subroutine int_str10 (fn,str)
c__________________________________________________TCG. 17.05.2001___71
c
c        convert interger 0<fn<1000 into string
c
c  Declare variable types
c
c     fn          Integer
c     str         string
c____________________________________________________________________71
	integer fn,i,temp, z1, z2, z3
	integer bib(0:9)
	Character str*8

	if (fn.gt.999) str='CRAC_xxx'

c load ANSII
	do i=0,9
	 bib(i)=48+i
	enddo

	temp=fn

c decompose integer
      z1 = temp/100
	temp = temp-(100*z1)
	z2 = temp/10
	temp = temp-(10*z2)
	z3 = temp
c bulid string
	str='CRAC_'//char(bib(z1))//char(bib(z2))//char(bib(z3))
      end



      subroutine int_str_ctet10 (fn,str)
c__________________________________________________TCG. 17.05.2001___71
c
c        convert interger 0<fn<1000 into string
c
c  Declare variable types
c
c     fn          Integer
c     str         string
c____________________________________________________________________71
	integer fn,i,temp, z1, z2, z3
	integer bib(0:9)
	Character str*8

	if (fn.gt.999) str='CTET_xxx'

c load ANSII
	do i=0,9
	 bib(i)=48+i
	enddo

	temp=fn

c decompose integer
      z1 = temp/100
	temp = temp-(100*z1)
	z2 = temp/10
	temp = temp-(10*z2)
	z3 = temp
c bulid string
	str='CTET_'//char(bib(z1))//char(bib(z2))//char(bib(z3))
      end


      subroutine int_str_disp11 (fn,str)
c__________________________________________________TCG. 15.07.2003___71
c
c        convert interger 0<fn<1000 into string
c
c  Declare variable types
c
c     fn          Integer
c     str         string
c____________________________________________________________________71
	integer fn,i,temp, z1, z2, z3
	integer bib(0:9)
	Character str*8

	if (fn.gt.999) str='DISP_xxx'

c load ANSII
	do i=0,9
	 bib(i)=48+i
	enddo

	temp=fn

c decompose integer
      z1 = temp/100
	temp = temp-(100*z1)
	z2 = temp/10
	temp = temp-(10*z2)
	z3 = temp
c bulid string
	str='DISP_'//char(bib(z1))//char(bib(z2))//char(bib(z3))
      end


      subroutine int_str_stre11 (fn,str)
c__________________________________________________TCG. 15.07.2003___71
c
c        convert interger 0<fn<1000 into string
c
c  Declare variable types
c
c     fn          Integer
c     str         string
c____________________________________________________________________71
	integer fn,i,temp, z1, z2, z3
	integer bib(0:9)
	Character str*8

	if (fn.gt.999) str='STRE_xxx'

c load ANSII
	do i=0,9
	 bib(i)=48+i
	enddo

	temp=fn

c decompose integer
      z1 = temp/100
	temp = temp-(100*z1)
	z2 = temp/10
	temp = temp-(10*z2)
	z3 = temp
c bulid string
	str='STRE_'//char(bib(z1))//char(bib(z2))//char(bib(z3))
      end


      subroutine int_str_hist11 (fn,str)
c__________________________________________________TCG. 15.07.2003___71
c
c        convert interger 0<fn<1000 into string
c
c  Declare variable types
c
c     fn          Integer
c     str         string
c____________________________________________________________________71
	integer fn,i,temp, z1, z2, z3
	integer bib(0:9)
	Character str*8

	if (fn.gt.999) str='HIST_xxx'

c load ANSII
	do i=0,9
	 bib(i)=48+i
	enddo

	temp=fn

c decompose integer
      z1 = temp/100
	temp = temp-(100*z1)
	z2 = temp/10
	temp = temp-(10*z2)
	z3 = temp
c bulid string
	str='HIST_'//char(bib(z1))//char(bib(z2))//char(bib(z3))
      end


                             

      subroutine write_crack(ivisu,ndm,numel,nen,nen1,x,xlen,ix,
     1                       ixlen,hist,hlen,n1,n3)
c__________________________________________________TCG. 11.09.2002___71
c
c      write x+u into VISU_xxx
c
c  Declare variable types
c
c     ivisu       File ID
c     ndm         Number of dimension
c     numel       Number of elements
c     numnp       Number of nodal points
c     x           Ref. coordinate of nodes
c     xlen        Length of x-vector
c     ix          Element-connectivity
c     ixlen       Length of ix-vector
c     hist        history vector
c     hlen        Length of history vector
c     h1len       Length of nh1 (and nh2) history vector per elment
c     h3len       Length of nh3 history vector per elment
c____________________________________________________________________71
      implicit none

	integer xlen, ivisu, hlen
	integer i, k, j,nen, ndm, ixlen, numel,n1,n3
      integer ix(ixlen), nodenr,  nn, nen1
	real*8 x(xlen), xbar(3), hist(hlen)
	real*8 N0(3), Ae, uh, Ve

c pointer of nh3 array
      nn = 2*n1 + n3
      
c   compute center of the cracked elments
      do i= 1, numel
       if (hist(1+(i-1)*nn+2*n1+48).gt.0.5d0) then 
        xbar=0.0d0
	  do j= 1,nen
	   nodenr = ix((i-1)*nen1+j)
	   do k= 1,ndm
	     xbar(k) = xbar(k) + 0.25d0*x((nodenr-1)*ndm + k)
	   enddo
	  enddo
c   retrieve normal
       do j= 1,3 
        N0(j) = Hist(1+(i-1)*nn+2*n1+43+j)
	 enddo
c   retrieve Ae
        Ae = Hist(2*n1+1+(i-1)*nn+42)
c   retrieve Ve
        Ve = Hist(2*n1+1+(i-1)*nn+43)
c   retrieve uh
        uh = Hist((i-1)*nn+1)
c   write line into CRACK.txt
	  write(ivisu,1000)  (xbar(k),k=1,3), (N0(k),k=1,3), Ae, uh,Ve
       endif
	enddo


1000  format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5)

	end



      subroutine write_ctet10(ivisu,ndm,numel,nen,nen1,x,xlen,ix,
     1                       ixlen,hist,hlen,n1,n3)
c__________________________________________________TCG. 11.09.2002___71
c
c      write x+u into VISU_xxx
c
c  Declare variable types
c
c     ivisu       File ID
c     ndm         Number of dimension
c     numel       Number of elements
c     numnp       Number of nodal points
c     x           Ref. coordinate of nodes
c     xlen        Length of x-vector
c     ix          Element-connectivity
c     ixlen       Length of ix-vector
c     hist        history vector
c     hlen        Length of history vector
c     h1len       Length of nh1 (and nh2) history vector per elment
c     h3len       Length of nh3 history vector per elment
c____________________________________________________________________71
      implicit none
      
	logical node_flag(4)
	integer xlen, ivisu, hlen
	integer i, k, j,nen, ndm, ixlen, numel,n1,n3
      integer ix(ixlen), nodenr,  nn, nen1
	real*8 x(xlen), xn(3,4), hist(hlen)
	real*8 uh, n_flag

c pointer of nh3 array
      nn = 2*n1 + n3
      
c   write element coordinates 
      do i= 1, numel
       if (hist(1+(i-1)*nn+2*n1+6).eq.1) then 
	  do j= 1,nen
	   nodenr = ix((i-1)*nen1+j)
	   do k= 1,ndm
	     xn(k,j) = x((nodenr-1)*ndm + k)
	   enddo
	  enddo
c   retrieve uh
        uh = Hist((i-1)*nn+1)
c   retrieve node_flag
        n_flag = hist(1+(i-1)*nn+2*n1+47)
	  
       call flag_convert10(int(n_flag),node_flag)


c   write line into CTET.txt
	  write(ivisu,1000)  
	1      i, ((xn(k,j),k=1,3),j=1,4), (node_flag(j),j=1,4)
       endif
	enddo


1000  format(i5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5
     1      ,e13.5,e13.5,e13.5,e13.5,L2,L2,L2,L2)

	end



      subroutine write_inerface_cracks11(ivisu,ndm,numel,nen,nen1,
     1                   x,xlen,ix,ixlen,hist,hlen,n1,n3)
c__________________________________________________TCG. 13.04.2006___71
c
c      write cracked interface elements into opened file
c
c  Declare variable types
c
c     ivisu       File ID
c     ndm         Number of dimension
c     numel       Number of elements
c     numnp       Number of nodal points
c     x           Ref. coordinate of nodes
c     xlen        Length of x-vector
c     ix          Element-connectivity
c     ixlen       Length of ix-vector
c     hist        history vector
c     hlen        Length of history vector
c     h1len       Length of nh1 (and nh2) history vector per elment
c     h3len       Length of nh3 history vector per elment
c____________________________________________________________________71
      implicit none

	include 'iface.h'
      
	logical node_flag(4)
	integer xlen, ivisu, hlen, nninter, nvinter
	integer i, k, j,nen, ndm, ixlen, numel,n1,n3
      integer ix(ixlen), nodenr,  nn, nen1
	real*8 x(xlen), xn(3,3), hist(hlen)
	real*8 uh, n_flag

c pointer of nh3 array
      nn = 2*n1 + n3

	nninter = (fielem-1)*nn

c	nvinter = 2*4 + 9
	nvinter = 2*1 + 9
      
c   write element coordinates 
      do i= fielem, lielem
c   retrieve uh
      uh = Hist(nninter+(i-fielem)*nvinter + 1)
       if (uh.gt.1.0d0/7.77207d0) then 
	  do j= 1,3
	   nodenr = ix((i-1)*nen1+j)
	   do k= 1,ndm
	     xn(k,j) = x((nodenr-1)*ndm + k)
	   enddo
	  enddo


c   write line into crack_xxx
	  write(ivisu,1000) 3, ((xn(k,j),k=1,3),j=1,3), 0.0,0.0,0.0
       endif
	enddo


1000  format(i6,  e13.5,e13.5,e13.5, 
     1            e13.5,e13.5,e13.5,
     2            e13.5,e13.5,e13.5,
     2            e13.5,e13.5,e13.5)

	end




      subroutine plot_ctip11
	integer i,j
	include 'tdata.h'
      include  'cTip.h'
 
      OPEN(UNIT=44,FILE='ctip.tex')
      do i= 1,nctip
	  write(44,1000) (ctip(i,j),j=1,6), (cTip_point(i,j),j=1,3)
	enddo

	CLOSE(44)

1000  format(i6,i6,i6,i6,i6,i6,e14.5,e14.5,e14.5)
      end
 
 
      subroutine plot_ctets11
	integer i,j
	include 'tdata.h'
      include  'cTets.h'
 
      OPEN(UNIT=45,FILE='ctets.tex')
      do i= 1,ntets
	  write(45,1000) (ctets(i,j),j=1,9)
	enddo

	CLOSE(45)

1000  format(i6,i6,i6,i6,i6,i6,i6,i6,i6)
      end



	subroutine solution_11 (ul, al, dul, ndf, ndm, nel, 
	1                        ulc, ule, alc,ale, dule)
c__________________________________________________TCG. 11.09.2002___71
c
c      get soluition data f solution vector
c
c  Declare variable types
c__________________________________________________TCG. 11.09.2002___71

      implicit none

	integer ndm, nel, ndf, i, j

	real*8 ul(ndf,nel), al(ndf,nel), dul(ndf,nel)
	real*8 ulc(ndm,nel), ule(ndm,nel), dule(ndm,nel)
	real*8 alc(ndm,nel), ale(ndm,nel)


c get compatible and enhanced displacements
      do i= 1,nel
	 do j= 1,ndm
	   ulc(j,i) = ul(j,i)       !compatible displacements
	   ule(j,i) = ul(j+3,i)     !enhanced displacements
	 enddo
	enddo

c get compatible and enhanced accelerations
      do i= 1,nel
	 do j= 1,ndm
	   alc(j,i) = al(j,i)       !compatible accelerations
	   ale(j,i) = al(j+3,i)     !enhanced accelerations
	 enddo
	enddo



c get increment of enhanced displacements
      do i= 1,nel
	 do j= 1,ndm
	   dule(j,i) = dul(j+3,i)     
	 enddo
	enddo


      end



      subroutine update_11(numel,n1,n3,nen1,hist,hlen,ix,ixlen)
c__________________________________________________TCG. 05.06.2003___71
c
c
c  Declare variable types
c
c     numel       Number of elements
c     hist        history vector
c     hlen        Length of history vector
c     h1len       Length of nh1 (and nh2) history vector per elment
c     h3len       Length of nh3 history vector per elment
c____________________________________________________________________71
      implicit none

	include 'cTip.h'

	logical flag
      
	integer hlen, nn, ixlen, nen1, nmat
	integer i, j, numel,n1,n3
	integer ix(ixlen)
	real*8 hist(hlen)

c pointer of nh3 array
      nn = 2*n1 + n3


      do i= 1, numel
c       do i= 1, 10689

        nmat = int(ix((i-1)*nen1 + nen1))
c limit to material where cracks are allowed
        flag=.false.
        do j= 1, ncmat
	    if(abs(cmat(j)).eq.nmat) flag=.true.
	  enddo
        if (flag) then
	    do j= 1,n1
            Hist((i-1)*nn+j) = Hist((i-1)*nn+n1+j)
	    enddo
	  endif
	enddo
      
	end



      subroutine WHIST_11(ivisu,numel,hist,hlen,n1,n3)
c__________________________________________________TCG. 17.11.2003___71
c
c      write HIST_xxx
c
c  Declare variable types
c
c     numel       Number of elements
c     hist        history vector
c     hlen        Length of history vector
c     h1len       Length of nh1 (and nh2) history vector per elment
c     h3len       Length of nh3 history vector per elment
c____________________________________________________________________71
      implicit none
      
	integer hlen, nn, ivisu
	integer i, j, numel,n1,n3
	real*8 hist(hlen)

      include  'nodes_id.h'


c pointer of nh3 array
      nn = 2*n1 + n3


      do i= 1, numel
c	 if(Hist((i-1)*nn+2*n1+7).eq.1) then
	   write(ivisu,1000) i
	   write(ivisu,1001) 
	   do j= 1,n1            
	    write(ivisu,1010) j-1, Hist((i-1)*nn+j)
	   enddo
	   write(ivisu,1002) 
	   do j= 1,n1            
	    write(ivisu,1010) j-1, Hist((i-1)*nn+n1+j)
	   enddo
	   write(ivisu,1003) 
	   do j= 1,n3            
	    write(ivisu,1010) j-1, Hist((i-1)*nn+2*n1+j)
	   enddo
c	 endif
	enddo
      
1000  format( 5x, 'Element   ',i5/)
1001  format( 5x, 'History nh1')
1002  format( 5x, 'History nh2')
1003  format( 5x, 'History nh3')
1010  format( 5x, '     ', i3, '     ', e20.12)

	end






      subroutine write_nodes11(ivisu,ndm,numnp,x,xlen)
c__________________________________________________TCG. 15.07.2003___71
c
c      write x into Coord
c
c  Declare variable types
c
c     ivisu       File ID
c     ndm         Number of dimension
c     numnp       Number of nodal points
c     x           Ref. coordinate of nodes
c     xlen        Length of x-vector
c____________________________________________________________________71
      implicit none

	integer xlen, ivisu, numnp
	integer i, j, k, ndm
	real*8 x(xlen)

      k=0
 
      do i= 1,numnp
	  write(ivisu,*)  (x(k+j),j=1,ndm)
	  k=k+ndm 
	enddo


	end

      subroutine write_connect11(ivisu,nen1,nen,numel,ix,ixlen)
c__________________________________________________TCG. 15.07.2003___71
c
c        write element-connectivity into VISU_xxx
c
c  Declare variable types
c
c     ivisu       File ID
c     nen1        Length of line in ix-vector
c     nen         Number of nods per element
c     ix          Element-connectivity
c     numel       Number of elements
c     ixlen       Length of ix-vector
c____________________________________________________________________71
	integer ivisu, nen1, numel, ixlen, nen
      integer ix(ixlen)
	integer i, k, j



	k=1
      do i= 1,numel
	 write(ivisu,1000)  (ix(k+j-1),j=1,nen)
	 k=k+nen1
	enddo

1000  format(i7,i7,i7,i7,i7,i7,i7,i7,i7,i7)
	end

      subroutine write_displacements11(ivisu,ndm,ndf,numnp,u,ulen)
c__________________________________________________TCG. 15.07.2003___71
c
c      write x into Coord
c
c  Declare variable types
c
c     ivisu       File ID
c     ndf         Number of degrees of freedom
c     numnp       Number of nodal points
c     u           displacementsof nodes
c     ulen        Length of u-vector
c____________________________________________________________________71
      implicit none

      include  'nodes_id.h'

	integer ulen, ivisu, numnp
	integer i, j, k, ndf,ndm
	real*8 u(ulen), ue(ndm), uc(ndm)

      k=0
 
      do i= 1,numnp
	  ue = 0.0d0
	  uc = 0.0d0
	  do j= 1,ndm
	    uc(j) = u(k + j)
	    if (nid(i).eq.1) then
	      ue(j) = u(k + j + ndm)
	    endif
	  enddo
	  write(ivisu,1000)  (uc(j) + ue(j),j=1,ndm)
	  k=k+ndf 
	enddo


1000  format(e13.5,e13.5,e13.5)
	end


      subroutine write_stress11(ivisu,st,numnp)
c__________________________________________________TCG. 17.05.2001___71
c
c        write Stresses into VISU_xxx
c
c  Declare variable types
c
c     ivisu       File ID
c     st          Nodal stresses
c     numnp       Number of nodal points
c____________________________________________________________________71

	integer numnp,i,j, ivisu
	real*8 st(numnp,*)


      do i= 1, numnp
	  write(ivisu,1000)  (st(i,j),j=1,6)
	enddo

1000  format(e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,
     1       e13.5,e13.5,e13.5)
	end


      subroutine id_mideside_11(n,ndm,nen,nen1,x,xlen,
	1                  ix,ixlen, N0,fp)
c__________________________________________________TCG. 17.07.2003___71
c
c        Define id of midside nodes for quadratic elements
c  Declare variable types
c
c     n          global node no.
c     xl         coordinates of node
c     fp         front point
c     N0         Normal vector
c____________________________________________________________________71

      implicit none

      include  'nodes_id.h'

	integer xlen, ixlen, nen, n, nen1, ndm
	integer i, k, j, ix(ixlen), nodenr(10)
	real*8 x(xlen), xl(3,10)
	real*8 N0(3), fp(3), temp



c retrieve node no. and lagrangian coordinates
      do j = 1,nen
	  nodenr(j) = ix((n-1)*nen1 + j)
	  do k= 1,ndm
	    xl(k,j) = x((nodenr(j)-1)*ndm + k) 
	  enddo
      enddo


      do i= 5,10
	if(nid(nodenr(i)).eq.0) then ! define nid 
	  temp = N0(1)*(xl(1,i)-fp(1)) + N0(2)*(xl(2,i)-fp(2)) +
	1         N0(3)*(xl(3,i)-fp(3))
	  if (temp.gt.0.0d0) then
	    nid(nodenr(i)) = -1
	  else
	    nid(nodenr(i)) = +1
	  endif
	endif
	enddo


	end


      subroutine get_nodes_elem(n,nen,nen1,ix,ixlen, nel)
c__________________________________________________TCG. 20.07.2006___71
c
c        Define number of nodes of elem no. n
c  Declare variable types
c
c     nel        number of nodes
c____________________________________________________________________71
      implicit none

	integer ixlen, nen, nel, n, nodenr, j, nen1
	integer ix(ixlen)


c retrieve node no. 
      do j = 1,nen
	  nodenr = ix((n-1)*nen1 + j)
        if(nodenr.gt.0) then
          nel = j
        endif
      enddo

      end

      subroutine id_mod_mideside_11(numel,ndm,nen,nen1, ix,ixlen)
c__________________________________________________TCG. 19.07.2003___71
c
c        Modify id of midside nodes for quadratic elements
c        in view of crackclosure at the tip
c  Declare variable types
c
c     n          global node no.
c     xl         coordinates of node
c     fp         front point
c     N0         Normal vector
c____________________________________________________________________71

      implicit none

      include  'nodes_id.h'
c      include  'iofile.h'

	integer ixlen, nen, numel, nen1, ndm
	integer i, j, ix(ixlen), nodenr(10)


c retrieve node no. and lagrangian coordinates
      do i= 1, numel
        do j = 1,nen
	    nodenr(j) = ix((i-1)*nen1 + j)
        enddo
c	  write(iow,1000) i, (nodenr(j),j=1,10)
c	  write(*,1000) i, (nodenr(j),j=1,10)
c modify nid if neccessary
	  if((nid(nodenr(1)).eq.0).and.((nid(nodenr(2)).eq.0))) then
	    nid(nodenr(5)) = 0
	  endif
	  if((nid(nodenr(2)).eq.0).and.((nid(nodenr(3)).eq.0))) then
	    nid(nodenr(6)) = 0
	  endif
	  if((nid(nodenr(1)).eq.0).and.((nid(nodenr(3)).eq.0))) then
	    nid(nodenr(7)) = 0
	  endif
	  if((nid(nodenr(1)).eq.0).and.((nid(nodenr(4)).eq.0))) then
	    nid(nodenr(8)) = 0
	  endif
	  if((nid(nodenr(2)).eq.0).and.((nid(nodenr(4)).eq.0))) then
	    nid(nodenr(9)) = 0
	  endif
	  if((nid(nodenr(3)).eq.0).and.((nid(nodenr(4)).eq.0))) then
	    nid(nodenr(10)) = 0
	  endif
	enddo

1000  format(5x, i6, i6, i6, i6, i6, i6, i6, i6, i6, i6, i6/)
	end



      subroutine close_front_11(ix,ixlen,hist,hlen,n1,n3,pcini,nen1)
c__________________________________________________TCG. 19.07.2003___71
c
c  constrain enhanced degree of freedom at the cracktip        
c
c
c  Declare variable types
c
c     ix         ix vector
c     ixlen      Length of ix vector
c     hist       history vector
c     hlen       length of history vector
c     n1         history pointer
c     n3         history pointer
c     pcini      pointer to crack initialization flag
c____________________________________________________________________71

      implicit none

      include  'nodes_id.h'
      include  'eldata.h'
      include  'cdata.h'
      include  'iofile.h'
      include  'cTip.h'

	logical flag
 
      integer hlen, ixlen, temp, nn, i, j, nen1, pcini
	integer n1, n3, ng(10), ix(ixlen)
	integer min, plus, count, nmat

	real*8 hist(hlen)

c pointer of n3 array
      nn = 2*n1 + n3


      count =0

      do i=1,numel
c get material no.
        nmat = ix((i-1)*nen1 + nen1)
c limit to materials where cracks are allowed
        flag=.false.
        do j= 1, ncmat
	    if(abs(cmat(j)).eq.nmat) flag=.true.
	  enddo
        if (flag) then 
          if(int(Hist(1+(i-1)*nn+2*n1+pcini)).ne.1) then ! displacement element
	      count = count +1
c retrieve element global nodes
            do j= 1,nel
	        ng(j) =  ix((i-1)*nen1 + j)
	      enddo
c	      write(iow,1001) i, (ng(j),j=1,nel), (nid(ng(j)),j=1,nel)
c check node id
            min=0
	      plus=0
	      do j= 1,nel
	        if(int(nid(ng(j))).eq.-1) min  = min  - 1
	        if(int(nid(ng(j))).eq.+1) plus = plus + 1
	      enddo
	      temp = min*plus
		  if(temp.lt.0) then !nodes are on the crack tip
	        write(*,1000) i
	        write(iow,1000) i
c redefine node_id
	        do j= 1,nel
	          nid(ng(j)) = 0
	        enddo      
		  endif	    
	    endif
	  endif
      enddo

	write(iow,*) 'checked', count-1, 'elements'

1000  format(5x, 'Inactivate enh. dof of elem.',i6)
1001  format(1x, i6, i5,i5,i5,i5, i3,i3,i3,i3)
      end



      subroutine asseble_ctets_11(d,dlen,n1,n3,nen1,pCini,
	1          pN0,pVe,pXbar,psig,pfi,pNelem,pT0,pNold,
     2          hist,hlen,ix,ixlen)

c__________________________________________________TCG. 26.11.2003___71
c
c  check crack initialization criterion using a non-local approach 
c  and store latent tets in ctets and n0 in history vector     
c
c
c  Declare variable types
c
c     nold       averaged normal vector of existing pu elements (not normalized)
c     sig        averaged Cauchy stress
c     sige       element Cauchy stress
c     fiv        averaged deformation gradient
c     five       element deformation gradient
c     sig0       cohesive strength
c     imat       bulk material id
c     pN0        Pointer to normal vector
c     pNelem     Pointer to neighbouring elements
c     psig       Pointer to stresses
c     pfi        Pointer to deformation
c     dlen       length of material property vector
c     d          material property vector
c     hlen       length of history vector
c     hist       history 
c     n1         history pointer
c     n3         history pointer
c     ix         ix vector
c     ixlen      Length of ix vector
c     ng         Global node no. of element corner nodes
c     nmat       no. of material
c     sb_r       strenght in radial direction
c     sb_t       strenght in circumferential direction
c     sb_z       strenght in axial direction
c     z          normal stress in radial direction
c     x          normal stress in circumferential direction
c     y          normal stress in axial direction
c____________________________________________________________________71

      implicit none

 
      include  'cdata.h'
      include  'cTets.h'
      include  'cTip.h'
      include  'eldata.h'
      include  'iofile.h'
      include  'vecdata.h'
      include  'nLength.h'
      include  'cGeom.h'
      include  'pvecdata.h'
      include  'Visco.h'

	logical cflag(2), flag

      integer ixlen, hlen, dlen, j, pSig, pFi, pNelem
	integer ii, jj, kk, ielem, imat, pN0, ng(4), pVe, pCini, itemp
	integer ix(ixlen), nen1, n1, n3, nn, k, cini, tn1
	integer nelem, nmat, pT0, pXbar, pNold
	integer icmat, ll, ccou, imatd

	real*8 sig0, hist(hlen), d(dlen), sig(6),fiv(9),sige(6),five(9)
	real*8 stress(3,3), n0(3), fi(3,3), lamb(3)
	real*8 Vei, sig_max, wVe, temp, xc(3), dist, tempv(3), R2
	real*8 sb_r, sb_t, sb_z
	real*8 t0, R, sig_root, sumN0, n0htemp(3)
      real*8 detH, detHe, n0h(3)

	real*8 xb_surf(3)
	real*8 z_coor, check(3)
      
      integer q, idv
      real*8 nc(3), normn, tr(3), sign, trs(3), sigs
      real*8 n0_case(3,4), sign_case(4), sigs_case(4)
      
      ! parameters
	sig0 = d(20)
	imat = int(d(70))
      imatd = int(d(72))

      ! pointer of n3 array
      tn1=2*n1
      nn = tn1 + n3

      ! clear ctets
      ctets=0
	ntets=0

      !ccou = 1
      !do ii= 1,numel 
      !  nmat = ix((ii-1)*nen1 + nen1)  
      !  !check = Hist(1+(ii-1)*nn+tn1+368)  
      !  if (nmat.eq.1) then
      !     check = Hist(1+(ccou-1)*nn+tn1+368)   
      !     ccou = ccou+1 
      !  endif    
      !enddo
      
      ccou=0
      do ii= 1,numel 
        ! get material no.
        nmat = ix((ii-1)*nen1 + nen1)
        ! limit to material where cracks are allowed
        flag=.false.
        do j= 1, ncmat
	    if(abs(cmat(j)).eq.nmat) then
		  flag=.true.
	      icmat=j
            !ccou = ccou+1 ! CJM (13/10/2022)=> Mate 1 is tissue, Mate 2 is clamp [in built lin elastic and so NO HISTORY], in ctet.dat cracking only allowed in Mate,1
            ! NEED THIS FOR CORRECT HISTORY ALLOCATION
	    endif
        enddo

c================ COMPUTE AVERAGE STRESS AND DEFORMATION GRADIENT ========
        
        if (flag) then
          ! get material properties
          imat = int(d((nmat-1)*401 + 70))
	    sig0 = d((nmat-1)*401 + 20) ! THE POINT (401) SEEMS TO BE WRONG IF MATERIAL ISNT NUMBER 1.
          sb_r = d((nmat-1)*401 + 25)
          sb_t = d((nmat-1)*401 + 26)
          sb_z = d((nmat-1)*401 + 27)
          ! get cini-flag
          !cini = int(Hist(1+(ii-1)*nn+tn1+pCini))
          !if (imat.ne.9) then ! CJM EDIT ()
             cini = int(Hist(1+(ii-1)*nn+tn1+pCini))       
          !else
          !   cini = int(Hist(1+(ccou-1)*nn+tn1+pCini))  
          !endif
          
	    if (cini.eq.0) then
            !if (imat.ne.9) then  
              ! get no. of neighbouring elements 
              nelem = int(Hist(1+(ii-1)*nn+tn1+pNelem))
              ! get z coordinate ** required for for bone model only ! 
              z_coor = int(Hist(1+(ii-1)*nn+tn1+pXbar+2))
              ! no average if required
	        if (nelem.eq.0) then
	          do k= 1,6
	            sig(k) = Hist(1+(ii-1)*nn+tn1+pSig-1+k)
	          enddo
	          do k= 1,9
	            Fiv(k) = Hist(1+(ii-1)*nn+tn1+pFi-1+k)
	          enddo
	          sige = sig
	          five = fiv
                if (imat.eq.9) then
                  detH = Hist(1+(ii-1)*nn+tn1+362)
                  detHe = detH
                  do k=1,3
                    N0h(k) = Hist(1+(ii-1)*nn+tn1+362+k)      
                  enddo    
                endif    
	        else
                ! averaging element quantities
	          do k= 1,6
	            sige(k) = Hist(1+(ii-1)*nn+tn1+pSig-1+k)
	          enddo
	          do k= 1,9
	            Five(k) = Hist(1+(ii-1)*nn+tn1+pFi-1+k)
	          enddo
                ! averaged quantities
                ! initialize quantities
	          sig  = 0.0d0
	          fiv  = 0.0d0
	          wVe  = 0.0d0
                if (imat.eq.9) then
                  detH = 0d0
                  n0h = 0d0
                endif  
                do jj= 1,nelem
                  ! get element no. of neighbouring element
	            ielem = int(Hist(1+(ii-1)*nn+tn1+pNelem+jj))
                  ! get element volume 
                  Vei = abs(Hist(1+(ielem-1)*nn+tn1+pVe))
                  ! summation
	            wVe = wVe + Vei
	            do k= 1,6
	              sig(k) = sig(k) + 
	1                   Vei*Hist(1+(ielem-1)*nn+2*n1+pSig-1+k)
	            enddo
	            do k= 1,9
	              fiv(k) = fiv(k) + 
	2                   Vei*Hist(1+(ielem-1)*nn+2*n1+pFi-1+k)
                  enddo
                  if (imat.eq.9) then
	              detH = detH + 
	2                   Vei*Hist(1+(ii-1)*nn+tn1+362)
                    do k=1,3
                      n0h(k) = n0h(k) + 
     2                     Vei*Hist(1+(ii-1)*nn+tn1+362+k)
                    enddo    
                  endif 
	          enddo  ! inner element loop
                ! normalization of the quantities
	          sig = sig/wVe
	          fiv = fiv/wVe
                if (imat.eq.9) then
                  detH = detH/wVe
                  n0h = n0h/wVe
                  temp = 0d0
                  do k=1,3 
                    temp = temp + (n0h(k))**2
                  enddo
                  n0h = n0h/sqrt(temp)   
                endif  
              endif
            !endif
          
c================ CHECK CRACK INITITIALIZATION CRITERION ========
          
	      cflag = .false.
            
            if (imatd.lt.8) then ! MAXIMAL PRINCIPAL STRESS CRITERIA ----------------------------- 
              ! compute stress criterion depending on imat
              call stress_matrix_11(sig,stress)
              
              call eigen_11(stress,lamb,itemp)
              call max_lamb_11 (lamb,k)
              sig_max = lamb(k)
              
              ! compute direction vector if crack initialization is satisfied
              t0=d(20)    
              if(sig_max.gt.t0) then    
                ! normal vector in spatial conf. (not normalized!)
                call get_normal_11_2(stress,lamb, n0)
                ! pull-back
                call defo_matrix_11(fiv,fi)
                call pushv10 (fi, n0)
                ! normalize
                temp = sqrt(N0(1)*N0(1) + N0(2)*N0(2) + N0(3)*N0(3))
                if (temp.ne.0.0d0) then
                  N0 = N0/temp 
                else 
                  write (*,*) 'ERROR: No normal available'
                  write (iow,*) 'ERROR: No normal available'
                endif
                ! set local cflag
                cflag(1) = .true.
              endif ! ini criterion 
              ! check root criterion
              if(sig_max.gt.t0) then
                cflag(2) =.true.
              endif
                
            elseif (imatd.eq.8) then ! GAUSS POINT LOCALISATION CRITERIA ----------------------------------------- 
              
              ! N0 is initialised as N0=[0,0,0] when simulation initially begins. Then
              ! when a gauss Point localises in an element, the normal to the discontinuity
              ! is written directly into the history in the kernel11 subroutine. So we now 
              ! check if the sum of all the components of N0 from the history are nonzero, and
              ! if so, we insert a discontinuity  
              
              if (detH.lt.0d0) then ! THERE IS A DISCONTINUITY IN THE ELEMENT
                
                ! Normal determined previously  
                  
                do k= 1,3
                  N0(k) = n0h(k)            
                enddo  
                  
                ! Correct format of the stress and deformation gradient (3x3)
                call stress_matrix_11(sig,stress)
                call defo_matrix_11(fiv,fi)
                
                ! Deformed discontinuity normal i.e. Push forward, (which is then normalised)
                do k= 1,3
                  do j= 1,3
                    nc(k)= nc(k)+fi(k,j)*n0(j)
                    !nc(k)= nc(k)+fi(j,k)*n0(j)
                  enddo
                enddo     
              
                normn = 0d0
                do j=1,3 
                  normn = normn + (nc(j))**2
                enddo
                
                if (normn.ne.0d0) then
                  nc = nc/sqrt(normn)  
                else  
                  write (*,*) 'ERROR: No normal available'
                  write (iow,*) 'ERROR: No normal available'
                endif                
            
                ! Cauchy traction vector and component of traction vector acting normal to discontinuity (sign)
                tr=0d0
                sign=0d0
                do k= 1,3
                  do j= 1,3
                    tr(k)= tr(k) + stress(k,j)*nc(j)
                    sign = sign + nc(k)*stress(k,j)*nc(j)  
                  enddo
                enddo                 
                
                ! Cauchy traction vector acting transversally
                do k= 1,3
                  trs(k)= tr(k)-sign*nc(k)
                enddo          
                
                ! Magnitude of this vector (sigs)
                sigs = 0d0
                do j=1,3 
                  sigs = sigs + (trs(j))**2
                enddo           
                sigs = sqrt(sigs)   
                 
                ! Store the peak normal and transverse tractions for the traction seperation law
                Hist(1+(ii-1)*nn+tn1+360) = sign
                Hist(1+(ii-1)*nn+tn1+361) = sigs  
                
                cflag(1) = .true.  
                cflag(2) = .true.
                
              endif                
              
            endif
       
c================ CHECK IF ELEMENT CAN BE A CRACK ROOT ========
       
c================ AUTOMATIC INITIALIZATION ONLY!!! ========
            if(cflag(2)) then
              ! initialize ctip (for automatic detection -1)
              ! if(cmat(icmat).lt.0) then
	        if(imo.eq.-1) then
                if(cmat(icmat).gt.0) then  ! initialize only if not yet initilized in that material
	            call auto_cini_11 (ii,N0,nmat)
                  ! inactivate automatic root definitin in that material
	            cmat(icmat) = -cmat(icmat)
	          endif
              ! initialize ctip (for automatic detection -2)
              elseif(imo.eq.-2) then
                ! compute nonlocal radius
                R = acrack_ini*Vei**0.33333333333d0
	          R2= R*R
                !R2 = acrack_ini**2*Vei**0.6666666667d0
                !get the center of the current element
	          do k= 1,3
	            xc(k) = Hist(1+(ii-1)*nn+tn1+pXbar-1+k)
	          enddo
                ! check distance to existing cracks
                dist=1.0d20
                ! search on existing crack surface
	          do k = 1,nsurf
	            if (matnc(k).eq.nmat) then  ! compare just within the same material
                    ! compute center of the surface element
	              xb_surf(1) = cnode(k,1) + cnode(k,4)+
     1                           cnode(k,7) + cnode(k,10)
	              xb_surf(2) = cnode(k,2) + cnode(k,5)+
     1                           cnode(k,8) + cnode(k,11)
	              xb_surf(3) = cnode(k,3) + cnode(k,6)+
     1                           cnode(k,9) + cnode(k,12)
	              xb_surf = xb_surf/ncorn(k)
	              tempv= xb_surf-xc
                    !temp = tempv(1)*N0(1)+tempv(2)*N0(2)+tempv(3)*N0(3)
 	              temp = tempv(1)
                    !temp=tempv(1)*tempv(1)+tempv(2)*tempv(2)+
                    !tempv(3)*tempv(3)
	              if(temp.lt.dist) dist=temp
	            endif
	          enddo
                ! search on current crack tip (includes automatic initilialized tets!)
	          do k = 1,ncTip
	            if (mTi(k).eq.nmat) then  ! compare just within the same material
                    ! get tip point
	              xb_surf(1) = cTip_point(k,1)
	              xb_surf(2) = cTip_point(k,2)
	              xb_surf(3) = cTip_point(k,3)
	              tempv= xb_surf-xc
                    !temp = tempv(1)*N0(1)+tempv(2)*N0(2)+tempv(3)*N0(3)
	              temp = tempv(1)
                    !temp=tempv(1)*tempv(1)+tempv(2)*tempv(2)+
                    !tempv(3)*tempv(3)
	              if(temp.lt.dist) dist=temp
	            endif
	          enddo
                !if (dist.gt.R2) then ! root element
	          if (dist.gt.R) then ! root element
	            call auto_cini_11 (ii,N0,nmat)
	          endif
c		     else
c	          write (*,*)
c                write(*,*)'==========================================='
c                write(*,*)'ERROR --> NO CRACK TIP INI.CRITERION CALLED'
c                write(*,*)'==========================================='
c	          write (*,*)
	        endif

c================ STORE DATA ==================================
c update history

              ! write n0 into history
              do k= 1,3
                Hist(1+(ii-1)*nn+tn1+pN0-1+k)= N0(k)                 
              enddo
              ! write t0 into history
              Hist(1+(ii-1)*nn+tn1+pT0)= t0 

              ! load ctets
              ! retrieve element global nodes (corner nodes only!)
              do k= 1,4
	          ng(k) =  ix((ii-1)*nen1 + k)
	        enddo
	        ntets = ntets+1
	        ctets(ntets,1) = ii
              ! write only the corner nodes
	        do k= 1,4
                ctets(ntets,k+1) = ng(k)
	        enddo	      
            endif
          endif ! cini
        endif
	enddo  ! outer element loop


	end



	subroutine stress_matrix_11(sig,stress)
      implicit none

	real*8 sig(6), stress(3,3)

	    stress(1,1) = sig(1) 
          stress(1,2) = sig(4) 
          stress(1,3) = sig(6) 
          stress(2,1) = sig(4) 
          stress(2,2) = sig(2) 
          stress(2,3) = sig(5) 
          stress(3,1) = sig(6) 
          stress(3,2) = sig(5) 
          stress(3,3) = sig(3) 
	end


	subroutine defo_matrix_11(fiv,fi)
      implicit none

	real*8 fiv(9), fi(3,3)

	    fi(1,1) = fiv(1) 
          fi(1,2) = fiv(2) 
          fi(1,3) = fiv(3) 
          fi(2,1) = fiv(4) 
          fi(2,2) = fiv(5) 
          fi(2,3) = fiv(6) 
          fi(3,1) = fiv(7) 
          fi(3,2) = fiv(8) 
          fi(3,3) = fiv(9) 
	end

      subroutine media_failure_11m(sig0,N0,stress,fi,cflag)

c__________________________________________________TCG. 02.06.2004___71
c
c  compute failure criterion of the media  
c
c  Declare variable types
c     stress     cauchy stress matrix
c     fi         deformation gradient
c     sb_r       strenght in radial direction
c     sb_t       strenght in circumferential direction
c     sb_z       strenght in axial direction
c     sig_r,     normal stress in radial direction
c     sig_t,     normal stress in circumferential direction
c     sig_z,     normal stress in axial direction
c     z      absolut value of the normal stress in radial direction
c     x      absolut value of the normal stress in circumferential direction
c     y      absolut value of the normal stress in axial direction
c     r0,rc  Radial direction in reference and current frame
c     t0,tc  Circumferential direction in reference and current frame
c     z0,zc  Axial direction in reference and current frame
c
c____________________________________________________________________71
      implicit none

      integer i, j

	logical cflag

      real*8 stress(3,3), sig0, fi(3,3), finv(3,3)
	real*8 nc(3),  ttemp(3), temp, phi, n0(3), sig_n

c invert deformation gradient          
      finv = fi
      call invert(finv,3,3)
c push-forward normal vectors
	nc = n0
      call pushv10 (finv, nc)
c normalize
      temp = sqrt(nc(1)*nc(1) + nc(2)*nc(2) + nc(3)*nc(3))
      nc = nc/temp 


c compute stress components criterion
      ttemp=0.0d0
      do i= 1,3
	  do j= 1,3
	    ttemp(i) = ttemp(i) + stress(i,j)*nc(j)
	  enddo
	enddo
	sig_n = 0.0d0
	do i= 1,3
	  sig_n = sig_n + ttemp(i)*nc(i)
	enddo

c maximum magnitude
      sig_n = dsqrt(ttemp(1)*ttemp(1) + ttemp(2)*ttemp(2) + 
	1              ttemp(3)*ttemp(3))

      phi = sig_n - sig0

c set local cflag
      if(phi.gt.0.0d0) cflag = .true.
              
	end






      subroutine media_failure_11(sb_r, sb_t, sb_z, elmtnr,stress,fi,
	1                         x,y,z,tc,zc,rc,phi)
c__________________________________________________TCG. 02.06.2004___71
c
c  compute failure criterion of the media  
c
c  Declare variable types
c     stress     cauchy stress matrix
c     fi         deformation gradient
c     sb_r       strenght in radial direction
c     sb_t       strenght in circumferential direction
c     sb_z       strenght in axial direction
c     sig_r,     normal stress in radial direction
c     sig_t,     normal stress in circumferential direction
c     sig_z,     normal stress in axial direction
c     z      absolut value of the normal stress in radial direction
c     x      absolut value of the normal stress in circumferential direction
c     y      absolut value of the normal stress in axial direction
c     r0,rc  Radial direction in reference and current frame
c     t0,tc  Circumferential direction in reference and current frame
c     z0,zc  Axial direction in reference and current frame
c
c____________________________________________________________________71
      implicit none
	include  'vecdata.h'

      integer i, elmtnr, j

      real*8 stress(3,3), sig_r, sig_t, sig_z, fi(3,3)
	real*8 ttemp(3), temp, phi
	real*8 a0(3), b0(3)
	real*8 ac(3), bc(3)
	real*8 rc(3),zc(3),tc(3)
	real*8 sb_r, sb_t, sb_z
	real*8  x,y,z


c get a0 vector
      do i= 1,3
	  a0(i) = rvec(1,i,elmtnr)
	enddo
c get b0 vector
      do i= 1,3
	  b0(i) = rvec(2,i,elmtnr)
	enddo

compute loacal coordinate system (rtz)
          
c push-forward fiber vectors
	ac = a0
	bc = b0
      call pushv10 (fi, ac)
      call pushv10 (fi, bc)
c normalize
      temp = sqrt(ac(1)*ac(1) + ac(2)*ac(2) + ac(3)*ac(3))
      ac = ac/temp 
c normalize
      temp = sqrt(bc(1)*bc(1) + bc(2)*bc(2) + bc(3)*bc(3))
      bc = bc/temp 

c compute theta and z direction
      tc = ac + bc
      zc = ac - bc
c normalize
      temp = sqrt(tc(1)*tc(1) + tc(2)*tc(2) + tc(3)*tc(3))
      tc = tc/temp 
c normalize
      temp = sqrt(zc(1)*zc(1) + zc(2)*zc(2) + zc(3)*zc(3))
      zc = zc/temp 
c compute r direction
      rc(1) =-tc(3)*zc(2) + tc(2)*zc(3)
	rc(2) = tc(3)*zc(1) - tc(1)*zc(3)
	rc(3) =-tc(2)*zc(1) + tc(1)*zc(2)

c compute stress components criterion
c r component
      ttemp=0.0d0
      do i= 1,3
	  do j= 1,3
	    ttemp(i) = ttemp(i) + stress(i,j)*rc(j)
	  enddo
	enddo
	sig_r = 0.0d0
	do i= 1,3
	  sig_r = sig_r + ttemp(i)*rc(i)
	enddo
c t component
      ttemp=0.0d0
      do i= 1,3
	  do j= 1,3
	    ttemp(i) = ttemp(i) + stress(i,j)*tc(j)
	  enddo
	enddo
	sig_t = 0.0d0
	do i= 1,3
	  sig_t = sig_t + ttemp(i)*tc(i)
	enddo
c z component
      ttemp=0.0d0
      do i= 1,3
	  do j= 1,3
	    ttemp(i) = ttemp(i) + stress(i,j)*zc(j)
	  enddo
	enddo
	sig_z = 0.0d0
	do i= 1,3
	  sig_z = sig_z + ttemp(i)*zc(i)
	enddo

c compute failure criterion plus pertubation
c define x,y,z
c      x =Dmax1(sig_t, teps)
c	y =Dmax1(sig_z, teps)
c	z =Dmax1(sig_r, teps)
c	x = x-eps
c	y = y-eps
c	z = z-eps
c compute failure criterion without pertubation
c define x,y,z
      x =Dmax1(sig_t, 0.0d0)
	y =Dmax1(sig_z, 0.0d0)
	z =Dmax1(sig_r, 0.0d0)

      phi=(x/sb_t)**2 + (y/sb_z)**2 + (x/sb_r)**2 - 1.0d0


c crit 2
      phi = Dmax1((x-sb_t),(y-sb_z),(z-sb_r))
c      phi = Dmax1((x-sb_t),-1.0d0)

	end

      subroutine normal_11(sb_r, sb_t, sb_z,fi,
	1                           x,y,z,tc,zc,rc,N0,t0)
c__________________________________________________TCG. 02.06.2004___71
c
c  compute normal to failure surface
c
c
c  Declare variable types
c
c     sb_r       strenght in radial direction
c     sb_t       strenght in circumferential direction
c     sb_z       strenght in axial direction
c     z          normal stress in radial direction
c     x          normal stress in circumferential direction
c     y          normal stress in axial direction
c     Ncb        Current Normal vector in local frame
c     N0         Referential Normal vector in global frame
c     t0         normal traction to the surface
c     nc         Spatial Normal vector 
c     r0         Radial direction in reference frame
c     t0         Circumferential direction in reference frame
c     z0         Axial direction in reference frame
c____________________________________________________________________71

      implicit none


	real*8 sb_t, sb_r, sb_z, x,y,z, N0(3), Ncb(3), nc(3)
	real*8 fi(3,3), alpha, beta, temp, t0, nalpha, zalpha
	real*8 rc(3),zc(3),tc(3)


c      zalpha = sb_r**2.0d0*(x**2*sb_z**4 + y**2*sb_t**4)
c      nalpha = x**2*sb_z**4*sb_t**2 + sb_z**2*(y**2-sb_z**2)*sb_t**4


      zalpha = x**2*(sb_z/sb_t)**2 + y**2*(sb_t/sb_z)**2
      nalpha = x**2*(sb_z/sb_r)**2 + (y**2-sb_z**2)*(sb_t/sb_r)**2

	alpha = sqrt(1.0d0- zalpha/nalpha)
      beta = sb_r/(sqrt(1.0d0-y**2/sb_z**2-x**2/(sb_t**2)))

	ncb(1) = x*beta/(alpha*sb_t**2)
	ncb(2) = y*beta/(alpha*sb_z**2)
	ncb(3) = 1.0d0/alpha

c compute t0
      t0 = sqrt(x*x + y*y + z*z)

c change to global coordinate system
c      nc(1) = ncb(1)*tc(1) + ncb(2)*tc(2) + ncb(3)*tc(3)
c      nc(2) = ncb(1)*zc(1) + ncb(2)*zc(2) + ncb(3)*zc(3)
c      nc(3) = ncb(1)*rc(1) + ncb(2)*rc(2) + ncb(3)*rc(3)


c change to global coordinate system
      nc(1) = ncb(1)*tc(1) + ncb(2)*zc(1) + ncb(3)*rc(1)
      nc(2) = ncb(1)*tc(2) + ncb(2)*zc(2) + ncb(3)*rc(2)
      nc(3) = ncb(1)*tc(3) + ncb(2)*zc(3) + ncb(3)*rc(3)


c pull back to reference configuration
      N0 = nc
      call pushv10 (fi, N0)

c normalize
      temp = sqrt(N0(1)*N0(1) + N0(2)*N0(2) + N0(3)*N0(3))
      N0 = N0/temp 



      end





      subroutine normal_para_11(sb_r, sb_t, sb_z,fi,
	1                           x,y,z,tc,zc,rc,N0,t0)
c__________________________________________________TCG. 15.06.2004___71
c
c  compute normal to failure surface in parameter description
c
c
c  Declare variable types
c
c     sb_r       strenght in radial direction
c     sb_t       strenght in circumferential direction
c     sb_z       strenght in axial direction
c     z          normal stress in radial direction
c     x          normal stress in circumferential direction
c     y          normal stress in axial direction
c     Ncb        Current Normal vector in local frame
c     N0         Referential Normal vector in global frame
c     t0         Magnitude of the traction 
c     nc         Spatial Normal vector 
c____________________________________________________________________71

      implicit none


	real*8 sb_t, sb_r, sb_z, x,y,z, N0(3), Ncb(3), nc(3)
	real*8 fi(3,3), alpha, beta, temp, t0, gamma
	real*8 rc(3),zc(3),tc(3), r, dzdb
	real*8 eps, ca, sa, reps

	data eps /1.0d-12/
c normalize eps
      temp = sb_t + sb_r + sb_z
	reps = temp*eps

c compute angles
      temp = dsqrt(x**2 + y**2) + reps
	alpha = dASin(y/temp)
	beta  = dATan(z/temp)
      ca = dcos(alpha)
	sa = dsin(alpha)

c compute radius
      temp= dcos(alpha)**2/sb_t**2 + dsin(alpha)**2/sb_z**2
      gamma= sqrt(dsin(beta)**2/sb_r**2 + dcos(beta)**2*temp) + reps

      r = 1.0d0/gamma

c compute t0
      t0 = r

c compute n0
      dzdb = dcos(beta)*(gamma**2 - 2.0d0*dsin(beta)**2*
	1      (1.0d0/sb_r**2 - temp))/gamma**3

	temp = dsqrt(1.0d0 + dzdb**2)

	ncb(1) = -dzdb/temp
	ncb(2) = 0.0d0
	ncb(3) = 1.0d0/temp

c change to global coordinate system
      nc(1) = rc(1)*ncb(3) + (ca*zc(1)-sa*tc(1))*ncb(2) +
	1       (ca*tc(1)-sa*zc(1))*ncb(1)
      nc(2) = rc(2)*ncb(3) + (ca*zc(2)-sa*tc(2))*ncb(2) +
	1       (ca*tc(2)-sa*zc(2))*ncb(1)
      nc(3) = rc(3)*ncb(3) + (ca*zc(3)-sa*tc(3))*ncb(2) + 
	1       (ca*tc(3)-sa*zc(3))*ncb(1)


c pull back to reference configuration
      N0 = nc
      call pushv10 (fi, N0)

c normalize
      temp = sqrt(N0(1)*N0(1) + N0(2)*N0(2) + N0(3)*N0(3))
      N0 = N0/temp 


      end


      subroutine normal_crit2_11(sb_r, sb_t, sb_z,fi,
	1                           x,y,z,tc,zc,rc,N0,t0)
c__________________________________________________TCG. 15.06.2004___71
c
c  compute normal to failure surface in parameter description
c
c
c  Declare variable types
c
c     sb_r       strenght in radial direction
c     sb_t       strenght in circumferential direction
c     sb_z       strenght in axial direction
c     z          normal stress in radial direction
c     x          normal stress in circumferential direction
c     y          normal stress in axial direction
c     Ncb        Current Normal vector in local frame
c     N0         Referential Normal vector in global frame
c     t0         Magnitude of the traction 
c     nc         Spatial Normal vector 
c____________________________________________________________________71

      implicit none


	real*8 sb_t, sb_r, sb_z, x,y,z, N0(3), Ncb(3), nc(3)
	real*8 fi(3,3), alpha, beta, temp, t0
	real*8 rc(3),zc(3),tc(3)


      
      If(x.gt.sb_t) then
	  ncb(1) = 1.0d0
	  ncb(2) = 0.0d0
	  ncb(3) = 0.0d0
	  alpha = dATan(y/x)
	  beta  = dATan(z/dsqrt(x**2+y**2))
	  t0 = sb_t/(dcos(alpha)*dcos(beta))
	elseif(y.gt.sb_z) then
	  ncb(1) = 0.0d0
	  ncb(2) = 1.0d0
	  ncb(3) = 0.0d0
	  alpha = dATan(x/y)
	  beta  = dATan(z/dsqrt(x**2+y**2))
	  t0 = sb_z/(dcos(alpha)*dcos(beta))
	elseif(z.gt.sb_r) then
	  ncb(1) = 0.0d0
	  ncb(2) = 0.0d0
	  ncb(3) = 1.0d0
	  temp = dsqrt(x**2+y**2)
	  if (temp.gt.1.0d-8)then
	    beta  = dATan(z/temp)
	  else
	    beta = 1.57079d0
	  endif
	  t0 = sb_r/dsin(alpha)
	elseif((x.gt.sb_t).and.(y.gt.sb_z)) then
	  ncb(1) = 0.707106781d0
	  ncb(2) = 0.707106781d0
	  ncb(3) = 0.0d0
	  temp = dsqrt(sb_t**2+sb_z**2)
	  beta  = dATan(z/temp)
	  t0 = temp/dcos(beta)
	elseif((x.gt.sb_t).and.(z.gt.sb_r)) then
	  ncb(1) = 0.707106781d0
	  ncb(2) = 0.0d0
	  ncb(3) = 0.707106781d0
	  alpha= datan(y/x)
	  t0 = dsqrt((sb_t/dcos(alpha))**2+sb_r**2)
	elseif((y.gt.sb_z).and.(z.gt.sb_r)) then
	  ncb(1) = 0.0d0
	  ncb(2) = 0.707106781d0
	  ncb(3) = 0.707106781d0
	  alpha= datan(x/y)
	  t0 = dsqrt((sb_z/dcos(alpha))**2+sb_r**2)
	elseif((x.gt.sb_t).and.(y.gt.sb_z).and.(z.gt.sb_r)) then
	  ncb(1) = 0.577350269d0
	  ncb(2) = 0.577350269d0
	  ncb(3) = 0.577350269d0
	  t0 = dsqrt(sb_t**2 + sb_z**2 + sb_r**2)
	else
	  write (*,*)
        write(*,*)'==========================================='
        write(*,*)'ERROR --> NO NORMAL and T0 COMPUTED'
        write(*,*)'==========================================='
	  write (*,*)
	endif



c change to global coordinate system
      nc(1) = ncb(1)*tc(1) + ncb(2)*zc(1) + ncb(3)*rc(1)
      nc(2) = ncb(1)*tc(2) + ncb(2)*zc(2) + ncb(3)*rc(2)
      nc(3) = ncb(1)*tc(3) + ncb(2)*zc(3) + ncb(3)*rc(3)


c pull back to reference configuration
      N0 = nc
      call pushv10 (fi, N0)

c normalize
      temp = sqrt(N0(1)*N0(1) + N0(2)*N0(2) + N0(3)*N0(3))
      N0 = N0/temp 


      end



      subroutine neighbour_11(d,dlen,n1,n3,nen1,
	1                       pXbar,pNelem,pVe, hist,hlen,ix, ixlen)
c__________________________________________________TCG. 11.06.2004___71
c
c  compute neighbouring elements and store in nh3 history
c
c
c    Type 1: use absolute measure (ichl=1)   chl <= chl
c    Type 2: use relative measure (ichl=2)   chl <= chl*Ve^(1/3)
c
c    Warning: works only proper if less that 101 elments are 
c             included in the neighborhood!!
c
c
c  Declare variable types
c
c     dlen       length of material property vector
c     d          material property vector
c     hlen       length of history vector
c     hist       history 
c     n1         history pointer
c     n3         history pointer
c     pXbar      Pointer to element center
c     pNelem     Pointer to neighbouring elements
c     xc         center of element
c     xc_curr    center of current element
c     chl        non local measure in terms of characteristic element length
c     chlb       non local measure 
c
c
c       history
c
c       hr(nh3 + pNelem)	       nelem     ... number of neighbouring elements
c       hr(nh3 + pNelem+1)	   no. elem1 
c       hr(nh3 + pNelem+2)	   no. elem2 
c        ....
c       hr(nh3 + pNelem+301)	   no. elem300
c____________________________________________________________________71
      implicit none

      include  'cdata.h'
      include  'iofile.h'
      include  'nLength.h'
      include  'cTip.h'

	logical flag1, flag2

      integer dlen, hlen, j, pXbar, pNelem, pVe, ixlen
	integer ii, jj, itemp, k, n1, nn, tn1, n3, nen1, count
	integer ix(ixlen),nmat1, nmat2

	real*8 hist(hlen), d(dlen), xc_curr(3), xc(3)
	real*8 r2, chl2, Vei, chlb

      if (ichl.eq.1) then  ! absolute measure
	  chl2 = chl*chl
	endif

c pointer of n3 array
      tn1=2*n1
      nn = tn1 + n3




      do ii= 1,numel

c use for orie! =============================================================
c       do ii= 1,2695 
c use for orie! =============================================================


c get material no.
        nmat1 = int(ix((ii-1)*nen1 + nen1))
c limit to material where cracks are allowed
        flag1=.false.
        do j= 1, ncmat
	    if(abs(cmat(j)).eq.nmat1) flag1=.true.
	  enddo

c use for orie! =============================================================
c	flag1=.true.
c use for orie! =============================================================

        if (flag1) then
	    if (ichl.eq.2) then ! relative measure
c get element volume 
            Vei = abs(Hist(1+(ii-1)*nn+2*n1+pVe))
c compute nonlocal length scale
            chlb = abs(chl)*(abs(Vei))**(1.0d0/3.0d0)
            chl2 = chlb**2
	    endif
c retrieve center of current element 
	    do k= 1,3
	      xc_curr(k) = Hist(1+(ii-1)*nn+2*n1+pXbar-1+k)
	    enddo
          count = 0
          do jj= 1,numel
c use for orie! =============================================================
c          do jj= 1, 2695
c use for orie! =============================================================
            nmat2 = int(ix((jj-1)*nen1 + nen1))
c limit to neighbours to be in the same material
            flag2=.false.
	      if(nmat1.eq.nmat2) flag2=.true.

c use for orie! =============================================================
c	flag2=.true.
c use for orie! =============================================================

            if (flag2) then

	        itemp = 1+(jj-1)*nn+tn1+pXbar
		    xc(1) = Hist(itemp)
              if (abs(xc(1)-xc_curr(1)).lt.chlb) then
	          xc(2) = Hist(itemp+1)
                if (abs(xc(2)-xc_curr(2)).lt.chlb) then
	            xc(3) = Hist(itemp+2)
                  if (abs(xc(3)-xc_curr(3)).lt.chlb) then

	              r2 = (xc(1)-xc_curr(1))*(xc(1)-xc_curr(1)) + 
     1                   (xc(2)-xc_curr(2))*(xc(2)-xc_curr(2)) + 
     1                   (xc(3)-xc_curr(3))*(xc(3)-xc_curr(3))
            
                    if(r2.lt.chl2) then 
	              count = count+1
c  write element no. into history if count .le. 300
                    if (count.le.300) then
                     Hist(1+(ii-1)*nn+2*n1+pNelem + count)= float(jj)
	              endif
	              endif
	            endif
	          endif
	        endif
	      endif
	    enddo  ! inner element loop
c  write no. of neighbouring elements
          if (count.le.300) then
            Hist(1+(ii-1)*nn+2*n1+pNelem)= float(count)
	    else
	      Hist(1+(ii-1)*nn+2*n1+pNelem)= 300.0d0
	    endif
c check if count gt 300
          if (count.gt.300) then
	      write (*,*)
	      write (*,1000) ii
	      write (*,*)
	      write (iow,*)
	      write (iow,1000) ii
	      write (iow,*)
	    endif
	  endif
	enddo  ! outer element loop


1000  format('==============================================='/
     1       'WARNING --> More than 300 neighbouring elements'/
     2       'of element no.', i5,                            /
     3       '===============================================')
	end


      subroutine output_11(fu,fm,pFi,pNelem,pVe,ndm,ndf,numel,nen,
     1                   nen1,n1,n3,hist,hlen,ix,ixlen,x,xlen,u,ulen)
c__________________________________________________TCG. 01.03.2004___71
c
c  compute neighbouring elements and store in nh3 history
c
c
c  Declare variable types
c
c     hlen       length of history vector
c     hist       history 
c     n1         history pointer
c     n3         history pointer
c     pFi        Pointer to inverse deformatipn gradient
c     xc         center of element
c
c
c       history
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
c____________________________________________________________________71
      implicit none

      logical flag
      integer hlen, i, j,  pFi, fu, nodenr, numel, pNelem, nelem
	integer ii,  k, n1, nn,  n3, nen1, ndm, nen, ixlen, xlen
	integer ulen, ielem, jj, io, fm, ndf

      integer ix(ixlen), imat, pVe, in, fb(100), nnf, des
	real*8 hist(hlen), x(xlen), u(ulen), xc(3), Fi(9), wVe, Vei
	real*8 xn(3,nen), xr(3,nen), xcr(3), un(3,nen)
	real*8 r0, x0, y0, r1, x1, y1, ex, ey
	real*8 alpha(100),  alp, theta(3), ff(3,3), finv(3,3)
	real*8 xi, eta(3), etat(3), thetat(3), rt(3)
	real*8 a0(3), b0(3), r, xx, yy, temp
	real*8 ttemp(4)

      include  'iofile.h'

c unit of data file
      in = 99


c pointer of n3 array
      nn = 2*n1 + n3

c current geometry


      OPEN(UNIT=in,FILE='Orie.dat',status='OLD',IOSTAT=io)
	  IF(io.eq.0) THEN
	    read (in,*) r0, x0, y0                 ! radius and center of small circle
	    read (in,*) r1, x1, y1                 ! radius and center of big circle
	    read (in,*) 
	    read (in,*) nnf
	    if (nnf.gt.100) then
	      write(*,*) 
     1		'CRITICAL WARNING ** nnf too large in output_11'
	      write(iow,*) 
     1		'CRITICAL WARNING ** nnf too large in output_11'
	    endif
	    read (in,*) 
	    do i= 1,nnf
	      read (in,*) fb(i), alpha(i)          ! No. of fibrous body, alpha angle
	    enddo
	    read (in,*) 
	    read (in,*) des

	  ELSE
          write(*,*) '** ERROR- No Orie.dat exists'
        ENDIF 

	Close(in)

	write(*,*) 'read Orie.dat'



      ex = x1-x0
	ey = y1-y0

      do ii= 1,numel
c retrieve center of current element 
c   write element coordinates 
	  xcr=0.0d0
	  xc=0.0d0
	  imat=ix((ii-1)*nen1+nen1)
c check if data is required
        flag=.false.
        do j= 1,nnf
	    if(imat.eq.fb(j)) then
		  flag=.true.
	      alp = alpha(j)
	    else
  		endif
	  enddo

c       flag=.true.

	  if(flag) then
	    do j= 1,4
	     nodenr = ix((ii-1)*nen1+j)

	ttemp(j) = nodenr

	     do k= 1,3
	       xr(k,j) = x((nodenr-1)*ndm + k)
	       un(k,j) = u((nodenr-1)*ndf + k)
	       xn(k,j) = xr(k,j) + un(k,j)          
	       xcr(k) = xcr(k) + xr(k,j)
	       xc(k)  = xc(k)  + xn(k,j)
	     enddo
	    enddo
	  if (nen.eq.0) then
	    write (*,105) 
	    write (iow,105)
	  else
	    xc = xc/4
	    xcr = xcr/4
	  endif

c get no. of neighbouring elements 
          nelem = int(Hist(1+(ii-1)*nn+2*n1+pNelem))

c initialize quantities
	    fi  = 0.0d0
	    wVe  = 0.0d0

          do jj= 1,nelem
c get element no. of neighbouring element
	      ielem = int(Hist(1+(ii-1)*nn+2*n1+pNelem+jj))
c get element volume 
            Vei = abs(Hist(1+(ielem-1)*nn+2*n1+pVe))
c summation
	      wVe = wVe + Vei
	      do k= 1,9
	        fi(k) = fi(k) + 
	2               Vei*Hist(1+(ielem-1)*nn+2*n1+pFi-1+k)
	      enddo

	    enddo  ! inner element loop
	    if (abs(wVe).lt.1.0d-12) then
	      write (*,106) 
	      write (iow,106)
	    else
	      fi = fi/wVe
	    endif

c no average if required
	    if (nelem.eq.0) then
	      do k= 1,9
	        Fi(k) = Hist(1+(ii-1)*nn+2*n1+pFi-1+k)
	      enddo
	    endif
c	    write(fu,100) imat, (xcr(j),j=1,3),(xc(j),j=1,3), (Fi(k),k=1,9)
c	    write(fu,99) imat, ii, (ttemp(j),j=1,4)
c compute fiber direction
          ff(1,1) = fi(1)
          ff(1,2) = fi(2)
          ff(1,3) = fi(3)
          ff(2,1) = fi(4)
          ff(2,2) = fi(5)
          ff(2,3) = fi(6)
          ff(3,1) = fi(7)
          ff(3,2) = fi(8)
          ff(3,3) = fi(9)
c invert ff
          finv = ff
          call invert(finv,3,3)
c compute circumferential direction vector 
	    xx = xc(1)
	    yy = xc(2)
          r= (r1*((xx - x0)*(x0 - x1) + (yy - y0)*(y0 - y1)) - 
     -    r0*((xx - x1)*(x0 - x1) + (yy - y1)*(y0 - y1)) + 
     -    Sqrt((r0 - r1)**2*(r1**2*((xx - x0)**2 + (yy - y0)**2)
     -       - 2*r0*r1*(xx**2 + x0*x1 -  xx*(x0 + x1) + (yy - y0)*
     -         (yy - y1)) + r0**2*((xx - x1)**2 + (yy - y1)**2)
     -       - (x0*yy -  x1*yy - xx*y0 + x1*y0 +xx*y1 - x0*y1) **2)))/
     -    ((r0 - r1)**2 - (x0 - x1)**2 - (y0 - y1)**2)          

	    if (abs(r1-r0).lt.1.0d-12) then
	      write (*,107) 
	      write (iow,107)
	    else
	      xi = (r-r0)/(r1-r0)
	    endif
	    theta(1) =   yy-y0-xi*ey
	    theta(2) = -(xx-x0-xi*ex)
	    theta(3) =   0.0d0
c pull back vector
          call pushv10_contra (finv, theta)
c normalize 
          temp= sqrt(theta(1)**2 + theta(2)**2 + theta(3)**2)
	    if (abs(temp).lt.1.0d-12) then
	      write (*,108) 
	      write (iow,108)
	    else
	      theta = theta / temp
	    endif

c compute axial direction vector 
	    eta(1) = 0.0d0
	    eta(2) = 0.0d0
	    eta(3) = 1.0d0
c pull back vector
          call pushv10_contra (finv, eta)
c normalize 
          temp= sqrt(eta(1)**2 + eta(2)**2 + eta(3)**2)
	    if (abs(temp).lt.1.0d-12) then
	      write (*,109) 
	      write (iow,109)
	    else
	      eta = eta / temp
	    endif

c define orthonormal system	  
	    if (des.eq.1) then        ! theta is ok
	      thetat= theta
	      call cross_prod_11(thetat,eta,rt)
            temp= sqrt(rt(1)**2 + rt(2)**2 + rt(3)**2)
	      if (abs(temp).lt.1.0d-12) then
	        write (*,110) 
	        write (iow,110)
	      else
	        rt = rt / temp
	      endif
	      call cross_prod_11(rt,thetat,etat)
		else				      ! eta is ok
	      etat = eta
	      call cross_prod_11(theta,etat,rt)
            temp= sqrt(rt(1)**2 + rt(2)**2 + rt(3)**2)
	      if (abs(temp).lt.1.0d-12) then
	        write (*,111) 
	        write (iow,111)
	      else
	        rt = rt / temp
	      endif
	      call cross_prod_11(etat,rt,thetat)
		endif 
c compute fiber vectors
          a0 = thetat + tan(alp)*etat
          temp= sqrt(a0(1)**2 + a0(2)**2 + a0(3)**2)
	    if (abs(temp).lt.1.0d-12) then
	      write (*,112) 
	      write (iow,112)
	    else
	      a0 = a0 / temp
	    endif
          b0 = thetat - tan(alp)*etat
          temp= sqrt(b0(1)**2 + b0(2)**2 + b0(3)**2)
	    if (abs(temp).lt.1.0d-12) then
	      write (*,113) 
	      write (iow,113)
	    else
	      b0 = b0 / temp
	    endif
	    write(fu,101) imat, ii, (a0(j),j=1,3), (b0(j),j=1,3), 
	1                            (xcr(j),j=1,3)
	    write(fm,102) ii,(a0(j),j=1,3),(b0(j),j=1,3),(rt(j),j=1,3)
	1                            
	  endif



	enddo  ! outer element loop


99    format(i5,i5, i5,i5,i5,i5)
100   format(i5, e13.5,',', e13.5,',',e13.5,',',e13.5,',',e13.5,
     1        ',',e13.5,',',e13.5,',',e13.5,',',e13.5,',',e13.5,',',
     2	    e13.5,',',e13.5,',',e13.5,',',e13.5,',',e13.5)
101   format(i5, i7, e13.5, e13.5, e13.5,
     1               e13.5, e13.5, e13.5,
     1               e13.5, e13.5, e13.5)
102   format(i5,  e13.5, e13.5, e13.5, e13.5, e13.5, e13.5, 
     1            e13.5, e13.5, e13.5)
105   format(
     1 5x, 'ERROR1 output_11: Division by zero'/)
106   format(
     1 5x, 'ERROR2 output_11: Division by zero'/)
107   format(
     1 5x, 'ERROR3 output_11: Division by zero'/)
108   format(
     1 5x, 'ERROR4 output_11: Division by zero'/)
109   format(
     1 5x, 'ERROR5 output_11: Division by zero'/)
110   format(
     1 5x, 'ERROR6 output_11: Division by zero'/)
111   format(
     1 5x, 'ERROR7 output_11: Division by zero'/)
112   format(
     1 5x, 'ERROR8 output_11: Division by zero'/)
113   format(
     1 5x, 'ERROR9 output_11: Division by zero'/)
	end


      subroutine cross_prod_11(a,b,c)
c__________________________________________________TCG. 01.03.2004___71
c
c  compute vector product c = a x b
c
c__________________________________________________TCG. 01.03.2004___71
      implicit none

	real*8 a(3), b(3), c(3)

	c(1) = -a(3)*b(2) + a(2)*b(3)
	c(2) =  a(3)*b(1) - a(1)*b(3)
	c(3) = -a(2)*b(1) + a(1)*b(2)

	end


      subroutine get_normal_11(nold, n, lamb, n0)
c__________________________________________________TCG. 06.08.2004___71
c
c  get normal for concrete model
c  regular and irregular cases
c
c  nold average normal of activated pu elements (not normalized)
c  n    eigenvectors
c  lamb eigenvalues
c  n0   normal
c
c__________________________________________________TCG. 06.08.2004___71

      implicit none

	logical abbr

	integer i, j

	real*8 nold(3), n(3,3), n0(3), lamb(3)
	real*8 temp, fact, t1,t2,t3

      fact = 0.2d0 

       
c sort values 
      
      abbr=.false.
      do while(not(abbr))
	  abbr =.true.
	  do i=1,2
	    if (lamb(i).lt.lamb(i+1)) then
	      temp      = lamb(i)
	      lamb(i)   = lamb(i+1)
	      lamb(i+1) = temp
	      do j= 1,3
	        temp     = n(j,i)
	        n(j,i)   = n(j,i+1)
	        n(j,i+1) = temp
	      enddo
	      abbr      =.false.
	    endif
	  enddo
	enddo

c regular case
	do i= 1,3
	  n0(i) = n(i,1)
	enddo
      

c close eigenvalues
	fact = abs(fact*lamb(1))
      
c case 1 
      t1=lamb(1)-lamb(2)
      if(t1.lt.fact) then
	  t1 = nold(1)*n(1,1)+nold(2)*n(2,1)+nold(3)*n(3,1)
	  t2 = nold(1)*n(1,2)+nold(2)*n(2,2)+nold(3)*n(3,2)
	  if (abs(t1).gt.abs(t2)) then
	    do i= 1,3
	      n0(i) = n(i,1)
	    enddo
	  else 
	    do i= 1,3
	      n0(i) = n(i,2)
	    enddo
	  endif
	endif

c case 2
      t1=lamb(1)-lamb(3)
      if(t1.lt.fact) then
	  t1 = nold(1)*n(1,1)+nold(2)*n(2,1)+nold(3)*n(3,1)
	  t3 = nold(1)*n(1,3)+nold(2)*n(2,3)+nold(3)*n(3,3)
	  if (abs(t1).gt.abs(t3)) then
	    do i= 1,3
	      n0(i) = n(i,1)
	    enddo
	  else 
	    do i= 1,3
	      n0(i) = n(i,3)
	    enddo
	  endif
	endif

c case 3
      t1=lamb(1)-lamb(2)
      t2=lamb(1)-lamb(3)
      if((t1.lt.fact).and.(t2.lt.fact)) then
	  t1 = nold(1)*n(1,1)+nold(2)*n(2,1)+nold(3)*n(3,1)
	  t2 = nold(1)*n(1,2)+nold(2)*n(2,2)+nold(3)*n(3,2)
	  t3 = nold(1)*n(1,3)+nold(2)*n(2,3)+nold(3)*n(3,3)
	  if ((abs(t1).gt.abs(t2)).and.(abs(t1).gt.abs(t3))) then
	    do i= 1,3
	      n0(i) = n(i,1)
	    enddo
	  endif 
	  if ((abs(t2).gt.abs(t1)).and.(abs(t2).gt.abs(t3))) then
	    do i= 1,3
	      n0(i) = n(i,2)
	    enddo
	  endif 
	  if ((abs(t3).gt.abs(t1)).and.(abs(t3).gt.abs(t2))) then
	    do i= 1,3
	      n0(i) = n(i,3)
	    enddo
	  endif 
	endif




c    
      t1=abs(lamb(1)-lamb(2))
      t2=abs(lamb(1)-lamb(3))
      if((t1.lt.fact).or.(t2.lt.fact)) then
	  n0 = nold
	endif
      end 




      subroutine get_normal_11_2(n, lamb, n0)
c__________________________________________________TCG. 06.08.2004___71
c
c  get normal for concrete model
c  regular and irregular cases
c
c  n    eigenvectors
c  lamb eigenvalues
c  n0   normal
c
c__________________________________________________TCG. 06.08.2004___71

      implicit none

      include  'iofile.h'
	logical abbr

	integer i, j

	real*8 n(3,3), n0(3), lamb(3)
	real*8 temp


       
c sort values 
      
      abbr=.false.
      do while(not(abbr))
	  abbr =.true.
	  do i=1,2
	    if (lamb(i).lt.lamb(i+1)) then
	      temp      = lamb(i)
	      lamb(i)   = lamb(i+1)
	      lamb(i+1) = temp
	      do j= 1,3
	        temp     = n(j,i)
	        n(j,i)   = n(j,i+1)
	        n(j,i+1) = temp
	      enddo
	      abbr      =.false.
	    endif
	  enddo
	enddo

	do i= 1,3
	  n0(i) = n(i,1)
	enddo


c normalize
      temp = sqrt(N0(1)*N0(1) + N0(2)*N0(2) + N0(3)*N0(3))
      N0 = N0/temp 

      end


      subroutine mod_N0_11(ndm, nen1, xlen, x, hlen, hist, ixlen, ix,  
	1                    el, ns, sflag, pN0, pcini, nn, tn1, N0)

c__________________________________________________TCG. 24.10.2004___71
c
c  correct normal if crack tries to run back into PU finite element
c
c
c  ndm   Number of dimension =3
c  xlen  Length of the x array
c  x     Nodal coordinate array
c  hlen  Length of the hist array
c  hist  history array
c  el    No. of current element
c  ns    No. of shared facets with the crack tip
c  sflag id of the shared facets with the crack tip
c  pN0   pointer to the discontinuity orientation
c  N0    (redefined) orientation of the discontinuity
c  nold  average normal of activated pu elements (not normalized)
c  noldb average of the old normals
c
c
c____________________________________________________________________71
      implicit none

      include  'cdata.h'
      include  'iofile.h'
      include  'cTip.h'

      integer ns, ixlen, nen1, nmat

	logical c,flagi(ns), flag


      integer xlen, hlen, pN0, nn, nbg, ng(4), tn1,  ix(ixlen)
	integer ndm, i,j,k,l, el, oel(ns), sflag(4), nodenr(3), cini
	integer pcini, nsa

	real*8 N0(3), Nold(3,ns), xtria(3,3), a(3), b(3), Nf(3,ns)
	real*8 x(xlen), Noldb(3)
	real*8 NfN0(3), NfNold(3), orie, hist(hlen), Nbar(3), norm

c init
      Nold=0.0d0


c step 1: find neighbor tets with embedded discontinuity (nsa)

      nsa=0
	do l=1,ns ! loop over all shared ctip facets
	  do j= 1,3
	    nodenr(j) = cTip(sflag(l),j)  ! nodes of the l-th ctip facets
	    do k= 1,ndm
	      xtria(k,j) = x((nodenr(j)-1)*ndm + k) ! coordinates of the nodes 
	    enddo
	  enddo


c find tets on the opposite side of the tipfacets
        do i=1, numel  ! loop over all tets
c get material no.
          nmat = int(ix((i-1)*nen1 + nen1))
c limit to material where cracks are allowed
          flag=.false.
          do j= 1, ncmat
	      if(abs(cmat(j)).eq.nmat) flag=.true.
	    enddo
	    if (flag) then
c get cini-flag
            cini = int(Hist(1+(i-1)*nn+tn1+pCini))
	      if (cini.eq.1) then
c retrieve element global nodes (only corner nodes!)
              do k= 1,4
	          ng(k) =  ix((i-1)*nen1 + k)
	        enddo
c compare current facet with the i-element
	        c=.false.
	        If(nodenr(1).eq.ng(1)) c=.true.
	        If(nodenr(1).eq.ng(2)) c=.true.
	        If(nodenr(1).eq.ng(3)) c=.true.
	        If(nodenr(1).eq.ng(4)) c=.true.
	        if (c) then
	          c=.false.
	          If(nodenr(2).eq.ng(1)) c=.true.
	          If(nodenr(2).eq.ng(2)) c=.true.
	          If(nodenr(2).eq.ng(3)) c=.true.
	          If(nodenr(2).eq.ng(4)) c=.true.
	          if (c) then
	            c=.false.
	            If(nodenr(3).eq.ng(1)) c=.true.
	            If(nodenr(3).eq.ng(2)) c=.true.
	            If(nodenr(3).eq.ng(3)) c=.true.
	            If(nodenr(3).eq.ng(4)) c=.true.
	            if (c) then  ! opposite tet found
c compte normal to facet
                    do j= 1,3
	                a(j) = xtria(j,2) - xtria(j,1)
	                b(j) = xtria(j,3) - xtria(j,1)
	              enddo
                    Nf(1,l) = -a(3)*b(2) + a(2)*b(3)
	              Nf(2,l) =  a(3)*b(1) - a(1)*b(3)
	              Nf(3,l) = -a(2)*b(1) + a(1)*b(2)
c store element number of opposite tet
                    oel(l) = i
c increase bad geometry counter
	              nsa=nsa+1
c read direction vectors of discontinuity in opposite tet
	              do k= 1,3
	                Nold(k,l) = Hist(1+(i-1)*nn+tn1+pN0-1+k)
	              enddo
  	            endif
	          endif
	        endif
	      endif	   	   
	    endif	   	   
	  enddo
	enddo

cstep2: modify normal if it is too far away from previous
c       non-uniqueness of the Rankine criterion!


      Noldb=0.0d0
      do i= 1, ns
	  do k= 1,3
	    Noldb(k) =  Noldb(k) + Nold(k,i)
	  enddo
	enddo

c normalize or define Noldb 
      write(*,*) 'compute norm'
	norm= Noldb(1)*Noldb(1) + Noldb(2)*Noldb(2) + 
     1            Noldb(3)*Noldb(3)
      write(*,*) 'norm=', norm

      norm = sqrt(Noldb(1)*Noldb(1) + Noldb(2)*Noldb(2) + 
     1            Noldb(3)*Noldb(3))
      write(*,*) 'norm=', norm
      if (norm.eq.0.0d0) then  
	  write (*,*) 'WARNING: No Nold available of element', el
        write (iow,*) 'WARNING: No Nold available of element', el
	else
      write(*,*) 'Normalize Noldb'
 	  Noldb = Noldb/norm
c check criterion
        write(*,*) 'check criterion'
        norm = N0(1)*Noldb(1) + N0(2)*Noldb(2) + N0(3)*Noldb(3)
      
        if (abs(norm).lt.0.707d0) then ! more than 45 degree change
	    N0 = Noldb
          write (*,*) 'Change N0 of elem. no. (step2 in mod_N0_11)',el
          write (iow,*) 
	1                'Change N0 of elem. no. (step2 in mod_N0_11)',el
          write (iow,*) 'N0.Nold=', norm
	  endif
	endif


cstep 3:  Check bad geometry criteria of the different facets
      flagi=.false.
	nbg=0
      do i=1,nsa
c compute x products
        do j= 1,3
          a(j) = Nf(j,i)
	    b(j) = N0(j)
	  enddo
        NfN0(1) = -a(3)*b(2) + a(2)*b(3)
	  NfN0(2) =  a(3)*b(1) - a(1)*b(3)
	  NfN0(3) = -a(2)*b(1) + a(1)*b(2)

        do j= 1,3
          a(j) = Nf(j,i)
	    b(j) = Nold(j,i)
	  enddo

        NfNold(1) = -a(3)*b(2) + a(2)*b(3)
	  NfNold(2) =  a(3)*b(1) - a(1)*b(3)
	  NfNold(3) = -a(2)*b(1) + a(1)*b(2)

c check orientation
        orie = NfN0(1)*NfNold(1)+NfN0(2)*NfNold(2)+NfN0(3)*NfNold(3)
	  if (orie.lt.0.0d0) then
	    flagi(i)=.true.
	    nbg = nbg+1
	  endif
	enddo



c redefine N if required

	if (nbg.gt.0) then
	  Nbar = 0.0d0
	  do j= 1,nbg
	    do i= 1,3
	      Nbar(i) = Nbar(i) + Nold(i,j)
	    enddo
	  enddo
c normalize
        norm = sqrt(Nbar(1)*Nbar(1) + Nbar(2)*Nbar(2) + 
	1              Nbar(3)*Nbar(3))
 	  Nbar = Nbar/norm

	  N0 = Nbar


        write (*,*) 'Change N0 of elem. no. (step3 in mod_N0_11)',el
        write (iow,*) 'Change N0 of elem. no. (step3 in mod_N0_11)',el
        write (iow,*) 'Number of facets', nbg

	endif

      end



      subroutine mod_princ_dir_11(ndm, ii, n1, n3, nen1, 
	1          pN0, pcini, pNelem, pNold,
     2          pAe, hlen, hist, ixlen, ix, N0)
c__________________________________________________TCG. 24.10.2004___71
c
c  modify normal due to non-uniqueness of the principal stress criterion
c
c
c     ndm        Number of dimension =3
c     ii         No. of current element
c     N0        (redefined) orientation of the discontinuity
c     pN0        Pointer to normal vector
c     pNelem     Pointer to neighbouring elements
c     pAe        Pointer to area
c     pNold      Pointer to old normal
c     hlen       length of history vector
c     hist       history 
c     n1         history pointer
c     n3         history pointer
c     ix         ix vector
c     ixlen      Length of ix vector
c     nmat       no. of material
c
c____________________________________________________________________71

      implicit none

      include  'cdata.h'
      include  'eldata.h'
      include  'iofile.h'
      include  'cTip.h'

       
      integer ndm, hlen, pN0, pcini, pNelem, pAe
	integer ii, jj, k, pNold, nelem, ielem
	integer cini, n1, n3, ixlen, nen1, tn1, nn
	integer ix(ixlen)

	real*8 hist(hlen), N0(3), Nold(3), Ae, temp


c pointer of n3 array
      tn1=2*n1
      nn = tn1 + n3

c get no. of neighbouring elements 
      nelem = int(Hist(1+(ii-1)*nn+tn1+pNelem))
c no average if required
	if (nelem.eq.0) then
c nothing to do yet
	else
c averaging
c initialize normal
	  nold = 0.0d0
        do jj= 1,nelem

c get element no. of neighbouring element
	   ielem = int(Hist(1+(ii-1)*nn+tn1+pNelem+jj))
c get cini-flag
         cini = int(Hist(1+(ielem-1)*nn+tn1+pCini))
	   if (cini.eq.1) then
c get discontinuity area 
            Ae = abs(Hist(1+(ielem-1)*nn+tn1+pAe))
c summation
	      do k= 1,3
	        Nold(k) = Nold(k) + 
	2        Ae*Hist(1+(ielem-1)*nn+tn1+pN0-1+k)
	      enddo
	    endif
	  enddo  ! inner element loop

c normalization 
        temp = sqrt(Nold(1)*Nold(1) + Nold(2)*Nold(2) + Nold(3)*Nold(3))

c        if (temp.eq.0.0d0) then  
cc define nold for cardif test
c	    nold(1) =0.707d0
c	    nold(2) =0.707d0
c	    nold(3) =0.0d0
c	    temp    =1.0d0
c	  endif


        if (temp.ne.0.0d0) then  
          Nold = Nold/temp      
          temp = nold(1)*n0(1) + nold(2)*n0(2) + nold(3)*n0(3)   
          if (abs(temp).lt.0.5d0) then
	      n0 = nold
	      write (*,*) 'Changed N0 to Nold of element', ii
            write (iow,*) 'Changed N0 to Nold of element', ii
	    endif
        else
	    write (*,*) 'WARNING: No Nold available of element', ii
          write (iow,*) 'WARNING: No Nold available of element', ii
        endif
      endif

c write nold into history
      do k= 1,3
	  Hist(1+(ielem-1)*nn+tn1+pNold-1+k)= Nold(k) 
      enddo
	end

      subroutine Change_ref_config_11(ndm,numnp,x,xlen,u,ulen)
c	subroutine Change_ref_config_11(nn,n1,numel,pFtot,pFres, 
c	1           hist,hlen)
c__________________________________________________TCG. 01.06.2006___71
c
c      updates reference configuration
c
c  Declare variable types
c
c     ndm         Number of dimension
c     numel       Number of elemnts
c     hist        History vector
c     hlen        Length of History vector
c     pFres       pointer to residual deformation
c     pFtot       pointer to total deformation
c____________________________________________________________________71
      implicit none

      include  'iofile.h'

	integer hlen, numnp, ndm
	integer xlen, ulen
	integer i, j, k, n1,nn
c	real*8 hist(hlen), pFtot, pFres
	real*8 x(xlen), u(ulen), uc(ndm), temp(xlen)


cc store Ftot into Fres
c      do i= 1,numel
cc      do i= 1,11411
c	  do k= 1,9
c	    Hist(1+(i-1)*nn+2*n1+pFres-1+k) = 
c	1    Hist(1+(i-1)*nn+2*n1+pFtot-1+k)
c	  enddo
c	enddo 


      
      k=0
      do i= 1,numnp
	  uc = 0.0d0
	  do j= 1,3
c	    uc(j) = u(numnp*3+k + j)
c	    x(k+j) = x(k+j) - uc(j)
	    uc(j) = u(k + j)
	    temp(k+j) = x(k+j) - uc(j)
	  enddo
	  k = k+3 
	enddo







c write to file
	write(73,*)  'COOR'
      k=0 
      do i= 1,numnp
	  write(73,1000)  i,(temp(k+j),j=1,3)
	  k=k+3 
	enddo




1000	format(1x,i5,' 0 ',f18.8,f18.8,f18.8)
	end


      subroutine multiply_F_11(f_res,f)
c__________________________________________________TCG. 30.06.2006___71
c
c      Set multiply deformation gradients
c
c  Declare variable types
c
c     f         deformation gradient
c     f_res     deformation gradient
c____________________________________________________________________71
      implicit none

	integer i, j, k
	real*8 f(3,3), f_res(3,3), temp(3,3)

      temp=0.0d0
      do i= 1,3
        do j= 1,3
          do k= 1,3
	      temp(i,j) = temp(i,j) + f_res(i,k)*f(k,j)
c	      temp(i,j) = temp(i,j) + f(i,k)*f_res(k,j)
	    enddo
	  enddo
	enddo
   
c store result into F
      f = temp

      end
      
      subroutine  reallo_hist_variables(hist,hlen,n1,n3)
c__________________________________________________TCG. 27.04.2021___71
c   transfers the histrory values from the displacement FEM element to the PUFEM element
c   (at the time point of inserting the cohesive zone) 
c   averages the history variables over all integration points of the 
c   displacement FEM element and then writes it into the history of the PUFEM element
c
c  Declare variable types
c     hist        history vector
c     hlen        Length of history vector per element
c     h1len       Length of n1 (and n2) history vector per elment
c     h3len       Length of n3 history vector per elment
c     l           Integration order
c     lintTET     # of integration points of the Q1 element (tet)
c     lintPU      # of integration points of the PUFEM element 
c     nhv         # of elements in the history vector per integration point
c____________________________________________________________________71
      implicit none
      include  'elHist.h'
      
	integer  n, hlen, n1,n3, nn, lintTET, lintPU
     	real*8   hist(hlen), value, sg(5,87)
      integer i,j,k, l, nhv  ! Moved sg(5,87) from this line to the above line CM (30/07/2022)

      nhv = nHistpGP   ! length of hist vector per Gauss point
      call tint3d_11(iorder,lintTET,sg)      
      lintPU = 2*(iorder*iorder*iorder)
     
!      if (nhv.gt.0) then
!        Do i= 1, nhv ! loop over all elements in the history vector per integration point
!          nn= i+ n1    ! read the history of t_n+1
!          nn= i+ n1    ! read the history of t_n+1
!!          nn= i    ! read the history of t_n
!          value= 0.0
!          Do j= 1, lintTET ! loop over all integration points of the tet
!            value = value + hist(nn)
!            nn = nn+ nhv
!          enddo
!          if (lintTET.gt.0) then 
!            value = value/lintTET   ! average all history values
!          else
!	      write (*,*)
!            write(*,*)'=================================='
!            write(*,*)'ERROR --> zero integration points'
!            write(*,*)'          failed to average'
!            write(*,*)'=================================='
!	      write (*,*)
!          endif
! ! write new value into history of PUFEM element         
!          nn= i
!          Do j= 1, lintPU ! loop over all integration points of the PUFEM element
!            hist(nn)     = value     ! history at t_n
!            hist(nn+ n1) = value     ! history at t_n+1
!            nn = nn+ nhv
!          enddo          
!        enddo
!      endif

      
      if (nhv.gt.0) then
        Do i= 1, nhv ! loop over all elements in the history vector per integration point
          nn= i+ n1    ! read the history of t_n+1
          !nn= n1 + (i-1)    ! read the history of t_n+1
!          nn= i    ! read the history of t_n
          value= 0.0
          Do j= 1, lintTET ! loop over all integration points of the tet
            value = value + hist(nn)
            nn = nn+ nhv
          enddo
          if (lintTET.gt.0) then 
            value = value/lintTET   ! average all history values
          else
	      write (*,*)
            write(*,*)'=================================='
            write(*,*)'ERROR --> zero integration points'
            write(*,*)'          failed to average'
            write(*,*)'=================================='
	      write (*,*)
          endif
 ! write new value into history of PUFEM element         
          nn= i
          Do j= 1, lintPU ! loop over all integration points of the PUFEM element
            hist(nn)     = value     ! history at t_n
            hist(nn+ n1) = value     ! history at t_n+1
            nn = nn+ nhv
          enddo          
        enddo
      endif      
      

      end


      subroutine str_ell(ddv,f)
c
c--------------------------------------------------------------------CM 2022
c
c     Compute Determinant of the Acoustic tensor for each of the
c     principal referential fiber directions. Take the minimum of
c     cases and store in the time independant history for use in the 
c     crack initialisation subroutine (asseble_ctets_11.f)     
c
c.... INPUT variables
c         ddv     Stiffness in Voigt notation
c         f       Deformation Gradient
c
c--------------------------------------------------------------------
      
      IMPLICIT NONE
      
      integer i, j, k, l, p, q, r, idx(6,2), idv
      real*8 F(3,3),Finv(3,3), detF, Ft(3,3)
      real*8 sigv(6), sig(3,3), FPK(3,3), m(3)
      real*8 H(3,3), eg(3)
      real*8 normn, normm
      real*8 ddv(6,6), dd(3,3,3,3), A(3,3,3,3), Am(9,9)
      real*8 detH, detH_case(4), check, detA, detD
      logical str_ell_flag
      real*8 fact(9,9), IPVT(9), det1, det2, det3
      real*8 che(3,3), f1, test(3)
      real*8 one3, jm13, m_out(3), Ftmp(3,3), Htmp(3,3)
      integer ipiv(9), info, str_ell_case
      data one3  /0.33333333333333333/
      
      !include 'Visco.h'
      include 'pvecdata.h'
      include 'eldata.h'
      include  'hdata.h'
      include  'comblk.h'
            
c.... Convert stiffness to index notation (3x3x3x3)     
            
      idx(1,1) = 1
      idx(2,1) = 2
      idx(3,1) = 3
      idx(4,1) = 1
      idx(5,1) = 2
      idx(6,1) = 1
      idx(1,2) = 1
      idx(2,2) = 2
      idx(3,2) = 3
      idx(4,2) = 2
      idx(5,2) = 3
      idx(6,2) = 3     
      
      dd = 0d0
      do i=1,6
          
        if (i.le.3) then   
          sig(idx(i,1),idx(i,2)) = sigv(i)     
        elseif (i.gt.3) then
          sig(idx(i,1),idx(i,2)) = sigv(i)  
          sig(idx(i,2),idx(i,1)) = sigv(i)  
        endif
            
        do j=1,6
          if ((i.le.3).and.(j.le.3)) then
            dd(idx(i,1),idx(i,2),idx(j,1),idx(j,2))=ddv(i,j)
          elseif ((i.gt.3).and.(j.le.3)) then
            dd(idx(i,1),idx(i,2),idx(j,1),idx(j,2))=ddv(i,j)
            dd(idx(i,2),idx(i,1),idx(j,1),idx(j,2))=ddv(i,j) 
          elseif ((i.le.3).and.(j.gt.3)) then
            dd(idx(i,1),idx(i,2),idx(j,1),idx(j,2))=ddv(i,j)
            dd(idx(i,1),idx(i,2),idx(j,2),idx(j,1))=ddv(i,j) 
          elseif ((i.gt.3).and.(j.gt.3)) then    
            dd(idx(i,1),idx(i,2),idx(j,1),idx(j,2))=ddv(i,j)
            dd(idx(i,2),idx(i,1),idx(j,1),idx(j,2))=ddv(i,j)
            dd(idx(i,1),idx(i,2),idx(j,2),idx(j,1))=ddv(i,j)
            dd(idx(i,2),idx(i,1),idx(j,2),idx(j,1))=ddv(i,j)
          endif  
        enddo
        
      enddo
      
c.... How many cases do we consider, based on the oritentation distribution function used
      
      if (ind.eq.4) then ! Bingham Distribution 
        ! Bingham Distribution has four possible directions 
        !idv = 4  ! In circumferential, axial and the two intermediate directions at +-45 degrees
        !idv = 2  ! In circumferential and axial 
        idv = 1  ! Just in circumferential
      elseif ((ind.eq.1).or.(ind.eq.2).or.(ind.eq.3)) then ! Any other distribution
        ! Uniaxial distribution, and Von mises distribution only have one fiber direction     
        idv=1
      endif
        
c.....We're using the deviatoric deformation gradient        
        
      Ftmp = F  
      call eig3(Ftmp,eg,i)
      detF = eg(1)*eg(2)*eg(3)  
        
      jm13 = detF**(-one3)
      Ft = jm13*F    
        
c.... Loop over each of the vectors to find acoustic vector determinant      
      
      !test(1) = 0.235702260395516
      !test(2) = 0.235702260395516
      !test(3) = 0.942809041582064
        
      !do q = 2
      do q=1,idv 
      
        m = 0.0d0
        do i= 1,3
          do j= 1,3
            !m(i)= m(i)+Ft(i,j)*qvec(n,j,q)
            m(i)= m(i)+Ft(j,i)*qvec(n,j,q)
          enddo
        enddo 
      
        normm = 0d0
        do j=1,3 
          normm = normm + (m(j))**2
        enddo
        m = m/sqrt(normm)     
      
c.... Determine the acoustic tensor  
      
        H = 0d0
        do j=1,3
          do k=1,3
            do i=1,3
              do l=1,3
                H(j,k) = H(j,k) + m(i)*dd(i,j,k,l)*m(l)
              enddo
            enddo
          enddo
        enddo 
         
c.... Determine the determinant of the acoustic tensor
      
        Htmp = H 
        call eig3(Htmp,eg,i)
        detH_case(q) = eg(1)*eg(2)*eg(3)
        
      enddo    
        
c.....Take the minimum of the determinant for each of the cases  
      
      q = 1
      detH = detH_case(1)
      do k = 1,idv  
        if (detH_case(k).lt.detH) then  
          detH = detH_case(k)
          q = k
        endif
      enddo

c.....Populate the minimum determinant in the time independant history               
      
      hr(nh3+362) = detH       
      
c.....Populate the corresponding referential normal in the time independant history       
      
      do i= 1,3
        hr(nh3+362+i) = qvec(n,i,q)
      enddo   
      
c.....Check if element alone meets crack initialisation criteria      
      
      if (detH.lt.0d0) then
         check = 1d0
      endif
      
      end
      
      
      
      
      
      
      
      
      
      
      
      
      
c.... Determine the intermediate stiffness     
      
     !! A = 0d0
     !! do i=1,3
     !!   do j=1,3
     !!     do k=1,3
     !!       do l=1,3
     !!         do p=1,3
     !!           do q=1,3
     !!             A(p,j,q,l) = A(p,j,q,l) + 
     !!1                            Finv(p,i)*Finv(q,k)*dd(i,j,k,l)
     !!           enddo    
     !!         enddo
     !!       enddo
     !!     enddo
     !!   enddo
     !! enddo   
     !! 
     !! A = A*detF
      
c.... Determine the First Piola Kirchoff stress       
      
      !FPK = 0d0
      !do i=1,3
      !  do j=1,3 
      !     do k=1,3
      !       FPK(i,j) = FPK(i,j) + sig(i,k)*Finv(j,k)  
      !     enddo
      !  enddo
      !enddo
      !
      !FPK = FPK*detF
      
      
      
      
      
      
      
      
      
      
      
c.... Convert the stiffness into the matrix representation (9x9) 
      
      !Am = 0d0
      !Dm = 0d0
      !do i=1,3
      !  do j=1,3
      !    do k=1,3
      !      do l=1,3
      !        !Am((i-1)*3+k,(j-1)*3+l) = A(i,j,k,l)   
      !        Dm((i-1)*3+k,(j-1)*3+l) = dd(i,j,k,l)   
      !      enddo
      !    enddo
      !  enddo
      !enddo

c.... Find the determinant of the intermediate stiffness
      
      !CALL dgetrf(9,9,Am,9,ipiv,info) 
      !
      !detA = 1d0
      !f1 = 1d0
      !do i=1,9
      !   detA = detA*Am(i,i) 
      !   if (ipiv(i).ne.i) then
      !       f1 = f1*-1d0
      !   endif
      !enddo
      !
      !!detA = detA*f1
      !
      !CALL dgetrf(9,9,Dm,9,ipiv,info) 
      !
      !detD = 1d0
      !f1 = 1d0
      !do i=1,9
      !   detD = detD*Dm(i,i) 
      !   if (ipiv(i).ne.i) then
      !       f1 = f1*-1d0
      !   endif
      !enddo      
      !
      !!detD = detD*f1   
      
c.... Determine the determinant of the stiffness
 
      !if (detD.lt.0d0) then
      !   check = 1d0    
      !   str_ell_flag = .true.
      !endif      
      !
      !if (detD.lt.0d0) then
      !   check = 1d0    
      !   str_ell_flag = .true.
      !endif      