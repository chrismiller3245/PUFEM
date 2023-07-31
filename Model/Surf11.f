      subroutine surf_11(d,xl,ul,ma,ndf,ndm,nel,mct,nst, p,s)
c
c--------------------------------------------------------------------71
c
c         Nodal force and tangent array for pressure loading
c            
c
c....  INPUT variables
c        d(*)       Material parameters
c        xl(3,*)    Nodal coordinates
c        ul(3,*)    Nodal displacements
c        ma         Option number (reference or current conf.)     
c        ndf        Number of DOF / node
c        ndm        Space dimension
c        nel        Number of nodes / Element
c        mct        Type of surface loading
c        nst        Dimension of residual vector
c    
c....  OUTPUT variables
c        p(nst)     Contribution to the residual
c        s(nst,nst) Contribution to the stiffness matrix
c
c
c....  PARAMATER set up:
c
c         ma   =  1  (for 3-d analysis in reference coordinates)
c              = -1  (for 3-d analysis in current coordinates;
c                     an unsymmetric tangent matrix also computed.)
c
c         nel  =  3  (for 3-node linear triangle)
c
c
c                        ^ eta
c                        | 
c  
c                        2                
c                        o
c                        | -     
c                        |   -    
c                        |     -  
c                        |       -
c                        |         -      
c                        |           -     
c                        |             -   
c         Nodes are:     o---------------o   ->xi
c                        3                1
c
c         mct  = 1   (for constant pressure on the face)
c                d(1) - value of constant pressure on the face
c
c              = 3   (for variable pressure on the face)
c                d(1) - value of pressure at node-1
c                d(2) - value of pressure at node-2
c                d(3) - value of pressure at node-3
c
c.... MACRO instruction for input in FEAP
c
c        sloa
c        iel, 3, 3, ma(1=ref. coord, 2=current coord)
c        node1,node2,node3,d(1),d(2),d(3)
c
c--------------------------------------------------------------------71
c
c....  Declare variable types
      integer lint, ma, nel, ndf, ndm, mct, nst, l, i, j, k, nn
      integer ii, jj, i1, j1
      real*8  dx(3,2), pn, pp, fact
c....  Declare array types
      real*8  sg, tg, wg, shp(3,3), xl(ndm,*), xu(3,3), d(*)
      real*8  ul(ndf,*), p(ndf,*), s(nst,*)
c
c....  Compute nodal coordinates in correct reference frame
c
	s(1,1)=0.0d0
c
      if(ma.gt.0) then
        fact = 0.d0
      else
	     fact = 1.d0
      endif
	do 100 nn = 1,3
		xu(1,nn) = xl(1,nn) + fact*ul(1,nn)
		xu(2,nn) = xl(2,nn) + fact*ul(2,nn)
		xu(3,nn) = xl(3,nn) + fact*ul(3,nn)
100   continue

c....  Get quadrature information
c
      lint  = 1
	sg    = 0.3333333333333333d0
	tg    = 0.3333333333333333d0
	wg    = 0.5d0
c
c....  First loop over quadrature points
c
      do 200 l = 1,lint
c
c....  Compute geometric factors
        call shaptriag (sg, tg, shp)
        call pzero (dx, 6)
        do 110 nn = 1,nel
          do 105 i = 1,3
            dx(i,1) = dx(i,1) + shp(1,nn)*xu(i,nn)
            dx(i,2) = dx(i,2) + shp(2,nn)*xu(i,nn)
105       continue
110     continue
c
c
c....  Compute pressure on face point
          if(mct.eq.4) then
            pn = 0.0d0
            do 120 nn = 1,nel
              pn = pn + shp(3,nn)*d(nn)
120         continue
          else
            pn = d(1)
          endif
          pn = pn*wg

c
c....  Compute nodal loads for pressures
          do 130 nn = 1,nel
            pp = shp(3,nn)*pn
            do 125 i = 1,3
              j = mod(i,3) + 1
              k = mod(j,3) + 1
              p(i,nn) = p(i,nn)
     1                + pp*(dx(j,1)*dx(k,2) - dx(k,1)*dx(j,2))
125         continue
130       continue
          

c....  Compute a tangent if necessary

          if(ma.lt.0) then
            i1 = 0
            do 150 ii = 1,nel
              pp = shp(3,ii)*pn
              j1 = 0
              do 140 jj = 1,nel
                do 135 i = 1,3
                  j = mod(i,3) + 1
                  k = mod(j,3) + 1
				s(i1+i,j1+j) = s(i1+i,j1+j) 
	1						- pp*(shp(1,jj)*dx(k,2) 
	2						- dx(k,1)*shp(2,jj))
				s(i1+i,j1+k) = s(i1+i,j1+k)
	1						- pp*(shp(2,jj)*dx(j,1) 
	2						- dx(j,2)*shp(1,jj))
135             continue	  
                j1 = j1 + ndf
140           continue
              i1 = i1 + ndf
150         continue
          endif

200   continue   !Gauss-loop
      end



      subroutine shaptriag (xi, eta, shp)
	implicit none

	real*8 xi, eta, shp(3,3)

c derivatives
	shp(1,1)=  1.0d0
	shp(1,2)=  0.0d0
	shp(1,3)= -1.0d0

	shp(2,1)=  0.0d0
	shp(2,2)=  1.0d0
	shp(2,3)= -1.0d0

c shapefunction
	shp(3,1)=  xi
	shp(3,2)=  eta
	shp(3,3)= 1.0d0 - xi- eta


      end