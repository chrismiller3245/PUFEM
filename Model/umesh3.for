      
      subroutine umesh3(prt)
	
c================================================================CJM.-24/04/2020========
c    
c            FIBER DISTRIBUTION SUBROUTINE for F E A P
c
c     PURPOSE 
c
c        Input routine to read in fiber distribution data and store in the header
c        file pvecdata
c
c     INPUT FORMAT
c      
c        LINE 1 => ind, elem_ind
c     
c           ind = 1          Uniaxial Distribution
c           ind = 2          Von Mises Distribution A    (b->0 give isotropy,b->inf gives out of plane isotropy)
c           ind = 3          Von Mises Distribution B    (b->0 give isotropy,b->inf gives in plane isotropy)      
c           ind = 4          Bingham Distribution
c      
c           elem_ind = 0     All elemtents have the same vectors defining the oreintation distribution
c           elem_ind = 1     Different vectors for different elements
c      
c        LINE 2 => Parameters
c
c           Uniaxial Distribution   -  Miss out this line
c      
c           Von Mises Distribution  - b
c      
c               b      Concentration Parameter
c               
c           Bingham Distribution - k1,k2,c
c      
c               k1     Concentration Parameter
c               k2     Concentration Parameter
c               c      Normalising Constant
c
c        LINE 3 => Vectors  
c      
c           elem_ind = 0
c            
c               Uniaxial Distribution - [n(1),n(2),n(3)]
c               
c               Von Mises Distribution - [n(1),n(2),n(3)]   (Mean direction vector)
c               
c               Bingham Distribution - [P1(1),P1(2),P1(3)],[P2(1),P2(2),P2(3)],[P3(1),P3(2),P3(3)]   (Principal direction vectors)
c                                      NOTE that P1, P2, P3 must be othagonal
c      
c           elem_ind = 1
c            
c               Uniaxial Distribution - [n(1),n(2),n(3)]
c                                       [n(1),n(2),n(3)]
c                                       [n(1),n(2),n(3)] 
c                                       etc....
c               
c               Von Mises Distribution - [n(1),n(2),n(3)]
c                                        [n(1),n(2),n(3)]
c                                        [n(1),n(2),n(3)]       
c                                        etc....     
c               
c               Bingham Distribution - [P1(1),P1(2),P1(3)],[P2(1),P2(2),P2(3)],[P3(1),P3(2),P3(3)]    
c                                      [P1(1),P1(2),P1(3)],[P2(1),P2(2),P2(3)],[P3(1),P3(2),P3(3)]    
c                                      [P1(1),P1(2),P1(3)],[P2(1),P2(2),P2(3)],[P3(1),P3(2),P3(3)]
c                                      etc....
c                                       
c======================================================================================================

  
      implicit none

      include  'umac1.h'
      include  'iofile.h'
	include  'pvecdata.h'
	include  'eldata.h'  
      include  'cdata.h'  

      logical   pcomp,prt,errck,pinput
      real*8    td(13),norm,vec(4,3), fact, b, dens, fx
      real*8    A(3,5), detA, check
	integer   i,j,k
      integer   elem_ind, test

      if(pcomp(uct,'mes3',4)) then      ! Usual    form
       uct = 'FDIStribution'            ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data
       
      elseif(urest.eq.2) then           ! Write restart data
       
      else                              ! Perform user operation 
          
c-----First Line---------------------------------------------------------------          
          
      ind  = 0
      elem_Ind = 0
      
      errck = pinput(td,2)
      if(errck) then
         write(iow,*) 'Cannot read the first line'
	   goto 10
      endif
      
      ind  = td(1)
      elem_ind  = td(2)
      
      if ((ind.eq.0).OR.(ind.eq.1)) elem_ind = 0

c-----Second Line--------------------------------------------------------------       
      
      param = 0d0

      if (ind.gt.1) then
      
         errck = pinput(td,4)
         if (errck) then
            write(iow,*) 'Cannot read the second line'
            goto 10
         endif          
      
         if ((Ind.eq.2).or.(Ind.eq.3)) then    
            param(1) = td(1)
            call erfi(sqrt(2*param(1)),fx)
            param(2) = fx
         elseif (Ind.eq.4) then          
            param(1) = td(1)
            param(2) = td(2)
            param(3) = td(3)          
         endif    
      
      endif
      
c-----Remaining Lines-----------------------------------------------------------      
      
      norm = 0d0
      
      if ((ind.eq.1).or.(ind.eq.2).or.(ind.eq.3)) then !UNIAXIAL OR VON MISES FISHER A OR VON MISES FISHER B
      
         if (elem_ind.eq.0) then !ALL ELEMENTS HAVE THE SAME MEAN DIRECTION
            
            errck = pinput(td,4)
            if(errck) then
               write(iow,*) 'Cannot read the element line'
               goto 10
            endif      

            do i=1,numel                
               
               do j=1,3 
                  pvec(i,j,1) = td(j+1)
                  norm = norm + td(j+1)**2
               enddo  

               do j=1,3 
                  pvec(i,j,1) = pvec(i,j,1)/sqrt(norm)
                  qvec(i,j,1) = pvec(i,j,1)
               enddo                
               
               norm = 0d0
               
            enddo             
            
         elseif (elem_Ind.ne.0) then !EACH ELEMENT HAS ITS OWN ASSOCIATED MEAN DIRECTION
          
            do i=1,numel
            
               errck = pinput(td,4)
               if(errck) then
                  write(iow,*) 'Cannot read the element line'
	            goto 10
               endif      
             
               do j=1,3 
                  pvec(i,j,1) = td(j+1)
                  norm = norm + td(j+1)**2
               enddo  

               do j=1,3 
                  pvec(i,j,1) = pvec(i,j,1)/sqrt(norm)
                  qvec(i,j,1) = pvec(i,j,1)
               enddo                
               
               norm = 0d0
               
            enddo    
             
         endif    
      
      elseif (ind.eq.4) then !BINGHAM DISTRIBUTION

         A = 0d0 
          
         if (elem_ind.eq.0) then !ALL ELEMENTS HAVE THE SAME PRINCIPAL DIRECTIONS

            errck = pinput(td,10)
            if(errck) then
               write(iow,*) 'Cannot read the element line'
               goto 10
            endif               
            
            do i=1,numel    
             
               A = 0d0  
                
               do k=1,3
                  norm = 0d0 
                  do j=1,3 
                     A(j,k) = td((k-1)*3+1+j)
                     norm = norm + td((k-1)*3+1+j)**2
                  enddo
                  do j=1,3 
                     A(j,k) = A(j,k)/sqrt(norm)
                     pvec(i,j,k) = A(j,k)
                  enddo  
               enddo
             
               detA = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
     1               -A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))
     2                +A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

               bvec(i,1,1) = (1d0/detA)*(A(2,2)*A(3,3) - A(3,2)*A(2,3))
               bvec(i,1,2) = (1d0/detA)*(A(1,3)*A(3,2) - A(3,3)*A(1,2)) 
               bvec(i,1,3) = (1d0/detA)*(A(1,2)*A(2,3) - A(2,2)*A(1,3)) 
               bvec(i,2,1) = (1d0/detA)*(A(2,3)*A(3,1) - A(3,3)*A(2,1)) 
               bvec(i,2,2) = (1d0/detA)*(A(1,1)*A(3,3) - A(3,1)*A(1,3)) 
               bvec(i,2,3) = (1d0/detA)*(A(1,3)*A(2,1) - A(2,3)*A(1,1)) 
               bvec(i,3,1) = (1d0/detA)*(A(2,1)*A(3,2) - A(3,1)*A(2,2)) 
               bvec(i,3,2) = (1d0/detA)*(A(1,2)*A(3,1) - A(3,2)*A(1,1)) 
               bvec(i,3,3) = (1d0/detA)*(A(1,1)*A(2,2) - A(2,1)*A(1,2)) 
               
               ! AVERAGE OF FIRST AND SECOND PRINCIPAL DIRECTIONS (45 DEGREES)
               norm = 0d0
               do j=1,3 
                  A(j,4) = 0.5d0*A(j,1) + 0.5d0*A(j,2)
                  norm = norm + A(j,4)**2
               enddo
               do j=1,3 
                  A(j,4) = A(j,4)/sqrt(norm)
                  !pvec(i,j,4) = A(j,4)
               enddo  
               
               ! VECTOR PERPENDICULAR TO THE 45 DEGREE VECTOR
               A(1,5) = A(2,3)*A(3,4) - A(3,3)*A(2,4)
               A(2,5) = A(3,3)*A(1,4) - A(1,3)*A(3,4)
               A(3,5) = A(1,3)*A(2,4) - A(2,3)*A(1,4)
               norm = A(1,5)**2 + A(2,5)**2 + A(3,5)**2
               
               do j=1,3 
                  A(j,5) = A(j,5)/sqrt(norm)
                  !pvec(i,j,5) = A(j,5)
               enddo             
               
               ! POPULATE THE QVEC FOR THE DIFFERENT LOCALISATION CASES
               do j=1,3
                  qvec(i,j,1) = A(j,1) 
                  qvec(i,j,2) = A(j,2) 
                  qvec(i,j,3) = A(j,4) 
                  qvec(i,j,4) = A(j,5)
                  !qvec(i,j,1,1) = A(j,1)
                  !qvec(i,j,1,2) = A(j,2)
                  !qvec(i,j,2,1) = A(j,2)
                  !qvec(i,j,2,2) = A(j,1)
                  !qvec(i,j,3,1) = A(j,4)
                  !qvec(i,j,3,2) = A(j,5)
                  !qvec(i,j,4,1) = A(j,5)
                  !qvec(i,j,4,2) = A(j,4)
               enddo  
               
               check = 1d0
               
            enddo    
          
         elseif (elem_Ind.ne.0) then !EACH ELEMENT HAS ITS OWN ASSOCIATED SET OF PRINCIPAL DIRECTIONS
          
            do i=1,numel 
            
               errck = pinput(td,10)
               if(errck) then
                  write(iow,*) 'Cannot read the element line'
	            goto 10
               endif      
             
               A = 0d0
               
               do k=1,3
                  do j=1,3 
                     A(j,k) = td((k-1)*3+1+j)
                     norm = norm + td((k-1)*3+1+j)**2
                  enddo
                  do j=1,3 
                     A(j,k) = A(j,k)/sqrt(norm)
                     pvec(i,j,k) = A(j,k)
                  enddo  
                  norm = 0d0
               enddo
               
               detA = A(1,1)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
     1               -A(1,2)*(A(2,1)*A(3,3)-A(3,1)*A(2,3))
     1                +A(1,3)*(A(2,1)*A(3,2)-A(3,1)*A(2,2))

               bvec(i,1,1) = (1d0/detA)*(A(2,2)*A(3,3)-A(3,2)*A(2,3))
               bvec(i,1,2) = (1d0/detA)*(A(1,3)*A(3,2)-A(3,3)*A(1,2)) 
               bvec(i,1,3) = (1d0/detA)*(A(1,2)*A(2,3)-A(2,2)*A(1,3)) 
               bvec(i,2,1) = (1d0/detA)*(A(2,3)*A(3,1)-A(3,3)*A(2,1)) 
               bvec(i,2,2) = (1d0/detA)*(A(1,1)*A(3,3)-A(3,1)*A(1,3)) 
               bvec(i,2,3) = (1d0/detA)*(A(1,3)*A(2,1)-A(2,3)*A(1,1)) 
               bvec(i,3,1) = (1d0/detA)*(A(2,1)*A(3,2)-A(3,1)*A(2,2)) 
               bvec(i,3,2) = (1d0/detA)*(A(1,2)*A(3,1)-A(3,2)*A(1,1)) 
               bvec(i,3,3) = (1d0/detA)*(A(1,1)*A(2,2)-A(2,1)*A(1,2)) 
               
               ! AVERAGE OF FIRST AND SECOND PRINCIPAL DIRECTIONS (45 DEGREES)
               norm = 0d0
               do j=1,3 
                  A(j,4) = 0.5d0*A(j,1) + 0.5d0*A(j,2)
                  norm = norm + A(j,4)**2
               enddo
               do j=1,3 
                  A(j,4) = A(j,4)/sqrt(norm)
                  !pvec(i,j,4) = A(j,4)
               enddo  
               
               ! VECTOR PERPENDICULAR TO THE 45 DEGREE VECTOR
               A(1,5) = A(2,3)*A(3,4) - A(3,3)*A(2,4)
               A(2,5) = A(3,3)*A(1,4) - A(1,3)*A(3,4)
               A(3,5) = A(1,3)*A(2,4) - A(2,3)*A(1,4)
               norm = A(1,5)**2 + A(2,5)**2 + A(3,5)**2
               
               do j=1,3 
                  A(j,5) = A(j,5)/sqrt(norm)
                  !pvec(i,j,5) = A(j,5)
               enddo 
               
               ! POPULATE THE CVEC FOR THE DIFFERENT LOCALISATION CASES
               do j=1,3
                  qvec(i,j,1) = A(j,1) 
                  qvec(i,j,2) = A(j,2) 
                  qvec(i,j,3) = A(j,4) 
                  qvec(i,j,4) = A(j,5)
                  !qvec(i,j,1,1) = A(j,1)
                  !qvec(i,j,1,2) = A(j,2)
                  !qvec(i,j,2,1) = A(j,2)
                  !qvec(i,j,2,2) = A(j,1)
                  !qvec(i,j,3,1) = A(j,4)
                  !qvec(i,j,3,2) = A(j,5)
                  !qvec(i,j,4,1) = A(j,5)
                  !qvec(i,j,4,2) = A(j,4)
               enddo  
               
               check = 1d0
               
            enddo    
             
         endif          
          
      endif    
      
c------------------------------------------------------------------------------      

      write(iow,*) 

      endif

500	format(1x,i6,'-',f9.4,f9.4,f9.4,f9.4,f8.4)
10    end

      subroutine erfi(x, fx)
c===========================================T.C.G.-14.08.2008==
c   computes imaginary error function by series about x=0 
c
c   Coefficients are taken from
c   Sloane, N. J. A. Sequences A000079/M1129, A001147/M3002, 
c   and A084253 in "The On-Line Encyclopedia of Integer Sequences." 
c==============================================================
      implicit none
      
      real*8 x, fx
      
      fx = 0.5641895835477563d0*(
     &	   2.0d0*x 
     &   + 0.6666666666666667*x**3 
     &   + x**5/5.0d0 
     &   + x**7/21.0d0 
     &   + x**9/108.0d0  
     &   + x**11/660.0d0
     &   + x**13/4680.0d0 
     &   + x**15/37800.0d0 
     &   + x**17/342720.0d0 
     &   + x**19/3447360.0d0 
     &   + x**21/38102400.0d0 
     &   + x**23/459043200.0d0 
     &   + x**25/5987520000.0d0 
     &   + x**27/84064780800.0d0 
     &   + x**29/1264085222400.0d0 
     &   + x**31/20268952704000.0d0 
     &   + x**33/345226033152000.0d0 
     &   + x**35/6224529991680000.0d0
     &   + x**37/118443913555968000.0d0 
     &   + x**39/2372079457972224000.0d0)

      end
      