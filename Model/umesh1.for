      subroutine umesh1(prt)
	
c================================================================T.C.G.-05.09.1997==
c   input routine to read the direction vectors and store in 
c   commom block /vecdata/
c   
c    RECORD 1.  first element, last element, nfib, step
c    RECORD 2.  if (step.EQ.0)   ...element specific generation
c
c                  node, a1, a2, a3, b1, b2, b3, c1, c2 ...... 
c
c           	  else			   ...automatical generation
c
c                  a1, a2, a3, b1, b2, b3	,c1, c2, ......
c
c               endif
c
c==============================================================

  
      implicit none

      include  'umac1.h'
      include  'iofile.h'
	include  'vecdata.h'

      logical   pcomp,prt,errck,pinput
      real*8    td(13),norm,vec(4,3)
	integer   i,von, bis, schw, nfib ,j, k
	!character tx(*)*15


      if(pcomp(uct,'mes1',4)) then      ! Usual    form
       uct = 'FIBRE'                     ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data
       
      elseif(urest.eq.2) then           ! Write restart data
       
      else                              ! Perform user operation
	 von =0
	 bis =0
	 schw=0
10     errck = pinput(td,4)
       if(errck) then
        write(iow,*) 'cannot read number of fibre elements!'
	  goto 10
	 endif
	 von  = td(1)
	 bis  = td(2)
	 nfib = td(3)
	 schw = td(4)

	if (von.NE.0) then
	write(iow,*) 'Fibre-directionvectors in each element as follows'
	write(iow,*) '---------------------------------------------------'
	write(iow,*) 'Elem x   x1-fib 1   x1-fib 2	  x1-fib 3   x1-fib 4'
	write(iow,*) '         x2-fib 1   x2-fib 2	  x2-fib 3   x2-fib 4'
	write(iow,*) '         x3-fib 1   x3-fib 2	  x3-fib 3   x3-fib 4'
	write(iow,*) '---------------------------------------------------'
	 if (schw.EQ.0) then
        do i=von ,bis
	   errck = pinput(td,1+3*nfib)
	   if (errck) then
	     write(iow,*) 'cannot read componets of direction vectors!'
	     goto 1000
	   endif
	   do j=1,nfib
	     norm = 0.0d0
	     do k=1,3
	       rvec(j,k,i) =td(1+k+(j-1)*3)
	       norm = norm + rvec(j,k,i)*rvec(j,k,i)
	     enddo
c normalize
           norm = SQRT(norm)
           do k= 1,3
	       rvec(j,k,i) = rvec(j,k,i)/norm
		 enddo
	   enddo
	  write (iow,500)i,rvec(1,1,i),rvec(2,1,i),rvec(3,1,i),rvec(4,1,i)
	  write (iow,501)  rvec(1,2,i),rvec(2,2,i),rvec(3,2,i),rvec(4,2,i)
	  write (iow,501)  rvec(1,3,i),rvec(2,3,i),rvec(3,3,i),rvec(4,3,i)
        write(iow,*) 
	 enddo		    
	 else
	   errck = pinput(td,3*nfib)
	   if (errck) then
	     write(iow,*) 'cannot read componets of direction vectors!'
	     goto 1000
	   endif
	   do j=1,nfib
	     norm = 0.0d0
	     do k=1,3
	       vec(j,k) =td(k+(j-1)*3)
	       norm = norm + vec(j,k)*vec(j,k)
	     enddo
c normalize
           norm = SQRT(norm)
           do k= 1,3
	       vec(j,k) = vec(j,k)/norm
		 enddo
	   enddo
	   do i= von,bis,schw
	   do j=1,nfib
	     do k=1,3
	      rvec(j,k,i) =vec(j,k)
	    enddo
	  enddo
	  write (iow,500)i,rvec(1,1,i),rvec(2,1,i),rvec(3,1,i),rvec(4,1,i)
	  write (iow,501)  rvec(1,2,i),rvec(2,2,i),rvec(3,2,i),rvec(4,2,i)
	  write (iow,501)  rvec(1,3,i),rvec(2,3,i),rvec(3,3,i),rvec(4,3,i)
	  write(iow,*) 
	  enddo
	  endif
	 endif

      endif



500	format(1x,'Elem',i5,'-',f10.4,f10.4,f10.4,f10.4)
501   format(11x, f10.4, f10.4, f10.4, f10.4)
1000  end
