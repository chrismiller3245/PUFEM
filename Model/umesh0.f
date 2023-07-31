c$Id: umesh0.f,v 1.2 2002/02/19 21:40:41 rlt Exp $
      subroutine umesh0(prt)

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose: read non-local measure for quantity smoothing

c      Inputs:
c         prt    - Flag, output results if true

c      Outputs:
c               store non-local measure in slength in
c               common block non_local
c-----[--.----+----.----+----.-----------------------------------------]
      implicit  none

      include  'umac1.h'
      include  'nLength.h'
      include  'iofile.h'

      logical   prt,pcomp,errck,pinput
      real*8    td(13)
      !character tx(*)*15

c     Set command

      if(pcomp(uct,'mes0',4)) then      ! Usual    form
       uct = 'NLENGTH'                     ! Specify 'name'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation
10      errck = pinput(td,5)
          if(errck) then
            write(iow,*) 'cannot read non-local data!'
	      goto 10
	    endif
        ichl  = int(td(1))
        chl  = td(2)
        nfit  = int(td(3))
        chr  = td(4)
        acrack_ini= td(5)
	  if ((ichl.ne.1).and.(ichl.ne.2)) then
	    write(*,1001) 
	    write(iow,1001) 
	  stop
	  endif
	  if ((nfit.ne.1).and.(nfit.ne.2)) then
	    write(*,1002) 
	    write(iow,1002) 
	  stop
	  endif
	  if (acrack_ini.lt.0) then
	    write(*,1003) 
	    write(iow,1003) 
	  stop
	  endif
	  write(iow,1000) ichl, chl, nfit, chr, acrack_ini
      endif

1000  format(
     1 5x, '------------------------------------------------------'/
     2  10x,'Non-local measures for quantity smoothing'/
     2  10x,'Tye:                    ',i5/
     3  10x,'1=abs. length,  2=rel. length'/
     4  10x,'chl:                    ',e12.5/
     4  10x,'                         '/
     2  10x,'Order of local surface  ',i5/
     3  10x,'1=linear,  2=quadratic'/
     4  10x,'chr:                    ',e12.5/
     4  10x,'acrack_ini:             ',e12.5/
     5 5x, '------------------------------------------------------')
1001  format('** ERROR: Unknown type for non-local measure!!')
1002  format('** ERROR: Unknown surface type !!')
1003  format('** ERROR: Negative non-local measure !!')
      end
