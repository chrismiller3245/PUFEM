      subroutine umacr5(lct,ctl)

c__________________________________________________TCG. 05.06.2003___71
c
c     updates history
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
c
c
c
c
c
c  Declare variable types
c     
c     IXpoint     Poiter of IX-Array (Element connections)
c     IXlen       Length of IX-Array 
c     IXpre       Precision of IX-Array 
c     IXflag      Flag of IX-Array 
c     Hpoint      Poiter of History
c     Hlen        Length of History 
c     Hpre        Precision of History 
c     Hflag       Flag of History 
c____________________________________________________________________71
      implicit  none

      include  'iofile.h'
      include  'umac1.h'
	include  'comblk.h'
	include  'pointer.h'
	include  'cdata.h'   
	include  'sdata.h' 
      include  'prstrs.h'
      include  'eldata.h'
	include  'hdatam.h'


      real*8             hist0(100000),hist1(100000)
      common /hist_tcg/  hist0, hist1



      logical   pcomp, prt
      character lct*15
      logical hflag
	integer hpoint, hlen, hpre
      logical ixflag
	integer ixpoint, ixlen, ixpre
      real*8    ctl(3)


	integer n1, n2, n3


      if(pcomp(uct,'mac5',4)) then      ! Usual    form
	  uct='UPDATE'
	  hist0=0.0d0
	  hist1=0.0d0
	write(*,*) 'history initialized'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation


c read history pointer from hdatam.h
       n1 = nhmax
	 n2 = nhmax
	 n3 = nh3max


c  get pointer and length of the History     
	  call pgetd ('H  ',hpoint,hlen,hpre,hflag)
c  get pointer and length of the IX     
	  call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)
	  write(*,2000)
	  write(iow,2000)
c  update 
	 if(hflag) then
	  call update_11(numel,n1,n3,nen1,
	1                 hr(hpoint),hlen,mr(ixpoint),ixlen)
       endif


       hist0=hist1


	endif


2000  format(' ** History vector updated')

1000  end
