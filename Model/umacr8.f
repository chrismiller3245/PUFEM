      subroutine umacr8(lct,ctl)

c      * * F E A P * * A Finite Element Analysis Program

c....  Copyright (c) 1984-2002: Regents of the University of California
c                               All rights reserved

c-----[--.----+----.----+----.-----------------------------------------]
c      Purpose:  User interface for adding solution command language
c                instructions.

c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true

c      Outputs:
c         N.B.  Users are responsible for command actions.  See
c               programmers manual for example.
c-----[--.----+----.----+----.-----------------------------------------]
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



	integer fn, io

      logical   pcomp, prt
      character lct*15
	character Fstring*8
      logical hflag
	integer hpoint, hlen, hpre
      logical xflag
	integer xpoint, xlen, xpre
      logical uflag
	integer upoint, ulen, upre
	integer nn
      real*8    ctl(3),pFres, pFtot


	integer n1, n2, n3

      save fn

      if(pcomp(uct,'mac8',4)) then      ! Usual    form
	  uct='CREF'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation



c read history pointer from hdatam.h
       n1 = nhmax
	 n2 = nhmax
	 n3 = nh3max

c pointer of n3 array
       nn = 2*n1 + n3



c  get pointer and length of the History     
	 call pgetd ('H  ',hpoint,hlen,hpre,hflag)
	 call pgetd ('X  ',xpoint,xlen,xpre,xflag)
	 call pgetd ('U  ',upoint,ulen,upre,uflag)



c      pFtot=44
      pFtot=38
      pFres=360

c pointer of n3 array
      nn = 2*n1 + n3

	if((hflag)) then
	  write(*,2000) 
c	  call Change_ref_config_11(nn,n1,numel,pFtot,pFres, 
c	1                               hr(hpoint),hlen)

        OPEN(UNIT=73,FILE='coor')



	  call Change_ref_config_11(ndm,numnp,
	1                 hr(xpoint),xlen,hr(upoint),ulen)

        close(73)

      else
	  write(*,1000) 
      endif


      endif

1000  format(' ** ERROR in umac8')
2000  format(' ** Change Reference Coordinates')
      end
