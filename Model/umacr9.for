      subroutine umacr9(lct,ctl)

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
	integer ivisu
      logical IXflag
	integer IXpoint, IXlen, IXpre
      logical Xflag
	integer Xpoint, Xlen, Xpre
      logical Uflag
	integer Upoint, Ulen
      real*8    ctl(3)


	integer n1, n2, n3

      save fn

c     Set file ID
        ivisu=99 

      if(pcomp(uct,'mac9',4)) then      ! Usual    form
	  uct='VISU'
c Set counter i
        OPEN(UNIT=ivisu,FILE='VISU_temp',status='OLD',IOSTAT=io)
	  IF(io.eq.0) THEN
	    read (ivisu,*) fn
	    CLOSE(ivisu)
	  ELSE
          fn=100
        ENDIF 
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation

c read history pointer from hdatam.h
       n1 = nhmax
	 n2 = nhmax
	 n3 = nh3max



c  Generating Output for Visualisation coord 
       Fstring = "Coord"
	 write (iow,2000) Fstring 
       if(ior.lt.0) write(*,2000)  Fstring
c  Open file 
       OPEN(UNIT=ivisu,FILE=Fstring)

c  get pointer and length of the X    
	 call pgetd ('X ',xpoint,xlen,xpre,xflag)
c  print x
	 if(xflag) then
	   call write_nodes11(ivisu,ndm,numnp,hr(xpoint),xlen)
	 endif

	  CLOSE(ivisu)

c  Generating Output for Visualisation connect
       Fstring = "Connect"
	 write (iow,2001) Fstring 
       if(ior.lt.0) write(*,2001)  Fstring
c  Open file 
       OPEN(UNIT=ivisu,FILE=Fstring)

c  get pointer and length of the IX     
	 call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)
c  print ix

	 if(ixflag) then
	  call write_connect11 (ivisu,nen1,nen,numel,mr(ixpoint),ixlen)
	 endif

	 CLOSE(ivisu)


c Determine Filename 
       call int_str_disp11 (fn,Fstring) 

c  Generating Output for Visualisation displacements 
	 write (iow,2002) Fstring 
       if(ior.lt.0) write(*,2002)  Fstring
c  Open file 
       OPEN(UNIT=ivisu,FILE=Fstring)

c  get pointer and length of the X    
	 call pgetd ('U ',upoint,ulen,xpre,uflag)
c  print u
	 if(xflag) then
	   call write_displacements11(ivisu,ndm,ndf,numnp,
	1                               hr(upoint),ulen)
	 endif
	 CLOSE(ivisu)


c Determine Filename 
       call int_str_stre11 (fn,Fstring) 


c  Generating Output for Visualisation stresses 
	 write (iow,2003) Fstring 
       if(ior.lt.0) write(*,2003)  Fstring
c  Open file 
       OPEN(UNIT=ivisu,FILE=Fstring)

c print fluxes
       call write_stress11 (ivisu,Hr(np(58)+numnp),numnp)

	 CLOSE(ivisu)
c increase file nb. counter
	 fn=fn+1

c write last file nb.+1 into VISU_temp
       OPEN(UNIT=ivisu,FILE='visu_temp')
	   write (ivisu,*) fn
	 CLOSE(ivisu)
      endif



2000  format(' ** Writing Reference coordinates in --> ', a8)
2001  format(' ** Writing Connectivity in --> ', a8)
2002  format(' ** Writing Displacements in --> ', a8)
2003  format(' ** Writing stress in --> ', a8)
      end
