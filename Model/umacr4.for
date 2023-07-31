      subroutine umacr4(lct,ctl)

c__________________________________________________TCG. 12.05.2003___71
c
c     Generates Outputfile for Visualisation purposes by calling the 
c     userdefined FEAP-command  CTET
c
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
c
c      Outputs:  serie of CTET_xxx files
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

	integer fn, io

      logical   pcomp, prt
      character lct*15
	character Fstring*8
	integer ivisu
      logical IXflag
	integer IXpoint, IXlen, IXpre
      logical Xflag
	integer Xpoint, Xlen, Xpre
      logical hflag
	integer hpoint, hlen, hpre
      real*8    ctl(3)


	integer n1, n2, n3

      save fn

c     Set file ID
        ivisu=103 

      if(pcomp(uct,'mac4',4)) then      ! Usual    form
	  uct='CTET'
c Set counter i
        OPEN(UNIT=ivisu,FILE='CTET_temp',status='OLD',IOSTAT=io)
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





c Determine Filename 

       call int_str_ctet10 (fn,Fstring) 

c  Generating Output for Visualisation --> CTET 
       if(.not.pcomp(uct,'mac4',4)) then
	  write (iow,2000) Fstring 
        if(ior.lt.0) write(*,2000)  Fstring
c  Open file visu (id=10)
        OPEN(UNIT=ivisu,FILE=Fstring)
c  get pointer and length of the X and U    
	  call pgetd ('X ',xpoint,xlen,xpre,xflag)
c  get pointer and length of the IX     
	  call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)
c  get pointer and length of the History     
	  call pgetd ('H  ',hpoint,hlen,hpre,hflag)


c  print 
	 if(xflag.and.ixflag.and.hflag) then
	  call write_ctet10(ivisu,ndm,numel,nen,nen1,hr(xpoint),xlen,
	1                  mr(ixpoint),ixlen, hr(hpoint),hlen,n1,n3)
       endif

	 CLOSE(ivisu)

c increase file nb. counter
	 fn=fn+1

c write last file nb.+1 into VISU_temp
       OPEN(UNIT=ivisu,FILE='CTET_temp')
	   write (ivisu,*) fn
	 CLOSE(ivisu)

	 endif
	endif


2000  format(' ** Generating Output for Visualisation --> ', a8)
2010  format('gmvinput ascii')
2020  format('gmvend')

1000  end
