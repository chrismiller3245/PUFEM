      subroutine umacr2(lct,ctl)

c__________________________________________________TCG. 12.05.2003___71
c
c     Generates Outputfile for Visualisation purposes by calling the 
c     userdefined FEAP-command  CRACk
c
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
c
c      Outputs:  serie of visu_xxx files
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
c     cTets       test with violated crack criterion
c     ntets       number of test with violated crack criterion 
c     c_Tip       Crack Tip (Collection of facects)
c     nc_Tip      No. of facets forming the Crack Tip
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

      include  'cGeom.h'
      include  'iface.h'


	integer fn, io

      logical   pcomp, prt
      character lct*15
	character Fstring*8
	integer ivisu
      real*8    ctl(3)
      logical IXflag
	integer IXpoint, IXlen, IXpre
      logical Xflag
	integer Xpoint, Xlen, Xpre
      logical Uflag
	integer Upoint, Ulen
      logical hflag
	integer hpoint, hlen, hpre


	integer n1, n2, n3, i,j

      save fn

c     Set file ID
        ivisu=102 

      if(pcomp(uct,'mac2',4)) then      ! Usual    form
	  uct='CRAC'
c Set counter i
        OPEN(UNIT=ivisu,FILE='CRACK_temp',status='OLD',IOSTAT=io)
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

       call int_str10 (fn,Fstring) 

c  Generating Output for Visualisation --> CRACK 
       if(.not.pcomp(uct,'mac2',4)) then
	  write (iow,2000) Fstring 
        if(ior.lt.0) write(*,2000)  Fstring
c  Open file visu (id=10)
        OPEN(UNIT=ivisu,FILE=Fstring)

c   write data into CRACK_xxx.txt
	   write(ivisu,100)  (ncorn(j),(cnode(j,i),i=1,12),j=1,nsurf) 

c add cracked interfaces

c  get pointer and length of the X and U    
	  call pgetd ('X ',xpoint,xlen,xpre,xflag)
c  get pointer and length of the IX     
	  call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)
c  get pointer and length of the History     
	  call pgetd ('H  ',hpoint,hlen,hpre,hflag)


c  print 
	 if(xflag.and.ixflag.and.hflag) then
c	  call write_inerface_cracks11(ivisu,ndm,numel,nen,nen1,
c	1   hr(xpoint), xlen, mr(ixpoint),ixlen, hr(hpoint),hlen,n1,n3)
       endif

	  CLOSE(ivisu)

c increase file nb. counter
	 fn=fn+1

c write last file nb.+1 into VISU_temp
       OPEN(UNIT=ivisu,FILE='Crack_temp')
	   write (ivisu,*) fn
	 CLOSE(ivisu)

	 endif


      endif

100   format(i6,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5,e13.5
     1	   ,e13.5,e13.5,e13.5,e13.5,e13.5)
2000  format(' ** Generating Output for Visualisation --> ', a8)
2010  format('gmvinput ascii')
2020  format('gmvend')

1000  end
