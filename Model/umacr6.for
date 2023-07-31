      subroutine umacr6(lct,ctl)

c__________________________________________________TCG. 17.11.2003___71
c
c     write history to file
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
c
c      Outputs:  serie of CHIST_xxx files
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



	integer fn, io, ivisu

      logical   pcomp, prt
      character lct*15
	character Fstring*8
      logical hflag
	integer hpoint, hlen, hpre
      real*8    ctl(3)


	integer n1, n2, n3

      save fn

c     Set file ID
        ivisu=98

      if(pcomp(uct,'mac6',4)) then      ! Usual    form
	  uct='WHIST'
c Set counter i
        OPEN(UNIT=ivisu,FILE='HIST_temp',status='OLD',IOSTAT=io)
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

        call int_str_hist11(fn,Fstring) 

c  Generating Output --> HIST 
	  write (iow,2000) Fstring 
        if(ior.lt.0) write(*,2000)  Fstring
c  Open file visu (id=10)
        OPEN(UNIT=ivisu,FILE=Fstring)

c  get pointer and length of the History     
	  call pgetd ('H  ',hpoint,hlen,hpre,hflag)
c  update 
	  if(hflag) then
	    call WHIST_11(ivisu,numel,hr(hpoint),hlen,n1,n3)
        endif

	  CLOSE(ivisu)


c increase file nb. counter
	  fn=fn+1

c write last file nb.+1 into VISU_temp
        OPEN(UNIT=ivisu,FILE='HIST_temp')
	   write (ivisu,*) fn
	  CLOSE(ivisu)


	endif


2000  format(' ** Generating Output of History --> ', a8)

1000  end
