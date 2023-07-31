      subroutine umacr7(lct,ctl)

c__________________________________________________TCG. 01.03.2004___71
c
c     Generates Outputfile for fiber orientation
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
      include  'nLength.h'

	integer pNelem, pFi, pVe, pXbar, imesh

      logical   pcomp, prt
      character lct*15
	integer ivisu
      logical hflag
	integer hpoint, hlen, hpre
      real*8    ctl(3)
      logical IXflag
	integer IXpoint, IXlen, IXpre
      logical Xflag
	integer Xpoint, Xlen, Xpre
      logical Uflag
	integer Upoint, Ulen, Upre
      logical dflag
	integer dpoint, dlen, dpre


	integer n1, n2, n3


c     Set file ID
        ivisu=103 

      if(pcomp(uct,'mac7',4)) then      ! Usual    form
	  uct='ORIE'
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation



       ivisu= 78
       imesh= 79


c read history pointer from hdatam.h
       n1 = nhmax
	 n2 = nhmax
	 n3 = nh3max

c Define pointer to history data
	pVe        = 1
	pXbar      = 26  
      pFi        = 44  
	pNelem     = 53  

c  get pointer and length of the History     
	  call pgetd ('H  ',hpoint,hlen,hpre,hflag)
c  get pointer and length of the X and U    
	  call pgetd ('X ',xpoint,xlen,xpre,xflag)
c  get pointer and length of the IX     
	  call pgetd ('IX ',ixpoint,ixlen,ixpre,ixflag)
c  get pointer and length of U    
	 call pgetd ('U ',Upoint,Ulen,Upre,Uflag)
c  get pointer and length of the D     
	 call pgetd ('D ',dpoint,dlen,dpre,dflag)

c compute neighbouring elements
        if (chl.ne.0.0d0) then
	    write(*,2001) 
	    write(iow,2001) 
	    write(*,*) 'neighbour_11 nees to be edit for orie!!!!!'
	    write(iow,*) 'neighbour_11 nees to be edit for orie!!!!!'

          call neighbour_11(hr(dpoint),dlen,n1,n3,nen1,pXbar,pNelem,
	1                  pVe,hr(hpoint),hlen,mr(ixpoint),ixlen)
	  endif

        OPEN(UNIT=ivisu,FILE='Orientation')
        OPEN(UNIT=imesh,FILE='fibers')
        
	  write(*,2000)
	  write(iow,2000)

	  call output_11(ivisu,imesh,pFi,pNelem,pVe,ndm,ndf,numel,nen,
	1               nen1,n1,n3,hr(hpoint),hlen,mr(ixpoint),ixlen,
     2                 hr(xpoint),xlen,hr(upoint),ulen)


	  CLOSE(imesh)
	  CLOSE(ivisu)


	endif


2000  format(' ** Generating Output for Fiberorientation')
2001  format(' ** Compute neighbour elements......')

1000  end
