      subroutine umacr1(lct,ctl)

c__________________________________________________TCG. 12.05.2003___71
c
c     Generates Outputfile for plotting curves
c
c      Inputs:
c         lct       - Command character parameters
c         ctl(3)    - Command numerical parameters
c         prt       - Flag, output if true
c
c      Outputs:  plot.txt file
c
c
c   NOTE:  Flux vector Hr(np(58)+numnp) must be updated by calling a 
c          stressplot-command befor calling VISU  
c
c
c  Declare variable types
c     
c     nr          Nodes of needed reactions
c     nd          Nodes of needed displ
c     rflag       Propper reading flag
c     icurv       reading UNIT number
c     nreac       Number of nodes	of needed reactions
c     ndispl      Number of nodes	of needed displ
c     Upoint      Poiter of U-Array (Displacements)
c     Ulen        Length of U-Array 
c     Upre        Precision of U-Array 
c     Uflag       Flag of U-Array 
c     DRpoint     Poiter of DR-Array (Displacements)
c     DRlen       Length of DR-Array 
c     DRpre       Precision of DR-Array 
c     DRflag      Flag of DR-Array 
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

	integer i

      character lct*15
      logical pcomp, prt
	integer icurv
	integer nr(1000), nd(5), nreac, ndispl
	logical rflag
      logical Uflag
	integer Upoint, Ulen, Upre
      logical DRflag
	integer DRpoint, DRlen, DRpre
      real*8    ctl(3)

      save  nr, nd, nreac, ndispl

c     Set file ID
        icurv=101 

      if(pcomp(uct,'mac1',4)) then      ! Usual    form
	  uct='CURV'
c load plott data
	  call rcurv (icurv, nr,  nreac, ndispl, nd, rflag)
	  IF(rflag) THEN
	    write(*,*) 'PLOT DATA loaded'
          OPEN(UNIT=icurv,FILE='CURVE.inf')
	      write(icurv,*) '----------------------------'
	    write(icurv,*) 'Plot Reaction Force of nodes'
	    do I= 1, nreac
	      write(icurv,*) nr(i)
	    enddo
	    write(icurv,*) 'versus Displacement of nodes'
	    do I= 1, ndispl
	      write(icurv,*) nd(i)
	    enddo
	      write(icurv,*) '----------------------------'
	      write(icurv,*) 
	      write(icurv,*) '!! Call REAC bevore calling CURVe !!'
	      write(icurv,*) 
	    Close(icurv)
	  ELSE
          Write (*,*) 'WARNING -  NO PLOT DATA loaded!!'
        ENDIF 
      elseif(urest.eq.1) then           ! Read  restart data

      elseif(urest.eq.2) then           ! Write restart data

      else                              ! Perform user operation


c     Set command word

c  get pointer and length of U    
	 call pgetd ('U ',Upoint,Ulen,Upre,Uflag)
c  get pointer and length of DR    
	 call pgetd ('DR ',DRpoint,DRlen,DRpre,DRflag)


       OPEN(UNIT=icurv,FILE='CURVE.txt')
c	call printdata (icurv, numnp, ndm, Ulen,hr(Upoint),DRlen,
c	1                hr(DRpoint), nreac,nr,nd)

	 call printdata (icurv, numnp, ndf, Ulen,hr(Upoint),DRlen,
	1                hr(DRpoint), nreac,nr, ndispl, nd)

c	 CLOSE(icurv)
      endif

1000  end
