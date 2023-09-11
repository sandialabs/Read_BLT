!*---------------------------------------------------------------------------------*
!* Read_BLT.DLL  - Copyright (2005) Sandia Corporation. Under the                  *
!* terms of Contract DE-AC04-94AL85000, there is a non-exclusive license for use of* 
!* this work by or on behalf of the U.S. Government. Export of this program may    *
!* require a license from the United States Government.                            * 
!*                                                                                 *
!***********************************************************************************!   
!*  DLL software for the integration of the BLT-MS (Breach Leach Transport 
!*  -Multi Species)code with GoldSim. 
!*
!*  The program is intended as a DLL to read a standard BLT-MS input file
!*  and allocate parameters to memory for subsequent use in GoldSim and in
!*  writing out input files for the iterative process of performing Monte
!*  Carlo analyses to quantify uncertainty.
!*
!*------------------------------------------------------------------------------*
!*
!*  Some routines ReadDeck() and MakeDeck() this code are taken from the 
!*  preprocessor for BLT-MS developed by:
!*  Terry Sullivan
!*  BNL Department of Advanced Technology
!*  
!*  This routine was written by 
!*  Patrick Mattie
!*  Sandia National Laboratories 
!*
!*  with assistance from
!*  Bob Knowlton
!*  Sandia National Laboratories 
!*     
!*  
!*------------------------------------------------------------------------------*
!*
!*
!*****************************************************************************
! Start of subroutines for GoldSim
!*****************************************************************************
      subroutine read_blt(method, state, in, out)
!
!
! attibute statements needed for the dll, per GoldSim specifications
!
!DEC$ ATTRIBUTES dllexport,c :: read_blt
!DEC$ ATTRIBUTES alais : "read_blt" : "READ_BLT"
!DEC$ ATTRIBUTES value :: method
!DEC$ ATTRIBUTES reference :: state
!DEC$ ATTRIBUTES reference :: in
!DEC$ ATTRIBUTES reference :: out



! ******************** Declare Statements ***********************************


integer(4) method, state, file_index,i,j,k,l, h, array_index, max_in
real(8) in(*), out(*)
integer,dimension (35) :: uncert_switch, array_size


!      PROGRAM read_blt
!*
      character title*79, cnam*10, inchar, in_file*16
!*
      common /data1/ nprob, title
      common /data2/ kvi, kss
      common /data3/ niso, ncon, nctype, nwtype, idiff, iact
      common /data4/ ilump, imid, iwet, ioptim, ivml
      common /data5/ nti, ndtchg
      common /data6/ kmesh, nnp, nel, nmat, ncm
      common /data7A/ ndnp(10), ndpr(10), nddp(10)
      common /data7B/ ncnp(10), nces(10), ncpr(10), ncdp(10)
      common /data7C/ nnnp(10), nnes(10), nnpr(10), nndp(10)
      common /data8/ nsel(10),nspr(10),nsdp(10),nwnp(10),nwpr(10)&
              ,nwdp(10)
      common /data9/  cnam(10), decayc(10), csat(10), atm(10)
      common /data9a/ ichain, lchain(10), idchain(10,10)
      common /data9b/ bfr(10,10)
      common /data10/ delt, chng, delmax, tmax, w, wv,tdtch(100)
      common /data12/ ni(2,600),nseq(2,600),nad(2,600),&
             xni(600),xad(600),xrd(600),&
             zni(600),zad(600),zrd(600),&
             nelx,nelz,nreg,&
             imin(50),jmin(50),imax(50),jmax(50),&
             x1(4,50),x2(4,50),mi(550),ie(550,5),&
             nlay(550), modl(550),widfac
      common /data13/ prop(8,5),prop1(8,5,10)
      common /data14/ mi14(600),nseq14(600),mad(600),mityp(600),&
             mtypad(600)
      common /data15/ ni15(10,600),nnode(10,600),nad15(10,600),&
             rni(10,600),rad(10,600),rrd(10,600)
      common /data16/ cvbf(10,10,6),tcvbf(10,10,6),mi16(10,600),&
             nseq16(10,600),miad16(10,600),mityp16(10,600),&
             mtypad16(10,600),npcb(600,10),mi_v16(10,600),&
             nseq_v16(10,600),m_v16(10,600),is1_v16(10,600),&
             is2_v16(10,600),miad_v16(10,600),mad_v16(10,600),&
             is1ad(10,600),is2ad(10,600)
      common /data16n/ cvbfn(10,10,6),tcvbfn(10,10,6),mi16n(10,600),&
             nseq16n(10,600),miad16n(10,600),mityp16n(10,600),&
             mtpad16n(10,600), npnb(600,10),&
             mi_v16n(10,600),nseq_16n(10,600),m_v16n(10,600),&
             is1_v16n(10,600),is2_v16n(10,600),miad_16n(10,600),&
             mad_v16n(10,600),is1adn(10,600),is2adn(10,600)
      common /data17/ cdbf(10,10,6),tcdbf(10,10,6),ni17(10,600),&
             nseq17(10,600),niad17(10,600),nityp17(10,600),&
             ntypad17(10,600),npdb(600,10)
      common /data18/ sosf(10,10,6),tsosf(10,10,6),ni18(10,600),&
             nseq18(10,600),niad18(10,600),nityp18(10,600),&
             ntypad18(10,600),les(600,10)
      common /data19/ wssf(10,10,6),twssf(10,10,6),ni19(10,600),&
             nseq19(10,600),niad19(10,600),nityp19(10,600),&
             ntypad19(10,600),npw(600,10)
      common /data20/ thick(20),pitn(20),pitk(20),area(20),ascale(20),&
             pits(20),grate(20),clay(20),sph(20),iaer(20),&
             nelcon(100),mi20(600),nseq20(600),mad20(600),&
             mityp20(600),mtypad20(600)
      common /data21/ sfract(20,10),pfract(20,10),bfract(20,10),&
             deff(20,10),disol(20,10),partk0(20,10),&
             partki(20,10),partd(20,10),porel(20),volwf(20),&
             vratio(20),mi21(600),nseq21(600),mad21(600),&
             mityp21(600),mtypad21(600),wtinit(100,10)
      common /data22/ kpr0, kpr(1000),ntrc,ntrf,nstptr,lctrc(20),&
             lftrc(20)
      common /data23/ ni23(600),nseq23(600),nad23(600),vxni23(600),&
             vzni23(600),vxad23(600),vzad23(600),ni23b(600),&
             nseq23b(600),nad23b(600),thni23(600),&
             thniad23(600)
!*
!* Initialize default values
      cnam(1) = 'xxxxx'
      nstptr  = 1
      niso = 1
      nmat = 1
      nnp = 1
      nel = 1
      w = 1.0
      wv = 1.0
      MI(1) = 1
      title = '  (none)'
      in_file = 'bltmsin.inp'

! ************* Goldsim DLL interface commands **************************

if (method.eq.0) then				!Initialize
	continue

elseif (method.eq.2) then			!Report Version
	out(1) = 1.0



elseif (method.eq.3) then   ! Confirm input/outputs to Goldsim
	max_in =2249
	out(1) =max_in          ! Pass to Goldsim the maximum expected in () array size
	out(2) =5

! *************  BLT-MS Input Files **************************

elseif (method.eq.1) then	
		
!                        ***********************
!                        * Begining of Main    *
!                        ***********************
!--------------------------------------------------------------------------------------

!  open the input file with a generic name - bltmsin.inp
      open (2,FILE=in_file,ACCESS='sequential',FORM='formatted',&
          STATUS='old')

!***************************************************************************
!  read in the BLT-MS input file and store parameters/variables in memory  *
!  for writing to input files for the probabilistic analysis               *
!***************************************************************************

         call ReadDeck()

		 close(2)

!**************************************************************************
!* Insert Uncertain Samples from GoldSim                                  *
!*                                                                        *
!* This subroutine is designed to insert uncertain paramters              *
!* in an already existing input deck                                      *
!*------------------------------------------------------------------------* 
!* Patrick D. Mattie - Sandia National Labortories                        *
!* October 2005                                                           *
!**************************************************************************
      open (3,FILE='RBLTMS_debug.txt',ACCESS='sequential',FORM='formatted',&
           STATUS='REPLACE')
	do h=1,max_in
		write (3,*) in(h)
	enddo
!	write (3, *) prop1(1,k,j), k, j, i, h, in(h)
!	close(3)     
	

!**************************************************************************
!*  DLL Inputs from GoldSim:                                              *
!*                                                                        *
!*    in(1) = Realization Number [file_index]                             *
!*    in(2)-(36) = Uncertain Array Switches : tells the DLL               *
!*               which of the 35 parameters in the input file you want to *
!*               replace with uncertain values. [uncert_switch]           *
!*    in(37)-(72) =Uncertain Array Sizes : tells the DLL the maximum size *
!*               of each of the 35 uncertain parameters passed by Goldsim *
!*               [array_size]                                             *
!*    array_index = the starting point in the Goldsim array [in(*)] at    *
!*                which the uncertain values begin. Used throughout the   *
!*                code to index the proper location for the uncert values.*
!**************************************************************************


    file_index=in(1)
	
	h=2	
	do i=1, 35
	  uncert_switch(i)=in(h)
	  write (3,*) uncert_switch(i)
	  h=h+1
	enddo
    
	do i=1, 35
	  array_size(i)=in(h)
	  write (3,*) array_size(i)
	  h=h+1
	enddo

    array_index = 72   !72 is the beginning of the uncertain data in the input array
		i=1
!************************************************************
!*                **** DATA 9:****                          *
!* CSAT(iso) : Solubility limit for contaminant iso (g/cm3) *
!************************************************************
	i=1
	    write (3,*) '------------- DATA 9: CSAT(iso) -----------------------'
		write (3, *) i, uncert_switch(i)
	if(uncert_switch(i).eq.1) then
			j=1
			h=array_index
	      do j=1, NISO
			csat(j)=in(h)
			h=h+1			
		  enddo
	endif
    
	array_index = array_index+array_size(i)
    write (3, *) array_index
!************************************************************
!*                **** DATA 13:****                         *
!* prop(1,i) : Molecular Diff. coef. (cm2/s) of material i  *
!*                                                          *
!**                    **********                          **
!* prop(2,i) : bulk density of soil (g/cm3)  of material i  *
!*                                                          *
!**                    **********                          **
!* prop(3,i) : longitudinal dispersivity (cm) of material i *
!*                                                          *
!**                    **********                          **
!* prop(4,i) : transverse dispersivity (cm)  of material i  *
!*                                                          *
!**                    **********                          **
!* prop(5,i) : porosity of material i	                    *
!************************************************************

	do k=1, 5	
		i=i+1
		write (3,*) '------------- DATA 13: prop(',k,',i) -----------------------'
		write (3, *) i, uncert_switch(i)
	if(uncert_switch(i).eq.1) then
			j=1
			h=array_index
	      do j=1, NMAT
			prop(k,j)=in(h)
			h=h+1			
		  enddo	
	endif
	array_index = array_index+array_size(i)
	enddo
    write (3, *) array_index
    
!************************************************************
!*                 **** DATA 13: ****                       *
!* prop1(1,k,i):dist. coef.(cm3/g) of material k & contam. i*
!************************************************************
	i=i+1
	write (3,*) '------------- DATA 13: prop1(1,k,i) -----------------------'
	write (3, *) i, uncert_switch(i)
	h=array_index
	if(uncert_switch(i).eq.1) then
		do k=1, NMAT	
			j=1		
	      do j=1, NISO
			prop1(1,k,j)=in(h)
			h=h+1
			write(3,*) prop1(1,k,j)			
		  enddo	
		h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index
!************************************************************
!*   **** DATA SET 16: Cauchy & Neuman Flux Values	*****   *
!* QCBF(j,i,iso) & QNBF(j,i,iso):flux(mass/cm2/yr) of jth   *
!* data point in the ith profile, for contaminant iso. For  *
!* this application, all points in grid are the same.       *
!************************************************************


!*----- Cauchy (total flux) profiles
		 i=i+1
		 h=array_index
       write (3,*) '------------- DATA 16A: CAUCHY -----------------------'
	   write (3, *) i, uncert_switch(i)
	   if(uncert_switch(i).eq.1) then
         do k=1,niso
		    write (3,*) k, ' of ',niso,' isotopes'
		    do l=1,ncpr(k)
			   write (3,*) l,' of ',ncpr(k),' profiles ', ncdp(k),' data points'
               do j=1,ncdp(k)
			      write(3,*) k,j,l,h,in(h)
				  cvbf(k,j,l)= in(h)
                  write (3,*) cvbf(k,j,l)
               enddo 
			enddo
			h=h+1 
	     enddo
	    endif
		array_index=array_index+array_size(i)
        write (3, *) array_index 
        
!*----- Neumann (dispersive flux) profiles
		 i=i+1
		 h=array_index
       write (3,*) '-------------- DATA 16B: NEUMANN ----------------------'
	   write (3, *) i, uncert_switch(i)
	   if(uncert_switch(i).eq.1) then
         do k=1,niso
		    write (3,*) k, ' of ',niso,' isotopes'
		    do l=1,nnpr(k)
			   write (3,*) l,' of ',nnpr(k),' profiles ', ncdp(k),' data points'
               do j=1,ncdp(k)
			      write(3,*) k,j,l,h,in(h)
				  cvbfn(k,j,l)= in(h)
                  write (3,*) cvbfn(k,j,l)
               enddo 
			enddo
			h=h+1 
	     enddo
	    endif
		array_index=array_index+array_size(i)
        write (3, *) array_index 

!************************************************************
!*     ****** DATA SET 17: Dirichlet Conc Value	*******     *
!* CDBF(j,i,iso):flux(mass/cm3) of jth data point in the    *
!*			     ith profile, for contaminant iso. For this *
!*		       application, all points in grid are the same *
!************************************************************

!*----- Dirichlet(concentrations) profiles
		 i=i+1
		 h=array_index
       write (3,*) '-------------- DATA 17: DIRICHLET CONCENTRATIONS --------------'
	   write (3, *) i, uncert_switch(i)
	   if(uncert_switch(i).eq.1) then
         do k=1,niso
		    write (3,*) k, ' of ',niso,' isotopes'
		    do l=1,ndpr(k)
			   write (3,*) l,' of ',ndpr(k),' profiles ', ncdp(k),' data points'
               do j=1,ncdp(k)
			      write(3,*) k,j,l,h,in(h)
				  cdbf(k,j,l)= in(h)
                  write (3,*) cdbf(k,j,l)
               enddo 
			enddo
			h=h+1 
	     enddo
	    endif
		array_index=array_index+array_size(i)
        write (3, *) array_index 

!************************************************************
!*  **** DATA SET 20: Container Degradation Parameters **** *
!* thick(j) :thickness(cm) of jth container                 *
!* pitn(j)  :pitting parameter n for jth container          *
!* pitk(j)  :pitting parameter k for jth container          *
!* area(j)  :area of jth container(cm2)                     *
!* ascale(j):area scaling exponent for jth container        *
!* pits(j)  :number of penetrating pits in jth container    *
!* grate(j) :general corosion rate (cm/s) in jth container  *
!* clay(j)  :clay fraction in the soil around jth container *
!* sph(j)   :soil pH around jth container                   *
!* iaer(j)  :aeration index for soil around jth container   *
!************************************************************


	    write (3,*) '------------- DATA 20: Container Degradation Parameters ------------------'
		
!*----- thick(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- thick(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			thick(j)=in(h)
			write (3, *) thick(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- pitn(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- pitn(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			pitn(j)=in(h)
			write (3, *) pitn(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- pitk(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- pitk(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			pitk(j)=in(h)
			write (3, *) pitk(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- area(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- area(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			area(j)=in(h)
			write (3, *) area(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- ascale(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- ascale(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			ascale(j)=in(h)
			write (3, *) ascale(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- pits(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- pits(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			pits(j)=in(h)
			write (3, *) pits(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- grate(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- grate(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			grate(j)=in(h)
			write (3, *) grate(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- clay(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- clay(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			clay(j)=in(h)
			write (3, *) clay(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- sph(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- sph(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			sph(j)=in(h)
			write (3, *) sph(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- iaer(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- iaer(j) ----'
		write (3, *) i, uncert_switch(i), NCTYPE
		h=array_index
		do j=1, NCTYPE
			iaer(j)=in(h)
			write (3, *) iaer(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!************************************************************
!*  **** DATA SET 21: Waste Form Leaching Parameters   **** *
!* sfract(j,iso) :fraction of mass in the jth waste form    *
!*    available for wash-off in rinse release for iso       *
!* pfract(j,iso) :fraction of mass in the jth waste form    *
!*    available for diffusion controlled release for iso    *
!* bfract(j,iso) :fraction of mass in the jth waste form    *
!*    available for dissolution controlled release for iso  *
!* deff(j,iso)   :efffective difffusion coefficient (cm2/s) *
!*    of jth waste formfor contaminant iso.                 *
!* disol(j,iso)  :fractional release rate (1/yr) for the jth*
!*    waste form for contaminant iso.                       *
!* partko(j,iso)  :partition coefficient (cm3/g) for the jth*
!*    waste form for contaminant iso at time t=0.           *
!* partki(j,iso)  :partition coefficient (cm3/g) for the jth*
!*    waste form for contaminant iso. time=infinity.        *
!* partd(j,iso)  :partition coefficient degradation rate    *
!*    constant(1/yr) for the jth waste form, contaminant iso*
!* porel(j)  :radius or half-length of jth waste form (cm)  *
!* volwf(j)  :volume of jth waste form (cm3)                *
!* vartion(j)  :ratio of the volume of jth waste form to the*
!*     finite element in which it is located.               *
!* wtinit(j,iso)  : initial inventory of iso in jth waste   *
!*     form (if IACT=0, grams or if IACT=1, curies)         *
!************************************************************

write (3,*) '------------- DATA 21: Waste Form Leaching Parameters ------------------'

!*-------sfract(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- sfract(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			sfract(k,j)=in(h)
			write (3, *) sfract(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------pfract(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- pfract(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			pfract(k,j)=in(h)
			write (3, *) pfract(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------bfract(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- bfract(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			bfract(k,j)=in(h)
			write (3, *) bfract(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------deff(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- deff(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			deff(k,j)=in(h)
			write (3, *) deff(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------disol(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- disol(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			disol(k,j)=in(h)
			write (3, *) disol(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------partko(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- partk0(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			partk0(k,j)=in(h)
			write (3, *) partk0(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------partki(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- partki(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			partki(k,j)=in(h)
			write (3, *) partki(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*-------partd(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- partd(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			partd(k,j)=in(h)
			write (3, *) partd(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- porel(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- porel(j) ----'
		write (3, *) i, uncert_switch(i), nwtype
		h=array_index
		do j=1, nwtype
			porel(j)=in(h)
			write (3, *) porel(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- volwf(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- volwf(j) ----'
		write (3, *) i, uncert_switch(i), nwtype
		h=array_index
		do j=1, nwtype
			volwf(j)=in(h)
			write (3, *) volwf(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*----- vratio(j) ----------
	i=i+1			
	if(uncert_switch(i).eq.1) then
		write (3, *) '---- vratio(j) ----'
		write (3, *) i, uncert_switch(i), nwtype
		h=array_index
		do j=1, nwtype
			vratio(j)=in(h)
			write (3, *) vratio(j)
			h=h+1			
		  enddo
	endif    
	array_index = array_index+array_size(i)
    write (3, *) array_index

!*------- wtinit(j,iso)-----------
	i=i+1
	h=array_index
	write (3, *) '---- wtinit(j,iso) ----'
	write (3, *) i, uncert_switch(i), niso, ' isotopes'
	if(uncert_switch(i).eq.1) then
		do k=1, nwtype	
			j=1
			write(3,*) 'waste form',nwtype  
	      do j=1, NISO
			wtinit(k,j)=in(h)
			write (3, *) wtinit(k,j)
			h=h+1			
		  enddo	
			h=array_index+(10*k)
		enddo	
	endif
	array_index = array_index+array_size(i)
    write (3, *) array_index

!************************************************************
!* ** DATA SET 23:Velocity & Moisture Content Parameters ** *
!*                                                          *
!* vxni23(j):darcy velocity for node i in the x-direction   *
!*           for element j.
!* vzni23(j):darcy velocity for node i in the y-direction   *
!*           for element j.                                 *
!* thni23(j):moisture content at nodal point IQ				
!************************************************************

!*------- vxni23(j)-----------
	i=i+1
	h=array_index
	write (3, *) '---- vxni23(j) ----'
	write (3, *) i, uncert_switch(i)
	if(uncert_switch(i).eq.1) then
         j = 0
 1550      j = j + 1
         if (ni23(j) .eq. 0) goto 1551
		 vxni23(j)=vxni23(j)*in(h)
		 write (3, *) ni23(j),' Nodes', vxni23(j)
!		 h=h+1										!comment out element specific multiplier term (decreased in(*) array by 599)
		 if (ni23(j) .ne. 0) goto 1550
         
    endif
1551	array_index = array_index+array_size(i)
	write (3, *) array_index
!*------- vzni23(j)-----------
	i=i+1
	h=array_index
	write (3, *) '---- vzni23(j) ----'
	write (3, *) i, uncert_switch(i)
	if(uncert_switch(i).eq.1) then
         j = 0
 1555      j = j + 1
		 if (ni23(j) .eq. 0) goto 1556
         vzni23(j)=vzni23(j)*in(h)
		 write (3, *) ni23(j),' Nodes', vzni23(j)
!		 h=h+1										!comment out element specific multiplier term (decreased in(*) array by 599)
		 if (ni23(j) .ne. 0) goto 1555
    endif
 1556	array_index = array_index+array_size(i)
	write (3, *) array_index
!*------- thni23(j)-----------
	i=i+1
	h=array_index
	write (3, *) '---- thni23(j) ----'
	write (3, *) i, uncert_switch(i)
	if(uncert_switch(i).eq.1) then
         j = 0
 1600      j = j + 1
		   if (ni23(j) .eq. 0) goto 1601
           thni23(j)=thni23(j)*in(h)
		   write (3, *) ni23b(j),' Nodes', thni23(j)
!		   h=h+1									!comment out element specific multiplier term (decreased in(*) array by 599)
         if (ni23b(j) .ne. 0) goto 1600
	endif 

 1601   array_index = array_index+array_size(i)
		write (3, *) array_index
!---------------------------------------------------------------------------------*
!*		                                                                          *
!* write the data back out to an input file with a generic name - goldblt_01.inp  *
!*                                                                                *
!*--------------------------------------------------------------------------------*
      
	   
      WRITE (in_file,1011) file_index
1011  FORMAT('BLTMS',I5.5,'.inp')

 
      open (1,FILE=in_file,ACCESS='sequential',FORM='formatted',&
           STATUS='REPLACE')
         call MakeDeck()

!----------------------------------------------------------------------------------
!                        ***********************
!                        *     End of Main     *
!                        ***********************
!----------------------------------------------------------------------------------

!*!*************Pass Out Array back to Goldsim ************************


	  
	  out(1)=niso		!Pass Back the number of Isotopes in the model run
	  out(2)=nti		!Pass Back the number of Time Steps
	  out(3)=NSTPTR		!Pass Back the number of steps between traces
	  out(4)=NTRC		!Pass Back the number of nodes with concentration traces
	  out(5)=NCON		!Pass Back the number of containers

	write(3,*)'----------Output back to Goldsim-------------'
	write(3,*) 'NISO=',out(1)
	write(3,*) 'NTI=',out(2)
	write(3,*) 'NSTPTR=', out(3)
	write(3,*) 'NTRC=', out(4)
	write(3,*) 'NCON=',out(5)
	

elseif (method.eq.99) then			!Cleanup
!  close the two files

  	  close(3)
      close(1)   
else
	write(2,*) 'DLL was called with an invald argument'

endif

return

end subroutine read_blt

!*------------------------------------------------------------------------------*
!*                 Begining of Subroutine ReadDeck
!*
!*------------------------------------------------------------------------------*

      SUBROUTINE ReadDeck()
!*-----
!*-----This subroutine is designed to read in an already existing input deck
!*-----
      character label*79, title*79, cnam*10
      character*79 blank
      integer index, rep
!*
      common /data1/ nprob, title
      common /data2/ kvi, kss
      common /data3/ niso, ncon, nctype, nwtype, idiff, iact
      common /data4/ ilump, imid, iwet, ioptim, ivml
      common /data5/ nti, ndtchg
      common /data6/ kmesh, nnp, nel, nmat, ncm
      common /data7A/ ndnp(10), ndpr(10), nddp(10)
      common /data7B/ ncnp(10), nces(10), ncpr(10), ncdp(10)
      common /data7C/ nnnp(10), nnes(10), nnpr(10), nndp(10)
      common /data8/ nsel(10),nspr(10),nsdp(10),nwnp(10),nwpr(10)&
              ,nwdp(10)
      common /data9/  cnam(10), decayc(10), csat(10), atm(10)
      common /data9a/ ichain, lchain(10), idchain(10,10)
      common /data9b/ bfr(10,10)
      common /data10/ delt, chng, delmax, tmax, w, wv,tdtch(100)
      common /data12/ ni(2,600),nseq(2,600),nad(2,600),&
              xni(600),xad(600),xrd(600),&
              zni(600),zad(600),zrd(600),&
              nelx,nelz,nreg,&
              imin(50),jmin(50),imax(50),jmax(50),&
              x1(4,50),x2(4,50),mi(550),ie(550,5),&
              nlay(550), modl(550),widfac
      common /data13/ prop(8,5),prop1(8,5,10)
      common /data14/ mi14(600),nseq14(600),mad(600),mityp(600),&
              mtypad(600)
      common /data15/ ni15(10,600),nnode(10,600),nad15(10,600),&
              rni(10,600),rad(10,600),rrd(10,600)
      common /data16/ cvbf(10,10,6),tcvbf(10,10,6),mi16(10,600),&
              nseq16(10,600),miad16(10,600),mityp16(10,600),&
              mtypad16(10,600),npcb(600,10),mi_v16(10,600),&
              nseq_v16(10,600),m_v16(10,600),is1_v16(10,600),&
              is2_v16(10,600),miad_v16(10,600),mad_v16(10,600),&
              is1ad(10,600),is2ad(10,600)
      common /data16n/ cvbfn(10,10,6),tcvbfn(10,10,6),mi16n(10,600),&
              nseq16n(10,600),miad16n(10,600),mityp16n(10,600),&
              mtpad16n(10,600), npnb(600,10),&
              mi_v16n(10,600),nseq_16n(10,600),m_v16n(10,600),&
              is1_v16n(10,600),is2_v16n(10,600),miad_16n(10,600),&
              mad_v16n(10,600),is1adn(10,600),is2adn(10,600)
      common /data17/ cdbf(10,10,6),tcdbf(10,10,6),ni17(10,600),&
              nseq17(10,600),niad17(10,600),nityp17(10,600),&
              ntypad17(10,600),npdb(600,10)
      common /data18/ sosf(10,10,6),tsosf(10,10,6),ni18(10,600),&
              nseq18(10,600),niad18(10,600),nityp18(10,600),&
              ntypad18(10,600),les(600,10)
      common /data19/ wssf(10,10,6),twssf(10,10,6),ni19(10,600),&
              nseq19(10,600),niad19(10,600),nityp19(10,600),&
              ntypad19(10,600),npw(600,10)
      common /data20/ thick(20),pitn(20),pitk(20),area(20),ascale(20),&
              pits(20),grate(20),clay(20),sph(20),iaer(20),&
              nelcon(100),mi20(600),nseq20(600),mad20(600),&
              mityp20(600),mtypad20(600)
      common /data21/ sfract(20,10),pfract(20,10),bfract(20,10),&
              deff(20,10),disol(20,10),partk0(20,10),&
              partki(20,10),partd(20,10),porel(20),volwf(20),&
              vratio(20),mi21(600),nseq21(600),mad21(600),&
              mityp21(600),mtypad21(600),wtinit(100,10)
      common /data22/ kpr0, kpr(1000),ntrc,ntrf,nstptr,lctrc(20),&
              lftrc(20)
      common /data23/ ni23(600),nseq23(600),nad23(600),vxni23(600),&
              vzni23(600),vxad23(600),vzad23(600),ni23b(600),&
              nseq23b(600),nad23b(600),thni23(600),&
              thniad23(600)
!*
!*
!*-----Data Set 1:Problem Identification and Description
!*
      read(2,10) label
      read(2,20) nprob,title
!*
!*-----data set 2: Integer Velocity Input and Steady-state control
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) kvi,kss
!*
!*-----Data Set 3: Integer Parameters for Containers & Waste
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) niso,ncon,nctype,nwtype,idiff,iact
!*
!*-----Data Set 4: Integration Integer Parameters
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) ilump,imid,iwet,ioptim,ivml
!*
!*-----Data Set 5: Integer Solution Control Parameters
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) nti,ndtchg
!*
!*-----Data Set 6: Grid And Element Parameters
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) nnp,nel,nmat,ncm,kmesh
!*
!*-----Data Set 7: Dirichlet Integer Parameters for Boundary Conditions
!*-----a) Dirichlet b.c. Data
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do iso = 1,niso
      read(2,*) ndnp(iso),ndpr(iso),nddp(iso)
      enddo
!*
!*-----b) Cauchy b.c. data
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do iso = 1,niso
      read(2,*) ncnp(iso),nces(iso),ncpr(iso),ncdp(iso)
      enddo

!*-----c) Neuman b.c. data
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do iso = 1,niso
      read(2,*) nnnp(iso),nnes(iso),nnpr(iso),nndp(iso)
      enddo
!*
!*-----Data Set 8: Integer Parameters for Sources
!*-----a) Element Source/Sink Data
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do iso = 1,niso
      read(2,*) nsel(iso),nspr(iso),nsdp(iso)
      enddo
!*
!*-----b) Point Source/Sink Data
!*
      read(2,10) label
      read(2,10) label
      do iso = 1,niso
      read(2,*) nwnp(iso),nwpr(iso),nwdp(iso)
      enddo
!*
!*-----Data Set 9: Chemical Component Information
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do 100 index=1,niso
         read(2,*) cnam(index),ATM(index),decayc(index),csat(index)
  100 continue
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) ichain
      if(ichain .ge. 1) then
      read(2,10) label
      read(2,10) label
      read(2,*) (lchain(j),j=1,ichain)
      read(2,10) label
      read(2,10) label
      do i = 1,ichain
      read(2,*) (idchain(j,i),j=1,lchain(i))
      enddo
      read(2,10) label
      read(2,10) label
      do i = 1,ichain
      read(2,*) (bfr(j,i),j=1,lchain(i)-1)
      enddo
      endif
!*
!*-----Data Set 10: Time Integration Control Parameters
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) delt,chng,delmax,tmax,w,wv
      delt0 = delt
      delt0 = delt0
      if(ndtchg .gt. 0) then
      read(2,10) label
      read(2,10) label
      read(2,*) (tdtch(i),i=1,ndtchg )
      endif
!*
!*-----Data Set 11: Real Parameters for Nonlinear Solve
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
!*
!*---- Data Set 12a: Node Coordinate via Card if KVI is less than 0
!*
!*   --read a blank line, title of data set 12 and a note that
!*   --data sets 12a and b are not Used
!*
         read(2,10) label
         read(2,10) label
         read(2,10) label
         read(2,*) widfac
         read(2,10) label
            read(2,10) label
      if (kvi .le. 0) then
         if (kmesh .gt. 0) then
            read(2,10) label
            call rnode()
         else
!*--------readr
            i = 0
  200         i = i + 1
              read(2,*) ni(1,i),nseq(1,i), nad(1,i),xni(i),xad(i),xrd(i)
            if (ni(1,i) .ne. 0) goto 200
!*
            read(2,10) label
            read(2,10) label
!*--------readr
            i = 0
  300         i = i + 1
              read(2,*) ni(2,i),nseq(2,i),nad(2,i),zni(i),zad(i),zrd(i)
            if (ni(2,i) .ne. 0) goto 300
         endif
!*
!*------ Data Set 12b: Element Indices if kmesh <= 0
!*
         mj = 0
         if (kmesh .le. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) label
              rep = 0
              mj = 0
!*------ Code taken from bltec code. (forgive the insanity)
  400         rep = rep + 1
              read(2,*) mi(rep),ai,bi,ci,di,ei
              if (mi(rep) .gt. 0) then
                 ie(mi(rep),1) = ai
                 ie(mi(rep),2) = bi
                 ie(mi(rep),3) = ci
                 ie(mi(rep),4) = di
                 ie(mi(rep),5) = ei
              endif
              read(2,10) label
              read(2,10) label
              read(2,*) modl(rep),nlay(rep)
  450         mj = mj + 1
              if (nel .eq. 0) goto 500
              if (mj .lt. mi(rep)) goto 450
              if (mj .eq. nel) goto 500
              if (modl(rep) .le. 0) goto 400
              mj = mj + (nlay(rep) * modl(rep)) - 1
              if (mj .lt. nel) goto 400
  500         continue
         endif
      endif
!*
!*---- Data Set 13: Porous Media Properties
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do 600 i=1,nmat
         read(2,*) (prop(j,i),j=1,5)
  600 continue
      read(2,10) label
      read(2,10) label
      read(2,10) label
      do k =1,nmat
      do i =1,niso
      read(2,*) (prop1(k,j,i),j=1,3)
      enddo
      enddo
!*
!*---- Data Set 14: Material Type Correction
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      if (ncm .gt. 0) then
         i = 0
 650       i = i + 1
           read(2,*) mi14(i),nseq14(i),mad(i),mityp(i),mtypad(i)
         if (mi14(i) .ne. 0) goto 650
      endif
!*
!*---- Data Set 15: Inital or Pre-Initial Conditions
!*
         read(2,10) label
         read(2,10) label
         read(2,10) label
         do 725 i=1,niso
            read(2,10) label
            j = 0
  700         j = j + 1
              read(2,*) ni15(i,j),nnode(i,j),nad15(i,j),rni(i,j),&
               rad(i,j),rrd(i,j)
            if (ni15(i,j) .ne. 0) goto 700
  725    continue
!*
!*---- Data Set 16: Cauchy Boundary Conditions
!*
      icheck = 0
      do iso=1,niso
      icheck = icheck + ncnp(iso)
      enddo
         read(2,10) label
      if (icheck .le. 0) then
         read(2,10) label
         read(2,10) label
      else
!*----- Read Cauchy (total flux) profiles
         do 800 k=1,niso
         if(nces(k) .eq. 0 .or. ncnp(k) .eq. 0) go to 800
         read(2,10) label
            read(2,10) label
            read(2,10) label
            do 750 i=1,ncpr(k)
               read(2,*) (tcvbf(k,j,i),cvbf(k,j,i),j=1,ncdp(k))
  750       continue
!*----- Read global node information
         read(2,10) label
         read(2,10) label
         read (2,*)  (NPCB(i,k),i=1,NCNP(k))
!*----- Read incoming flux type assigned to each of Cauchy sides
            read(2,10) label
            read(2,10) label
            j = 0
  775         j = j + 1
              read(2,*) mi16(k,j),nseq16(k,j),miad16(k,j),mityp16(k,j),&
              mtypad16(k,j)
            if (mi16(k,j) .ne. 0) goto 775
!*----- Read Cauchy boundary element sides
         read(2,10) label
         read(2,10) label
         j = 0
  797      j = j + 1
         read(2,*) mi_v16(k,j),nseq_v16(k,j),m_v16(k,j),is1_v16(k,j),&
         is2_v16(k,j),miad_v16(k,j),mad_v16(k,j),is1ad(k,j),is2ad(k,j)
         if (mi_v16(k,j) .ne. 0) goto 797
  800    continue
         endif
!  Neumann boundary conditions
!*
      icheck = 0
      do iso=1,niso
      icheck = icheck + nnnp(iso)
      enddo
      if (icheck .le. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) label
      else
         read(2,10) label
!*----- Read Neumann (diffusive flux) profiles
         do 900 k=1,niso
         if(nnnp(k) .eq. 0 .or. nnes(k) .eq. 0) go to 900
            read(2,10) label
            read(2,10) label
            read(2,10) label
            do i=1,nnpr(k)
               read(2,*) (tcvbfn(k,j,i),cvbfn(k,j,i),j=1,nndp(k))
            enddo
!*----- Read global node information
         read(2,10) label
         read(2,10) label
         read (2,*)  (NPNB(i,k),i=1,NNNP(k))

!*----- Read incoming profile type assigned to each of Neumann sides
            read(2,10) label
            read(2,10) label
            j = 0
  875         j = j + 1
              read(2,*) mi16n(k,j),nseq16n(k,j),miad16n(k,j),&
                     mityp16n(k,j), mtpad16n(k,j)
      if (mi16n(k,j) .ne. 0) goto 875
!*----- Read Neumann boundary element sides
         read(2,10) label
         read(2,10) label
         j = 0
 890       j = j + 1
       read(2,*) mi_v16n(k,j),nseq_16n(k,j),m_v16n(k,j),is1_v16n(k,j),&
        is2_v16n(k,j),miad_16n(k,j),mad_v16n(k,j),is1adn(k,j),is2adn(k,j)
         if (mi_v16n(k,j) .ne. 0) goto 890
  900  continue
      endif
!*
!*---- Data Set 17: Dirichlet Boundary Conditions
!*
      icheck = 0
      do iso = 1,niso
      icheck =icheck + ndnp(iso)
      enddo
      if (icheck .eq. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) label
      else
         read(2,10) label
!*----- Read concentration profiles
         do 1000 k=1,niso
            if(ndnp(k) .ne. 0) then
            read(2,10) label
            read(2,10) label
            read(2,10) label
            do 950 i=1,ndpr(k)
               read(2,*) (tcdbf(k,j,i),cdbf(k,j,i),j=1,nddp(k))
  950       continue
!*----- Read global node information
         read(2,10) label
         read(2,10) label
           read(2,*)  (NPDB(j,k),j=1,ndnp(k))
!*----- Read incoming concentration type assigned to each of Dirichlet point
            read(2,10) label
            read(2,10) label
            j = 0
  975         j = j + 1
              read(2,*) ni17(k,j),nseq17(k,j),niad17(k,j),nityp17(k,j),&
              ntypad17(k,j)
            if (ni17(k,j) .ne. 0) goto 975
            endif
 1000    continue
      endif
!*
!*---- Data Set 18: External Volumetric Sources
!*
      icheck = 0
      do iso = 1,niso
      icheck =icheck + nsel(iso)
      enddo
      if (icheck .eq. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) blank
      else
      read(2,10) label
      read(2,10) label
      read(2,10) label
!*----- Read source strength profiles
         do 1001 k=1,niso
             if(nsel(k) .eq. 0) go to 1001
!*----- Read source type assigned to each element
      read(2,10) label
             read(2,*)  (LES(j,k),j=1,nsel(k))
             read(2,10) blank
             read(2,10) label
             do 951 i=1,nspr(k)
               read(2,*) (tsosf(k,j,i),sosf(k,j,i),j=1,nsdp(k))
             read(2,10) blank
  951     continue
!*----- Read volumetric source element assigned to a specific profile
             read(2,10) label
             j = 0
  976         j = j + 1
              read(2,*) ni18(k,j),nseq18(k,j),niad18(k,j),nityp18(k,j),&
              ntypad18(k,j)
             if (ni18(k,j) .ne. 0) goto 976
      if(k.ne.niso) then
      read(2,10) blank
      endif
 1001    continue
      endif
!*
!*-----Data Set 19: External Point (Well) Sources
!*
      icheck = 0
      do iso = 1,niso
      icheck =icheck + nwnp(iso)
      enddo
      if (icheck .eq. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) blank
      else
      read(2,10) label
      read(2,10) label
      read(2,10) label
!*----- Read source strength profiles
         do 1002 k=1,niso
             if(nwnp(k) .eq. 0) go to 1002
!*----- Read source type assigned to each element
      read(2,10) label
             read(2,*)  (NPW(j,k),j=1,nwnp(k))
             read(2,10) blank
      read(2,10) label
             do 952 i=1,nwpr(k)
               read(2,*) (twssf(k,j,i),wssf(k,j,i),j=1,nwdp(k))
             read(2,10) blank
  952     continue
!*----- Read point source element assigned to a specific profile
             read(2,10) label
             j = 0
  977         j = j + 1
              read(2,*) ni19(k,j),nseq19(k,j),niad19(k,j),nityp19(k,j),&
              ntypad19(k,j)
             if (ni19(k,j) .ne. 0) goto 977
      if(k.ne.niso) then
      read(2,10) blank
      endif
 1002    continue
      endif
!*
!*-----Data Set 20: Waste Container Parameters
!*
      if (ncon .eq. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) label
      else
         read(2,10) label
         do 1200 j=1,nctype
            read(2,10) label
      read(2,10) label
            read(2,*) thick(j),pitn(j),pitk(j),area(j),ascale(j),pits(j)
            read(2,10) label
            read(2,10) label
            read(2,*) grate(j),clay(j),sph(j),iaer(j)
 1200    continue
!*-----Read global element no. and the cont. type assignment
         read(2,10) label
         read(2,10) label
         read(2,*) (nelcon(i),i=1,ncon)
         read(2,10) label
         read(2,10) label
         j = 0
 1225      j = j + 1
           read(2,*) mi20(j),nseq20(j),mad20(j),mityp20(j),mtypad20(j)
         if (mi20(j) .ne. 0) goto 1225
      endif
!*
!*-----Data Set 21a: Waste Form and Leaching Information
!*
      if (ncon .eq. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) label
      else
         read(2,10) label
         do 1350 j=1,nwtype
            read(2,10) label
            read(2,10) label
            do 1300 i=1,niso
               read(2,10) label
               read(2,*) sfract(j,i),pfract(j,i),bfract(j,i),deff(j,i),&
               disol(j,i)
               read(2,10) label
               read(2,10) label
               read(2,*) partk0(j,i),partki(j,i),partd(j,i)
               read(2,10) label
 1300       continue
            read(2,10) label
            read(2,*) porel(j),volwf(j),vratio(j)
 1350    continue
!*
!*-----Data Set 21b: Read Initial Masses in each Waste Element
!*
         read(2,10) label
         do 1400 i=1,niso
            read(2,10) label
            read(2,10) label
            read(2,10) label
            read(2,*) (wtinit(j,i),j=1,ncon)
 1400    continue
!*
!*-----Data Set 21c: Read Waste Type Assignments
!*
         read(2,10) label
         read(2,10) label
         read(2,10) label
         j = 0
 1450      j = j + 1
           read(2,*) mi21(j),nseq21(j),mad21(j),mityp21(j),mtypad21(j)
         if (mi21(j) .ne. 0) goto 1450
      endif
!*
!*-----Data Set 22a: Printer and Auxiliary Storage Control
!*
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,10) label
      nti1 = nti
      if(nti .gt. 1000) nti1 = 1000
      read(2,30) kpr0,(kpr(i),i=1,nti1)
!*
!*-----Data Set 22b: Chemical Output and Chemical Property Type Indicator
      read(2,10) label
      read(2,10) label
      read(2,10) label
!*-----Data Set 22c: Trace File Indicators.
      read(2,10) label
      read(2,10) label
      read(2,10) label
      read(2,*) ntrc,ntrf,nstptr
      if (ntrc .gt. 0) then
      read(2,10) label
      read(2,10) label
      read(2,*) (lctrc(i),i=1,ntrc)
      endif
      if (ntrf .gt. 0) then
      read(2,10) label
      read(2,10) label
      read(2,*) (lftrc(i),i=1,ntrf)
      endif
!*
!*-----Data Set 23: Hydrological Data
!*
      if (kvi .gt. 0) then
         read(2,10) label
         read(2,10) label
         read(2,10) label
      else
         read(2,10) label
         read(2,10) label
         read(2,10) label
         j = 0
 1550      j = j + 1
           read(2,*) ni23(j),nseq23(j),nad23(j),vxni23(j),vzni23(j),&
           vxad23(j),vzad23(j)
         if (ni23(j) .ne. 0) goto 1550
         read(2,10) label
         read(2,10) label
         j = 0
 1600      j = j + 1
           read(2,*) ni23b(j),nseq23b(j),nad23b(j),thni23(j),thniad23(j)
         if (ni23b(j) .ne. 0) goto 1600
      endif
!*
      label = label
      return
!*
!*----- Format Section
!*
   10 FORMAT(A79)
   11 format(1x,a79)
   20 FORMAT(I5,A79)
   30 FORMAT(80(I1))
   40 FORMAT(A13)
   50 format(5i10)
      END
!* End of Subroutine ReadDeck
!*------------------------------------------------------------------------------*
!*
!*------------------------------------------------------------------------------*
!*


!*------------------------------------------------------------------------------*
!* Begining of Subroutine MakeDeck
!*
!*
!*------------------------------------------------------------------------------*
      SUBROUTINE MakeDeck()
!*
      character label*79, vlabl*79, title*79, cnam*10, ans
      character vlabl2*79, nlabl*79
      integer rep
!*
      common /data1/ nprob, title
      common /data2/ kvi, kss
      common /data3/ niso, ncon, nctype, nwtype, idiff, iact
      common /data4/ ilump, imid, iwet, ioptim, ivml
      common /data5/ nti, ndtchg
      common /data6/ kmesh, nnp, nel, nmat, ncm
      common /data7A/ ndnp(10), ndpr(10), nddp(10)
      common /data7B/ ncnp(10), nces(10), ncpr(10), ncdp(10)
      common /data7C/ nnnp(10), nnes(10), nnpr(10), nndp(10)
      common /data8/ nsel(10),nspr(10),nsdp(10),nwnp(10),nwpr(10)&
              ,nwdp(10)
      common /data9/  cnam(10), decayc(10), csat(10), atm(10)
      common /data9a/ ichain, lchain(10), idchain(10,10)
      common /data9b/ bfr(10,10)
      common /data10/ delt, chng, delmax, tmax, w, wv,tdtch(100)
      common /data12/ ni(2,600),nseq(2,600),nad(2,600),&
              xni(600),xad(600),xrd(600),&
              zni(600),zad(600),zrd(600),&
              nelx,nelz,nreg,&
              imin(50),jmin(50),imax(50),jmax(50),&
              x1(4,50),x2(4,50),mi(550),ie(550,5),&
              nlay(550), modl(550),widfac
      common /data13/ prop(8,5),prop1(8,5,10)
      common /data14/ mi14(600),nseq14(600),mad(600),mityp(600),&
              mtypad(600)
      common /data15/ ni15(10,600),nnode(10,600),nad15(10,600),&
              rni(10,600),rad(10,600),rrd(10,600)
      common /data16/ cvbf(10,10,6),tcvbf(10,10,6),mi16(10,600),&
              nseq16(10,600),miad16(10,600),mityp16(10,600),&
              mtypad16(10,600),npcb(600,10),mi_v16(10,600),&
              nseq_v16(10,600),m_v16(10,600),is1_v16(10,600),&
              is2_v16(10,600),miad_v16(10,600),mad_v16(10,600),&
              is1ad(10,600),is2ad(10,600)
      common /data16n/ cvbfn(10,10,6),tcvbfn(10,10,6),mi16n(10,600),&
              nseq16n(10,600),miad16n(10,600),mityp16n(10,600),&
              mtpad16n(10,600), npnb(600,10),&
              mi_v16n(10,600),nseq_16n(10,600),m_v16n(10,600),&
              is1_v16n(10,600),is2_v16n(10,600),miad_16n(10,600),&
              mad_v16n(10,600),is1adn(10,600),is2adn(10,600)
      common /data17/ cdbf(10,10,6),tcdbf(10,10,6),ni17(10,600),&
              nseq17(10,600),niad17(10,600),nityp17(10,600),&
              ntypad17(10,600),npdb(600,10)
      common /data18/ sosf(10,10,6),tsosf(10,10,6),ni18(10,600),&
              nseq18(10,600),niad18(10,600),nityp18(10,600),&
              ntypad18(10,600),les(600,10)
      common /data19/ wssf(10,10,6),twssf(10,10,6),ni19(10,600),&
              nseq19(10,600),niad19(10,600),nityp19(10,600),&
              ntypad19(10,600),npw(600,10)
      common /data20/ thick(20),pitn(20),pitk(20),area(20),ascale(20),&
              pits(20),grate(20),clay(20),sph(20),iaer(20),&
              nelcon(100),mi20(600),nseq20(600),mad20(600),&
              mityp20(600),mtypad20(600)
      common /data21/ sfract(20,10),pfract(20,10),bfract(20,10),&
              deff(20,10),disol(20,10),partk0(20,10),&
              partki(20,10),partd(20,10),porel(20),volwf(20),&
              vratio(20),mi21(600),nseq21(600),mad21(600),&
              mityp21(600),mtypad21(600),wtinit(100,10)
      common /data22/ kpr0, kpr(1000),ntrc,ntrf,nstptr,lctrc(20),&
              lftrc(20)
      common /data23/ ni23(600),nseq23(600),nad23(600),vxni23(600),&
              vzni23(600),vxad23(600),vzad23(600),ni23b(600),&
              nseq23b(600),nad23b(600),thni23(600),&
              thniad23(600)
!*
!*
      nlabl ='            (DATA SET NOT USED)'
!*
      label ='*** DATA SET 1: BLT-MS PROBLEM IDENTIFICATION'
      write (1,2) label
      write (1,6) nprob,title
      write (1,2)
!*
      label ='*** DATA SET 2: INTEGER PARAMETERS FOR VELOCITY INPUT AND&
     & STEADY STATE'
      write (1,2) label
      vlabl ='     KVI         KSS'
      write (1,2) vlabl
      write (1,10) kvi,kss
      write (1,2)
!*
      label ='*** DATA SET 3: PARAMETERS FOR CONTAINERS AND WASTE'
      vlabl ='    NISO        NCON      NCTYPE      NWTYPE       IDIFF&
     &      IACT'
      write (1,2) label
      write (1,2) vlabl
      write (1,10) niso,ncon,nctype,nwtype,idiff,iact
      write (1,2)
!*
      label ='*** DATA SET 4: INTEGRATION INTEGER PARAMETERS'
      vlabl ='   ILUMP        IMID        IWET      IOPTIM        IVML'
      write (1,2) label
      write (1,2) vlabl
      write (1,10) ilump,imid,iwet,ioptim,ivml
      write (1,2)
!*
      label ='*** DATA SET 5: INTEGER SOLUTION CONTROL PARAMETERS'
      vlabl ='    NTI      NDTCHG'
      write (1,2) label
      write (1,2) vlabl
      write (1,10) nti,ndtchg
      write (1,2)
!*
      label ='*** DATA SET 6: GRID AND ELEMENT PARAMETERS'
      vlabl ='    NNP         NEL        NMAT         NCM        KMSH'
      write (1,2) label
      write (1,2) vlabl
      write (1,10) nnp,nel,nmat,ncm,kmesh
      write (1,2)
!*
      label ='*** DATA SET 7a: INTEGER PARAMETERS FOR DIRICHLET BOUNDARY&
     & COND.'
      vlabl ='  NDNP(I)     NDPR(I)     NDDP(I)'
      write (1,2) label
      write (1,2) vlabl
      do iso = 1,niso
      write (1,11) ndnp(iso),ndpr(iso),nddp(iso)
      enddo
      write (1,2)
      label ='*** DATA SET 7b: INTEGER PARAMETERS FOR CAUCHY BOUNDARY '&
     &      //'COND.'
      write(1,2) label
      vlabl ='  NCNP(I)     NCES(I)     NCPR(I)     NCDP(I)   CONTAMINAN&
     &T'
      write (1,2) vlabl
      do iso = 1,niso
      write (1,12) ncnp(iso),nces(iso),ncpr(iso),ncdp(iso),cnam(iso)
      enddo
      write (1,2)
      label ='*** DATA SET 7c: INTEGER PARAMETERS FOR NEUMANN BOUNDARY '&
     &     //'COND.'
      write(1,2) label
      vlabl ='  NNNP(I)     NNES(I)     NNPR(I)     NNDP(I)   CONTAMINAN&
     &T'
      write (1,2) vlabl
      do iso = 1,niso
      write (1,12) nnnp(iso),nnes(iso),nnpr(iso),nndp(iso),cnam(ISO)
      enddo
      write (1,2)
!*
      label ='*** DATA SET 8: INTEGER PARAMETERS FOR SOURCE/SINKS'
      vlabl ='  NSEL(I)     NSPR(I)     NSDP(I)  VOLUME/ELEMENT SOURCES'
      write (1,2) label
      write (1,2) vlabl
      do iso = 1,niso
      write (1,11) nsel(iso),nspr(iso),nsdp(iso),cnam(iso)
      enddo
      write (1,2)
      vlabl ='  NWNP(I)     NWPR(I)     NWDP(I)   WELL/POINT SOURCES'
      write (1,2) vlabl
      do iso = 1,niso
      write (1,11) nwnp(iso),nwpr(iso),nwdp(iso),cnam(iso)
      enddo
      write (1,2)
!*
      label ='*** DATA SET 9: CHEMICAL COMPONENT INFORMATION AND DECAY P&
     &ROPERTIES'
      write(1,2) label
      vlabl ='  RADNUC(J)      ATM(J)       THALF(J)   CSAT(J) '
      write (1,2) vlabl
      do 400 j = 1,NISO
         write (1,20) cnam(j),atm(j),decayc(j),csat(j)
 400  continue
       write (1,2)
       label ='     DECAY CHAIN INFORMATION:'
       vlabl = '     NUMBER OF CHAINS '
       write(1,2) label
       write(1,2) vlabl
       write(1,10) ichain
       write(1,2)
       if(ichain .ge. 1) then
       vlabl = '     LENGTH OF THE CHAIN(S) '
       write(1,2) vlabl
       write(1,10) (lchain(i),i =1,ichain)
       write(1,2)
       vlabl = '     COMPONENTS OF THE CHAIN '
       write(1,2) vlabl
       do i=1,ichain
       write(1,10) (idchain(j,i),j=1,lchain(i))
       enddo
       write(1,2)
       vlabl = '     BRANCHING FRACTIONS OF THE CHAIN '
       write(1,2) vlabl
       do i=1,ichain
       write(1,16) (bfr(j,i),j=1,lchain(i)-1)
       enddo
       write(1,2)
       endif
!*
      label ='*** DATA SET 10: TIME INTEGRATION CONTROL PARAMETERS'
      vlabl ='    DELT        CHNG        DELMAX      TMAX        W&
     &      WV'
      write (1,2) label
      write (1,2) vlabl
      write (1,16) delt,chng,delmaX,tmaX,w,wv
      write (1,2)
      if(ndtchg .gt. 0) then
      vlabl = 'SET 10A:  DELTA T CHANGES AT YEAR(S) :  '
      write(1,2) vlabl
      write(1,16) (tdtch(i),i=1,ndtchg)
      write(1,2)
      endif
!*
      label ='*** DATA SET 11: REAL PARAMETERS FOR NONLINEAR SOLVE'
      write (1,2) label
      vlabl = ' DUMMY DATA SET: NOT USED IN BLT-MS '
      write (1,2) vlabl
      write (1,2)
!*
      label ='*** DATA SET 12a: NODE COORDINATE VIA CARD IF KVI <= 0'
!*   --read a blank line, title of data set 12 and a note that
!*   --data sets 12a and b are not Used
!*
      if (kvi .gt. 0) then
         write(1,2) label
         vlabl = 'WIDTH OF THE WASTE CONTAINING REGION (CM) '
         write(1,2) vlabl
         write(1,18) widfac
         write(1,2)
         vlabl ='    GENERATION OF MESH THROUGH AUXILIARY INPUT FILE'
         write(1,2) vlabl
      else
         label ='*** DATA SET 12a: NODAL POINT COORDINATE DATA '
         write(1,2) label
         vlabl = 'WIDTH OF THE WASTE CONTAINING REGION (CM) '
         write(1,2) vlabl
         write(1,18) widfac
         write(1,2)
         if (kmesh .gt. 0) then
            vlabl = '     AUTOMATIC MESH GENERATION '
            write(1,2) vlabl
            vlabl ='    NELX        NELZ        NREG'
            write(1,2) vlabl
            write(1,10) nelx,nelz,nreg
            vlabl ='REGION DEFINITION (NREG)'
            write(1,2) vlabl
            do 450 i=1,nreg
      jj=i/10
      if(jj.eq.0) jj=-16
      l=mod(i,10)
      vlabl = 'MIN AND MAX NODE NUMBERS FOR REGION: '//char(jj+48)//char&
     &(l+48)
               write(1,2) vlabl
               write(1,10) imin(i),jmin(i),imax(i),jmax(i)
       vlabl = 'VALUES FOR THE CORNERS OF REGION:   '//char(jj+48)//char&
     &(l+48)
               write(1,2) vlabl
               write(1,46) (x1(j,i),x2(j,i),j=1,4)
  450       continue
         else
!*---readr
      vlabl ='     NI     NSEQ      NAD    XNI        XAD        XRD'
            write(1,2) vlabl
            i = 0
  500         i = i + 1
              write(1,40) ni(1,i),nseq(1,i), nad(1,i),xni(i),xad(i),&
             xrd(i)
            if (ni(1,i) .ne. 0) goto 500
            vlabl = 'Z COORDINATES '
            write(1,2) vlabl
!*---readr
      vlabl ='     NI     NSEQ      NAD    ZNI        ZAD        ZRD'
            write(1,2) vlabl
            i = 0
  600         i = i + 1
              write(1,40) ni(2,i),nseq(2,i),nad(2,i),zni(i),zad(i),&
              zrd(i)
            if (ni(2,i) .ne. 0) goto 600
            write(1,2)
         endif
!*
!*------ DATA SET 12b: ELEMENT INDICES IF KMESH <= 0
!*
         label ='*** DATA SET 12b: ELEMENT INDICES'
      vlabl ='      MI         IE1         IE2         IE3         IE4'&
     &//'         IE5'
         if (kmesh .le. 0) then
            rep = 0
            mj = 0
!*------ Code taken from bltec code (forgive the insanity)
            write(1,2) label
            write(1,2) vlabl
  700       rep = rep + 1
            if (mi(rep) .eq. 0) then
!c fix in case MI not defined in BLTMSIN. Prevent search for IE(0,1) which
!c is out of bounds.
            ii = 1
            write(1,10) mi(rep),ii,ii,ii,ii,ii
            else
            write(1,10) mi(rep),(ie(mi(rep),i),i=1,5)
            endif
            write(1,2)
            vlabl2 ='    MODL        NLAY'
            write(1,2) vlabl2
            write(1,10) modl(rep),nlay(rep)
  750       mj = mj + 1
            if (nel .eq. 0) goto 800
            if (mj .lt. mi(rep)) goto 750
            if (mj .eq. nel) goto 800
            if (modl(rep) .le. 0) goto 700
            mj = mj + (nlay(rep) * modl(rep)) - 1
            if (mj .lt. nel) goto 700
  800       continue
         endif
      endif
      write (1,2)
!*
!*
      label ='*** DATA SET 13: POROUS MEDIA PROPERTIES'
      VLABL =' PROPERTIES DEPENDENT ON MATERIAL TYPE '
      write(1,2) label
      write(1,2) vlabl
      VLABL = 'MOL-DIF-CO   DENSITY LONG-DISP TRANS-DSP  POROSITY  MATER&
     &IAL'
      write(1,2) vlabl
      do 850 i=1,nmat
         write(1,45) (prop(j,i),j=1,5),i
  850 continue
      write(1,2)
      VLABL = 'PROPERTIES DEPENDENT ON MATERIAL AND CONTAMINANT '
      write(1,2) vlabl
      VLABL = '          KD    LIQ DEGR    SOL DEGR  CONTAMINANT '&
     &      //' MATERIAL'
      write(1,2) vlabl
      do i = 1,nmat
         do j = 1,niso
         write(1,47) (prop1(i,k,j),k=1,3),cnam(j),i
         enddo
         write(1,2)
      enddo
!*
      label ='*** DATA SET 14: MATERIAL TYPE CORRECTION'
      vlabl ='     NI        NSEQ         NAD       NITYP      NTYPAD'
      write(1,2) label
      if (ncm .le. 0) then
         write(1,2) nlabl
      else
      write(1,2) vlabl
         i = 0
  900      i = i + 1
           write(1,10) mi14(i),nseq14(i),mad(i),mityp(i),mtypad(i)
         if (mi14(i) .ne. 0) goto 900
      endif
      write(1,2)
!*
      label ='*** DATA SET 15: INITIAL CONDITIONS'
      write(1,2) label
      vlabl ='          NI    NNODE      NAD    RNI        RAD        RR&
     &D'
       write(1,2) vlabl
       do i = 1,niso
       vlabl = cnam(i)
       write(1,2) vlabl
       j = 0
  925         j = j + 1
              write(1,41) ni15(i,j),nnode(i,j),nad15(i,j),rni(i,j),&
                         rad(i,j),rrd(i,j)
            if (ni15(i,j) .ne. 0) goto 925
        enddo
        write(1,2)
!*
      label ='*** DATA SET 16: CAUCHY BOUNDARY CONDITIONS'
      write (1,2) label
      icheck = 0
      do iso=1,niso
          icheck = icheck + ncnp(iso)
      enddo
      if (icheck .le. 0) then
         write(1,2) nlabl
         write(1,2)
      else
         do 1050 k=1,niso
            if (ncnp(k) .eq. 0 .or. nces(k) .eq. 0) go to 1050
      LABEL = CNAM(K)//' CAUCHY BOUNDARY CONDITIONS: TIME UNIT IN YEARS'
            write(1,2) label
            WRITE(1,50) ('TIME','FLUX',I=1,3)
            do 1000 i=1,ncpr(k)
               write(1,60) (tcvbf(k,j,i),cvbf(k,j,i),j=1,ncdp(k))
 1000       continue
            write(1,2)
         VLABL ='    CAUCHY GLOBAL NODE NUMBERS '
         write(1,2) vlabl
         write(1,48) (NPCB(i,k),i=1,ncnp(k))
        write(1,2)
            vlabl ='      NI        NSEQ        NIAD        NITYP'//&
     &'     NTYPAD - CAUCHY NODE TYPES'
            write(1,2) vlabl
            j = 0
 1025         j = j + 1
              write(1,10) mi16(k,j),nseq16(k,j),miad16(k,j),&
              mityp16(k,j),mtypad16(k,j)
            if (mi16(k,j) .ne. 0) goto 1025
            write(1,2)
      vlabl = '   MI    NSEQ       M     IS1     IS2    MIAD     MAD   I&
     &S1AD   IS2AD'
       write(1,2) vlabl
         j = 0
 1200      j = j + 1
        write(1,15) mi_v16(k,j),nseq_v16(k,j),m_v16(k,j),is1_v16(k,j),&
       is2_v16(k,j),miad_v16(k,j),mad_v16(k,j),is1ad(k,j),is2ad(k,j)
         if (mi_v16(k,j) .ne. 0) goto 1200
      write(1,2)
 1050    continue
      endif
      label ='*** DATA SET 16N: NEUMANN BOUNDARY CONDITIONS'
      write (1,2) label
      icheck = 0
      do iso=1,niso
          icheck = icheck + nnnp(iso)
      enddo
      if (icheck .le. 0) then
         write(1,2) nlabl
         write(1,2)
      else
         do 1055 k=1,niso
            if (nnnp(k) .eq. 0 .or. nnes(k) .eq. 0) go to 1055
      LABEL = CNAM(K)//' NEUMANN BOUNDARY CONDITIONS: TIME UNIT IN &
     &YEARS'
            write(1,2) label
            write(1,50) ('TIME','FLUX',i=1,3)
            do 1052 i=1,nnpr(k)
               write(1,60) (tcvbfn(k,j,i),cvbfn(k,j,i),j=1,nndp(k))
 1052       continue
            write(1,2)
         VLABL ='     NEUMANN GLOBAL NODE NUMBERS '
         write(1,2) vlabl
         write(1,48) (NPNB(i,k),i=1,nnnp(k))
         write(1,2)
            vlabl ='      NI        NSEQ         NAD       NITYP'//&
     &'     NTYPAD - NEUMANN NODE TYPES'
            write(1,2) vlabl
            j = 0
 1041         j = j + 1
              write(1,10) mi16n(k,j),nseq16n(k,j),miad16n(k,j),&
              mityp16n(k,j),mtpad16n(k,j)
            if (mi16n(k,j) .ne. 0) goto 1041
            write(1,2)
      vlabl = '   MI    NSEQ       M     IS1     IS2    MIAD     MAD   I&
     &S1AD   IS2AD'
        write(1,2) vlabl
         j = 0
 1241      j = j + 1
           write(1,15) mi_v16n(k,j),nseq_16n(k,j),m_v16n(k,j),&
           is1_v16n(k,j),is2_v16n(k,j),miad_16n(k,j),mad_v16n(k,j),&
           is1adn(k,j),is2adn(k,j)
         if (mi_v16n(k,j) .ne. 0) goto 1241
         write(1,2)
 1055    continue
      endif
!*
      label ='*** DATA SET 17: DIRICHLET BOUNDARY CONDITIONS'
      write (1,2) label
      icheck = 0
      do iso = 1,niso
          icheck = icheck + ndnp(iso)
      end do
      if (icheck .eq. 0) then
         write(1,2) nlabl
         write(1,2)
      else
         do 1350 k=1,niso
            if ( ndnp(k) .ne. 0) then
            VLABL=CNAM(K)//'DIRICHLET BOUNDARY CONDITION TABLE'
            write(1,2) vlabl
            write(1,50) ('TIME','CONC',i=1,3)
            do 1300 i=1,ndpr(k)
               write(1,60) (tcdbf(k,j,i),cdbf(k,j,i),j=1,nddp(k))
 1300       continue
            write(1,2)
         VLABL ='           ** GLOBAL DIRICHLET NODES **'
         write(1,2) vlabl
         write(1,55) (npdb(j,k),j=1,ndnp(k))
         write(1,2)
            vlabl ='      NI        NSEQ        NIAD       NITYP        NTYPA&
     &D'
            write(1,2) vlabl
            j = 0
 1325         j = j + 1
              write(1,10) ni17(k,j),nseq17(k,j),niad17(k,j),&
              nityp17(k,j),ntypad17(k,j)
            if (ni17(k,j) .ne. 0) goto 1325
            write(1,2)
         endif
 1350    continue
      endif
!*
!*---- Data Set 18: External Volumetric Sources
      label ='*** DATA SET 18: EXTERNAL VOLUMETRIC SOURCES'
!*
      icheck = 0
      do iso = 1,niso
      icheck =icheck + nsel(iso)
      enddo
      if (icheck .eq. 0) then
      write (1,2) label
      write (1,2) nlabl
      write (1,2)
      else
      write (1,2) label
      vlabl = ' '
      write (1,2) vlabl
!*----- Write source strength profiles
         do 1001 k=1,niso
             if(nsel(k) .eq. 0) go to 1001
!*----- Write source type assigned to each element
      vlabl = 'GLOBAL ELEMENT NUMBER(S)'//'  '//cnam(k)
      write (1,2) vlabl
             write (1,51)  (LES(j,k),j=1,nsel(k))
             write (1,2)
      vlabl = 'SOURCE STRENGTH VS. TIME PROFILE(S)'//'  '//cnam(k)
             write (1,2) vlabl
             do 951 i=1,nspr(k)
               write (1,65) (tsosf(k,j,i),sosf(k,j,i),j=1,nsdp(k))
             write (1,2)
  951     continue
!*----- Write volumetric source element assigned to a specific profile
      vlabl = '   NI NSEQ  NAD NITYP NTYPAD'//'  '//cnam(k)
             write (1,2) vlabl
             j = 0
  976         j = j + 1
          write (1,55) ni18(k,j),nseq18(k,j),niad18(k,j),nityp18(k,j),&
             ntypad18(k,j)
             if (ni18(k,j) .ne. 0) goto 976
!c     if(k.ne.niso) then
      write (1,2)
!c     endif
 1001    continue
      endif
!*
!*-----Data Set 19: External Point (Well) Sources
      label ='*** DATA SET 19: EXTERNAL POINT (WELL) SOURCES'
!*
      icheck = 0
      do iso = 1,niso
      icheck =icheck + nwnp(iso)
      enddo
      if (icheck .eq. 0) then
      write (1,2) label
      write (1,2) nlabl
      write (1,2)
      else
      write (1,2) label
      vlabl = '  '
      write (1,2) vlabl
!*----- Write source strength profiles
         do 1002 k=1,niso
             if(nwnp(k) .eq. 0) go to 1002
!*----- Write source type assigned to each element
      vlabl = 'GLOBAL ELEMENT NUMBER(S)'//'  '//cnam(k)
      write (1,2) vlabl
             write (1,51)  (NPW(j,k),j=1,nwnp(k))
             write (1,2)
      vlabl = 'SOURCE STRENGTH VS. TIME PROFILE(S)'//'  '//cnam(k)
      write (1,2) vlabl
             do 952 i=1,nwpr(k)
               write (1,65) (twssf(k,j,i),wssf(k,j,i),j=1,nwdp(k))
             write (1,2)
  952     continue
!*----- Write point source element assigned to a specific profile
      vlabl = '   NI NSEQ  NAD NITYP NTYPAD'//'  '//cnam(k)
             write (1,2) vlabl
             j = 0
  977         j = j + 1
         write (1,55) ni19(k,j),nseq19(k,j),niad19(k,j),nityp19(k,j),&
             ntypad19(k,j)
             if (ni19(k,j) .ne. 0) goto 977
!c     if(k.ne.niso) then
      write (1,2)
!c     endif
 1002    continue
      endif
      label ='*** DATA SET 20: WASTE CONTAINER PARAMETERS'
      write(1,2) label
      if (ncon .eq. 0) then
         write(1,2) nlabl
      else
         do 1400 j=1,nctype
        vlabl ='    THICK       PITN        PITK        AREA        ASC'&
     &//'ALE      PITS'
            write(1,2) vlabl
            write(1,16) thick(j),pitn(j),pitk(j),area(j),ascale(j),&
            pits(j)
            write(1,2)
            vlabl ='    GRATE        CLAY        SPH              IAER'
            write(1,2) vlabl
            write(1,17) grate(j),clay(j),sph(j),iaer(j)
         write(1,2)
 1400    continue
         vlabl ='    GLOBAL ELEMENT NUMBERS FOR CONTAINER LOCATIONS  '
         write(1,2) vlabl
         write(1,48) (nelcon(j),j=1,ncon)
         write(1,2)
      vlabl ='     NI        NSEQ         NAD       NITYP      NTYPAD'
         write(1,2) vlabl
         j = 0
 1450      j = j + 1
           write(1,10) mi20(j),nseq20(j),mad20(j),mityp20(j),mtypad20(j)
         if (mi20(j) .ne. 0) goto 1450
      endif
      write(1,2)
!*
      label ='*** DATA SET 21: WASTE FORM AND LEACHING INFORMATION'
      write(1,2) label
      if (ncon .eq. 0) then
         write(1,2) nlabl
         write(1,2)
      else
         do 1550 j=1,nwtype
      jj=j/10
      if(jj.eq.0) jj=-16
      l=mod(j,10)
            VLABL2 = ' WASTE FORM TYPE '//CHAR(JJ+48)//CHAR(L+48)
            write(1,2) vlabl2
            do 1500 i=1,niso
      vlabl ='    SFRACT      PFRACT      BFRACT      DEFF        DISOL'&
     &//'     '//cnam(i)
               write(1,2) vlabl
               write(1,16) sfract(j,i),pfract(j,i),bfract(j,i),&
               deff(j,i),disol(j,i)
               write(1,2)
        vlabl2 ='    PARTK0      PARTKI      PARTD'//'     '//cnam(i)
               write(1,2) vlabl2
               write(1,16) partk0(j,i),partki(j,i),partd(j,i)
               write(1,2)
 1500       continue
            vlabl2 ='    POREL       VOLWF      VRATIO'
            write(1,2) vlabl2
            write(1,16) porel(j),volwf(j),vratio(j)
            write(1,2)
 1550    continue
!*
         label ='*** DATA SET 21b: INITIAL MASSES OF EACH WASTE ELEMENT'
         write(1,2) label
         VLABL='               INITIAL INVENTORY IN EACH CONTAINER'
         do 1600 j=1,niso
            write(1,70) cnam(j)
            write(1,2) vlabl
            write(1,65) (wtinit(i,j),i=1,ncon)
            write(1,2)
 1600    continue
!*
         label ='*** DATA SET 21c: READ WASTE TYPE ASSIGNMENTS'
      vlabl ='     NI        NSEQ         NAD       NITYP      NTYPAD'
         write(1,2) label
         write(1,2) vlabl
         j = 0
 1650      j = j + 1
           write(1,10) mi21(j),nseq21(j),mad21(j),mityp21(j),&
           mtypad21(j)
         if (mi21(j) .ne. 0) goto 1650
         write(1,2)
      endif
!*
      label ='*** DATA SET 22a: PRINTER AND TRACE VARIABLE CONTROLS'
      write(1,2) label
      LABEL ='        PRINT INDICATOR FLAGS'
      write(1,2) label
      vlabl ='KPR0,KPR(I); I = 1 to NTI'
      write(1,2) vlabl
      nti1 = nti
      if(nti .gt. 1000) nti1 = 1000
      write(1,80) kpr0,(kpr(i),i=1,nti1)
      write(1,2)
!*
      label ='*** DATA SET 22b: CHEMICAL OUTPUT AND CHEMICAL PROPERTY &
     &INDICATOR'
        write(1,2) label
        vlabl = '     DATA SET NOT USED IN BLT-MS '
      write(1,2) vlabl
      write(1,2)
!*
      label = '*** DATA SET 22c: CONCENTRATION AND FLUX TRACE FILE &
     &INDICATORS'
      write(1,2) label
      vlabl = '    NTRC        NTRF      NSTPTR'
      write(1,2) vlabl
      write(1,10) ntrc,ntrf,nstptr
      write(1,2)
      if(ntrc .gt. 0) then
      VLABL = '     GLOBAL NODE NUMBERS FOR CONCENTRATION TRACES'
      write(1,2) vlabl
      write(1,10) (lctrc(i),i=1,ntrc)
      write(1,2)
      endif
      if(ntrf .gt. 0) then
      VLABL = '     GLOBAL NODE NUMBERS FOR FLUX TRACES '
      write(1,2) vlabl
      write(1,10) (lftrc(i),i=1,ntrf)
      write(1,2)
      endif
!*
      label ='*** DATA SET 23: VELOCITY AND MOISTURE CONTENT'
      vlabl ='     NI     NSEQ      NAD    VXNI       VZNI       VXAD'&
     &//'       VZAD'
      write(1,2) label
      if (kvi .gt. 0) then
         write(1,2) nlabl
      else
         write(1,2) vlabl
         j = 0
 1750      j = j + 1
           write(1,40)  ni23(j),nseq23(j),nad23(j),vxni23(j),vzni23(j),&
           vxad23(j),vzad23(j)
         if (ni23(j) .ne. 0) goto 1750
         write(1,2)
         vlabl ='     NI     NSEQ      NAD    THNI       THNIAD   '
         write(1,2) vlabl
         j = 0
 1800      j = j + 1
           write(1,40)  ni23b(j),nseq23b(j),nad23b(j),thni23(j),&
           thniad23(j)
         if (ni23b(j) .ne. 0) goto 1800
      endif
       write(1,2)
       write(1,2)
       write(1,2)
       write(1,2)
       write(1,2)
       write(1,2)
      return
!*
!*------ Format Section
!*
    2 FORMAT (A79)
    6 FORMAT (I5,A70)
   10 FORMAT (6(3X,I4,5X))
   11 FORMAT (1X,3(3X,I4,5X),5X,A10)
   12 FORMAT (1X,4(3X,I4,5X),2X,A10)
   15 FORMAT (9(I5,3X))
   16 FORMAT (2X,6(1X,1PE10.3,1X))
   17 FORMAT (2x,3(1X,1PE10.3,1X),11X,I1)
   18 FORMAT (2X,1PE10.3)
   20 FORMAT (4X,A10,2X,1PE10.3,2x,1PE10.3,2X,1PE10.3)
   25 FORMAT (2X,I3)
   30 FORMAT (A13)
   35 FORMAT (A)
   40 FORMAT (1X,3(3X,I3,3X),4(1PE10.3,1X))
   41 FORMAT (6X,3(3X,I3,3X),4(1PE10.3,1X))
   45 FORMAT (5(1PE10.3),2x,i5)
   46 FORMAT (8(1PE9.2,1x))
   47 FORMAT (3(2x,1pe10.3),5x,a10,i5)
   48 FORMAT (14I5)
   50 FORMAT (' ',3x,3(A4,7X,A4,7x))
   51 FORMAT (9(I5,4X))
   55 FORMAT (14I5)
   60 FORMAT (6(1PD10.3,1x))
   65 FORMAT (1p,8e10.3)
   70 FORMAT (10X,a10)
   80 FORMAT (80(I1))
   85 FORMAT (600(I3,1X))
   95 FORMAT (1X,20(1PE10.3,1X))
!
!*
      END
!*
!* End of Subroutine MakeDeck
!*----------------------------------------------------------------------------*
!*
!*----------------------------------------------------------------------------*
!* Begining of Subroutine rnode()
!*
      SUBROUTINE rnode()
!*
      character*79 label
      common /data12/ ni(2,600),nseq(2,600),nad(2,600),&
              xni(600),xad(600),xrd(600),&
              zni(600),zad(600),zrd(600),&
              nelx,nelz,nreg,&
              imin(50),jmin(50),imax(50),jmax(50),&
              x1(4,50),x2(4,50),mi(550),ie(550,5),&
              nlay(550), modl(550),widfac
!
      data KIND2 /4/
!*
      read(2,*) nelx,nelz,nreg
         read(2,10) label
      do 100 i=1,nreg
         read(2,10) label
         read(2,*) imin(i),jmin(i),imax(i),jmax(i)
         read(2,10) label
         read(2,*) (x1(j,i),x2(j,i), j=1,KIND2)
 100  continue
      label = label
   10 FORMAT(A79)
!*
      return
!*
!*-----Format Section
!*
   5  FORMAT(A79)
      END
!*
!* End of Subroutine rnode()
!*-----------------------------------------------------------------------------*
!*
!*-----------------------------------------------------------------------------*
