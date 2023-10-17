!Last Changed by Prateek and Farhan
!Tue Jan 9 16:16:15 IST 2015
program boundarylayer
!
! This Program does the 3D simulation of Turbulent Flows in a Rectangular Cuboid
! This is the parallel version of Program 
! Developed and Written by Fuaad P A 
! Aligarh Muslim University 
! For DST SERC Project 'Controlling Near Wall turbulence using Buoyancy Forces'
! Under the Guidance of Prof Mirza Faisal Baig 
! This is the 2D-parallel code.

! 
  use mpi

  implicit none
 
! number types
!
INTEGER, parameter :: nstart = 50
!
INTEGER             :: nstep                       ! current number of time integrations 
INTEGER             :: maxstep      =30000000       ! maximum number of time integrations 
!INTEGER             :: snapshot     =50000
DOUBLE PRECISION    :: errorstop    =1.0d-6        ! Convergence Criteria of Transport equations
!INTEGER             :: bckp         =10000           ! display frequency 
INTEGER             :: disp         =1           ! display frequency 
INTEGER             :: disp1        =1            ! display frequency of divergence and massflow.
INTEGER             :: istat        =1            ! stat frequency 
INTEGER             :: nfile
INTEGER             :: nsnapshot 
INTEGER             :: snapshot     =5                     ! display frequency 
DOUBLE PRECISION    :: time                        ! current non-dimensional time
DOUBLE PRECISION    :: dt          =1.0d-4         ! timestep
DOUBLE PRECISION    :: dXi,dx,dy,invdx,invdy,invdXi,onetwelfth
DOUBLE PRECISION, dimension(3) :: divgrad  ! Grad of Ustar components
DOUBLE PRECISION, dimension(:,:,:,:), allocatable    :: Ume,Umedist
DOUBLE PRECISION, dimension(:,:,:,:), allocatable    :: Ufluc
DOUBLE PRECISION    :: ke1,ke
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: Ut,Q,conv,Umrf  ! 3d Flow Field (t:i:j:k)
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: Udist      ! disturbances in  Flow Field (i,j,k,direction)
DOUBLE PRECISION, dimension(:,:,:), allocatable   :: PCorr,Press
DOUBLE PRECISION, dimension(:,:), allocatable   :: Utmean,fnlUtmean
INTEGER             :: kcnt,jt
!INTEGER  :: alpha=10, beta=5 , gama=1
!******************************************************************************************

DOUBLE PRECISION               :: Re  =150.92d0                                ! Reynolds Number 
DOUBLE PRECISION               :: Pran=0.71d0                                 ! Prandtl Number
DOUBLE PRECISION               :: Ri=0.0d0                                   ! Richardson Number 
DOUBLE PRECISION, dimension(3) :: forc                          ! forcing in 3 dimensions
DOUBLE PRECISION :: Pprod, Preq, Punc!,CF
DOUBLE PRECISION :: Pprodpercent, Preqpercent,Pnetpercent, h1,h2
DOUBLE PRECISION, PARAMETER :: UbN=15.730d0
DOUBLE PRECISION, parameter :: delta = 1.95d0
DOUBLE PRECISION, parameter :: delprmtr = 5.0d0
! Mesh Parameters
!
INTEGER  , dimension(3) :: n3d                                           ! grid points in 3 dimensions
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Aw,An,Ae,As,Ap,At,Ab,rAp  ! Stencil Parameters for Pressure
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Aww,Aee,Ass,Att,Abb   ! Stencil Parameters for Pressure
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Anw,Ane,Ans,Anp,Ant,Anb,rAnp
DOUBLE PRECISION, dimension(:,:,:), allocatable :: Avw,Ann,Avn,Ave,Avs,Avp,Avt,Avb! Stencil Parameters for Pressure
DOUBLE PRECISION, dimension(:),allocatable      :: x1,x2,x3,Jac,Jac2,Xi,App,Umr,Umr1    ! Physical Geometery (Cartesian)
DOUBLE PRECISION,allocatable,dimension(:)      :: tauxz,tmp2
DOUBLE PRECISION,allocatable,dimension(:)      :: fnltauxz
DOUBLE PRECISION,allocatable,dimension(:)      :: sumxz
DOUBLE PRECISION,allocatable,dimension(:)    :: bltkns
DOUBLE PRECISION,allocatable,dimension(:)      :: bltk,fnlbltk
!
! gmres internal Parameters
!
DOUBLE PRECISION,allocatable,dimension(:,:,:,:):: v
DOUBLE PRECISION,allocatable,dimension(:,:,:)  :: r,rhat
DOUBLE PRECISION,allocatable,dimension(:,:)    :: h
DOUBLE PRECISION,allocatable,dimension(:)      :: g,c,s,y,mag
DOUBLE PRECISION,allocatable, dimension(:) :: cf  
DOUBLE PRECISION  :: flux
DOUBLE PRECISION                               :: rho,hik,hipk,nu,gk,gkp,w,errtol
DOUBLE PRECISION                               :: loc_rho,tempgmres
INTEGER                            :: kmax,kit,mmax,m,sgn
INTEGER                            :: my_rank,proc,source,ierr,loc_nr,loc_nt,my
INTEGER                            :: status(mpi_status_size),bn,bnt,en,ent,oneplane,oneplanet,twoplanes,fourvar,sevenvars,oneplan,oneplaned

DOUBLE PRECISION,dimension(3)                  :: df,db,d2b,d2f        ! update locals       ! update locals
DOUBLE PRECISION                               :: Pi,t,div,adv,fder,fder1,fder2
DOUBLE PRECISION                               :: sumav,sumnor,res,tmp,tmp1,diver
DOUBLE PRECISION                               :: tf,Pe,term1,term2,Nuss,Nuss1,Nuss2,b
DOUBLE PRECISION, dimension(3)                 :: divuu    ! Grad of Ustar components
DOUBLE PRECISION, dimension(3)                 :: gradp    ! Grad of Ustar components
INTEGER :: i, j, k,it,istart,iend,npln,jstart,jend                     ! loop variables
INTEGER :: jcnt,im1,im2,ip1,ip2,jm1,jm2,jp1,jp2,im3,ip3,jm3,jp3
INTEGER :: cnt,sbn,loc_nx,loc_nl,itemp,loc_n
DOUBLE PRECISION,dimension(5)    :: UbDNS,Ubp1                   ! Random 5 points where mass flow has to calculated
integer  ,dimension(5)           :: plane =(/5,20,40,60,85/)! These are those points
DOUBLE PRECISION                 :: dAv,Av,a1
DOUBLE PRECISION                 :: A,dA,su,sv,sw,sT,sp
DOUBLE PRECISION,allocatable,dimension(:,:) :: uusum,vvsum,wwsum,ttsum,uvsum,uwsum,vwsum,usum,vsum,wsum,tsum,twsum,Psum  
DOUBLE PRECISION,allocatable,dimension(:,:) :: uusum1,vvsum1,wwsum1,ttsum1,uvsum1,uwsum1,vwsum1,usum1,vsum1,wsum1,tsum1,twsum1,Psum1
DOUBLE PRECISION,allocatable,dimension(:,:) :: uusum2,vvsum2,wwsum2,ttsum2,uvsum2,uwsum2,vwsum2,usum2,vsum2,wsum2,tsum2,twsum2,Psum2
DOUBLE PRECISION,allocatable,dimension(:,:) :: Ruu , Rvv , Rww ,Rtt , Ruv , Ruw , Rvw ,Rtw 
DOUBLE PRECISION :: tmpusum,tmpvsum,tmpwsum,tmpTsum,tmpuusum,tmpvvsum,tmpwwsum,tmpttsum,tmpuvsum,tmpuwsum,tmpvwsum,tmptwsum,tmpPsum
  CHARACTER(LEN=25)                :: filname
integer :: ds,ids,nstat,procz,procy,myr,myp
!DOUBLE PRECISION,dimension(2) :: tauxz                  ! Tau container
DOUBLE PRECISION,allocatable, dimension(:,:) :: mfb1,mfb2   ! local vars

DOUBLE PRECISION :: N,nstatinv
!DOUBLE PRECISION,dimension(2,4)  ::Ufluc
DOUBLE PRECISION,allocatable,dimension(:,:,:,:) ::RUUcorr
DOUBLE PRECISION,allocatable,dimension(:,:) ::Umean
DOUBLE PRECISION        :: Pbar1,Pbar,invPi,ftime 
DOUBLE PRECISION, PARAMETER :: Ue=19.798d0!23.15d0
integer :: counter,const
integer :: cum=5
integer :: cont
!***********************************************************************************************************************************
!DOUBLE PRECISION,allocatable,dimension(:) :: utau,lamda
DOUBLE PRECISION :: beta=10.0d0, gama=10.0d0
DOUBLE PRECISION        :: lamda,utau1,utau2,kp,kk,eps,lamda03,lamda1 
DOUBLE PRECISION, dimension(:,:,:,:), allocatable    :: Uflucinn,Uflucout
DOUBLE PRECISION, dimension(:,:), allocatable    :: Utmeaninn,Vtmeaninn,Wtmeaninn,Utmeanout,Vtmeanout,Wtmeanout,Wt,Ume1
DOUBLE PRECISION, dimension(:), allocatable      :: theta,deltastar,Ht,zrecy,eta,Reth
INTEGER, DIMENSION(:), ALLOCATABLE      	:: new
DOUBLE PRECISION, dimension(:,:,:,:), allocatable :: Utold  ! 3d Flow Field (t:i:j:k)
integer :: threevar,ksep,cnt2,prompt

!***********************************************************************************************************************************
Pi=4.d0*datan(1.d0)
onetwelfth =1.0d0/12.0d0
invPi = 1.0d0/(4.0d0*Pi)

forc(1) =0.0d0
forc(2) =0.0d0
forc(3) =0.0d0

errtol=1.0d-5
kmax=50
mmax=1000
w = 1.835

!WRITE(*,*),'Enter 1 to restart from 3D.tp or Enter 2 to start from t=0'
!READ*,prompt
prompt=1
SELECT CASE(prompt)
	CASE(1)
		OPEN (UNIT=12,FILE='cnt',STATUS='unknown')
		READ (12,*)n3d(3),cnt,proc
		CLOSE(12)

		OPEN (UNIT=23,FILE='cnt2',STATUS='unknown')
		READ (23,*)cnt2
		CLOSE(23)

		!Reading Tecplot ASCII FORMAT Directly
		OPEN (UNIT=12,FILE='3D.tp',STATUS='unknown')
		READ (12,01)time,nstep
		READ (12,*)
		READ (12,02)n3d(1),n3d(2),n3d(3)
	CASE(2)
		cnt=0
		cnt2=0
		nstep=0
		time=0.0
		n3d(1)=90
		n3d(2)=60
		n3d(3)=50
	CASE DEFAULT
		WRITE(*,*),'Invalid Input'
END SELECT

     allocate( RUUcorr(2, 4, n3d(3),n3d(3)) )
     allocate( Umean(n3d(3),4))
     allocate( Ume(n3d(1),n3d(2),n3d(3),3))
     allocate( Umrf(n3d(1),n3d(2),n3d(3),3))
     allocate( Umedist(n3d(1),n3d(2),n3d(3),3))
     allocate( Ufluc(n3d(1),n3d(2),n3d(3),3))
     allocate( Ut(n3d(1), n3d(2), n3d(3),4) )
     allocate( Udist(n3d(1), n3d(2), n3d(3),3) )
     allocate( x1(n3d(1)))
     allocate( x2(n3d(2))) 
     allocate( x3(n3d(3)))
     allocate( Umr(n3d(3)))
     allocate( Umr1(n3d(3)))
     allocate( Xi(n3d(3)))
     allocate(Jac(n3d(3)) )
     allocate(Jac2(n3d(3)) )
     allocate(App(n3d(3)) )
     allocate( Conv( n3d(1), n3d(2), n3d(3),4) )
     allocate( Q ( n3d(1), n3d(2), n3d(3),8) )
     allocate( PCorr(n3d(1), n3d(2), n3d(3)) )
     allocate( Press(n3d(1), n3d(2), n3d(3)) )
     allocate( cf(n3d(1) ) )
     allocate( tauxz(n3d(1) ) )   
     allocate( tmp2(n3d(3) ) )   
     allocate( fnltauxz(n3d(1) ) )   
     allocate( sumxz(n3d(1) ) )   
     allocate( bltkns(n3d(1)))   
     allocate( bltk(n3d(1)) )  
     allocate( fnlbltk(n3d(1) ))   
     allocate( Utmean(n3d(1),n3d(3) )) 
     allocate( fnlUtmean(n3d(1),n3d(3) ))  
!*********************************************************************
     allocate( Uflucinn(n3d(1),n3d(2),n3d(3),3))
     allocate( Uflucout(n3d(1),n3d(2),n3d(3),3))
     allocate( Utold(n3d(1), n3d(2), n3d(3),3) )
     allocate( Utmeaninn(n3d(1),n3d(3) )) 
     allocate( Vtmeaninn(n3d(1),n3d(3) )) 
     allocate( Wtmeaninn(n3d(1),n3d(3) )) 
     allocate( Utmeanout(n3d(1),n3d(3) )) 
     allocate( Vtmeanout(n3d(1),n3d(3) )) 
     allocate( Wtmeanout(n3d(1),n3d(3) )) 
     allocate( Ume1(n3d(1),n3d(3) )) 
     allocate( Wt(n3d(1),n3d(3) )) 
     allocate( zrecy(n3d(3) )) 
     allocate( new(n3d(3) )) 
     allocate( eta(n3d(3) )) 
     allocate( theta(n3d(1) )) 
     allocate( Reth(n3d(1) )) 
     allocate( Ht(n3d(1) )) 
     allocate( deltastar(n3d(1) )) 

!*********************************************************************
     allocate(Anw(n3d(1), n3d(2), n3d(3)) )
     allocate(Ann(n3d(1), n3d(2), n3d(3)) )
     allocate(Ane(n3d(1), n3d(2), n3d(3)) )
     allocate(Ans(n3d(1), n3d(2), n3d(3)) )
     allocate(Anb(n3d(1), n3d(2), n3d(3)) )
     allocate(Ant(n3d(1), n3d(2), n3d(3)) )
     allocate(Anp(n3d(1), n3d(2), n3d(3)) )
     allocate( rAnp(n3d(1),n3d(2),n3d(3)) )

     allocate(Aw(n3d(1), n3d(2), n3d(3)) )
     allocate(An(n3d(1), n3d(2), n3d(3)) )
     allocate(Ae(n3d(1), n3d(2), n3d(3)) )
     allocate(As(n3d(1), n3d(2), n3d(3)) )
     allocate(Ab(n3d(1), n3d(2), n3d(3)) )
     allocate(At(n3d(1), n3d(2), n3d(3)) )
     allocate(Ap(n3d(1), n3d(2), n3d(3)) )
     allocate(RAp(n3d(1), n3d(2), n3d(3)) )

     allocate(Avw(n3d(1), n3d(2), n3d(3)) )
     allocate(Avn(n3d(1), n3d(2), n3d(3)) )
     allocate(Ave(n3d(1), n3d(2), n3d(3)) )
     allocate(Avs(n3d(1), n3d(2), n3d(3)) )
     allocate(Avb(n3d(1), n3d(2), n3d(3)) )
     allocate(Avt(n3d(1), n3d(2), n3d(3)) )
     allocate(Avp(n3d(1), n3d(2), n3d(3)) )
     allocate(Aee(n3d(1), n3d(2), n3d(3)) )
     allocate(Aww(n3d(1), n3d(2), n3d(3)) )
     allocate(Att(n3d(1), n3d(2), n3d(3)) )
     allocate(Abb(n3d(1), n3d(2), n3d(3)) )
     allocate(Ass(n3d(1), n3d(2), n3d(3)) )

     allocate( mfb1(n3d(1),n3d(2)) )
     allocate( mfb2(n3d(1),n3d(2)) )

!These variables are internal to GMRES 

        allocate(v(n3d(1),n3d(2),n3d(3),nstart+1))
        allocate(r(n3d(1),n3d(2),n3d(3)))
        allocate(rhat(n3d(1),n3d(2),n3d(3)))
        allocate(h(nstart+1,nstart+1))
        allocate(g(nstart+1))
        allocate(c(nstart+1))
        allocate(s(nstart+1))
        allocate(y(nstart+1))
        allocate(mag(nstart+1))

!These variables are internal to statistics 
     allocate( uusum(n3d(1),n3d(3)))
     allocate( vvsum(n3d(1),n3d(3)))
     allocate( wwsum(n3d(1),n3d(3)))
     allocate( ttsum(n3d(1),n3d(3)))
     allocate( uvsum(n3d(1),n3d(3)))
     allocate( uwsum(n3d(1),n3d(3)))
     allocate( vwsum(n3d(1),n3d(3)))
     allocate( usum(n3d(1),n3d(3)))
     allocate( vsum(n3d(1),n3d(3)))
     allocate( wsum(n3d(1),n3d(3)))
     allocate( tsum(n3d(1),n3d(3)))
     allocate( Psum(n3d(1),n3d(3)))
     allocate( twsum(n3d(1),n3d(3)))

     allocate( uusum1(n3d(1),n3d(3)))
     allocate( vvsum1(n3d(1),n3d(3)))
     allocate( wwsum1(n3d(1),n3d(3)))
     allocate( ttsum1(n3d(1),n3d(3)))
     allocate( uvsum1(n3d(1),n3d(3)))
     allocate( uwsum1(n3d(1),n3d(3)))
     allocate( vwsum1(n3d(1),n3d(3)))
     allocate( usum1(n3d(1),n3d(3)))
     allocate( vsum1(n3d(1),n3d(3)))
     allocate( wsum1(n3d(1),n3d(3)))
     allocate( tsum1(n3d(1),n3d(3)))
     allocate( Psum1(n3d(1),n3d(3)))
     allocate( twsum1(n3d(1),n3d(3)))
     
     allocate( uusum2(n3d(1),n3d(3)))
     allocate( vvsum2(n3d(1),n3d(3)))
     allocate( wwsum2(n3d(1),n3d(3)))
     allocate( ttsum2(n3d(1),n3d(3)))
     allocate( uvsum2(n3d(1),n3d(3)))
     allocate( uwsum2(n3d(1),n3d(3)))
     allocate( vwsum2(n3d(1),n3d(3)))
     allocate( usum2(n3d(1),n3d(3)))
     allocate( vsum2(n3d(1),n3d(3)))
     allocate( wsum2(n3d(1),n3d(3)))
     allocate( tsum2(n3d(1),n3d(3)))
     allocate( Psum2(n3d(1),n3d(3)))
     allocate( twsum2(n3d(1),n3d(3)))

     allocate(Ruu(n3d(1),n3d(3)))
     allocate(Rvv(n3d(1),n3d(3)))
     allocate(Rww(n3d(1),n3d(3)))
     allocate(Rtt(n3d(1),n3d(3)))
     allocate(Ruv(n3d(1),n3d(3)))
     allocate(Ruw(n3d(1),n3d(3)))
     allocate(Rvw(n3d(1),n3d(3)))
     allocate(Rtw(n3d(1),n3d(3)))

Ume(1:n3d(1),1:n3d(2),1:n3d(3),1:3)= 0.0d0
OPEN (UNIT=15,FILE='U_innerdns.dat') !Reading Inflow Velocity from File
DO k=1,37
	READ(15,98) Ume(1,1,k,1)
END DO
CLOSE(15)

98 FORMAT (24x,E18.8)

DO k=38,n3d(3)
	Ume(1,1,k,1) = Ume(1,1,37,1) !free-stream velocity to all higher points
END DO

DO i=1,n3d(1)
	DO j=1,n3d(2)
		DO k=1,n3d(3)
			Ume(i,j,k,1) = Ume(1,1,k,1) !Copy the array to all points
		END DO
	END DO
END DO

DO k=1,n3d(3)
	Ume1(1,k)=Ume(1,1,k,1)
END DO


SELECT CASE(prompt)
	CASE(1)
		DO i=1,n3d(1)
			READ(12,*)
			DO j=1,n3d(2)
				READ(12,*)
				DO k=1,n3d(3)
					read(12,97)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Press(i,j,k),usum2(i,k),vsum2(i,k),wsum2(i,k),Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Wt(i,k),new(k)
				END DO
			END DO
		END DO
		CLOSE(12)
		
		do k=1,n3d(3)
			do j=1,n3d(2)
				do i=1,n3d(1)
					Ut(i,j,k,1)= Ume(i,j,k,1)+(0.05d0*Ue*dsin(beta*x2(j))*dcos(gama*x3(k)))
					Ut(i,j,k,2)= 			(0.05d0*Ue*dsin(beta*x2(j))*dcos(gama*x3(k)))
					Ut(i,j,k,3)= 			(0.05d0*Ue*dsin(beta*x2(j))*dcos(gama*x3(k)))
				enddo
			enddo
		enddo
		
		Ufluc(i,j,k,1)=0.0d0
		Ufluc(i,j,k,2)=0.0d0
		Ufluc(i,j,k,3)=0.0d0
		
	CASE(2)
				!*****************************Initial Flow Field**************************
		do jcnt=1,3
			do k=1,n3d(3)
				do j=1,n3d(2)
					do i=1,n3d(1)
						Ut(i,j,k,1)= Ume(i,j,k,1)+(0.02d0*Ue*dsin(beta*x2(j))*dcos(gama*x3(k)))
						Ut(i,j,k,2)=              (0.02d0*Ue*dsin(beta*x2(j))*dcos(gama*x3(k)))
						Ut(i,j,k,3)=              (0.02d0*Ue*dsin(beta*x2(j))*dcos(gama*x3(k)))
					enddo
			 	enddo
		 	enddo
		enddo
		!*************************************************************************
END SELECT

!---------------------------Mesh Setup Starts---------------------------
x1(1) = 0.0d0
x2(1) = 0.0d0
Xi(1) =0.0d0
x3(1) =0.0d0
x1(n3d(1)) = 12.0d0
x2(n3d(2)) =  2.0d0
x3(n3d(3)) =  3.4d0

dx = x1(n3d(1))/dble(n3d(1)-1)
dy = x2(n3d(2))/dble(n3d(2)-1)
dXi= x3(n3d(3))/dble(n3d(3)-1)

invdx = 1.0d0/dx
invdy = 1.0d0/dy
invdXi= 1.0d0/dXi

DO i=1,n3d(1)-1
	x1(i+1)=x1(i)+dx
end do

DO i=1,n3d(2)-1
	x2(i+1)=x2(i)+dy
end do

DO k=1,n3d(3)-1
	Xi(k+1) =Xi(k) + dXi
END DO

!Mapping of Non-Uniform Points from Uniform Points
!Gives Symmetric spacing with respect to 1
DO k=1,n3d(3)
	x3(k) = 3.4d0 + 3.4d0*dTanh(delta*0.50*(Xi(k)-3.4))/(dTanh(1.7d0*delta))
	Jac(k) = 1.7d0*delta/(dTanh(1.7*delta)*(dCosh(0.5 *delta *(Xi(k)-3.4)))**2)
	Jac2(k) = -1.7d0*delta**2*dTanh(0.5*delta*(Xi(k)-3.4))/(dTanh(1.7*delta)*(dCosh(0.5*delta*(Xi(k)-3.4)))**2)
END DO
!---------------------------Mesh Setup Ends-------------------------------

!x0=32.0d0
!DO i= 2,n3d(1)
!	Rex(i)=(x1(i)+x0)*Re*Ue
!END DO

!Initial guess.
rhat = 0.0d0
h = 0.0d0
v= 0.0d0
c= 0.0d0
s= 0.0d0
g = 0.0d0
y = 0.0d0
cont=0
!nstep=0
!ftime =0
lamda=1.0d0

      dA = dy
      Av = x3(n3d(3))*(x2(n3d(2)-1)-x2(2)) ! Total vertical area 

!      A=(x2(n3d(2))-x2(2))*(x1(n3d(1))-x1(2)) ! Total planar area
       A=x2(n3d(2))-x2(2)


        call mpi_init(ierr)
        call mpi_comm_rank(mpi_comm_world,my_rank,ierr)
        call mpi_comm_size(mpi_comm_world,proc,ierr) 


        call mpi_type_contiguous(  n3d(1)*n3d(2),mpi_double_precision,oneplane,ierr)
        call mpi_type_contiguous(2*n3d(1)*n3d(2),mpi_double_precision,twoplanes,ierr)
        call MPI_TYPE_VECTOR(3,2*n3d(1)*n3d(2),n3d(1)*n3d(2)*n3d(3),mpi_double_precision,threevar,ierr)
                           !(variables,2planes,totalsizeOfVar,type,newtype,ierr)
        call mpi_type_commit(oneplane   ,ierr)
        call mpi_type_commit(twoplanes  ,ierr)
        call mpi_type_commit(threevar,ierr)



   loc_n = (n3d(3)-1)/proc
         bn= 2+(my_rank)*loc_n 
         en=bn+loc_n-1
        if (my_rank.eq.proc-1) en=n3d(3)-1
                               loc_n=en-bn+1
                               loc_nx=loc_n
        if (my_rank.eq.proc-1) loc_nx=n3d(3)-bn+1
        call mpi_barrier(mpi_comm_world,ierr)

        if (my_rank.eq.0)      CALL MPI_Recv(loc_nx, 1, MPI_INTEGER, proc-1, 50, MPI_COMM_WORLD, status, ierr)
        if (my_rank.eq.proc-1) CALL MPI_Send(loc_nx , 1, MPI_INTEGER, 0     , 50, MPI_COMM_WORLD, ierr)


   istart =bn
   iend   =en
if(my_rank.eq.0)        istart =bn-1
if (my_rank.eq.proc-1)  iend =en+1


                  db(1)=x3(n3d(3))-x3(n3d(3)-1)
                  db(2)=x3(n3d(3))-x3(n3d(3)-2)
                  db(3)=x3(n3d(3))-x3(n3d(3)-3)

                  df(1)=x3(2)-x3(1)
                  df(2)=x3(3)-x3(1)
                  df(3)=x3(4)-x3(1)


! Set up equations
        do k=bn,en
                do j=1,n3d(2)
                        do i=1,n3d(1)
      Anw(i,j,k)= 1.0d0/(dx*dx)
      Ann(i,j,k)= 1.0d0/(dy*dy)
      Ane(i,j,k)= 1.0d0/(dx*dx)
      Ans(i,j,k)= 1.0d0/(dy*dy)
      Anb(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))+(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Ant(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))-(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Anp(i,j,k)=-2.0d0*((1.0d0/(dx*dx)) +(1.0d0/(dy*dy))+(1.0/(dXi*dXi*Jac(k)*Jac(k))))
     rAnp(i,j,k)= 1.0d0/Anp(i,j,k)

                        enddo
                enddo
        enddo

    ! Set up equations
        do k=bn,en
                do j=1,n3d(2)
                        do i=1,n3d(1)
      Avw(i,j,k)=1.0d0/(dx*dx)
      Avn(i,j,k)=1.0d0/(dy*dy)
      Ave(i,j,k)=1.0d0/(dx*dx)
      Avs(i,j,k)=1.0d0/(dy*dy)
      Avb(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))+(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Avt(i,j,k)=(1.0d0/(dXi*dXi*Jac(k)*Jac(k)))-(Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k)))
      Avp(i,j,k)=-2.0d0*((1.0d0/(dx*dx)) +(1.0d0/(dy*dy))+(1.0/(dXi*dXi*Jac(k)*Jac(k))))

                        enddo
                enddo
        enddo    

! Boundary Conditions at the bottom Wall 
             if (my_rank.eq.0)      Ut(2:n3d(1)-1,2:n3d(2)-1,bn-1,4) =0.0d0 ! Bottom wall cold
             if (my_rank.eq.proc-1) Ut(2:n3d(1)-1,2:n3d(2)-1,en+1,4) =0.0d0 ! Top wall heated      

        do k=istart,iend
  Jac(k)= 1.0d0/Jac(k) ! This is made to make all divisions a multiplication
        enddo

  do while( nstep < maxstep)

 time = time+dt 
 nstep = nstep +1 
ftime = ftime+dt


!*************************************giving disturbances to inflow filed******************************!!

!write(*,*)proc,'proc'
!if(my_rank.eq.6)write(*,*)lamda,'out'
Utold(:,:,:,:)=Ut(:,:,:,:)

!eps=cnt2*(5.0d-5)

!if(eps.gt.1.0d0)then
!eps=1.0d0
!endif

if(nstep.gt.5) then

do k=bn,en
if(usum2(68,k).lt.(0.99d0*Ue)) then

!if(lamda.le.1.20d0)then
Utmeaninn(1,k)=Ume1(1,k)
!else
!Utmeaninn(1,k)=0.9d0*lamda*Ume(68,j,k,1)
!endif
!Utmeaninn(1,k)= eps*(usum2(68,new(k)))+(1-eps)*Ume1(1,k)
Vtmeaninn(1,k)=vsum2(68,new(k))
Wtmeaninn(1,k)=wsum2(68,new(k))

!if(lamda.le.1.20d0)then
Utmeanout(1,k)=Ume1(1,k)
!Utmeanout(1,k)=lamda*Ume(68,j,k,1)
!else
!Utmeanout(1,k)=0.9d0*lamda*Ume(68,j,k,1)
!endif
!Utmeanout(1,k)=eps*(usum2(68,new(k)))+(1-eps)*Ume1(1,k)
Vtmeanout(1,k)=vsum2(68,new(k))
Wtmeanout(1,k)=wsum2(68,new(k))

do j=2,n3d(2)-1
Uflucinn(1,j,k,1)=lamda*Ufluc(68,j,new(k),1)
Uflucinn(1,j,k,2)=lamda*Ufluc(68,j,new(k),2)
Uflucinn(1,j,k,3)=lamda*Ufluc(68,j,new(k),3)
enddo

do j=2,n3d(2)-1
Uflucout(1,j,k,1)=lamda*Ufluc(68,j,new(k),1)
Uflucout(1,j,k,2)=lamda*Ufluc(68,j,new(k),2)
Uflucout(1,j,k,3)=lamda*Ufluc(68,j,new(k),3)
enddo

!Wt(i,k)=0.5d0*(1+(dtanh((4.0d0*((x3(k)/bltkns(i))-0.2d0))/((1-2.0d0*0.2d0)*(x3(k)/bltkns(i))+0.2d0))/dtanh(4.0d0)))
do i=1,n3d(1)-1
do j=2,n3d(2)-1
Ut(1,j,k,1)=(Utmeaninn(1,k)+Uflucinn(1,j,k,1))*(1.0d0-Wt(i,k))+(Utmeanout(1,k)+Uflucout(1,j,k,1))*Wt(i,k)
Ut(1,j,k,2)=(Vtmeaninn(1,k)+Uflucinn(1,j,k,2))*(1.0d0-Wt(i,k))+(Vtmeanout(1,k)+Uflucout(1,j,k,2))*Wt(i,k)
Ut(1,j,k,3)=(Wtmeaninn(1,k)+Uflucinn(1,j,k,3))*(1.0d0-Wt(i,k))+(Wtmeanout(1,k)+Uflucout(1,j,k,3))*Wt(i,k)
enddo
enddo

else

do j=2,n3d(2)-1
!Ut(1,j,k,1)=Ut(68,j,k,1)
!Ut(1,j,k,2)=Ut(68,j,k,2)
!Ut(1,j,k,3)=Ut(68,j,k,3)
Ut(1,j,k,1)=usum2(68,k)
Ut(1,j,k,2)=vsum2(68,k)
Ut(1,j,k,3)=wsum2(68,k)
enddo


endif
enddo
!write(*,*)eps,'eps'

endif


    do jcnt=1,3
         do k=bn,en
               do j=2,n3d(2)-1
                        do i=2,n3d(1)-1
     
                       tf = Re 
              if(jcnt.eq.4) tf = Re*Pran
!***********************************************************************
                im1 =i-1
   !           if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
   !           if (i.eq.3) im2 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
!***********************************************************************

!      Find  Cell Peclet Number         
               sgn = sign(1,ceiling(Ut(i,j,k,1)))  
              Pe = sgn*tf*Ut(i,j,k,1)*dx ! taking only the positive quantity
 if(i.le.2.or.i.ge.(n3d(1)-1)) then
        divuu(1) = Ut(i,j,k,1)*(Ut(ip1,j,k,jcnt)-Ut(im1,j,k,jcnt))/(2.0d0*dx)
      elseif(Pe.gt.2.0d0 .and. Pe.lt.50.0d0) then
term1 = (Abs(Ut(i,j,k,1)))*(Ut(ip2,j,k,jcnt)-(4*(Ut(ip1,j,k,jcnt)))+(6*(Ut(i,j,k,jcnt)))-(4*(Ut(im1,j,k,jcnt)))+Ut(im2,j,k,jcnt))*(0.25d0*invdx)    
term2 =     (Ut(i,j,k,1))*(-Ut(ip2,j,k,jcnt)+(8*(Ut(ip1,j,k,jcnt)))-(8*(Ut(im1,j,k,jcnt)))+Ut(im2,j,k,jcnt))*(onetwelfth*invdx)    
 adv = term1+term2
divuu(1) = adv
         else
  adv=Ut(i,j,k,1)*((Ut(im2,j,k,jcnt)-8*Ut(im1,j,k,jcnt)+8*Ut(ip1,j,k,jcnt)-Ut(ip2,j,k,jcnt))*(onetwelfth*invdx))

 divuu(1) = adv
        endif
!The skew symmetric  form ! See Morinshi et al JCP 143,90 Fully conservative finite difference schemes for incompressible flows 
!***********************************************************************
             sgn = sign(1,ceiling(Ut(i,j,k,2)))  
             Pe = sgn*tf*Ut(i,j,k,2)*dy

if((j.eq.2).or.(j.eq.(n3d(2)-1)))then
a1=max(Ut(i,j,k,2),0.0d0)
b=min(Ut(i,j,k,2),0.0d0)

 term1=a1*((Ut(i,j,k,jcnt)-Ut(i,jm1,k,jcnt))*(invdy))
 term2=b*((Ut(i,jp1,k,jcnt)-Ut(i,j,k,jcnt))*(invdy))
adv=term1+term2
divuu(2)=adv

          elseif(Pe.gt.2.0d0 .and. Pe.lt.50.0d0) then
term1 = (Abs(Ut(i,j,k,2)))*(Ut(i,jp2,k,jcnt)-(4*(Ut(i,jp1,k,jcnt)))+(6*(Ut(i,j,k,jcnt)))-(4*(Ut(i,jm1,k,jcnt)))+Ut(i,jm2,k,jcnt))*(0.25d0*invdy)    
term2 =     (Ut(i,j,k,2))*(-Ut(i,jp2,k,jcnt)+(8*(Ut(i,jp1,k,jcnt)))-(8*(Ut(i,jm1,k,jcnt)))+Ut(i,jm2,k,jcnt))*(onetwelfth*invdy)    
  adv = term1+term2
 divuu(2) =adv
         else
  adv=Ut(i,j,k,2)*((Ut(i,jm2,k,jcnt)-8*Ut(i,jm1,k,jcnt)+8*Ut(i,jp1,k,jcnt)-Ut(i,jp2,k,jcnt))*(onetwelfth*invdy))
!!endif
 divuu(2) =adv
endif
!***********************************************************************
               sgn = sign(1,ceiling(Ut(i,j,k,3)))  
!!              Pe = (sgn*tf*Ut(i,j,k,3)*dXi)*Jac(k)
              Pe = (sgn*tf*Ut(i,j,k,3)*(x3(k+1)-x3(k)))
          if( (k.le.2).or.(k.ge.(n3d(3)-1)) ) then
!        div = ((Ut(i,j,k+1,jcnt)*(Ut(i,j,k+1,3)))-(Ut(i,j,k-1,jcnt)*(Ut(i,j,k-1,3))))/(2.0*dXi)
!!        adv = (Ut(i,j,k,3)*(Ut(i,j,k+1,jcnt)-Ut(i,j,k-1,jcnt))*(0.5d0*invdXi))
        divuu(3) = (Ut(i,j,k,3)*(Ut(i,j,k+1,jcnt)-Ut(i,j,k-1,jcnt))*(0.5d0*invdXi))
!!        divuu(3) =adv
!!else


       elseif(Pe.gt.2.0d0 .and. Pe.lt.5.0d0) then
term1 = (Abs(Ut(i,j,k,3)))*(Ut(i,j,k+2,jcnt)-(4*(Ut(i,j,k+1,jcnt)))+(6*(Ut(i,j,k,jcnt)))-(4*(Ut(i,j,k-1,jcnt)))+Ut(i,j,k-2,jcnt))*(0.25d0*invdXi)   
term2 =     (Ut(i,j,k,3))*(-Ut(i,j,k+2,jcnt)+(8*(Ut(i,j,k+1,jcnt)))-(8*(Ut(i,j,k-1,jcnt)))+Ut(i,j,k-2,jcnt))*(onetwelfth*invdXi)    
     adv = (term1+term2)
 divuu(3) =adv
!        endif
! divuu(3) =adv
else
  adv=Ut(i,j,k,3)*((Ut(i,j,k-2,jcnt)-8*Ut(i,j,k-1,jcnt)+8*Ut(i,j,k+1,jcnt)-Ut(i,j,k+2,jcnt))*(onetwelfth*invdXi))
divuu(3)=adv
endif
           conv(i,j,k,jcnt) = divuu(1)+divuu(2)+(divuu(3)*Jac(k))
                                
            end do
           end do  
         end do
       end do
       
      do k=bn,en
               do j=2,n3d(2)-1
                        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
           !   if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
           !   if (i.eq.3) im2 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2

!***********************************************************************

      gradp(1)=0.5d0*(Press(ip1,j,k)-Press(im1,j,k))*(invdx)
!      gradp(1)=(-Press(ip2,j,k)+8*Press(ip1,j,k)-8*Press(im1,j,k)+Press(im2,j,k))*(onetwelfth*invdx)
      gradp(2)=0.5d0*(Press(i,jp1,k)-Press(i,jm1,k))*(invdy)
!      gradp(2)=(-Press(i,jp2,k)+8*Press(i,jp1,k)-8*Press(i,jm1,k)+Press(i,jm2,k))*(onetwelfth*invdy)
!       if(k.lt.4.or.k.gt.(n3d(3)-3)) then
      gradp(3)=0.5d0*(Press(i,j,k+1)-Press(i,j,k-1))*(invdXi*Jac(k))
!       else
!      gradp(3)=(-Press(i,j,k+2)+8*Press(i,j,k+1)-8*Press(i,j,k-1)+Press(i,j,k-2))*(onetwelfth*invdXi*Jac(k))
!       endif
               kcnt=kcnt+1
      Q(i,j,k,5)=Ut(i,j,k,1)+(dt*(forc(1)-conv(i,j,k,1)))
      Q(i,j,k,1)= Q(i,j,k,5)-(dt*gradp(1))

      Q(i,j,k,6)=Ut(i,j,k,2)+(dt*(forc(2)-conv(i,j,k,2)))
      Q(i,j,k,2)= Q(i,j,k,6)-(dt*gradp(2))

      Q(i,j,k,7)=Ut(i,j,k,3)+(dt*(forc(3)-conv(i,j,k,3)+(Ri*Ut(i,j,k,4))))
      Q(i,j,k,3)= Q(i,j,k,7)-(dt*gradp(3))

      Q(i,j,k,4)=Ut(i,j,k,4)-(dt*(conv(i,j,k,4)))

       end do
       end do  
       end do
!

      do jcnt =1,3 
             t=-dt/(Re)
             if(jcnt.eq.4) t =-dt/(Re*Pran)

!Gauss siedel iteration

      do it=1,mmax

       tmp=0.0d0

if (mod(it,2).eq.1) then  ! check for residual every 2 times

      do k=bn,en
               do j=2,n3d(2)-1
                        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
!              if (i.eq.3) im2 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
!***********************************************************************
!   if((k.eq.2).or.(k.eq.(n3d(3)-1)).or.(j.eq.2).or.(j.eq.(n3d(2)-1))) then
res=Q(i,j,k,jcnt)-                       &
     &(1.0d0+ (t*Avp(i,j,k)))*Ut(i,j,k,jcnt)-    &
     &(t*Ave(i,j,k))*Ut(ip1,j,k,jcnt)-           &
     &(t*Avw(i,j,k))*Ut(im1,j,k,jcnt)-           &
     &(t*Avn(i,j,k))*Ut(i,jp1,k,jcnt)-           &
     &(t*Avs(i,j,k))*Ut(i,jm1,k,jcnt)-           &
     &(t*Avt(i,j,k))*Ut(i,j,k+1,jcnt)-           &
     &(t*Avb(i,j,k))*Ut(i,j,k-1,jcnt)


       tmp=tmp+abs(res)

                     end do
             end do
      end do
 
       call mpi_allreduce(tmp,diver,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
       if(it.eq.1)sumnor=1.0d0/diver
                  sumav=diver*sumnor 
       if(sumav.lt.errorstop)goto 191
 endif
 
      do k=bn,en
               do j=2,n3d(2)-1
                        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
!              if (i.eq.3) im2 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
!***********************************************************************
! This is a hardwired patch at the boundaries.. Please Generalize .. not recommended 
!   if((k.eq.2).or.(k.eq.(n3d(3)-1)).or.(j.eq.2).or.(j.eq.(n3d(2)-1))) then
            Ut(i,j,k,jcnt)=(Q(i,j,k,jcnt)-                &
     &           (t*Ave(i,j,k))*Ut(ip1,j,k,jcnt)-   &
     &           (t*Avw(i,j,k))*Ut(im1,j,k,jcnt)-   &
     &           (t*Avn(i,j,k))*Ut(i,jp1,k,jcnt)-   &
     &           (t*Avs(i,j,k))*Ut(i,jm1,k,jcnt)-   &
     &           (t*Avt(i,j,k))*Ut(i,j,k+1,jcnt)-   &
     &           (t*Avb(i,j,k))*Ut(i,j,k-1,jcnt))/  &
     &      ( 1.0d0 +(t*Avp(i,j,k) ))


                     end do
             end do
      end do

!if(my_rank.eq.7)write(*,*)Ut(10,10,49,1),'Ut'


  if (my_rank.eq.0) then
               call mpi_recv(Ut(1,1,en+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,en-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
           end if
           if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1)) then
               call mpi_send(Ut(1,1,en-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,en+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,bn  ,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,bn-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
           end if
           if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0)) then
               call mpi_recv(Ut(1,1,bn-2,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,bn  ,jcnt),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,en+1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,en-1,jcnt),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
           end if 
           if (my_rank.eq.proc-1) then
               call mpi_send(Ut(1,1,bn  ,jcnt),1,twoplanes,my_rank-1,50, mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,bn-2,jcnt),1,twoplanes,my_rank-1,50, mpi_comm_world,status,ierr)
           end if 
      end do 

191    CONTINUE
      END DO


!if(my_rank.eq.6)write(*,*)Ut(10,10,en-1,1),'1'
!if(my_rank.eq.7)write(*,*)Ut(10,10,bn-2,1),'2'



! First, exchange  Updated velocity vector field between processors for 2 planes.         
  do k=bn,en
               do j=2,n3d(2)-1
                        do i=2,n3d(1)-1

!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************

               divuu(1)= (((Ut(ip1,j,k,1) + Ut(i,j,k,1))*0.5d0) -    (((Ut(i,j,k,1) + Ut(im1,j,k,1))*0.5d0)))*invdx
 
 
               divuu(2)=(((Ut(i,jp1,k,2) + Ut(i,j,k,2))*0.5d0) - (((Ut(i,j,k,2) + Ut(i,jm1,k,2))*0.5d0)))*invdy
 
 
               divuu(3)=2.0*(((Ut(i,j,k+1,3) + Ut(i,j,k,3))*0.5d0) - (((Ut(i,j,k,3) + Ut(i,j,k-1,3))*0.5d0)))/(x3(k+1)-x3(k-1))
         
         Q(i,j,k,8)=(divuu(1)+divuu(2)+divuu(3))/dt
           
          end do
        end do
      end do

!***************************************************************************************************************************
!*******initializing Pcorr before solution of pressure-poisson equation by GMRES*************************************
!on inflow and outflow faces in streamwise directions
                               Pcorr (      1, 1:n3d(2) ,istart:iend) =0.0d0
                               Pcorr ( n3d(1), 1:n3d(2) ,istart:iend) =0.0d0
!on both the faces along spanwise directions   
                               Pcorr ( 1:n3d(1),       1,  istart:iend) =0.0d0
                               Pcorr ( 1:n3d(1),  n3d(2),  istart:iend) =0.0d0

!   Top Wall and bottom wall zero 
!if (my_rank.lt.procy)      PCorr(1:n3d(1),jstart:jend,bn-1 ) =0.0d0
! if ((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1))) PCorr(1:n3d(1),jstart:jend,en+1 ) =0.0d0
     if (my_rank.eq.0)         PCorr(1:n3d(1),1:n3d(2),bn-1 ) =0.0d0
     if (my_rank.eq.proc-1)    PCorr(1:n3d(1),1:n3d(2),en+1 ) =0.0d0

        sumav=1.0d0
        m = 0

!  Begin restart loop.
        do while ((sumav>errtol).and.(m<mmax))
        m = m+1              
        h = 0.0d0
        v= 0.0d0
        c= 0.0d0
        s= 0.0d0
        g = 0.0d0
        y = 0.0d0
! Matrix vector product for the initial residual.      
   if (my_rank.eq.0) then
     call mpi_recv(PCorr(1,1,en+1),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
     call mpi_send(PCorr(1,1,en  ),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
        end if
        if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1)) then

     call mpi_send(PCorr(1,1,en  ),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)

     call mpi_recv(PCorr(1,1,en+1),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
     call mpi_send(PCorr(1,1,bn  ),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
     call mpi_recv(PCorr(1,1,bn-1),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)

         end if

         if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0)) then

      call mpi_recv(PCorr(1,1,bn-1),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_send(PCorr(1,1,bn  ),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(PCorr(1,1,en+1),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(PCorr(1,1,en  ),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)

         end if 

         if (my_rank.eq.proc-1) then
      call mpi_send(PCorr(1,1,bn  ),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(PCorr(1,1,bn-1),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
         end if  

!if(my_rank.eq.6)write(*,*)Ut(10,10,en+1,1),'1'
!if(my_rank.eq.7)write(*,*)Ut(10,10,bn,1),'2'
    do k=bn,en
             do j=2,n3d(2)-1
!               do j=bnt,ent
                        do i=2,n3d(1)-1
!**********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************

 r(i,j,k)=Q(i,j,k,8)-( Anw(i,j,k)*PCorr(im1,j  ,k)  + Ane(i,j,k)*PCorr(ip1,j  ,k)+     &
                        Ans(i,j,k)*PCorr(i  ,jm1,k)  + Ann(i,j,k)*PCorr(i  ,jp1,k)+     &
                        Anb(i,j,k)*PCorr(i  ,j  ,k-1)+ Ant(i,j,k)*PCorr(i  ,j  ,k+1)+   &
                        Anp(i,j,k)*PCorr(i  ,j  ,k))  
      
           enddo
           enddo
           enddo
!       This preconditioner changes with the number of processors!
    do k=bn,en
             do j=2,n3d(2)-1
!               do j=bnt,ent
                        do i=2,n3d(1)-1
!            rhat(i,j,k) = -w*(r(i,j,k)+Aw(i,j,k)*rhat(i-1,j,k)+As(i,j,k)*rhat(i,j-1,k)+Ab(i,j,k)*rhat(i,j,k-1))*rAp(i,j,k)
             rhat(i,j,k) = -w*(r(i,j,k)+Anw(i,j,k)*rhat(i-1,j,k)+Ans(i,j,k)*rhat(i,j-1,k)+Anb(i,j,k)*rhat(i,j,k-1))*rAnp(i,j,k)
           end do
          end do
        end do

!        rhat(2:n3d(1)-1,bnt:ent,bn:en) =((2-w)/w)*Ap(2:n3d(1)-1,bnt:ent,bn:en)*rhat(2:n3d(1)-1,bnt:ent,bn:en)
         rhat(2:n3d(1)-1,2:n3d(2)-1,bn:en) =((2-w)/w)*Anp(2:n3d(1)-1,2:n3d(2)-1,bn:en)*rhat(2:n3d(1)-1,2:n3d(2)-1,bn:en)

        do k= en,bn,-1
          do j = n3d(2)-1,2,-1
             do i = n3d(1)-1,2,-1
!                rhat(i,j,k) =-w*(rhat(i,j,k)+Ae(i,j,k)*rhat(i+1,j,k)+An(i,j,k)*rhat(i,j+1,k)+At(i,j,k)*rhat(i,j,k+1))*rAp(i,j,k)
                 rhat(i,j,k) =-w*(rhat(i,j,k)+Ane(i,j,k)*rhat(i+1,j,k)+Ann(i,j,k)*rhat(i,j+1,k)+Ant(i,j,k)*rhat(i,j,k+1))*rAnp(i,j,k)
             end do
          end do
        end do

        r(2:n3d(1)-1,2:n3d(2)-1,bn:en) = rhat(2:n3d(1)-1,2:n3d(2)-1,bn:en)                                  

        loc_rho=(sum(r(2:n3d(1)-1,2:n3d(2)-1,bn:en)*r(2:n3d(1)-1,2:n3d(2)-1,bn:en)))

        call mpi_allreduce(loc_rho,rho,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

        rho = sqrt(rho)
        if(m.eq.1)sumnor=rho
                  sumav =rho/sumnor 

        g(1) =rho
!        v(2:n3d(1)-1,bnt:ent,bn:en,1)=r(2:n3d(1)-1,bnt:ent,bn:en)/rho
        v(2:n3d(1)-1,2:n3d(2)-1,bn:en,1)=r(2:n3d(1)-1,2:n3d(2)-1,bn:en)/rho
        kit=0
! Begin gmres loop.
        do while((sumav > errtol).and.(kit < kmax))          
                   kit=kit+1                  
!**********Exchange information in phi b/w processors *****************
   if (my_rank.eq.0) then
               call mpi_recv(v(1,1,en+1,kit),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(v(1,1,en  ,kit),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
           end if
           if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1)) then
               call mpi_send(v(1,1,en  ,kit),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
               call mpi_recv(v(1,1,en+1,kit),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(v(1,1,bn  ,kit),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
               call mpi_recv(v(1,1,bn-1,kit),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
           end if

           if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0)) then
               call mpi_recv(v(1,1,bn-1,kit),1,oneplane,my_rank-1,50,mpi_comm_world,status,ierr)
               call mpi_send(v(1,1,bn  ,kit),1,oneplane,my_rank-1,50,mpi_comm_world,ierr)
               call mpi_recv(v(1,1,en+1,kit),1,oneplane,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(v(1,1,en  ,kit),1,oneplane,my_rank+1,50,mpi_comm_world,ierr)
           end if 

           if (my_rank.eq.proc-1) then
               call mpi_send(v(1,1,bn  ,kit),1,oneplane,my_rank-1,50, mpi_comm_world,ierr)
               call mpi_recv(v(1,1,bn-1,kit),1,oneplane,my_rank-1,50, mpi_comm_world,status,ierr)
           end if 


             do k=bn,en
              do j=2,n3d(2)-1
!               do j=bnt,ent
                        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************


 v(i,j,k,kit+1)=  Anw(i,j,k)*v(im1,j  ,k,kit) +Ane(i,j,k)*v(ip1,j  ,k,kit)&
                    +Ans(i,j,k)*v(i  ,jm1,k,kit) +Ann(i,j,k)*v(i  ,jp1,k,kit)&
                    +Anb(i,j,k)*v(i  ,j  ,k-1,kit) +Ant(i,j,k)*v(i  ,j  ,k+1,kit)&
                    +Anp(i,j,k)*v(i  ,j  ,k,kit) 
        
           enddo
           enddo
           enddo

!       This preconditioner changes with the number of processors!
  do k=bn,en
             do j=2,n3d(2)-1
!               do j=bnt,ent
                        do i=2,n3d(1)-1
             rhat(i,j,k) =-w*(v(i,j,k,kit+1) + Anw(i,j,k)*rhat(i-1,j,k) + Ans(i,j,k)*rhat(i,j-1,k) + Anb(i,j,k)*rhat(i,j,k-1))*rAnp(i,j,k)

           end do
           end do
           end do

!         rhat(2:n3d(1)-1,bnt:ent,bn:en) =  ((2-w)/w)*Anp(2:n3d(1)-1,bnt:ent,bn:en)*rhat(2:n3d(1)-1,bnt:ent,bn:en)
       rhat(2:n3d(1)-1,2:n3d(2)-1,bn:en) =  ((2-w)/w)*Anp(2:n3d(1)-1,2:n3d(2)-1,bn:en)*rhat(2:n3d(1)-1,2:n3d(2)-1,bn:en)
          do k= en,bn,-1
!          do j = ent,bnt,-1
           do j = n3d(2)-1,2,-1
          do i = n3d(1)-1,2,-1
!                rhat(i,j,k) =-w*(rhat(i,j,k)+Ae(i,j,k)*rhat(i+1,j,k) + An(i,j,k)*rhat(i,j+1,k)+ At(i,j,k)*rhat(i,j,k+1))*rAp(i,j,k)
                 rhat(i,j,k) =-w*(rhat(i,j,k)+Ane(i,j,k)*rhat(i+1,j,k) + Ann(i,j,k)*rhat(i,j+1,k)+ Ant(i,j,k)*rhat(i,j,k+1))*rAnp(i,j,k)
          end do
          end do
          end do
!        v(2:n3d(1)-1,bnt:ent,bn:en,kit+1) = rhat(2:n3d(1)-1,bnt:ent,bn:en)
         v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1) = rhat(2:n3d(1)-1,2:n3d(2)-1,bn:en)
                                                 
! Begin modified GS. May need to reorthogonalize.   
                do k=1,kit                               

!       tempgmres=sum(v(2:n3d(1)-1,bnt:ent,bn:en,k)*v(2:n3d(1)-1,bnt:ent,bn:en,kit+1))
        tempgmres=sum(v(2:n3d(1)-1,2:n3d(2)-1,bn:en,k)*v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1))

        call mpi_allreduce(tempgmres,h(k,kit),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

!        v(2:n3d(1)-1,bnt:ent,bn:en,kit+1)=v(2:n3d(1)-1,bnt:ent,bn:en,kit+1)-h(k,kit)*v(2:n3d(1)-1,bnt:ent,bn:en,k)
        v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)=v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)-h(k,kit)*v(2:n3d(1)-1,2:n3d(2)-1,bn:en,k)
        
                end do

!                       tempgmres=(sum(v(2:n3d(1)-1,bnt:ent,bn:en,kit+1)*v(2:n3d(1)-1,bnt:ent,bn:en,kit+1)))
                        tempgmres=(sum(v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)*v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)))

                    call mpi_allreduce(tempgmres,h(kit+1,kit),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

                       h(kit+1,kit) = sqrt(h(kit+1,kit))

                       if (h(kit+1,kit).gt.0.0.or.h(kit+1,kit).lt.0.0) then

 !                          v(2:n3d(1)-1,bnt:ent,bn:en,kit+1)=v(2:n3d(1)-1,bnt:ent,bn:en,kit+1)/h(kit+1,kit)
                            v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)=v(2:n3d(1)-1,2:n3d(2)-1,bn:en,kit+1)/h(kit+1,kit)
                       end if

                       if (kit>1)  then                           
! Apply old Givens rotations to h(1:kit,kit).
                               do i=1,kit-1
                                 hik    =c(i)*h(i,kit)-s(i)*h(i+1,kit)
                                 hipk   =s(i)*h(i,kit)+c(i)*h(i+1,kit)
                                 h(i,kit) = hik
                                 h(i+1,kit) = hipk
                              end do
                    end if
                    nu=sqrt(h(kit,kit)**2 + h(kit+1,kit)**2)                
! May need better Givens implementation.
! Define and Apply new Givens rotations to h(kit:kit+1,kit).  
                    if (nu.gt.0.0) then
                        c(kit)=h(kit,kit)/nu
                        s(kit)=-h(kit+1,kit)/nu
                        h(kit,kit)=c(kit)*h(kit,kit)-s(kit)*h(kit+1,kit)
                        h(kit+1,kit)=0
                        gk    =c(kit)*g(kit) -s(kit)*g(kit+1)
                        gkp  =s(kit)*g(kit) +c(kit)*g(kit+1)
                        g(kit) = gk
                        g(kit+1) = gkp
                end if
                rho=abs(g(kit+1))
                  sumav =rho/sumnor 
                    mag(kit) = rho
! End of gmres loop.
        end do
! h(1:kit,1:kit) is upper triangular matrix in QR.                                        
        y(kit) = g(kit)/h(kit,kit)
        do i = kit-1,1,-1
                y(i) = g(i)
                do k = i+1,kit
                        y(i) = y(i) -h(i,k)*y(k)
                end do
                y(i) = y(i)/h(i,i)
        end do
! Form linear combination.
         do i = 1,kit                                   
!           PCorr(2:n3d(1)-1,bnt:ent,bn:en) = PCorr(2:n3d(1)-1,bnt:ent,bn:en) + v(2:n3d(1)-1,bnt:ent,bn:en,i)*y(i)
            PCorr(2:n3d(1)-1,2:n3d(2)-1,bn:en) = PCorr(2:n3d(1)-1,2:n3d(2)-1,bn:en) + v(2:n3d(1)-1,2:n3d(2)-1,bn:en,i)*y(i)
         end do
! End restart loop.         
        
        end do


        if (my_rank.eq.0) then
     call mpi_recv(PCorr(1,1,en+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
     call mpi_send(PCorr(1,1,en-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
        end if
        if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1)) then

     call mpi_send(PCorr(1,1,en-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
     call mpi_recv(PCorr(1,1,en+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
     call mpi_send(PCorr(1,1,bn  ),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
     call mpi_recv(PCorr(1,1,bn-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)

         end if

         if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0)) then

      call mpi_recv(PCorr(1,1,bn-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_send(PCorr(1,1,bn  ),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(PCorr(1,1,en+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(PCorr(1,1,en-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)

         end if 

         if (my_rank.eq.proc-1) then
      call mpi_send(PCorr(1,1,bn  ),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(PCorr(1,1,bn-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
         end if  





!         Press(2:n3d(1)-1,bnt:ent,bn:en)=Press(2:n3d(1)-1,bnt:ent,bn:en)+PCorr(2:n3d(1)-1,bnt:ent,bn:en)+Pdist(2:n3d(1)-1,bnt:ent,bn:en)


   if (my_rank.eq.0) PCorr(1:n3d(1),2:n3d(2)-1,bn-1) = ( (4.0d0*PCorr(1:n3d(1),2:n3d(2)-1,bn)) -PCorr(1:n3d(1),2:n3d(2)-1,bn+1) )/3.0d0
!                        Pcorr (1:n3d(1),2:n3d(2)-1,1) = Pcorr (1:n3d(1),2:n3d(2)-1,2)

!if((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1)))  PCorr(1:n3d(1),bnt:ent,en+1)= (5.0*PCorr(1:n3d(1),bnt:ent,en)+              &
!          &                  (-4.0*PCorr(1:n3d(1),bnt:ent,en-1)+ PCorr(1:n3d(1),bnt:ent,en-2) ))/2.0d0 
! if (my_rank.eq.proc-1) Pcorr (1:n3d(1),2:n3d(2)-1,en+1) = 0.0d0
if(my_rank.eq.proc-1)  PCorr(1:n3d(1),2:n3d(2)-1,en+1)= (5.0*PCorr(1:n3d(1),2:n3d(2)-1,en)+              &
          &                  (-4.0*PCorr(1:n3d(1),2:n3d(2)-1,en-1)+ PCorr(1:n3d(1),2:n3d(2)-1,en-2) ))/2.0d0 
!!if((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1)))  PCorr(1:n3d(1),bnt:ent,en+1)= 0.0d0
! ensure periodicity by mapping in spanwise direction for the correction pressure (Pcorr).
                   Pcorr (1:n3d(1),     1,istart:iend) = Pcorr(1:n3d(1),n3d(2)-1,istart:iend)
                   Pcorr (1:n3d(1),n3d(2),istart:iend) = Pcorr(1:n3d(1),       2,istart:iend) 

!*********************************************************************************************
                PCorr(1,2:n3d(2)-1,istart:iend) =  (4.0d0*PCorr(2,2:n3d(2)-1,istart:iend) -PCorr(3,2:n3d(2)-1,istart:iend) )/3.0d0

 PCorr(n3d(1),2:n3d(2)-1,istart:iend)= (5.0*PCorr(n3d(1)-1,2:n3d(2)-1,istart:iend)+              &
          &                  (-4.0*PCorr(n3d(1)-2,2:n3d(2)-1,istart:iend)+ PCorr(n3d(1)-3,2:n3d(2)-1,istart:iend) ))/2.0d0 

!          PCorr(n3d(1),2:n3d(2)-1,istart:iend)= 0.0d0
!******************************correcting Pressure values******************************************************************


     Press(2:n3d(1)-1,2:n3d(2)-1,bn:en)=Press(2:n3d(1)-1,2:n3d(2)-1,bn:en)+PCorr(2:n3d(1)-1,2:n3d(2)-1,bn:en)


  do k=bn,en
             do j=2,n3d(2)-1
!               do j=bnt,ent
                        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
!              if (i.eq.3) im2 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
!************************************************************************
!************************************************************************
      gradp(1)=0.5d0*(PCorr(ip1,j,k)-PCorr(im1,j,k))*(invdx)
!      gradp(1)=(-Press(ip2,j,k)+8*Press(ip1,j,k)-8*Press(im1,j,k)+Press(im2,j,k))*(onetwelfth*invdx)
      gradp(2)=0.5d0*(PCorr(i,jp1,k)-PCorr(i,jm1,k))*(invdy)
!      gradp(2)=(-Press(i,jp2,k)+8*Press(i,jp1,k)-8*Press(i,jm1,k)+Press(i,jm2,k))*(onetwelfth*invdy)
!       if(k.lt.4.or.k.gt.(n3d(3)-3)) then
      gradp(3)=0.5d0*(PCorr(i,j,k+1)-PCorr(i,j,k-1))*(invdXi*Jac(k))

!************************************************************************

 
!************************************************************************
!*                  V^n+1   =V^n      - dt*grad( P')
                 Ut(i,j,k,1)= Ut(i,j,k,1)-(dt*gradp(1))
                 Ut(i,j,k,2)= Ut(i,j,k,2)-(dt*gradp(2))
                 Ut(i,j,k,3)= Ut(i,j,k,3)-(dt*gradp(3))
          end do
        end do
      end do

!***********************************************************************************************************
                    if(my_rank.eq.0)  Ut(1:n3d(1),2:n3d(2)-1,bn-1,1) = 0.0d0    ! u',v',w' are zero on the plate. 
                    if(my_rank.eq.0)  Ut(1:n3d(1),2:n3d(2)-1,bn-1,2) = 0.0d0    !    
                    if(my_rank.eq.0)  Ut(1:n3d(1),2:n3d(2)-1,bn-1,3) = 0.0d0    !



  if (my_rank.eq.proc-1)    Ut(1:n3d(1),2:n3d(2)-1,en+1,1)=(4.0*Ut(1:n3d(1),2:n3d(2)-1,en,1)-Ut(1:n3d(1),2:n3d(2)-1,en-1,1))/3.0d0
!2.0d0! first derivative of u' zero at top!
!   if ((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1)))  Ut(1:n3d(1),bnt:ent,en+1,1) = Ue


!  if ((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1)))  Ut(1:n3d(1),bnt:ent,en+1,2)=(5.0*Ut(1:n3d(1),bnt:ent,en,2)-4.0*Ut(1:n3d(1),bnt:ent,en-1,2)+    &
!                &         Ut(1:n3d(1),bnt:ent,en-2,2))/2.0d0! second derivative of v' zero at top!

    if (my_rank.eq.proc-1)  Ut(1:n3d(1),2:n3d(2)-1,en+1,2)=(4.0*Ut(1:n3d(1),2:n3d(2)-1,en,2)-Ut(1:n3d(1),2:n3d(2)-1,en-1,2))/3.0d0

    if (my_rank.eq.proc-1)  Ut(1:n3d(1),2:n3d(2)-1,en+1,3)=(4.0*Ut(1:n3d(1),2:n3d(2)-1,en,3)-Ut(1:n3d(1),2:n3d(2)-1,en-1,3))/3.0d0 
! if (my_rank.eq.proc-1)  Ut(1:n3d(1),2:n3d(2)-1,en+1,3)=(5.0*Ut(1:n3d(1),2:n3d(2)-1,en,3)-4.0*Ut(1:n3d(1),2:n3d(2)-1,en-1,3)+    &
!                &         Ut(1:n3d(1),2:n3d(2)-1,en-2,3))/2.0d0! second derivative of v' zero at top!
!***********************************************************************************************************

! Outflow BCs at the streamwise 
!!Ut(n3d(1),bnt:ent,istart:iend,1:3) =(5.0*Ut(n3d(1)-1,bnt:ent,istart:iend,1:3)-4.0*Ut(n3d(1)-2,bnt:ent,istart:iend,1:3)      &
!!               &           +Ut(n3d(1)-3,bnt:ent,istart:iend,1:3))/2.0d0  ! second derivative of u,v,w zero at front

Ut(n3d(1),2:n3d(2)-1,istart:iend,1:3) =  (4.0d0*Ut(n3d(1)-1,2:n3d(2)-1,istart:iend,1:3) -Ut(n3d(1)-2,2:n3d(2)-1,istart:iend,1:3) )/3.0d0
!Ut(n3d(1),2:n3d(2)-1,istart:iend,1) =(5.0*Ut(n3d(1)-1,2:n3d(2)-1,istart:iend,1)-4.0*Ut(n3d(1)-2,2:n3d(2)-1,istart:iend,1)      &
!               &           +Ut(n3d(1)-3,2:n3d(2)-1,istart:iend,1))/2.0d0  ! second derivative of u,v,w zero at front

!Ut(n3d(1),2:n3d(2)-1,istart:iend,2) =  (4.0d0*Ut(n3d(1)-1,2:n3d(2)-1,istart:iend,2) -Ut(n3d(1)-2,2:n3d(2)-1,istart:iend,2) )/3.0d0
!Ut(n3d(1),2:n3d(2)-1,istart:iend,3) =  (4.0d0*Ut(n3d(1)-1,2:n3d(2)-1,istart:iend,3) -Ut(n3d(1)-2,2:n3d(2)-1,istart:iend,3) )/3.0d0
!************************************************************************************************************      

! Update Pressure at boundaries using momentum equation
 !!! Bottom Boundary 
   if(my_rank.eq.0)then

          do j=2,n3d(2)-1
    do i=2,n3d(1)-1

!    fder =(-11.0*Ut(i,j,1,3) + 18.0*Ut(i,j,2,3) - 9.0*Ut(i,j,3,3) + 2.0*Ut(i,j,4,3))/(6.0*dXi)
    fder=(Ut(i,j,2,3)-Ut(i,j,1,3))/df(1)
!    divuu(3)=Ut(i,j,1,3)*(fder*Jac(1))
    divuu(3)=Ut(i,j,1,3)*fder


divgrad(1)=(Ut(i-1,j,1,3)+Ut(i+1,j,1,3)-2*Ut(i,j,1,3))/(dx*dx)
divgrad(2)=(Ut(i,j-1,1,3)+Ut(i,j+1,1,3)-2*Ut(i,j,1,3))/(dy*dy)
divgrad(3)=(2.d0*df(1)*Ut(i,j,3,3)-2.d0*df(2)*Ut(i,j,2,3)+2.d0*(df(2)-df(1))*Ut(i,j,1,3))/(df(1)*df(2)**2-(df(2)*df(1)**2))
!divgrad(3)=(2*Ut(i,j,1,3) - 5*Ut(i,j,2,3) + 4*Ut(i,j,3,3) - Ut(i,j,4,3))/(dXi*dXi)

!divgrad(3)  = (divgrad(3)*(Jac(1)*Jac(1))) - ((Jac2(1)*fder)*(Jac(1)*Jac(1)*Jac(1)))

! Finding The H and Pressure at the Top Wall          

Flux= -divuu(3) +((divgrad(1)+divgrad(2)+divgrad(3))/Re) 

!Press(i,j,1)=((-6.0*dXi*Flux/Jac(1))+ 18.0*Press(i,j,2)-9.0*Press(i,j,3)+2.d0*Press(i,j,4))/11.0d0
Press(i,j,1)=Press(i,j,2)-Flux*df(1)
            enddo
          enddo
endif

 !!! Top Boundary 


!if ((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1)))then
if(my_rank.eq.proc-1) then
           npln =n3d(3)

!          do j=bnt,ent
          do j=2,n3d(2)-1
          do i=2,n3d(1)-1

fder =(11.0*Ut(i,j,npln,3) - 18.0*Ut(i,j,npln-1,3) + 9.0*Ut(i,j,npln-2,3) - 2.0*Ut(i,j,npln-3,3))/(6.0*dXi)
divuu(3)=(fder*Jac(npln))

!!fder =(Ut(i+1,j,npln,3)-Ut(i-1,j,npln,3))/2.0d0*dx
!!divuu(1)=Ut(i,j,npln,1)*fder

!!fder1 =(Ut(i,j+1,npln,3)-Ut(i,j-1,npln,3))/2.0d0*dy
!!divuu(2)=Ut(i,j,npln,2)*fder1


!!divuu(3)=(Ut(i,j,npln,3)-Utold(i,j,npln,3))/dt

!! Press(i,j,npln)=((4.0*Press(i,j,npln-1)-Press(i,j,npln-2))-((divuu(1)+divuu(2)+divuu(3))*2.0d0*dXi*Jac(npln)))/3.0d0
Press(i,j,npln)=(divuu(3)/Re)!-0.25d0*(((Ut(i,j,npln,1))**2)+((Ut(i,j,npln,2))**2) + ((Ut(i,j,npln,3))**2))*(1.0d0-dtanh((Ut(i,j,npln,1)/dble(delprmtr))))

            enddo
          enddo
endif
!! if ((my_rank.lt.proc).and.(my_rank.gt.(proc-procy-1)))  Press(1:n3d(1),bnt:ent,n3d(3))=(4.0*Press(1:n3d(1),bnt:ent,n3d(3)-1)-Press(1:n3d(1),bnt:ent,n3d(3)-2))/3.0d0

!!Press(n3d(1),bnt:ent,istart:iend) = (4.0d0*Press(n3d(1)-1,bnt:ent,istart:iend)-Press(n3d(1)-2,bnt:ent,istart:iend))/3.0d0 ! first derivative  of p zero  at left
do k=bn,en
!do j=bnt,ent
do j=2,n3d(2)-1
!Press(n3d(1),j,k)=((3.0d0*Ut(n3d(1),j,k,1)-4.0d0*Ut(n3d(1)-1,j,k,1) + Ut(n3d(1)-2,j,k,1))/(2.0d0*dx*Re)) &
!               & -0.25d0*(((Ut(n3d(1),j,k,1))**2)+((Ut(n3d(1),j,k,2))**2) + ((Ut(n3d(1),j,k,3))**2))*(1.0d0-dtanh((Ut(n3d(1),j,k,1)/dble(delprmtr))))

Press(n3d(1),j,k)=((3.0d0*Ut(n3d(1),j,k,1)-4.0d0*Ut(n3d(1)-1,j,k,1) + Ut(n3d(1)-2,j,k,1))/(2.0d0*dx*Re)) 
!               & -0.25d0*(((Ut(n3d(1),j,k,1))**2)+((Ut(n3d(1),j,k,2))**2) + ((Ut(n3d(1),j,k,3))**2))*(1.0d0-dtanh((Ut(n3d(1),j,k,1)/dble(delprmtr))))

!!fder =(Ut(n3d(1),j+1,k,1)-Ut(n3d(1),j-1,k,1))/2.0d0*dy
!!divuu(1)=Ut(n3d(1),j,k,2)*fder

!!fder1 =(Ut(n3d(1),j,k+1,1)-Ut(n3d(1),j,k-1,1))/2.0d0*dXi
!!divuu(2)=Ut(n3d(1),j,k,3)*fder1*Jac(k)

!!divuu(3)=(Ut(n3d(1),j,k,1)-Utold(n3d(1),j,k,1))/dt

!! Press(n3d(1),j,k)=((4.0*Press(n3d(1)-1,j,k)-Press(n3d(1)-2,j,k))-((divuu(1)+divuu(2)+divuu(3))*2.0d0*dx))/3.0d0
enddo
enddo   

do k=bn,en
!do j=bnt,ent
do j=2,n3d(2)-1

         !divuu(1) = Ut(1,j,k,1)*(((4*Ut(2,j,k,1)-Ut(3,j,k,1)-3*Ut(1,j,k,1)))/(2.0d0*dx))
         divuu(1) = Ut(1,j,k,1)*((Ut(2,j,k,1)-Ut(1,j,k,1))/(dx))
         divuu(2) = Ut(1,j,k,2)*((Ut(1,j+1,k,1)-Ut(1,j-1,k,1))/(2.0d0*dy))
         divuu(3) = Ut(1,j,k,3)*((Ut(1,j,k+1,1)-Ut(1,j,k-1,1))/(2.0d0*dXi))

       divgrad(1) =(Ut(1,j,k,1)+Ut(3,j,k,1)-2*Ut(2,j,k,1))/(dx*dx)
       divgrad(2) =(Ut(1,j+1,k,1)+Ut(1,j-1,k,1)-2*Ut(1,j,k,1))/(dy*dy)
if (k.le.2) then
        fder =     (Ut(1,j,k+1,1)-Ut(1,j,k,1))/dXi
divgrad(3) = (Ut(1,j,k+2,1) -2.0d0*Ut(1,j,k+1,1) + Ut(1,j,k,1))/(dXi*dXi)
 elseif (k.ge.n3d(3)-1) then
        fder =     (Ut(1,j,k,1)-Ut(1,j,k-1,1))/dXi
divgrad(3) = (Ut(1,j,k-2,1) -2.0d0*Ut(1,j,k-1,1) + Ut(1,j,k,1))/(dXi*dXi)
 else
        fder =     (Ut(1,j,k+1,1)-Ut(1,j,k-1,1))/(2.0d0*dXi)
divgrad(3) = (-2*Ut(1,j,k,1) + Ut(1,j,k+1,1) + Ut(1,j,k-1,1))/(dXi*dXi)
endif
!divgrad(3) = (-2*Ut(1,j,k,1) + Ut(1,j,k+1,1) + Ut(1,j,k-1,1))/(dXi*dXi)
divgrad(3)  = (divgrad(3)*(Jac(k)*Jac(k))) - ((Jac2(k)*fder)*(Jac(k)*Jac(k)*Jac(k)))

fder2=(Ut(1,j,k,1)-Utold(1,j,k,1))/dt

Press(1,j,k) =((2*dx*((fder2+divuu(1)+divuu(2)+divuu(3)*Jac(k))-((divgrad(1)+divgrad(2)+divgrad(3))/Re)))-forc(1)-Press(3,j,k)+(4*Press(2,j,k))) /3.0d0
enddo
enddo        
! Exchange  Updated velocity vector field between processors for 2 planes.         
 

                Ut(:,1,istart:iend,1:3)= Ut(:,n3d(2)-1,istart:iend,1:3)         ! ensure periodicity of u,v,w by mapping in spanwise direction
                Ut(:,n3d(2),istart:iend,1:3)= Ut(:,2,istart:iend,1:3)           ! ensure periodicity of u,v,w by mapping in spanwise direction

                Press(:,     1,istart:iend)= Press(:,n3d(2)-1,istart:iend)      ! ensure periodicity by mapping in spanwise direction
                Press(:,n3d(2),istart:iend)= Press(:,       2,istart:iend)      ! ensure periodicity by mapping in spanwise direction


  !-----------------------------------------------------------------------------

! Exchange  Updated velocity vector field between processors for 2 planes.         

           if (my_rank.eq.0) then
               call mpi_recv(Ut(1,1,en+1,1),1,threevar,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,en-1,1),1,threevar,my_rank+1,50,mpi_comm_world,ierr)
           end if
           if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1)) then
               call mpi_send(Ut(1,1,en-1,1),1,threevar,my_rank+1,50,mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,en+1,1),1,threevar,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,bn  ,1),1,threevar,my_rank-1,50,mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,bn-2,1),1,threevar,my_rank-1,50,mpi_comm_world,status,ierr)
           end if

           if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0)) then
               call mpi_recv(Ut(1,1,bn-2,1),1,threevar,my_rank-1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,bn  ,1),1,threevar,my_rank-1,50,mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,en+1,1),1,threevar,my_rank+1,50,mpi_comm_world,status,ierr)
               call mpi_send(Ut(1,1,en-1,1),1,threevar,my_rank+1,50,mpi_comm_world,ierr)
           end if 

           if (my_rank.eq.proc-1) then
               call mpi_send(Ut(1,1,bn  ,1),1,threevar,my_rank-1,50, mpi_comm_world,ierr)
               call mpi_recv(Ut(1,1,bn-2,1),1,threevar,my_rank-1,50, mpi_comm_world,status,ierr)
           end if 

!Exchange updated pressure field between processors for plane after interface.         
        if (my_rank.eq.0) then
     call mpi_recv(Press(1,1,en+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
     call mpi_send(Press(1,1,en-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
        end if
        if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.1)) then

     call mpi_send(Press(1,1,en-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
     call mpi_recv(Press(1,1,en+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
     call mpi_send(Press(1,1,bn  ),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
     call mpi_recv(Press(1,1,bn-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)

         end if

         if ((my_rank.gt.0).and.(my_rank.lt.proc-1).and.(mod(my_rank,2).eq.0)) then
      call mpi_recv(Press(1,1,bn-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
      call mpi_send(Press(1,1,bn  ),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Press(1,1,en+1),1,twoplanes,my_rank+1,50,mpi_comm_world,status,ierr)
      call mpi_send(Press(1,1,en-1),1,twoplanes,my_rank+1,50,mpi_comm_world,ierr)
         end if 

         if (my_rank.eq.proc-1) then
      call mpi_send(Press(1,1,bn  ),1,twoplanes,my_rank-1,50,mpi_comm_world,ierr)
      call mpi_recv(Press(1,1,bn-2),1,twoplanes,my_rank-1,50,mpi_comm_world,status,ierr)
         end if  


 if(mod(nstep,disp1)==0) then
! Calculate Divergence 
          tmp =0.0d0
!************************************************************************
        do k=bn,en
             do j=2,n3d(2)-1
!               do j=bnt,ent
                        do i=2,n3d(1)-1
!***********************************************************************
                im1 =i-1
!              if (i.eq.2) im1 =n3d(1)-1
                im2 =im1-1
!              if (i.eq.3) im2 =n3d(1)-1
!***********************************************************************
                jm1 =j-1
              if (j.eq.2) jm1 =n3d(2)-1
                jm2 =jm1-1
              if (j.eq.3) jm2 =n3d(2)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
                ip2 =ip1+1
!              if (i.eq.n3d(1)-2) ip2 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
                jp2 =jp1+1
              if (j.eq.n3d(2)-2) jp2 =2
!***********************************************************************

 if((i.eq.2).or.(i.eq.(n3d(1)-1)))then
  divuu(1)= (Ut(ip1,j,k,1)-Ut(im1,j,k,1))*(0.5d0*invdx)
else
  divuu(1)= ((Ut(im2,j,k,1)-8*Ut(im1,j,k,1)+8*Ut(ip1,j,k,1)-Ut(ip2,j,k,1))*(onetwelfth*invdx))
endif
 if((j.eq.2).or.(j.eq.(n3d(2)-1)))then
  divuu(2)= (Ut(i,jp1,k,2)-Ut(i,jm1,k,2))*(0.5d0*invdy)
else 
 divuu(2)= ((Ut(i,jm2,k,2)-8*Ut(i,jm1,k,2)+8*Ut(i,jp1,k,2)-Ut(i,jp2,k,2))*(onetwelfth*invdy))
endif 
if(k.le.2.or.k.ge.(n3d(3)-1)) then
  divuu(3)=((Ut(i,j,k+1,3)-Ut(i,j,k-1,3))*(0.5d0*invdXi))
else
  divuu(3)= ((Ut(i,j,k-2,3)-8*Ut(i,j,k-1,3)+8*Ut(i,j,k+1,3)-Ut(i,j,k+2,3))*(onetwelfth*invdXi))
endif
          tmp =tmp + (divuu(1)+divuu(2)+(divuu(3)*Jac(k)))
           
          end do
        end do
      end do
      call MPI_Reduce(tmp, diver,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
    if (my_rank.eq.0) write(*,21)time,nstep,(abs(diver)/((n3d(2)-2)*(n3d(1)-2)*(n3d(3)-2)))
        do i=1,5
         UbDNS(i) =0.0d0
     itemp= iend
         if(my_rank.eq.proc-1) itemp= iend -1
        do k=istart,itemp
          do j=2,n3d(2)-1
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************
    dAv = dy*(x3(k+1)-x3(k))
   UbDNS(i) =UbDNS(i)+(0.25d0*(Ut(plane(i),j,k,1)+Ut(plane(i),jp1,k,1)+Ut(plane(i),j,k+1,1)+Ut(plane(i),jp1,k+1,1))*(dAv))
          enddo
      enddo
      call MPI_Reduce(UbDNS(i),tmp,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
      UbDNS(i)=tmp/Av
      enddo


                   Nuss =0.0d0

 do j=2,n3d(2)-1
 do i=2,n3d(1)-1
if (my_rank.eq.proc-1) Nuss= Nuss -(((Ut(i,j,n3d(3),4)- Ut(i,j,n3d(3)-2,4))*db(1)**2+ (Ut(i,j,n3d(3)-1,4)-Ut(i,j,n3d(3),4))*db(2)**2)/(db(2)*db(1)**2-db(2)**2*db(1)))
if (my_rank.eq.0)      Nuss= Nuss  +(((Ut(i,j,3,4)-Ut(i,j,1,4))*df(1)**2+(Ut(i,j,1,4)-Ut(i,j,2,4))*df(2)**2)/(df(2)*df(1)**2-df(2)**2*df(1)))
 enddo
 enddo
  if ((my_rank.eq.0).or.(my_rank.eq.proc-1)) Nuss = Nuss/((n3d(1)-2)*(n3d(2)-2))

         if (my_rank.eq.proc-1) call mpi_send(Nuss,1,mpi_double_precision,0,50,mpi_comm_world,ierr)
         if (my_rank.eq.0)      call mpi_recv(tmp,1,mpi_double_precision,proc-1,50,mpi_comm_world,status,ierr)


 if (my_rank.eq.0) then

  open (unit=12,file='massflow.xy',status='unknown',access='append')
  write(12,99)time,UbDNS(1),UbDNS(2),UbDNS(3),UbDNS(4),UbDNS(5)!,Nuss   !,tmp,Pprod,Punc,Preq,Pprodpercent,Preqpercent,Pnetpercent,CF,Pbar
  close(12)
 endif
 endif


!***************************************************subroutine for wall-shear-stress********************starts here**************************************
       ! tauxz(:)=0.0d0

 if(my_rank.eq.0)then

        tauxz(:)=0.0d0

   df =x3(2)-x3(1)
   d2f =x3(3)-x3(1)

                  df(1)=x3(2)-x3(1)
                  df(2)=x3(3)-x3(1)
                  df(3)=x3(4)-x3(1)



        do i=2,n3d(1)-1
        !do j=bnt,ent
        do j=2,n3d(2)-1

                im1 =i-1
                im2 =i-2
                ip1= i+1
                ip2 =i+2
!              if (i.eq.2) im2 =n3d(1)-2
!              if (i.eq.n3d(1)-1) ip2 =3
                im3 =im2-1
                ip3 =ip2+1
!              if (i.eq.3) im3 =n3d(1)-2
!              if (i.eq.n3d(1)-2) ip3 =3

 if (i.le.3.or.i.ge.n3d(1)-2) then
 mfb1(i,j) = ( Ut(ip1,j,1,3)-Ut(im1,j,1,3) )/(2.0d0*dx)
else
 mfb1(i,j) = (-Ut(im3,j,1,3)+(9.0*Ut(im2,j,1,3))-(45.0*Ut(im1,j,1,3))+(45.0*Ut(ip1,j,1,3))-(9.0*Ut(ip2,j,1,3))+Ut(ip3,j,1,3))/(60.0d0*dx)
endif
 fder = ( ((Ut(i,j,3,1)-Ut(i,j,1,1))*(df(1)**2)) + ((Ut(i,j,1,1)-Ut(i,j,2,1))*(df(2)**2)) )/( (df(2)*(df(1)**2))-((df(2)**2)*df(1)) )
!  fder=(Ut(i,j,2,1)-Ut(i,j,1,1))/df(1)
   mfb2(i,j) = fder

      tauxz(i)=tauxz(i)+((mfb1(i,j)+mfb2(i,j))/Re)

        enddo
        enddo


      do i=2,n3d(1)-1
          sumxz(i)= tauxz(i)+sumxz(i)
      enddo

if(mod(nstep,500).eq.0) then

sumxz(:) = sumxz(:)/500.0d0

     open(2,file='numrcltauwall.xy',status='unknown',access='append')
        do i=2,n3d(1)-1
         write(2,94)time, x1(i), sumxz(i)/dble((n3d(2)-2)),2.0d0*(sumxz(i)/(Ue*Ue*dble(n3d(2)-2)))
        enddo
      close(2)

 sumxz(:)=0.0d0

endif
 
endif
!***********************************************unformatted file write*******************************************!!

IF(mod(nstep,1000)==0) THEN
		WRITE(filname,'(a,i3.3,a)')'3D',my_rank+1,'.tp'
		!Writing tecplot ascii FORMAT directly
		OPEN (UNIT=12,FILE=filname,STATUS='unknown')
		IF(my_rank.eq.0) WRITE (12,03)'TITLE ="',time,nstep,'"'
		IF(my_rank.eq.0) WRITE (12,*)'Variables = "x","y","z","U","V","W","P","Um","Vm","Wm","Uf","Vf","Wf","Wt","Nk"'
		IF(my_rank.eq.0) WRITE (12,04)'ZONE k=',n3d(1),',j=',n3d(2)-2,',i=',n3d(3),',DATAPACKING="POINT"'
		DO i=1,n3d(1)
			WRITE(12,*)
			DO j=2,n3d(2)-1
				WRITE(12,*)
				DO k=istart,iend
					WRITE(12,97)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Press(i,j,k),usum2(i,k),vsum2(i,k),wsum2(i,k),Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Wt(i,k),new(k)
				END DO
			END DO
		END DO
		CLOSE(12)
	END IF

!*******************************************************************************************************************
if(nstep.gt.5)then

ke1 =0.0d0
ke =0.0d0

    do k=bn,en
!     do j=bnt,ent
     do j=2,n3d(2)-1
      do i=1,n3d(1)
  Ufluc(i,j,k,1) = Ut(i,j,k,1)-usum2(i,k)
  Ufluc(i,j,k,2) = Ut(i,j,k,2)-vsum2(i,k)
  Ufluc(i,j,k,3) = Ut(i,j,k,3)-wsum2(i,k)
  
      enddo
     enddo
    enddo

    do k=bn,en
!      do j=bnt,ent
     do j=2,n3d(2)-1
        do i=1,n3d(1)
    ke1 = ke1 + Ufluc(i,j,k,1)*Ufluc(i,j,k,1) + Ufluc(i,j,k,2)*Ufluc(i,j,k,2) + Ufluc(i,j,k,3)*Ufluc(i,j,k,3)
        enddo
       enddo
      enddo

      call MPI_Reduce(ke1,ke,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

 if (my_rank.eq.0) then

ke = ke/dble(n3d(1)*(n3d(2)-2)*(n3d(3)-2))

 open (unit=19,file='kineticenergy.dat',status='unknown',access='append')
    write(19,255) time, ke
 close(19)
255 format(2(2x,E20.10))   

 endif

endif

!************************************************************************************************************
      do i=1,n3d(1)
       do j=2,n3d(2)-1
        do k=bn,en


if (Ufluc(i,j,k,1).gt.25.0d0) write(*,*) nstep,time,i,j,k,Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Press(i,j,k)
if (Ufluc(i,j,k,2).gt.25.0d0) write(*,*) nstep,time,i,j,k,Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Press(i,j,k)
if (Ufluc(i,j,k,3).gt.25.0d0) write(*,*) nstep,time,i,j,k,Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Press(i,j,k)


if (Ufluc(i,j,k,1).gt.25.0d0) then
print*,"uprime>25"
stop 
endif

if (Ufluc(i,j,k,2).gt.25.0d0) then
print*,"vprime>25"
stop 
endif

if (Ufluc(i,j,k,3).gt.25.0d0) then
print*,"wprime>25"
stop 
endif
         enddo
          enddo
           enddo
!****************************************************************************************************************

!!if(nstep.gt.15)then
!! if (my_rank.eq.0) then


!!        do j=2,n3d(2)-1
!!utau1= sqrt ((((Ume1(1,2)-Ume1(1,1))/(x3(2)-x3(1)))+lamda*((Ufluc(68,j,2,1)-Ufluc(68,j,1,1))/(x3(2)-x3(1))))/Re)
!!utau2= sqrt (((Ut(68,j,2,1)-Ut(68,j,1,1))/(x3(2)-x3(1)))/Re)
!!        end do

!!lamda=(utau1)/(utau2)
!lamda=1.0d0
!utau1=sqrt  (((usum2(1,2)-usum2(1,1))/(x3(2)-x3(1)))/Re)
!utau2=sqrt  (((usum2(68,2)-usum2(68,1))/(x3(2)-x3(1)))/Re)
!lamda1=lamda

! write(*,*) lamda,'lamda'
!!if(mod(nstep,15).eq.0) then
!! write(*,*) lamda,'lamda'
!!endif

!!if(mod(nstep,15).eq.0) then
!!      open (unit=19,file='lamda.dat',status='unknown',access='append')
!!    write(19,333) time,lamda
!! close(19)
!endif
!!333 format(2(2x,F20.10))

!!endif
!!endif


!! call  MPI_Bcast(lamda,1,mpi_double_precision,0,mpi_comm_world,ierr)

!!endif

   
!****************************************************************************************************************
!This Subroutine is used for Calculating velocity correlations
        nstat=nstat+1

        nstatinv=1.0d0/dble(nstat)
        N=1.0d0/dble((n3d(1)-2)*(n3d(2)-2))

        do k=istart,iend
        do j=2,n3d(2)-1
        do i=2,n3d(1)-1
        do jcnt=1,3
        Umean(k,jcnt)=Umean(k,jcnt)+Ut(i,j,k,jcnt)*N
        end do
        end do
        end do
        end do


   do k=istart,iend
   do ds=1,n3d(1)/2+1
!   do j=bnt,ent
   do j=2,n3d(2)-1
   do i=2,n3d(1)-1

             ids=ds
             if(i+ds.ge.n3d(1))ids=ds-n3d(1)+2
   do jcnt=1,3

! Finding Mean Dynamically 
!        Ufluc(1,jcnt)=Ut(i,j,k,jcnt)-Umean(k,jcnt)*nstatinv
!        Ufluc(2,jcnt)=Ut(i+ids,j,k,jcnt)-Umean(k,jcnt)*nstatinv

!        RUUcorr(1,jcnt,ds,k)=RUUcorr(1,jcnt,ds,k)+Ufluc(1,jcnt)*Ufluc(2,jcnt)*N

   end do
   end do
   end do
   end do
   end do


!   do k=istart,iend
!   do ds=1,n3d(2)/2+1
!   do j=bnt,ent
!   do i=2,n3d(1)-1

!             ids=ds
!             if(j+ds.ge.n3d(2))ids=ds-n3d(2)+2
!   do jcnt=1,4

!        Ufluc(1,jcnt)=Ut(i,j,k,jcnt)-Umean(k,jcnt)*nstatinv
!        Ufluc(2,jcnt)=Ut(i,j+ids,k,jcnt)-Umean(k,jcnt)*nstatinv

!        RUUcorr(2,jcnt,ds,k)=RUUcorr(2,jcnt,ds,k)+ Ufluc(1,jcnt)*Ufluc(2,jcnt)*N
!
!   end do
!   end do
!   end do
!   end do
!   end do

! if(mod(nstat,disp)==0) then
!  write(filname,'(a,i3.3,a)')'corr',my_rank+1,'.tp'
!  open(unit=20,file=filname,status='unknown')
!if (my_rank.eq.0) write(20,*)nstat,n3d(1),n3d(2),n3d(3)
!  do k=istart,iend
!  do ds=1,n3d(1)/2+1
!     write(20,12)x3(k),RUUcorr(1,1,ds,k)*nstatinv,RUUcorr(1,2,ds,k)*nstatinv,RUUcorr(1,3,ds,k)*nstatinv,RUUcorr(1,4,ds,k)*nstatinv
!  end do
!  end do

!  do k=istart,iend
!  do ds=1,n3d(2)/2+1
!     write(20,12)x3(k),RUUcorr(2,1,ds,k)*nstatinv,RUUcorr(2,2,ds,k)*nstatinv,RUUcorr(2,3,ds,k)*nstatinv,RUUcorr(2,4,ds,k)*nstatinv
!  end do
!  end do
!  close(20)
! endif
!if(my_rank.eq.0)write(*,*)Ut(2,2,en,1),'Ut8'
!if(my_rank.eq.1)write(*,*)Ut(2,2,bn-1,1),'Ut8'

!12  format(5(1x,F15.7))
!13  format(3(1X,I5))

!write(*,*)proc,'proc'


 cont=cont+1
if(mod(nstep,istat)==0) then

!         usum(:,:) =0.0d0
!         vsum(:,:) =0.0d0
!         wsum(:,:) =0.0d0
!         tsum(:,:) =0.0d0

!        uusum(:,:) =0.0d0
!        vvsum(:,:) =0.0d0
!        wwsum(:,:) =0.0d0
!        ttsum(:,:) =0.0d0

!        uvsum(:,:) =0.0d0
!        uwsum(:,:) =0.0d0
!        vwsum(:,:) =0.0d0
!        twsum(:,:) =0.0d0

  do i=1,n3d(1)
    do k=istart,iend
      

        tmpusum=0.0d0                          
        tmpvsum=0.0d0                          
        tmpwsum=0.0d0                          
        tmpTsum=0.0d0                          
!        tmpPsum=0.0d0                          
                                               
        tmpuusum=0.0d0                         
        tmpvvsum=0.0d0                         
        tmpwwsum=0.0d0
        tmpttsum=0.0d0

        tmpuvsum=0.0d0
        tmpuwsum=0.0d0
        tmpvwsum=0.0d0
        tmptwsum=0.0d0

   do j=2,n3d(2)-1
!     do i=2,n3d(1)-1
!***********************************************************************
                ip1 =i+1
!              if (i.eq.n3d(1)-1) ip1 =2
!***********************************************************************
                jp1 =j+1
              if (j.eq.n3d(2)-1) jp1 =2
!***********************************************************************
!* Averaging on the four points to find phi@ the centroid 
       su=0.50d0*(Ut(i,j,k,1)+ Ut(i,jp1,k,1))
       sv=0.50d0*(Ut(i,j,k,2)+ Ut(i,jp1,k,2))
       sw=0.50d0*(Ut(i,j,k,3)+ Ut(i,jp1,k,3))
       sT=0.50d0*(Ut(i,j,k,4)+ Ut(i,jp1,k,4))
!       sp=0.50d0*(Press(i,j,k)+ Press(i,jp1,k))
!* First order stats integrating over a plane 
        tmpusum=tmpusum+(su*dA)
        tmpvsum=tmpvsum+(sv*dA)
        tmpwsum=tmpwsum+(sw*dA)
        tmptsum=tmptsum+(sT*dA)
!        tmpPsum=tmpPsum+(sp*dA)
!* Second order stats integrating over a plane 
        tmpuusum=tmpuusum+(su*su*dA)
        tmpvvsum=tmpvvsum+(sv*sv*dA)
        tmpwwsum=tmpwwsum+(sw*sw*dA)
        tmpttsum=tmpttsum+(sT*sT*dA)

        tmpuvsum=tmpuvsum+(su*sv*dA)
        tmpuwsum=tmpuwsum+(su*sw*dA)
        tmpvwsum=tmpvwsum+(sv*sw*dA)
        tmptwsum=tmptwsum+(sT*sw*dA)

!!     enddo
   enddo

        usum(i,k) = usum(i,k) +(tmpusum/A)
        vsum(i,k) = vsum(i,k) +(tmpvsum/A)
        wsum(i,k) = wsum(i,k) +(tmpwsum/A)
        tsum(i,k) = tsum(i,k) +(tmptsum/A)
!        Psum(i,k,my_rank) = Psum(i,k,my_rank) +(tmpPsum/A/procy)

        uusum(i,k)      = uusum(i,k)  +(tmpuusum/A)
        vvsum(i,k)      = vvsum(i,k)  +(tmpvvsum/A)
        wwsum(i,k)      = wwsum(i,k)  +(tmpwwsum/A)
        ttsum(i,k)      = ttsum(i,k)  +(tmpttsum/A)

        uvsum(i,k) =uvsum(i,k) +(tmpuvsum/A)
        uwsum(i,k) =uwsum(i,k) +(tmpuwsum/A)
        vwsum(i,k) =vwsum(i,k) +(tmpvwsum/A)
        twsum(i,k) =twsum(i,k) +(tmptwsum/A)
enddo
enddo


endif

if(mod(cont,cum)==0)then

!do i=1,n3d(1)
   if (my_rank.eq.0) then
             do source = 1,proc-1 
                    sbn = 2+(source)*loc_n
                                  loc_nl=loc_n
      if(source.eq.proc-1) loc_nl=loc_nx
             call mpi_recv(usum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(vsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(wsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(tsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)

             call mpi_recv(uusum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(vvsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(wwsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(ttsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)

             call mpi_recv(uvsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(uwsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(vwsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             call mpi_recv(twsum(1,sbn),n3d(1)*loc_nl,mpi_double_precision,source,50,mpi_comm_world,status,ierr)
             end do
         else
              call mpi_send(usum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(vsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(wsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(tsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   

              call mpi_send(uusum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(vvsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(wwsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(ttsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   

              call mpi_send(uvsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(uwsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(vwsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   
              call mpi_send(twsum(1,bn),n3d(1)*loc_nx,mpi_double_precision,0,50,mpi_comm_world,ierr)   

         end if

!enddo
!if(my_rank.eq.0)write(*,*)usum(2,40),'0'
!if(my_rank.eq.6)write(*,*)usum(2,40),'in'


!********************************************************************************************************************************************************************
! Send the local solutions to processor zero.         

 if (my_rank.eq.0) then
      open(unit=2,file='means',status='unknown')
!      write(2,104,advance='no')time
         do i=1,n3d(1)
         do k=1,n3d(3)
      write(2,103,advance='no')(usum(i,k)/cum),(vsum(i,k)/cum),(wsum(i,k)/cum),(tsum(i,k)/cum)
        enddo
        enddo
      write (2, *)
      close(2)
!* Store averaged values of u^2,v^2 and w^2
     open(unit=2,file='medians',status='unknown')
!     write(2,104,advance='no')time
        do i=1,n3d(1)
        do k=1,n3d(3)
     write(2,103,advance='no')(uusum(i,k)/cum),(vvsum(i,k)/cum),(wwsum(i,k)/cum),(ttsum(i,k)/cum)
       enddo
       enddo
    write (2, *)
    close(2)
!*Store averaged values of uv,uw and vw
     open(unit=2,file='taustress',status='unknown')
!     write(2,104,advance='no')time
        do i=1,n3d(1)
        do k=1,n3d(3)
     write(2,103,advance='no')(uvsum(i,k)/cum),(uwsum(i,k)/cum),(vwsum(i,k)/cum),(twsum(i,k)/cum)
       enddo
       enddo
     write (2,*)
     close(2)
     cnt=cnt+1
     open (unit=12,file='cnt',status='unknown')
     write (12,*)n3d(3),cnt,proc
     close(12)
!endif

        do i=1,n3d(1)
        do k=1,n3d(3)


usum2(i,k)=(usum(i,k)/cum)
vsum2(i,k)=(vsum(i,k)/cum)
wsum2(i,k)=(wsum(i,k)/cum)
tsum2(i,k)=(tsum(i,k)/cum)
!Psum2(i,k)=(Psum1(i,k)/procy/cum)


uusum2(i,k)=(uusum(i,k)/cum)
vvsum2(i,k)=(vvsum(i,k)/cum)
wwsum2(i,k)=(wwsum(i,k)/cum)
ttsum2(i,k)=(ttsum(i,k)/cum)


uvsum2(i,k)=(uvsum(i,k)/cum)
uwsum2(i,k)=(uwsum(i,k)/cum)
vwsum2(i,k)=(vwsum(i,k)/cum)
twsum2(i,k)=(twsum(i,k)/cum)

        enddo
        enddo

      do i=1,n3d(1)
      do k=1,n3d(3)
       Ruu(i,k) = uusum2(i,k)-(usum2(i,k)*usum2(i,k))
       Rvv(i,k) = vvsum2(i,k)-(vsum2(i,k)*vsum2(i,k))
       Rww(i,k) = wwsum2(i,k)-(wsum2(i,k)*wsum2(i,k))
       Rtt(i,k) = ttsum2(i,k)-(tsum2(i,k)*tsum2(i,k))
       Ruv(i,k) = uvsum2(i,k)-(usum2(i,k)*vsum2(i,k))
       Ruw(i,k) = uwsum2(i,k)-(usum2(i,k)*wsum2(i,k))
       Rvw(i,k) = vwsum2(i,k)-(vsum2(i,k)*wsum2(i,k))
       Rtw(i,k) = twsum2(i,k)-(tsum2(i,k)*tsum2(i,k))
      enddo
      enddo

!    write (*,*) 'Averaging for Non-Dimensional time length ',time(cnt)-time(start) ,'starting from ' ,time(start)
      open (unit=19,file='Umeanvalues.dat',status='unknown')
       do i=1,n3d(1)
         do k=1,n3d(3)
    write(19,213) x1(i), x3(k), usum2(i,k), vsum2(i,k), wsum2(i,k),x3(k)*Re!, Psum2(i,k)
          enddo
        enddo
close(19)

      open (unit=12,file='Reyfull',status='unknown')
            i=77
            do k=1,n3d(3)
      write (12,108)x3(k),usum2(i,k),vsum2(i,k),wsum2(i,k),sqrt(Ruu(i,k)), sqrt(Rvv(i,k)),sqrt(Rww(i,k)) &
     & ,Ruv(i,k), -(Ruw(i,k)),Rvw(i,k),x3(k)*Re!,tau(k),tau(k)-Ruw(k)
            enddo
close(12)

      open (unit=12,file='Reyfull10',status='unknown')
            i=10
            do k=1,n3d(3)
      write (12,119)x3(k),usum2(i,k),vsum2(i,k),wsum2(i,k),sqrt(Ruu(i,k)), sqrt(Rvv(i,k)),sqrt(Rww(i,k)) &
     & ,Ruv(i,k), Ruw(i,k),Rvw(i,k)!,tau(k),tau(k)-Ruw(k)
            enddo
close(12)

      open (unit=12,file='Reyfull30',status='unknown')
            i=30
            do k=1,n3d(3)
      write (12,110)x3(k),usum2(i,k),vsum2(i,k),wsum2(i,k),sqrt(Ruu(i,k)), sqrt(Rvv(i,k)),sqrt(Rww(i,k)) &
     & ,Ruv(i,k), Ruw(i,k),Rvw(i,k)!,tau(k),tau(k)-Ruw(k)
            enddo
close(12)

!utau1=sqrt(( ( ((usum2(1,3)-usum2(1,1))*(df(1)**2)) + ((usum2(1,1)-usum2(1,2))*(df(2)**2)) )/( (df(2)*(df(1)**2))-((df(2)**2)*df(1)) ))/Re)
utau1=sqrt  (((usum2(1,2)-usum2(1,1))/(x3(2)-x3(1)))/Re)
utau2=sqrt  (((usum2(68,2)-usum2(68,1))/(x3(2)-x3(1)))/Re)
!utau2=sqrt(( ( ((usum2(68,3)-usum2(68,1))*(df(1)**2)) + ((usum2(68,1)-usum2(68,2))*(df(2)**2)) )/( (df(2)*(df(1)**2))-((df(2)**2)*df(1)) ))/Re)

lamda=(utau1)/(utau2)
!lamda=1.0d0
!lamda1=lamda

 write(*,*) lamda,'lamda'

!if(mod(nstep,15).eq.0) then
!      open (unit=19,file='lamda.dat',status='unknown',access='append')
!    write(19,333) time,lamda
! close(19)
!endif
!333 format(2(2x,F20.10))
!if(lamda1.gt.1.0d0)then
!lamda1=1.0d0
!endif


do i=1,n3d(1)
do k=1,n3d(3)
    if(usum2(i,k).ge.(0.99d0*Ue)) then
     bltkns(i) = x3(k)
   go to 109 
endif
enddo
109 continue 
enddo


lamda03=bltkns(68)/bltkns(1)

!if(lamda03.gt.1.3d0)then
!lamda03=1.3d0
!endif

!if(mod(nstep,15).eq.0) then
      open (unit=19,file='lamda.dat',status='unknown',access='append')
    write(19,333) time,lamda,lamda03
 close(19)
!endif
333 format(3(2x,F20.10))

do i=1,n3d(1)
do k=1,n3d(3)
Wt(i,k)=0.5d0*(1+(dtanh((4.0d0*((x3(k)/bltkns(i))-0.2d0))/((1-2.0d0*0.2d0)*(x3(k)/bltkns(i))+0.2d0))/dtanh(4.0d0)))
enddo
enddo

      open (unit=19,file='W.dat',status='unknown')
       do i=1,n3d(1)
         do k=1,n3d(3)
    write(19,314) x1(i), x3(k), Wt(i,k)!, vsum2(i,k), wsum2(i,k), Psum2(i,k)
          enddo
        enddo
close(19)



new(:)=n3d(3)-1

!do i=1,n3d(1)
do k=1,n3d(3)-1
if(Wt(68,k).le.0.99d0) then
ksep=k
eta(k)= lamda
!eta(k)= 1.0d0
else
eta(k)= lamda03
!eta(k)=1.3d0
endif
enddo
!enddo


do k=1,n3d(3)-1
zrecy(k)=eta(k)*x3(k)
enddo

counter=1
do k=1,n3d(3)-1
if(counter.le.(ksep)) then
if(x3(k).ge.zrecy(counter))then
new(counter)=k
counter= counter+1
endif
endif
enddo


counter=ksep+1
do k=ksep+1,n3d(3)-1
if(counter.le.(n3d(3)-1)) then
if(x3(k).ge.zrecy(counter))then
new(counter)=k
counter= counter+1
endif
endif
enddo


do i=1,n3d(1)
do k=1,n3d(3)-1
deltastar(i)=deltastar(i)+((1.0d0-(usum2(i,k)/Ue))*(x3(k+1)-x3(k)))
if(x3(k).eq.bltkns(i)) go to 144
enddo
144 continue
enddo

do i=1,n3d(1)
do k=1,n3d(3)-1
theta(i)=theta(i)+((usum2(i,k)/Ue)*(1.0d0-(usum2(i,k)/Ue))*(x3(k+1)-x3(k)))
if(x3(k).ge.bltkns(i)) go to 155
enddo
155 continue
enddo

do i= 1,n3d(1)
Ht(i)=deltastar(i)/theta(i)
enddo

do i= 1,n3d(1)
Reth(i)=theta(i)*Re*Ue
enddo


!if(mod(nstep,100).eq.0) then
 open (unit=19,file='boundarylthickness.dat',status='unknown')
       do i=1,n3d(1)
    write(19,433) time,x1(i), bltkns(i), deltastar(i), theta(i), Ht(i), Reth(i)
       enddo
 close(19)
!endif

433 format(7(2x,F20.10))

               
            theta(:)=0.0d0 
        deltastar(:)=0.0d0   
!         usum2(:,:) =0.0d0
!         vsum2(:,:) =0.0d0
!         wsum2(:,:) =0.0d0
!         tsum2(:,:) =0.0d0

!        uusum2(:,:) =0.0d0
!        vvsum2(:,:) =0.0d0
!        wwsum2(:,:) =0.0d0
!        ttsum2(:,:) =0.0d0

!        uvsum2(:,:) =0.0d0
!        uwsum2(:,:) =0.0d0
!        vwsum2(:,:) =0.0d0
!        twsum2(:,:) =0.0d0

if(mod(nstep,100).eq.0) then
 open (unit=19,file='k.dat',status='unknown',access='append')
do k=1,n3d(3)
    write(19,733) k, new(k),eta(k) 
 enddo
733 format(I5,2(F20.10))
endif

if(mod(nstep,100).eq.0) then
      open (unit=19,file='meanvalues.dat',status='unknown',access='append')
       do i=1,n3d(1)
         do k=1,n3d(3)
    write(19,313) x1(i), x3(k), usum2(i,k), vsum2(i,k), wsum2(i,k)!, Psum2(i,k)
          enddo
        enddo
close(19)
endif


endif

! call  MPI_Bcast(lamda,1,oneplan,0,mpi_comm_world,ierr)

  call  MPI_Bcast(lamda,1,mpi_double_precision,0,mpi_comm_world,ierr)
! call  MPI_Bcast(lamda1,1,mpi_double_precision,0,mpi_comm_world,ierr)
 call  MPI_Bcast(new,n3d(3),mpi_double_precision,0,mpi_comm_world,ierr)
 call  MPI_Bcast(bltkns,n3d(1),mpi_double_precision,0,mpi_comm_world,ierr)
 call  MPI_Bcast(usum2,(n3d(1)-1)*n3d(3),mpi_double_precision,0,mpi_comm_world,ierr)
 call  MPI_Bcast(vsum2,(n3d(1)-1)*n3d(3),mpi_double_precision,0,mpi_comm_world,ierr)
 call  MPI_Bcast(wsum2,(n3d(1)-1)*n3d(3),mpi_double_precision,0,mpi_comm_world,ierr)
 call  MPI_Bcast(Wt,(n3d(1)-1)*n3d(3),mpi_double_precision,0,mpi_comm_world,ierr)


 do k=bn,en
!     do j=bnt,ent
      do j=2,n3d(2)-1
  Ufluc(68,j,k,1) = Ut(68,j,k,1)-usum2(68,k)
  Ufluc(68,j,k,2) = Ut(68,j,k,2)-vsum2(68,k)
  Ufluc(68,j,k,3) = Ut(68,j,k,3)-wsum2(68,k)
  
      enddo
     enddo

         usum(:,:) =0.0d0
         vsum(:,:) =0.0d0
         wsum(:,:) =0.0d0
         tsum(:,:) =0.0d0


        uusum(:,:) =0.0d0
        vvsum(:,:) =0.0d0
        wwsum(:,:) =0.0d0
        ttsum(:,:) =0.0d0

        uvsum(:,:) =0.0d0
        uwsum(:,:) =0.0d0
        vwsum(:,:) =0.0d0
        twsum(:,:) =0.0d0


endif

      cnt2= cnt2+1
     open (unit=23,file='cnt2',status='unknown')
     write (23,*)cnt2
     close(23)
!********************************************************************************************************************************************************************
 end do ! while loop

01    format(8X,F20.10,I9,1X)
02    format(7X,I4,3X,I4,3X,I6)
03    format(A,F20.10,I9,A)
04    format(A,I4,A,I4,A,I4,A)
97    FORMAT(3(1X,F10.6),11(1X,E22.15),(1X,I4))
103   format(4(1X,F15.7))
213   format(6(1X,F15.7))
313   format(5(1X,F15.7))
314   format(3(1X,F15.7))
!104   format(1(1X,F15.7))
99    format(6(1X,F15.7))
21    format((1X,F14.7,1X,I9,1X,F14.7))
94    format(4(1X,F15.7))
93    format(4(1X,F15.7))
108   format(11(1X,F15.7))
119   format(10(1X,F15.7))
110   format(10(1X,F15.7))
208   format(13(1X,F15.7))

    deallocate(rAp)
    deallocate(v)
    deallocate(r)
    deallocate(rhat)
    deallocate(h)
    deallocate(g)
    deallocate(c)
    deallocate(s)
    deallocate(y)
    deallocate(mag)
    deallocate(Aw)
    deallocate(An)
    deallocate(Ae)
    deallocate(As)
    deallocate(Ab)
    deallocate(At)
    deallocate(Ap)
    deallocate(Ut)
    deallocate(Umrf)
    deallocate(Ut)
    deallocate(Utold)
    deallocate(PCorr)
    deallocate(Press)
    deallocate(Q)
    deallocate(x1)
    deallocate(x2)
    deallocate(x3)
    deallocate(Conv)
    deallocate(Jac)
    deallocate(Xi)
    deallocate(Jac2)
    deallocate(tmp2)
    deallocate(uusum)
    deallocate(vvsum)
    deallocate(wwsum)
    deallocate(ttsum)
    deallocate(uvsum)
    deallocate(uwsum)
    deallocate(vwsum)
    deallocate(usum)
    deallocate(vsum)
    deallocate(wsum)
    deallocate(tsum)
    deallocate(Psum)
    deallocate(Ume1)
    deallocate(twsum)
    deallocate(mfb1)
    deallocate(mfb2)
    deallocate(cf)
    deallocate(eta)
    deallocate(Umr)
    deallocate(Umr1)
    deallocate(Reth)
    deallocate(Udist)

    call mpi_type_free(oneplane,ierr)
    call mpi_type_free(twoplanes,ierr)
!    call mpi_type_free(threevar,ierr)
    call mpi_finalize(ierr)
!

end
