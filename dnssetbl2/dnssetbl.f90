PROGRAM revamptbl
!=======================================================================
!Student: Parvez Ahmad, Course: M.Tech
!3D Parallel DNS of Turbulent Boundary Layer over a Flat Plate
!CFD Lab,ZHCET,Aligarh Muslim University
!Under the Guidance of Prof Mirza Faisal S. Baig
!Last Modified on 02/10/16
!=======================================================================
USE mpi
IMPLICIT NONE
!=======================================================================
INTEGER, PARAMETER			::	prmpt		=1,&		!1=>Restart from file,2=>From t=0
								nX			=90,&
								nY			=80,&
								nZ			=50,&
								nstart		=51,&
								maxstep		=30000000,&	!Max no of time iterations
								freqDiv		=1,&   		!Freq of divergence & massflow disp.
								freqStat	=1,&      	!Freq of stat calc
								freqTimAvg	=500,&		!Freq of time averaging
								freqWrite	=500,&		!Freq of Soln file write
								freqLog		=500,&		!Freq of log file write
								kmax		=50,&
								mmax		=1000
						
DOUBLE PRECISION, PARAMETER	::	errorstop	=1.0d-6,& 	!Conver Crit of Transport Eqns
								dt			=1.0d-4,&   !timestep
								Re 			=150.92d0,&	!Based on utau
								delta		=1.95d0,&	!Grid Condensation Param
								Ue			=19.798d0,&	!Freestream Vel
								beta		=10.0d0,&
								gama		=10.0d0,&
								errtol		=1.0d-5,&
								w			=1.835
								
CHARACTER(LEN=25)           ::filname
LOGICAL						::isFile
INTEGER						::nstep, counter, i, j, k, new(nZ), it, threevar, ksep
INTEGER 					::p, nproc, id, ierr, it1, it2, cnt
INTEGER                     ::STATUS(mpi_status_size), oneplane, twoplanes, kit, m
DOUBLE PRECISION			::df(2), gradp(3), massLoc,mass(5),x2(nY)
DOUBLE PRECISION			::rt1, rt2, rt3, rt4, rt5, rt6, rt11, rt12
DOUBLE PRECISION			::time, dXi, dx, dy, invdx, invdy, invdXi, onetwelfth
DOUBLE PRECISION			::kinEnLoc, kinEn, Pe
DOUBLE PRECISION			::lamda, utau1, utau2, eps, lamda03, dAv, A, su, sv, sw
DOUBLE PRECISION           	::rhoLoc, tempgmres, rho, hik, hipk, nu, gk, gkp
DOUBLE PRECISION          	::sumav, sumnor, resLoc, tmp, res, v(nX, nY, nZ, nstart), h(nstart, nstart)
DOUBLE PRECISION 			::tmpusum, tmpvsum, tmpwsum, tmpuusum, tmpvvsum, tmpwwsum, tmpuvsum, tmpuwsum, tmpvwsum
DOUBLE PRECISION 			::Aw, Ae, An, As
INTEGER, DIMENSION(5)     	::plane=[5, 20, 40, 60, 85]								
INTEGER, DIMENSION(nY)		::jp1, jp2, jm1, jm2
DOUBLE PRECISION, DIMENSION(nX)				::x1, theta, deltastar, Ht, Reth, Rex, bltkns
DOUBLE PRECISION, DIMENSION(nZ)				::Ap, At, Ab, rAp
DOUBLE PRECISION, DIMENSION(nZ)				::x3, Jac, JacInv, Jac2, Xi, Ume, zrecy, eta
DOUBLE PRECISION, DIMENSION(nstart) 		:: g, c, s, y, mag
DOUBLE PRECISION, DIMENSION(nX, nZ)			::Utmeaninn, Vtmeaninn, Wtmeaninn, Utmeanout, Vtmeanout, Wtmeanout, Wt
DOUBLE PRECISION, DIMENSION(nX, nZ)			::uusum, vvsum, wwsum, uvsum, uwsum, vwsum, usum, vsum, wsum
DOUBLE PRECISION, DIMENSION(nX, nZ)			::uusum2, vvsum2, wwsum2, uvsum2, uwsum2, vwsum2, usum2, vsum2, wsum2
DOUBLE PRECISION, DIMENSION(nX, nZ)			::Ruu, Rvv, Rww, Ruv, Ruw, Rvw
DOUBLE PRECISION, DIMENSION(nX, nY, nZ)		::pRhs, PCorr, Press, r, rhat
DOUBLE PRECISION, DIMENSION(nX, nY, nZ, 3)	::Ut, Utold, Ufluc, Uflucinn, Uflucout, conv, Q
INTEGER, ALLOCATABLE, DIMENSION(:)			::b,e,bn,en,siz

!==============================MPI====================================== 
CALL mpi_init(ierr)
CALL mpi_comm_rank(mpi_comm_world,id,ierr)
CALL mpi_comm_size(mpi_comm_world,nproc,ierr)

CALL mpi_type_contiguous(  nX*nY,mpi_double_precision,oneplane,ierr)
CALL mpi_type_contiguous(2*nX*nY,mpi_double_precision,twoplanes,ierr)
CALL MPI_TYPE_VECTOR(3,2*nX*nY,nX*nY*nZ,mpi_double_precision,threevar,ierr)
CALL mpi_type_commit(oneplane   ,ierr)
CALL mpi_type_commit(twoplanes  ,ierr)
CALL mpi_type_commit(threevar,ierr)

ALLOCATE(b(0:nproc-1),e(0:nproc-1),bn(0:nproc-1),en(0:nproc-1),siz(0:nproc-1))

it1=nZ/nproc
it2=nproc-mod(nZ,nproc)

b(0)=1
do i=0,nproc-1
	if(i==it2) it1=it1+1
	e(i)=b(i)+it1-1
	b(i+1)=e(i)+1
	siz(i)=e(i)-b(i)+1
enddo

bn(:)=b(:)
en(:)=e(:)
bn(0)=b(0)+1
en(nproc-1)=e(nproc-1)-1

CALL mpi_barrier(mpi_comm_world,ierr)
!-----------------------------MPI---------------------------------------
!=======================================================================
if(id==0) then
INQUIRE(directory='Output', EXIST=isFile)
IF (.not.(isFile)) CALL SYSTEM('mkdir Output')

INQUIRE(directory='Output/backup', EXIST=isFile)
IF (.not.(isFile)) CALL SYSTEM('mkdir Output/backup')
endif

!=========================Read Initial Profile==========================
OPEN (UNIT=15,FILE='Input/initProf.dat',status='old') 
DO k=1,37
	READ(15,103) Ume(k)
END DO
CLOSE(15)
103 FORMAT (24x,E18.8)

DO k=38,nZ
	Ume(k) = Ume(37)
END DO
!-----------------------------------------------------------------------
!==============================IC Setup=================================
SELECT CASE(prmpt)
	CASE(1)
		OPEN (UNIT=12,FILE='Input/3D.dat',STATUS='old')
		READ (12,101)time,nstep
		READ (12,*)
		READ (12,*)
		DO i=1,nX
			READ(12,*)
			DO j=1,nY
				READ(12,*)
				DO k=1,nZ
					READ(12,104)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Press(i,j,k)
				END DO
			END DO
		END DO
		CLOSE(12)
		101 FORMAT(8X,F20.10,I9,1X)
		104 FORMAT(3(1X,F10.6),4(1X,E22.15))
		
		OPEN (UNIT=20,FILE='Input/statData.dat')
		DO k=1,nZ
			READ(20,115) usum2(68,k),vsum2(68,k),wsum2(68,k),Wt(1,k),new(k)
		END DO
		CLOSE(20)
		115 FORMAT(4(1X,F15.7),I5.3)
		
		DO k=bn(id),en(id)
			DO j=2,nY-1
				Ufluc(68,j,k,1) = Ut(68,j,k,1)-usum2(68,k)
				Ufluc(68,j,k,2) = Ut(68,j,k,2)-vsum2(68,k)
				Ufluc(68,j,k,3) = Ut(68,j,k,3)-wsum2(68,k)
			END DO
		END DO
	CASE(2)
		nstep=0
		time=0.0
		!=========================Initial Flow Field====================
		do p=1,3
			do k=1,nZ
				do j=1,nY
					do i=1,nX
						Ut(i,j,k,1)= Ume(k)+(0.01d0*Ue*sin(beta*x2(j))*cos(gama*x3(k)))
						Ut(i,j,k,2)=        (0.01d0*Ue*sin(beta*x2(j))*cos(gama*x3(k)))
						Ut(i,j,k,3)=        (0.01d0*Ue*sin(beta*x2(j))*cos(gama*x3(k)))
					enddo
			 	enddo
		 	enddo
		enddo
		!---------------------------------------------------------------
END SELECT
!-----------------------------------------------------------------------
if(id==0) then
OPEN(UNIT=25,FILE='Output/simLog.dat',POSITION='append')
WRITE(25,*)
WRITE(25,'(2(a))')'Simulation date and Time :',dateTime()
WRITE(25,'(a,i2.2)') 'Restart Status(1=>Restart,2=>t=0) :',prmpt
WRITE(25,'(a,f10.6,3x,i10.10))') 'Time & No of Iter :',time,nstep
WRITE(25,'(a,i4.3)') 'No of processes started :',nproc
WRITE(25,'(a,3(i4.3,3x))') 'Mesh Dimensions :',nX,nY,nZ
WRITE(25,'(a,2(f10.6,3x))') 'Reynolds No. & Timestep :',Re,dt
WRITE(25,'(a,i5.5,a)') 'Period of time average :',freqTimAvg,' iterations'
WRITE(25,*)
CLOSE(25)
endif
!===========================Mesh Setup Starts===========================
x1(1) = 0.0d0
x2(1) = 0.0d0
Xi(1) =0.0d0
x3(1) =0.0d0
x1(nX) = 12.0d0
x2(nY) =  4.0d0 !2.0d0
x3(nZ) =  3.4d0

dx = x1(nX)/dble(nX-1)
dy = x2(nY)/dble(nY-1)
dXi= x3(nZ)/dble(nZ-1)

invdx = 1.0d0/dx
invdy = 1.0d0/dy
invdXi= 1.0d0/dXi

DO i=1,nX-1
	x1(i+1)=x1(i)+dx
end do

DO i=1,nY-1
	x2(i+1)=x2(i)+dy
end do

DO k=1,nZ-1
	Xi(k+1) =Xi(k) + dXi
END DO

!Mapping of Non-Uniform Points from Uniform Points
DO k=1,nZ
	x3(k) = 3.4d0 + 3.4d0 * tanh(delta*0.50*(Xi(k)-3.4)) / (tanh(1.7d0*delta))
	Jac(k) = 1.7d0*delta / ( tanh(1.7*delta) * (cosh(0.5*delta*(Xi(k)-3.4)))**2 )
	Jac2(k) = -1.7d0*delta**2 * tanh(0.5*delta*(Xi(k)-3.4)) / (tanh(1.7*delta) * (cosh(0.5*delta*(Xi(k)-3.4)))**2)
END DO
!---------------------------Mesh Setup Ends-----------------------------
!========================Invariant Expression===========================

df(1)=x3(2)-x3(1)
df(2)=x3(3)-x3(1)

Aw= 1.0d0/(dx*dx)
An= 1.0d0/(dy*dy)
Ae= 1.0d0/(dx*dx)
As= 1.0d0/(dy*dy)
DO k=bn(id),en(id)
	Ab(k)=1.0d0/(dXi*dXi*Jac(k)*Jac(k)) + Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k))
	At(k)=1.0d0/(dXi*dXi*Jac(k)*Jac(k)) - Jac2(k)/(2.0*dXi*Jac(k)*Jac(k)*Jac(k))
	Ap(k)=-2.0d0*((1.0d0/(dx*dx)) +(1.0d0/(dy*dy))+(1.0/(dXi*dXi*Jac(k)*Jac(k))))
	rAp(k)= 1.0d0/Ap(k)
END DO

DO j=2,nY-1	
	jm1(j) =j-1
	IF (j.eq.2) jm1(j) =nY-1
	jm2(j) =jm1(j)-1
	IF (j.eq.3) jm2(j) =nY-1
	jp1(j) =j+1
	IF (j.eq.nY-1) jp1(j) =2
	jp2(j) =jp1(j)+1
	IF (j.eq.nY-2) jp2(j) =2
ENDDO

DO k=b(id),e(id)
	JacInv(k)= 1.0d0/Jac(k)
END DO

onetwelfth =1.0d0/12.0d0
rhat = 0.0d0
cnt=0
!-----------------------------------------------------------------------
DO WHILE(nstep < maxstep)
	time = time+dt
	nstep = nstep+1
	cnt=cnt+1
	Utold(:,:,:,:)=Ut(:,:,:,:)
	
	eps=nstep*(5.0d-5)
	IF(eps.gt.1.0d0)THEN
		eps=1.0d0
	END IF

	IF(nstep.gt.freqTimAvg) THEN
!=========================Inflow BC=====================================	
		DO k=bn(id),en(id)
			IF(usum2(68,k).lt.(0.99d0*Ue)) THEN
				Utmeaninn(1,k)= eps*(usum2(68,new(k)))+(1-eps)*Ume(k)
				!Utmeaninn(1,k)= eps*lamda*(usum2(68,new(k)))+(1-eps)*Ume(k)
				!Utmeaninn(1,k)= Ume(k)
				Vtmeaninn(1,k)=lamda*vsum2(68,new(k))
				Wtmeaninn(1,k)=wsum2(68,new(k))
				
				Utmeanout(1,k)=eps*(usum2(68,new(k)))+(1-eps)*Ume(k)
				!Utmeanout(1,k)=eps*lamda*(usum2(68,new(k)))+(1-eps)*Ume(k)
				!Utmeanout(1,k)=Ume(k)
				Vtmeanout(1,k)=lamda*vsum2(68,new(k))
				Wtmeanout(1,k)=wsum2(68,new(k))
				
				DO j=2,nY-1
					Uflucinn(1,j,k,1)=lamda*Ufluc(68,j,new(k),1)
					Uflucinn(1,j,k,2)=lamda*Ufluc(68,j,new(k),2)
					Uflucinn(1,j,k,3)=lamda*Ufluc(68,j,new(k),3)
				END DO

				DO j=2,nY-1
					Uflucout(1,j,k,1)=lamda*Ufluc(68,j,new(k),1)
					Uflucout(1,j,k,2)=lamda*Ufluc(68,j,new(k),2)
					Uflucout(1,j,k,3)=lamda*Ufluc(68,j,new(k),3)
				END DO

				DO j=2,nY-1
					Ut(1,j,k,1)=(Utmeaninn(1,k)+Uflucinn(1,j,k,1))*(1.0d0-Wt(1,k))+(Utmeanout(1,k)+Uflucout(1,j,k,1))*Wt(1,k)
					Ut(1,j,k,2)=(Vtmeaninn(1,k)+Uflucinn(1,j,k,2))*(1.0d0-Wt(1,k))+(Vtmeanout(1,k)+Uflucout(1,j,k,2))*Wt(1,k)
					Ut(1,j,k,3)=(Wtmeaninn(1,k)+Uflucinn(1,j,k,3))*(1.0d0-Wt(1,k))+(Wtmeanout(1,k)+Uflucout(1,j,k,3))*Wt(1,k)
				END DO
			ELSE
				Ut(1,2:nY-1,k,1)=usum2(68,k)
				Ut(1,2:nY-1,k,2)=vsum2(68,k)
				Ut(1,2:nY-1,k,3)=wsum2(68,k)
			END IF
		END DO
	END IF
!-------------------------------Inflow BC-------------------------------	
!=======================Discretize Convective Terms=====================
	DO p=1,3
		DO k=bn(id),en(id)
			DO j=2,nY-1
				DO i=2,nX-1
				!-------------------------------------------------------
				Pe = abs(Re*Ut(i,j,k,1)*dx)
				
				IF(i.le.2.or.i.ge.(nX-1)) THEN
					rt1 = Ut(i,j,k,1)*(Ut(i+1,j,k,p)-Ut(i-1,j,k,p))*0.5d0*invdx
				ELSEIF(Pe.gt.2.0d0) THEN
					rt11 = Abs(Ut(i,j,k,1))*(Ut(i+2,j,k,p)-(4*(Ut(i+1,j,k,p)))+(6*(Ut(i,j,k,p)))-(4*(Ut(i-1,j,k,p)))+Ut(i-2,j,k,p))*(0.25d0*invdx)
					rt12 =     (Ut(i,j,k,1))*(-Ut(i+2,j,k,p)+(8*(Ut(i+1,j,k,p)))-(8*(Ut(i-1,j,k,p)))+Ut(i-2,j,k,p))*(onetwelfth*invdx)
					rt1= rt11+rt12
				ELSE
					rt1=Ut(i,j,k,1)*((Ut(i-2,j,k,p)-8*Ut(i-1,j,k,p)+8*Ut(i+1,j,k,p)-Ut(i+2,j,k,p))*(onetwelfth*invdx))
				END IF
				!-------------------------------------------------------
				Pe = abs(Re*Ut(i,j,k,2)*dy)

				IF((j.eq.2).or.(j.eq.(nY-1)))THEN
					rt2 = Ut(i,j,k,2)*(Ut(i,jp1(j),k,p)-Ut(i,jm1(j),k,p))*0.5d0*invdy
				ELSEIF(Pe.gt.2.0d0) THEN
					rt11 = (Abs(Ut(i,j,k,2)))*(Ut(i,jp2(j),k,p)-(4*(Ut(i,jp1(j),k,p)))+(6*(Ut(i,j,k,p)))-(4*(Ut(i,jm1(j),k,p)))+Ut(i,jm2(j),k,p))*(0.25d0*invdy)
					rt12 =     (Ut(i,j,k,2))*(-Ut(i,jp2(j),k,p)+(8*(Ut(i,jp1(j),k,p)))-(8*(Ut(i,jm1(j),k,p)))+Ut(i,jm2(j),k,p))*(onetwelfth*invdy)
					rt2 = rt11+rt12
				ELSE
					rt2=Ut(i,j,k,2)*((Ut(i,jm2(j),k,p)-8*Ut(i,jm1(j),k,p)+8*Ut(i,jp1(j),k,p)-Ut(i,jp2(j),k,p))*(onetwelfth*invdy))
				END IF
				!-------------------------------------------------------
				Pe = abs(Re*Ut(i,j,k,3)*(x3(k+1)-x3(k)))

				IF( (k.le.2).or.(k.ge.(nZ-1)) ) THEN
					rt3 = (Ut(i,j,k,3)*(Ut(i,j,k+1,p)-Ut(i,j,k-1,p))*(0.5d0*invdXi))
				ELSEIF(Pe.gt.2.0d0) THEN
					rt11 = (Abs(Ut(i,j,k,3)))*(Ut(i,j,k+2,p)-(4*(Ut(i,j,k+1,p)))+(6*(Ut(i,j,k,p)))-(4*(Ut(i,j,k-1,p)))+Ut(i,j,k-2,p))*(0.25d0*invdXi)
					rt12 =     (Ut(i,j,k,3))*(-Ut(i,j,k+2,p)+(8*(Ut(i,j,k+1,p)))-(8*(Ut(i,j,k-1,p)))+Ut(i,j,k-2,p))*(onetwelfth*invdXi)
					rt3 = (rt11+rt12)
				ELSE
					rt3=Ut(i,j,k,3)*((Ut(i,j,k-2,p)-8*Ut(i,j,k-1,p)+8*Ut(i,j,k+1,p)-Ut(i,j,k+2,p))*(onetwelfth*invdXi))
				END IF

				conv(i,j,k,p) = rt1+rt2+(rt3*JacInv(k))

				END DO
			END DO
		END DO
	END DO
!----------------------Discretize Convective Terms----------------------
!=======================Discretize Pressure Terms=======================
	DO k=bn(id),en(id)
		DO j=2,nY-1
			DO i=2,nX-1
			gradp(1)=0.5d0*(Press(i+1,j,k)-Press(i-1,j,k))*(invdx)
			gradp(2)=0.5d0*(Press(i,jp1(j),k)-Press(i,jm1(j),k))*(invdy)
			gradp(3)=0.5d0*(Press(i,j,k+1)-Press(i,j,k-1))*(invdXi*JacInv(k))
!-----------------------------------------------------------------------
!==========================Construct RHS================================
			Q(i,j,k,1)= Ut(i,j,k,1)-dt*(conv(i,j,k,1)+gradp(1))
			Q(i,j,k,2)= Ut(i,j,k,2)-dt*(conv(i,j,k,2)+gradp(2))
			Q(i,j,k,3)= Ut(i,j,k,3)-dt*(conv(i,j,k,3)+gradp(3))
!-----------------------------------------------------------------------			
			END DO
		END DO
	END DO

	DO p =1,3 !!GS iteration
		DO it=1,mmax
			tmp=0.0d0
			IF (mod(it,2).eq.1) THEN  ! check for residual every 2 times
				DO k=bn(id),en(id)
					DO j=2,nY-1
						DO i=2,nX-1
						
						rt1=Ae*Ut(i+1,j,k,p) + Aw*Ut(i-1,j,k,p) + An*Ut(i,jp1(j),k,p) + As*Ut(i,jm1(j),k,p)
						rt2=At(k)*Ut(i,j,k+1,p) + Ab(k)*Ut(i,j,k-1,p)
						
						resLoc=Q(i,j,k,p) - (1.0d0 - dt/Re*Ap(k))*Ut(i,j,k,p) + dt*(rt1+rt2)/Re
						tmp=tmp+abs(resLoc)
						END DO
					END DO
				END DO

				CALL mpi_allreduce(tmp,res,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
				
				IF(it.eq.1)sumnor=res
				sumav=res/sumnor
				IF(sumav.lt.errorstop .or. res.lt.1e-12) GOTO 191
			END IF

			DO k=bn(id),en(id)
				DO j=2,nY-1
					DO i=2,nX-1
					
					rt1=Ae*Ut(i+1,j,k,p) + Aw*Ut(i-1,j,k,p) + An*Ut(i,jp1(j),k,p) + As*Ut(i,jm1(j),k,p)
					rt2=At(k)*Ut(i,j,k+1,p) + Ab(k)*Ut(i,j,k-1,p)
					
					Ut(i,j,k,p)=(Q(i,j,k,p) + dt/Re*(rt1+rt2))/(1.0d0-dt/Re*Ap(k))

					END DO
				END DO
			END DO

			IF (id.eq.0) THEN
				CALL mpi_recv(Ut(1,1,en(id)+1,p),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(Ut(1,1,en(id)-1,p),1,twoplanes,id+1,50,mpi_comm_world,ierr)
			END IF
			
			IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.1)) THEN
				CALL mpi_send(Ut(1,1,en(id)-1,p),1,twoplanes,id+1,50,mpi_comm_world,ierr)
				CALL mpi_recv(Ut(1,1,en(id)+1,p),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(Ut(1,1,bn(id)  ,p),1,twoplanes,id-1,50,mpi_comm_world,ierr)
				CALL mpi_recv(Ut(1,1,bn(id)-2,p),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
			END IF
			
			IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.0)) THEN
				CALL mpi_recv(Ut(1,1,bn(id)-2,p),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(Ut(1,1,bn(id)  ,p),1,twoplanes,id-1,50,mpi_comm_world,ierr)
				CALL mpi_recv(Ut(1,1,en(id)+1,p),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(Ut(1,1,en(id)-1,p),1,twoplanes,id+1,50,mpi_comm_world,ierr)
			END IF
			
			IF (id.eq.nproc-1) THEN
				CALL mpi_send(Ut(1,1,bn(id)  ,p),1,twoplanes,id-1,50, mpi_comm_world,ierr)
				CALL mpi_recv(Ut(1,1,bn(id)-2,p),1,twoplanes,id-1,50, mpi_comm_world,STATUS,ierr)
			END IF		
		END DO
		191    CONTINUE
	END DO

	!First,exchange Updated velocity vector field between processors for 2 planes.
	!=======================Construct Divergence========================
	DO k=bn(id),en(id)
		DO j=2,nY-1
			DO i=2,nX-1
				rt1= (((Ut(i+1,j,k,1) + Ut(i,j,k,1))*0.5d0)-(((Ut(i,j,k,1) + Ut(i-1,j,k,1))*0.5d0)))*invdx
				rt2=(((Ut(i,jp1(j),k,2) + Ut(i,j,k,2))*0.5d0)-(((Ut(i,j,k,2) + Ut(i,jm1(j),k,2))*0.5d0)))*invdy
				rt3=2.0*(((Ut(i,j,k+1,3) + Ut(i,j,k,3))*0.5d0) - (((Ut(i,j,k,3) + Ut(i,j,k-1,3))*0.5d0)))/(x3(k+1)-x3(k-1))

				pRhs(i,j,k)=(rt1+rt2+rt3)/dt
			END DO
		END DO
	END DO
	!-------------------------------------------------------------------
	!====================Solution of PPE by GMRES=======================
	!Inflow & Outflow BC
	Pcorr (  1, 1:nY ,b(id):e(id)) =0.0d0
	Pcorr ( nX, 1:nY ,b(id):e(id)) =0.0d0
	
	!Spanwise BC
	Pcorr ( 1:nX,  1, b(id):e(id)) =0.0d0
	Pcorr ( 1:nX, nY, b(id):e(id)) =0.0d0

	!Bottom & Top BC
	IF (id.eq.0)         PCorr(1:nX,1:nY,bn(id)-1 ) =0.0d0
	IF (id.eq.nproc-1)    PCorr(1:nX,1:nY,en(id)+1 ) =0.0d0

	sumav	= 1.0d0
	m		= 0
	DO WHILE ((sumav>errtol).and.(m<mmax)) !Begin restart loop.
		m = m+1
		h = 0.0d0
		v = 0.0d0
		c = 0.0d0
		s = 0.0d0
		g = 0.0d0
		y = 0.0d0

		! Matrix vector product for the initial residual.
		IF (id.eq.0) THEN
			CALL mpi_recv(PCorr(1,1,en(id)+1),1,oneplane,id+1,50,mpi_comm_world,STATUS,ierr)
			CALL mpi_send(PCorr(1,1,en(id)  ),1,oneplane,id+1,50,mpi_comm_world,ierr)
		END IF
		
		IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.1)) THEN
			CALL mpi_send(PCorr(1,1,en(id)  ),1,oneplane,id+1,50,mpi_comm_world,ierr)
			CALL mpi_recv(PCorr(1,1,en(id)+1),1,oneplane,id+1,50,mpi_comm_world,STATUS,ierr)
			CALL mpi_send(PCorr(1,1,bn(id)  ),1,oneplane,id-1,50,mpi_comm_world,ierr)
			CALL mpi_recv(PCorr(1,1,bn(id)-1),1,oneplane,id-1,50,mpi_comm_world,STATUS,ierr)
		END IF

		IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.0)) THEN
			CALL mpi_recv(PCorr(1,1,bn(id)-1),1,oneplane,id-1,50,mpi_comm_world,STATUS,ierr)
			CALL mpi_send(PCorr(1,1,bn(id)  ),1,oneplane,id-1,50,mpi_comm_world,ierr)
			CALL mpi_recv(PCorr(1,1,en(id)+1),1,oneplane,id+1,50,mpi_comm_world,STATUS,ierr)
			CALL mpi_send(PCorr(1,1,en(id)  ),1,oneplane,id+1,50,mpi_comm_world,ierr)
		END IF

		IF (id.eq.nproc-1) THEN
			CALL mpi_send(PCorr(1,1,bn(id)  ),1,oneplane,id-1,50,mpi_comm_world,ierr)
			CALL mpi_recv(PCorr(1,1,bn(id)-1),1,oneplane,id-1,50,mpi_comm_world,STATUS,ierr)
		END IF

		DO k=bn(id),en(id)
			DO j=2,nY-1
				DO i=2,nX-1
					r(i,j,k)=	pRhs(i,j,k)-(Aw*PCorr(i-1,j,k) + Ae*PCorr(i+1,j,k) + &
								As*PCorr(i,jm1(j),k) + An*PCorr(i,jp1(j),k) + &
								Ab(k)*PCorr(i,j,k-1) + At(k)*PCorr(i,j,k+1) + &
								Ap(k)*PCorr(i,j,k))
				END DO
			END DO
		END DO
		!This preconditioner changes with the number of processors
		DO k=bn(id),en(id)
			DO j=2,nY-1
				DO i=2,nX-1
					rhat(i,j,k) = -w*(r(i,j,k)+Aw*rhat(i-1,j,k)+As*rhat(i,j-1,k)+Ab(k)*rhat(i,j,k-1))*rAp(k)
				END DO
			END DO
		END DO
		
		DO k=bn(id),en(id)
			DO j=2,nY-1
				DO i=2,nX-1
					rhat(i,j,k) =((2-w)/w)*Ap(k)*rhat(i,j,k)
				END DO
			END DO
		END DO
		

		DO k= en(id),bn(id),-1
			DO j = nY-1,2,-1
				DO i = nX-1,2,-1
					rhat(i,j,k) =-w*(rhat(i,j,k)+Ae*rhat(i+1,j,k)+An*rhat(i,j+1,k)+At(k)*rhat(i,j,k+1))*rAp(k)
				END DO
			END DO
		END DO

		r(2:nX-1,2:nY-1,bn(id):en(id)) = rhat(2:nX-1,2:nY-1,bn(id):en(id))

		rhoLoc=(sum(r(2:nX-1,2:nY-1,bn(id):en(id))*r(2:nX-1,2:nY-1,bn(id):en(id))))

		CALL mpi_allreduce(rhoLoc,rho,1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)

		rho = sqrt(rho)
		IF(m.eq.1)sumnor=rho
		sumav =rho/sumnor

		g(1) =rho
		v(2:nX-1,2:nY-1,bn(id):en(id),1)=r(2:nX-1,2:nY-1,bn(id):en(id))/rho
		kit=0
		DO WHILE((sumav > errtol).and.(kit < kmax)) ! Begin gmres loop.
			kit=kit+1
			!***********Exchange information in phi b/w processors *****************
			IF (id.eq.0) THEN
				CALL mpi_recv(v(1,1,en(id)+1,kit),1,oneplane,id+1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(v(1,1,en(id)  ,kit),1,oneplane,id+1,50,mpi_comm_world,ierr)
			END IF
			
			IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.1)) THEN
				CALL mpi_send(v(1,1,en(id)  ,kit),1,oneplane,id+1,50,mpi_comm_world,ierr)
				CALL mpi_recv(v(1,1,en(id)+1,kit),1,oneplane,id+1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(v(1,1,bn(id)  ,kit),1,oneplane,id-1,50,mpi_comm_world,ierr)
				CALL mpi_recv(v(1,1,bn(id)-1,kit),1,oneplane,id-1,50,mpi_comm_world,STATUS,ierr)
			END IF

			IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.0)) THEN
				CALL mpi_recv(v(1,1,bn(id)-1,kit),1,oneplane,id-1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(v(1,1,bn(id)  ,kit),1,oneplane,id-1,50,mpi_comm_world,ierr)
				CALL mpi_recv(v(1,1,en(id)+1,kit),1,oneplane,id+1,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_send(v(1,1,en(id)  ,kit),1,oneplane,id+1,50,mpi_comm_world,ierr)
			END IF

			IF (id.eq.nproc-1) THEN
				CALL mpi_send(v(1,1,bn(id)  ,kit),1,oneplane,id-1,50, mpi_comm_world,ierr)
				CALL mpi_recv(v(1,1,bn(id)-1,kit),1,oneplane,id-1,50, mpi_comm_world,STATUS,ierr)
			END IF

			DO k=bn(id),en(id)
				DO j=2,nY-1
					DO i=2,nX-1
						v(i,j,k,kit+1)=	Aw*v(i-1,j,k,kit) + Ae*v(i+1,j,k,kit) + &
										As*v(i,jm1(j),k,kit) + An*v(i,jp1(j),k,kit) + &
										Ab(k)*v(i,j,k-1,kit) + At(k)*v(i,j,k+1,kit) + &
										Ap(k)*v(i,j,k,kit)
					END DO
				END DO
			END DO

			!This preconditioner changes with the number of processors
			DO k=bn(id),en(id)
				DO j=2,nY-1
					DO i=2,nX-1
						rhat(i,j,k) =-w*(v(i,j,k,kit+1) + Aw*rhat(i-1,j,k) + As*rhat(i,j-1,k) + Ab(k)*rhat(i,j,k-1))*rAp(k)
					END DO
				END DO
			END DO

			DO k=bn(id),en(id)
				DO j=2,nY-1
					DO i=2,nX-1
						rhat(i,j,k) =  ((2-w)/w)*Ap(k)*rhat(i,j,k)
					END DO
				END DO
			END DO
			
			DO k= en(id),bn(id),-1
				DO j = nY-1,2,-1
					DO i = nX-1,2,-1
						rhat(i,j,k) =-w*(rhat(i,j,k)+Ae*rhat(i+1,j,k) + An*rhat(i,j+1,k)+ At(k)*rhat(i,j,k+1))*rAp(k)
					END DO
				END DO
			END DO
			
			v(2:nX-1,2:nY-1,bn(id):en(id),kit+1) = rhat(2:nX-1,2:nY-1,bn(id):en(id))

			! Begin modified GS. May need to reorthogonalize.
			DO k=1,kit
				tempgmres=sum(v(2:nX-1,2:nY-1,bn(id):en(id),k)*v(2:nX-1,2:nY-1,bn(id):en(id),kit+1))
				CALL mpi_allreduce(tempgmres,h(k,kit),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
				v(2:nX-1,2:nY-1,bn(id):en(id),kit+1)=v(2:nX-1,2:nY-1,bn(id):en(id),kit+1)-h(k,kit)*v(2:nX-1,2:nY-1,bn(id):en(id),k)
			END DO
			
			tempgmres=(sum(v(2:nX-1,2:nY-1,bn(id):en(id),kit+1)*v(2:nX-1,2:nY-1,bn(id):en(id),kit+1)))
			CALL mpi_allreduce(tempgmres,h(kit+1,kit),1,mpi_double_precision,mpi_sum,mpi_comm_world,ierr)
			h(kit+1,kit) = sqrt(h(kit+1,kit))

			IF (h(kit+1,kit).gt.0.0.or.h(kit+1,kit).lt.0.0) THEN
				v(2:nX-1,2:nY-1,bn(id):en(id),kit+1)=v(2:nX-1,2:nY-1,bn(id):en(id),kit+1)/h(kit+1,kit)
			END IF

			IF (kit>1)  THEN
				!Apply old Givens rotations to h(1:kit,kit).
				DO i=1,kit-1
					hik			=c(i)*h(i,kit)-s(i)*h(i+1,kit)
					hipk		=s(i)*h(i,kit)+c(i)*h(i+1,kit)
					h(i,kit)	=hik
					h(i+1,kit)	=hipk
				END DO
			END IF
			
			nu=sqrt(h(kit,kit)**2 + h(kit+1,kit)**2)
			!May need better Givens implementation.
			!Define and Apply new Givens rotations to h(kit:kit+1,kit).
			IF (nu.gt.0.0) THEN
				c(kit)		= h(kit,kit)/nu
				s(kit)		= -h(kit+1,kit)/nu
				h(kit,kit)	= c(kit)*h(kit,kit)-s(kit)*h(kit+1,kit)
				h(kit+1,kit)= 0
				gk    		= c(kit)*g(kit) -s(kit)*g(kit+1)
				gkp  		= s(kit)*g(kit) +c(kit)*g(kit+1)
				g(kit) 		= gk
				g(kit+1) 	= gkp
			END IF
			
			rho		= abs(g(kit+1))
			sumav	= rho/sumnor
			mag(kit)= rho
			
		END DO !End of gmres loop.
		!h(1:kit,1:kit) is upper triangular matrix in QR.
		y(kit) = g(kit)/h(kit,kit)
		
		DO i = kit-1,1,-1
		y(i) = g(i)
			DO k = i+1,kit
				y(i) = y(i) -h(i,k)*y(k)
			END DO
		y(i) = y(i)/h(i,i)
		END DO
		
		DO i = 1,kit ! Form linear combination.
			PCorr(2:nX-1,2:nY-1,bn(id):en(id)) = PCorr(2:nX-1,2:nY-1,bn(id):en(id)) + v(2:nX-1,2:nY-1,bn(id):en(id),i)*y(i)
		END DO

	END DO ! End restart loop.

	IF (id.eq.0) THEN
		CALL mpi_recv(PCorr(1,1,en(id)+1),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(PCorr(1,1,en(id)-1),1,twoplanes,id+1,50,mpi_comm_world,ierr)
	END IF
	
	IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.1)) THEN
		CALL mpi_send(PCorr(1,1,en(id)-1),1,twoplanes,id+1,50,mpi_comm_world,ierr)
		CALL mpi_recv(PCorr(1,1,en(id)+1),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(PCorr(1,1,bn(id)  ),1,twoplanes,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(PCorr(1,1,bn(id)-2),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
	END IF

	IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.0)) THEN
		CALL mpi_recv(PCorr(1,1,bn(id)-2),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(PCorr(1,1,bn(id)  ),1,twoplanes,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(PCorr(1,1,en(id)+1),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(PCorr(1,1,en(id)-1),1,twoplanes,id+1,50,mpi_comm_world,ierr)
	END IF

	IF (id.eq.nproc-1) THEN
		CALL mpi_send(PCorr(1,1,bn(id)  ),1,twoplanes,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(PCorr(1,1,bn(id)-2),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
	END IF

	IF (id.eq.0) PCorr(1:nX,2:nY-1,bn(id)-1) = ( (4.0d0*PCorr(1:nX,2:nY-1,bn(id))) -PCorr(1:nX,2:nY-1,bn(id)+1) )/3.0d0
	IF(id.eq.nproc-1)  PCorr(1:nX,2:nY-1,en(id)+1)= (5.0*PCorr(1:nX,2:nY-1,en(id))+              &
	&                  (-4.0*PCorr(1:nX,2:nY-1,en(id)-1)+ PCorr(1:nX,2:nY-1,en(id)-2) ))/2.0d0
	
	Pcorr (1:nX,     1,b(id):e(id)) = Pcorr(1:nX,nY-1,b(id):e(id))
	Pcorr (1:nX,nY,b(id):e(id)) = Pcorr(1:nX,       2,b(id):e(id))

	!***********************************************************************
	PCorr(1,2:nY-1,b(id):e(id)) =  (4.0d0*PCorr(2,2:nY-1,b(id):e(id)) -PCorr(3,2:nY-1,b(id):e(id)) )/3.0d0

	PCorr(nX,2:nY-1,b(id):e(id))= (5.0*PCorr(nX-1,2:nY-1,b(id):e(id))+              &
	&                  (-4.0*PCorr(nX-2,2:nY-1,b(id):e(id))+ PCorr(nX-3,2:nY-1,b(id):e(id)) ))/2.0d0

	!*******************Correcting Pressure values**************************
	Press(2:nX-1,2:nY-1,bn(id):en(id))=Press(2:nX-1,2:nY-1,bn(id):en(id))+PCorr(2:nX-1,2:nY-1,bn(id):en(id))

	DO k=bn(id),en(id)
		DO j=2,nY-1
			DO i=2,nX-1
				gradp(1)=0.5d0*(PCorr(i+1,j,k)-PCorr(i-1,j,k))*(invdx)
				gradp(2)=0.5d0*(PCorr(i,jp1(j),k)-PCorr(i,jm1(j),k))*(invdy)
				gradp(3)=0.5d0*(PCorr(i,j,k+1)-PCorr(i,j,k-1))*(invdXi*JacInv(k))
				!-------------------------------------------------------
				! V^n+1 = V^n - dt*grad(P')
				Ut(i,j,k,1)= Ut(i,j,k,1)-(dt*gradp(1))
				Ut(i,j,k,2)= Ut(i,j,k,2)-(dt*gradp(2))
				Ut(i,j,k,3)= Ut(i,j,k,3)-(dt*gradp(3))
			END DO
		END DO
	END DO
	
	!===================Boundary Condition Starts=======================
	IF(id.eq.0)	Ut(1:nX,2:nY-1,bn(id)-1,1:3) = 0.0d0   !Bottom Vel BC: V=0.

	IF (id.eq.nproc-1)	Ut(1:nX,2:nY-1,en(id)+1,1:3)=(4.0*Ut(1:nX,2:nY-1,en(id),1:3)-Ut(1:nX,2:nY-1,en(id)-1,1:3))/3.0d0 !Top BC: 2nd Ord V'=0

	Ut(nX,2:nY-1,b(id):e(id),1:3) =(4.0d0*Ut(nX-1,2:nY-1,b(id):e(id),1:3)-Ut(nX-2,2:nY-1,b(id):e(id),1:3))/3.0d0 !Outflow BC: 2nd Ord V'=0

	Ut(:,1,b(id):e(id),1:3)= Ut(:,nY-1,b(id):e(id),1:3) !Spanwise Vel BC
	Ut(:,nY,b(id):e(id),1:3)= Ut(:,2,b(id):e(id),1:3)          
	!---------------------Bottom Press BC-------------------------------
	IF(id.eq.0)THEN
		DO j=2,nY-1
			DO i=2,nX-1
				rt11=Ut(i,j,1,3)*(Ut(i,j,2,3)-Ut(i,j,1,3))/df(1)

				rt1=(Ut(i-1,j,1,3)-2.0d0*Ut(i,j,1,3)+Ut(i+1,j,1,3))/(dx*dx)
				rt2=(Ut(i,j-1,1,3)-2.0d0*Ut(i,j,1,3)+Ut(i,j+1,1,3))/(dy*dy)
				rt3=(2.d0*df(1)*Ut(i,j,3,3)-2.d0*df(2)*Ut(i,j,2,3)+2.d0*(df(2)-df(1))*Ut(i,j,1,3))/(df(1)*df(2)**2-(df(2)*df(1)**2))

				! Finding The H and Pressure at the Top Wall
				rt12= -rt11 +(rt1+rt2+rt3)/Re
				Press(i,j,1)=Press(i,j,2)-rt12*df(1)
			END DO
		END DO
	END IF
	!--------------------Top Press BC-----------------------------------
	IF(id.eq.nproc-1) THEN
		DO j=2,nY-1
			DO i=2,nX-1
				rt1 =(11.0*Ut(i,j,nZ,3) - 18.0*Ut(i,j,nZ-1,3) + 9.0*Ut(i,j,nZ-2,3) - 2.0*Ut(i,j,nZ-3,3))/(6.0*dXi)
				rt2=(rt1*JacInv(nZ))
				Press(i,j,nZ)=rt2/Re
			END DO
		END DO
	END IF
	!----------------------Outflow Press BC-------------------------------------
	DO k=bn(id),en(id)
		DO j=2,nY-1
			Press(nX,j,k)=((3.0d0*Ut(nX,j,k,1)-4.0d0*Ut(nX-1,j,k,1) + Ut(nX-2,j,k,1))/(2.0d0*dx*Re))
		END DO
	END DO
	!----------------------Inflow Press BC--------------------------------------
	DO k=bn(id),en(id)
		DO j=2,nY-1
			rt1 = Ut(1,j,k,1)*((Ut(2,j,k,1)-Ut(1,j,k,1))/(dx))
			rt2 = Ut(1,j,k,2)*((Ut(1,j+1,k,1)-Ut(1,j-1,k,1))/(2.0d0*dy))
			rt3 = Ut(1,j,k,3)*((Ut(1,j,k+1,1)-Ut(1,j,k-1,1))/(2.0d0*dXi))

			rt4 =(Ut(1,j,k,1)+Ut(3,j,k,1)-2*Ut(2,j,k,1))/(dx**2)
			rt5 =(Ut(1,j+1,k,1)+Ut(1,j-1,k,1)-2*Ut(1,j,k,1))/(dy**2)
			
			IF (k.le.2) THEN
				rt11 =     (Ut(1,j,k+1,1)-Ut(1,j,k,1))/dXi
				rt6 = (Ut(1,j,k+2,1) -2.0d0*Ut(1,j,k+1,1) + Ut(1,j,k,1))/(dXi**2)
			ELSEIF (k.ge.nZ-1) THEN
				rt11 =     (Ut(1,j,k,1)-Ut(1,j,k-1,1))/dXi
				rt6 = (Ut(1,j,k-2,1) -2.0d0*Ut(1,j,k-1,1) + Ut(1,j,k,1))/(dXi**2)
			ELSE
				rt11 =     (Ut(1,j,k+1,1)-Ut(1,j,k-1,1))/(2.0d0*dXi)
				rt6 = (-2*Ut(1,j,k,1) + Ut(1,j,k+1,1) + Ut(1,j,k-1,1))/(dXi**2)
			END IF
			
			rt6  = rt6*JacInv(k)**2 - Jac2(k)*rt11*JacInv(k)**3
			rt12=(Ut(1,j,k,1)-Utold(1,j,k,1))/dt
			Press(1,j,k) =((2*dx*((rt12+rt1+rt2+rt3*JacInv(k))-((rt4+rt5+rt6)/Re)))-Press(3,j,k)+(4*Press(2,j,k))) /3.0d0
		END DO
	END DO
	!---------------------------Spanwise Press BC-------------------------------
	Press(:,     1,b(id):e(id))= Press(:,nY-1,b(id):e(id))      !
	Press(:,nY,b(id):e(id))= Press(:,       2,b(id):e(id))      
	!======================Boundry Condition Ends=======================

	!Exchange  Updated velocity vector field between processors for 2 planes.
	IF (id.eq.0) THEN
		CALL mpi_recv(Ut(1,1,en(id)+1,1),1,threevar,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Ut(1,1,en(id)-1,1),1,threevar,id+1,50,mpi_comm_world,ierr)
	END IF
	
	IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.1)) THEN
		CALL mpi_send(Ut(1,1,en(id)-1,1),1,threevar,id+1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Ut(1,1,en(id)+1,1),1,threevar,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Ut(1,1,bn(id)  ,1),1,threevar,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Ut(1,1,bn(id)-2,1),1,threevar,id-1,50,mpi_comm_world,STATUS,ierr)
	END IF

	IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.0)) THEN
		CALL mpi_recv(Ut(1,1,bn(id)-2,1),1,threevar,id-1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Ut(1,1,bn(id)  ,1),1,threevar,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Ut(1,1,en(id)+1,1),1,threevar,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Ut(1,1,en(id)-1,1),1,threevar,id+1,50,mpi_comm_world,ierr)
	END IF

	IF (id.eq.nproc-1) THEN
		CALL mpi_send(Ut(1,1,bn(id)  ,1),1,threevar,id-1,50, mpi_comm_world,ierr)
		CALL mpi_recv(Ut(1,1,bn(id)-2,1),1,threevar,id-1,50, mpi_comm_world,STATUS,ierr)
	END IF

	!Exchange updated pressure field between processors for plane after interface.
	IF (id.eq.0) THEN
		CALL mpi_recv(Press(1,1,en(id)+1),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Press(1,1,en(id)-1),1,twoplanes,id+1,50,mpi_comm_world,ierr)
	END IF
	
	IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.1)) THEN
		CALL mpi_send(Press(1,1,en(id)-1),1,twoplanes,id+1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Press(1,1,en(id)+1),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Press(1,1,bn(id)  ),1,twoplanes,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Press(1,1,bn(id)-2),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
	END IF

	IF ((id.gt.0).and.(id.lt.nproc-1).and.(mod(id,2).eq.0)) THEN
		CALL mpi_recv(Press(1,1,bn(id)-2),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Press(1,1,bn(id)  ),1,twoplanes,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Press(1,1,en(id)+1),1,twoplanes,id+1,50,mpi_comm_world,STATUS,ierr)
		CALL mpi_send(Press(1,1,en(id)-1),1,twoplanes,id+1,50,mpi_comm_world,ierr)
	END IF

	IF (id.eq.nproc-1) THEN
		CALL mpi_send(Press(1,1,bn(id)  ),1,twoplanes,id-1,50,mpi_comm_world,ierr)
		CALL mpi_recv(Press(1,1,bn(id)-2),1,twoplanes,id-1,50,mpi_comm_world,STATUS,ierr)
	END IF

	!==========================Calculate Divergence=====================
	IF(mod(nstep,freqDiv)==0) THEN
		tmp =0.0d0
		!res=0.0d0
		DO k=bn(id),en(id)
			DO j=2,nY-1
				DO i=2,nX-1
					IF((i.eq.2) .or. (i.eq.(nX-1)))THEN
						rt1= (Ut(i+1,j,k,1)-Ut(i-1,j,k,1))*(0.5d0*invdx)
					ELSE
						rt1= ((Ut(i-2,j,k,1)-8*Ut(i-1,j,k,1)+8*Ut(i+1,j,k,1)-Ut(i+2,j,k,1))*(onetwelfth*invdx))
					END IF
					
					IF((j.eq.2) .or. (j.eq.(nY-1)))THEN
						rt2= (Ut(i,jp1(j),k,2)-Ut(i,jm1(j),k,2))*(0.5d0*invdy)
					ELSE
						rt2= ((Ut(i,jm2(j),k,2)-8*Ut(i,jm1(j),k,2)+8*Ut(i,jp1(j),k,2)-Ut(i,jp2(j),k,2))*(onetwelfth*invdy))
					END IF
					
					IF(k.le.2 .or. k.ge.(nZ-1)) THEN
						rt3=((Ut(i,j,k+1,3)-Ut(i,j,k-1,3))*(0.5d0*invdXi))
					ELSE
						rt3= ((Ut(i,j,k-2,3)-8*Ut(i,j,k-1,3)+8*Ut(i,j,k+1,3)-Ut(i,j,k+2,3))*(onetwelfth*invdXi))
					END IF
					
					tmp =tmp + (rt1+rt2+rt3*JacInv(k))
				END DO
			END DO
		END DO
		
		CALL MPI_Reduce(tmp, res,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
		IF (id.eq.0) WRITE(*,105) time,nstep,abs(res)/((nY-2)*(nX-2)*(nZ-2))
		105 FORMAT(1X,F14.7,1X,I9,1X,F14.7)
		!---------------------------------------------------------------
		!======================Calculate Massflow=======================
		DO i=1,5
			massLoc =0.0d0
			
			DO k=bn(id),en(id)
				DO j=2,nY-1
					dAv = dy*(x3(k+1)-x3(k))
					massLoc =massLoc+(0.25d0*(Ut(plane(i),j,k,1)+Ut(plane(i),jp1(j),k,1)&
										+Ut(plane(i),j,k+1,1)+Ut(plane(i),jp1(j),k+1,1))*(dAv))
				END DO
			END DO
			CALL MPI_Reduce(massLoc,mass(i),1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)
			mass(i)=mass(i)/((x3(nZ-1)-x3(2))*(x2(nY-1)-x2(2)))
		END DO

		IF (id.eq.0) THEN
			OPEN (UNIT=12,FILE='Output/massflow.dat',POSITION='append')
			WRITE(12,106)time,mass(1),mass(2),mass(3),mass(4),mass(5)
			CLOSE(12)
			106 FORMAT(6(1X,F15.7))
		END IF
	END IF
	!---------------------------------------------------------------------
	!====================Write Intermediate Solution File=================
	IF(mod(nstep,freqWrite)==0) THEN
	
		IF(mod(nstep/freqWrite,2)==0) THEN
		WRITE(filname,'(a,i3.3,a)')'Output/3D',id+1,'.dat'
		ELSE
		WRITE(filname,'(a,i3.3,a)')'Output/backup/3D',id+1,'.dat'
		ENDIF
		
		OPEN (UNIT=12,FILE=filname)
		IF(id.eq.0) WRITE (12,107)'TITLE ="',time,nstep,nproc,'"'
		IF(id.eq.0) WRITE (12,'(a)')'Variables = "x","y","z","U","V","W","P"'
		IF(id.eq.0) WRITE (12,108)'ZONE k=',nX,',j=',nY-2,',i=',nZ,',DATAPACKING="POINT"'
		DO i=1,nX
			WRITE(12,*)
			DO j=2,nY-1
				WRITE(12,*)
				DO k=b(id),e(id)
					WRITE(12,104)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Press(i,j,k)
				END DO
			END DO
		END DO
		CLOSE(12)
		107 FORMAT(A,F20.10,I9,I3,A)
		108 FORMAT(A,I4,A,I4,A,I4,A)
		!104 FORMAT(3(1X,F10.6),4(1X,E22.15))
	END IF
	!---------------------------------------------------------------------
	!=================Calculating velocity correlations===================
	IF(mod(nstep,freqStat)==0) THEN
		DO i=1,nX
			DO k=b(id),e(id)
				tmpusum=0.0d0
				tmpvsum=0.0d0
				tmpwsum=0.0d0

				tmpuusum=0.0d0
				tmpvvsum=0.0d0
				tmpwwsum=0.0d0

				tmpuvsum=0.0d0
				tmpuwsum=0.0d0
				tmpvwsum=0.0d0

				DO j=2,nY-1 !Avg along Spanwise Direction
					su=0.50d0*(Ut(i,j,k,1)+ Ut(i,jp1(j),k,1))
					sv=0.50d0*(Ut(i,j,k,2)+ Ut(i,jp1(j),k,2))
					sw=0.50d0*(Ut(i,j,k,3)+ Ut(i,jp1(j),k,3))

					!First order stats integrating over a plane
					tmpusum=tmpusum+su
					tmpvsum=tmpvsum+sv
					tmpwsum=tmpwsum+sw

					!Second order stats integrating over a plane
					tmpuusum=tmpuusum+su*su
					tmpvvsum=tmpvvsum+sv*sv
					tmpwwsum=tmpwwsum+sw*sw

					tmpuvsum=tmpuvsum+su*sv
					tmpuwsum=tmpuwsum+su*sw
					tmpvwsum=tmpvwsum+sv*sw
				END DO
				A=nY-2
				usum(i,k) = usum(i,k) +(tmpusum/A)
				vsum(i,k) = vsum(i,k) +(tmpvsum/A)
				wsum(i,k) = wsum(i,k) +(tmpwsum/A)

				uusum(i,k)= uusum(i,k)  +(tmpuusum/A)
				vvsum(i,k)= vvsum(i,k)  +(tmpvvsum/A)
				wwsum(i,k)= wwsum(i,k)  +(tmpwwsum/A)

				uvsum(i,k) =uvsum(i,k) +(tmpuvsum/A)
				uwsum(i,k) =uwsum(i,k) +(tmpuwsum/A)
				vwsum(i,k) =vwsum(i,k) +(tmpvwsum/A)
			END DO
		END DO
	END IF
	!====================Send Mean Statistics to proc Zero================
	IF(mod(nstep,freqTimAvg)==0)THEN
		IF (id.eq.0) THEN
			DO i = 1,nproc-1				
				CALL mpi_recv(usum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_recv(vsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_recv(wsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)

				CALL mpi_recv(uusum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_recv(vvsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_recv(wwsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)

				CALL mpi_recv(uvsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_recv(uwsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
				CALL mpi_recv(vwsum(1,bn(i)),nX*siz(i),mpi_double_precision,i,50,mpi_comm_world,STATUS,ierr)
			END DO
		ELSE
			CALL mpi_send(usum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
			CALL mpi_send(vsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
			CALL mpi_send(wsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)

			CALL mpi_send(uusum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
			CALL mpi_send(vvsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
			CALL mpi_send(wwsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)

			CALL mpi_send(uvsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
			CALL mpi_send(uwsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
			CALL mpi_send(vwsum(1,bn(id)),nX*siz(id),mpi_double_precision,0,50,mpi_comm_world,ierr)
		END IF
		!-------------------------------------------------------------------

		IF (id.eq.0) THEN
			
			DO i=1,nX
				DO k=1,nZ
					usum2(i,k)=(usum(i,k)/freqTimAvg)
					vsum2(i,k)=(vsum(i,k)/freqTimAvg)
					wsum2(i,k)=(wsum(i,k)/freqTimAvg)

					uusum2(i,k)=(uusum(i,k)/freqTimAvg)
					vvsum2(i,k)=(vvsum(i,k)/freqTimAvg)
					wwsum2(i,k)=(wwsum(i,k)/freqTimAvg)

					uvsum2(i,k)=(uvsum(i,k)/freqTimAvg)
					uwsum2(i,k)=(uwsum(i,k)/freqTimAvg)
					vwsum2(i,k)=(vwsum(i,k)/freqTimAvg)
				END DO
			END DO
				
			DO k=bn(id),en(id)
				DO j=2,nY-1
					DO i=1,nX-1
						Ufluc(i,j,k,1) = Ut(i,j,k,1)-usum2(i,k)
						Ufluc(i,j,k,2) = Ut(i,j,k,2)-vsum2(i,k)
						Ufluc(i,j,k,3) = Ut(i,j,k,3)-wsum2(i,k)
					ENDDO
				END DO
			END DO

			DO i=1,nX
				DO k=1,nZ
					Ruu(i,k) = uusum2(i,k)-(usum2(i,k)*usum2(i,k))
					Rvv(i,k) = vvsum2(i,k)-(vsum2(i,k)*vsum2(i,k))
					Rww(i,k) = wwsum2(i,k)-(wsum2(i,k)*wsum2(i,k))
					Ruv(i,k) = uvsum2(i,k)-(usum2(i,k)*vsum2(i,k))
					Ruw(i,k) = uwsum2(i,k)-(usum2(i,k)*wsum2(i,k))
					Rvw(i,k) = vwsum2(i,k)-(vsum2(i,k)*wsum2(i,k))
				END DO
			END DO

			utau1=sqrt(((usum2(1,2)-usum2(1,1))/(x3(2)-x3(1)))/Re)
			utau2=sqrt(((usum2(68,2)-usum2(68,1))/(x3(2)-x3(1)))/Re)

			lamda=(utau1)/(utau2)

			WRITE(*,*) lamda,'lamda'

			DO i=1,nX
				DO k=1,nZ
					IF(usum2(i,k).ge.(0.99d0*Ue)) THEN
						bltkns(i) = x3(k)
						GOTO 500
					END IF
				END DO
			500 CONTINUE
			END DO

			lamda03=bltkns(68)/bltkns(1)

			OPEN (UNIT=19,FILE='Output/lamda.dat',POSITION='append')
			!WRITE (19,*)'Variables = "time","lamda","lamda03"'
			WRITE(19,110) time,lamda,lamda03
			CLOSE(19)
			110 FORMAT(3(2x,F20.10))

			DO i=1,nX
				DO k=1,nZ
					Wt(i,k)=0.5d0*(1+(dtanh((4.0d0*((x3(k)/bltkns(i))-0.2d0))/((1-2.0d0*0.2d0)*(x3(k)/bltkns(i))+0.2d0))/dtanh(4.0d0)))
				END DO
			END DO

			OPEN (UNIT=19,FILE='Output/W.dat')
			WRITE (19,*)'Variables = "x1","x3","Wt"'
			DO i=1,nX
				DO k=1,nZ
					WRITE(19,111) x1(i), x3(k), Wt(i,k)
				END DO
			END DO
			CLOSE(19)
			111 FORMAT(3(1X,F15.7))

			new(:)=nZ-1

			DO k=1,nZ-1
				IF(Wt(68,k).le.0.99d0) THEN
					ksep=k
					eta(k)= lamda
				ELSE
					eta(k)= lamda03
				END IF
			END DO

			DO k=1,nZ-1
				zrecy(k)=eta(k)*x3(k)
			END DO

			counter=1
			DO k=1,nZ-1
				IF(counter.le.(ksep)) THEN
					IF(x3(k).ge.zrecy(counter))THEN
						new(counter)=k
						counter= counter+1
					END IF
				END IF
			END DO

			counter=ksep+1
			DO k=ksep+1,nZ-1
				IF(counter.le.(nZ-1)) THEN
					IF(x3(k).ge.zrecy(counter))THEN
						new(counter)=k
						counter= counter+1
					END IF
				END IF
			END DO

			DO i=1,nX
				DO k=1,nZ-1
					deltastar(i)=deltastar(i)+((1.0d0-(usum2(i,k)/Ue))*(x3(k+1)-x3(k)))
					IF(x3(k).eq.bltkns(i)) GOTO 144
				END DO
			144 CONTINUE
			END DO

			DO i=1,nX
				DO k=1,nZ-1
					theta(i)=theta(i)+((usum2(i,k)/Ue)*(1.0d0-(usum2(i,k)/Ue))*(x3(k+1)-x3(k)))
					IF(x3(k).ge.bltkns(i)) GOTO 155
				END DO
			155 CONTINUE
			END DO

			DO i= 1,nX
				Ht(i)=deltastar(i)/theta(i)
			END DO

			DO i= 1,nX
				Reth(i)=theta(i)*Re*Ue
			END DO

			OPEN (UNIT=19,FILE='Output/boundLyrThik.dat')
			WRITE (19,*)'Variables = "x1","BLthick","deltastar","theta","Ht","Reth","Rex"'
			DO i=1,nX
				WRITE(19,112) x1(i), bltkns(i), deltastar(i), theta(i), Ht(i), Reth(i),Rex(i)
			END DO
			CLOSE(19)
			112 FORMAT(8(2x,F20.10))

			theta(:)=0.0d0
			deltastar(:)=0.0d0

			OPEN (UNIT=19,FILE='Output/newK.dat')
			WRITE (19,*)'Variables = "k","newk","etak"'
				DO k=1,nZ
					WRITE(19,113) k, new(k),eta(k)
				END DO
			113 FORMAT(2(I5),F20.10)
			
			OPEN (UNIT=19,FILE='Output/RMS.dat')
			WRITE (19,*)'Variables = "x1","x3","Ruu","Rvv","Rww","Ruv","Ruw","Rvw"'
			DO i=1,nX
				DO k=1,nZ
					WRITE(19,116) x1(i), x3(k), Ruu(i,k), Rvv(i,k), Rww(i,k), Ruv(i,k), Ruw(i,k), Rvw(i,k)
				END DO
			END DO
			CLOSE(19)
			116 FORMAT(8(1X,F15.7))

			OPEN (UNIT=19,FILE='Output/meanVel.dat')
			WRITE (19,*)'Variables = "x1","x3","usum2","vsum2","wsum2"'
			DO i=1,nX
				DO k=1,nZ
					WRITE(19,114) x1(i), x3(k), usum2(i,k), vsum2(i,k), wsum2(i,k)
				END DO
			END DO
			CLOSE(19)
			114 FORMAT(5(1X,F15.7))
			
			OPEN (UNIT=20,FILE='Input/statData.dat')
			DO k=1,nZ
				WRITE(20,115) usum2(68,k),vsum2(68,k),wsum2(68,k),Wt(1,k),new(k)
			END DO
			CLOSE(20)
			!115 FORMAT(4(1X,F15.7),I5.3)
		END IF

		CALL  MPI_Bcast(lamda,1,mpi_double_precision,0,mpi_comm_world,ierr)
		CALL  MPI_Bcast(new,nZ,mpi_double_precision,0,mpi_comm_world,ierr)
		CALL  MPI_Bcast(bltkns,nX,mpi_double_precision,0,mpi_comm_world,ierr)
		CALL  MPI_Bcast(usum2,(nX-1)*nZ,mpi_double_precision,0,mpi_comm_world,ierr)
		CALL  MPI_Bcast(vsum2,(nX-1)*nZ,mpi_double_precision,0,mpi_comm_world,ierr)
		CALL  MPI_Bcast(wsum2,(nX-1)*nZ,mpi_double_precision,0,mpi_comm_world,ierr)
		CALL  MPI_Bcast(Wt,(nX-1)*nZ,mpi_double_precision,0,mpi_comm_world,ierr)

		usum(:,:) =0.0d0
		vsum(:,:) =0.0d0
		wsum(:,:) =0.0d0

		uusum(:,:) =0.0d0
		vvsum(:,:) =0.0d0
		wwsum(:,:) =0.0d0

		uvsum(:,:) =0.0d0
		uwsum(:,:) =0.0d0
		vwsum(:,:) =0.0d0
	END IF
	
		!==========================Calculate Kinetic Energy===================
	IF(cnt.gt.freqTimAvg)THEN
		kinEnLoc =0.0d0
		kinEn =0.0d0
		
		DO k=bn(id),en(id)
			DO j=2,nY-1
				DO i=2,nX-1
					Ufluc(i,j,k,1) = Ut(i,j,k,1)-usum2(i,k)
					Ufluc(i,j,k,2) = Ut(i,j,k,2)-vsum2(i,k)
					Ufluc(i,j,k,3) = Ut(i,j,k,3)-wsum2(i,k)
					kinEnLoc = kinEnLoc + Ufluc(i,j,k,1)**2 + Ufluc(i,j,k,2)**2 + Ufluc(i,j,k,3)**2
				END DO
			END DO
		END DO

		CALL MPI_Reduce(kinEnLoc,kinEn,1,mpi_double_precision,mpi_sum,0,mpi_comm_world,ierr)

		IF (id.eq.0) THEN
			kinEn = kinEn/dble(nX*(nY-2)*(nZ-2))
			OPEN (UNIT=19,FILE='Output/kinEn.dat',POSITION='append')
			!WRITE (19,*)'Variables = "time","nstep","kinetic energy"'
			WRITE(19,109) time,kinEn
			CLOSE(19)
			109 FORMAT(2(2x,E20.10))
		END IF
	
		!---------------------------------------------------------------------
		!=====================Check for large Fluctuations====================
		DO i=1,nX
			DO j=2,nY-1
				DO k=bn(id),en(id)
					IF (maxval(Ufluc(i,j,k,:)).gt.25.0d0) THEN
						write(*,*)"Very Large Fluctuations"
						WRITE(*,*) nstep,time,i,j,k,Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3)
						STOP
					ENDIF
				END DO
			END DO
		END DO
	END IF
	!---------------------------------------------------------------------
	IF(id==0 .and. mod(nstep,freqLog)==0) THEN
		OPEN(UNIT=25,FILE='Output/simLog.dat',POSITION='append')
		WRITE(25,'(i10.6,a,a)')nstep,' iterations completed on ',dateTime()
		CLOSE(25)
	ENDIF
END DO ! Time loop

CALL mpi_type_free(oneplane,ierr)
CALL mpi_type_free(twoplanes,ierr)
CALL mpi_finalize(ierr)

contains
function dateTime()

implicit none
character(len=30)::dateTime
character(len=3):: ampm
integer:: d,h,m,n,s,y,mm,values(8)
character(len=3), parameter, dimension(12) :: &
month=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']

call date_and_time(values=values)

y=values(1)
m=values(2)
d=values(3)
h=values(5)
n=values(6)
s=values(7)
mm=values(8)

if(h<12) then
	ampm='AM'
elseif(h==12) then
	if(n==0 .and. s==0) then
		ampm='Noon'
	else
		ampm='PM'
	endif
else
	h=h-12
	if(h<12) then
		ampm='PM'
	elseif(h==12) then
		if(n==0 .and. s==0) then
			ampm='Midnight'
		else
			ampm='AM'
		endif
	endif
endif

write(dateTime,'(i2,1x,a,1x,i4,2x,i2,a1,i2.2,a1,i2.2,a1,i3.3,1x,a)')&
d,trim(month(m)),y,h,':',n,':',s,'.',mm,trim(ampm)
end function dateTime

END PROGRAM revamptbl
