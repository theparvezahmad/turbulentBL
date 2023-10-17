!Last Changed 06/05/16 Parvez Ahmad
program restart
!This Program changes the parallel write to a restart file

implicit none

DOUBLE PRECISION,dimension(:,:,:,:),allocatable::Ut,Ufluc
DOUBLE PRECISION,dimension(:,:,:),allocatable::Press
DOUBLE PRECISION,dimension(:,:),allocatable::usum2,vsum2,wsum2,Wt
DOUBLE PRECISION,dimension(:),allocatable::x1,x2,x3	!Physical Geometery (Cartesian)
DOUBLE PRECISION::time
INTEGER, DIMENSION(:), ALLOCATABLE :: new
INTEGER::my_rank,proc,loc_n,nstep
INTEGER::i,j,k,istart,iend,bn,en,nX,nY,nZ,cnt
CHARACTER(LEN=25)::filname

open(unit=12,file='cnt',status='unknown')
read(12,*)nZ,cnt,proc
close(12)

do my_rank=0,proc-1
	write(filname,'(a,i3.3,a)')'3D',my_rank+1,'.tp'
	open(unit=12,file=filname,status='unknown')
	if(my_rank.eq.0)then
		read(12,01)time,nstep
		read(12,*)
		read(12,02)nX,nY,nZ
		nY=nY+2

		allocate(Ut(nX,nY,nZ,3))
		allocate(Ufluc(nX,nY,nZ,3))
		allocate(Press(nX,nY,nZ))
		allocate(usum2(nX,nZ))
		allocate(vsum2(nX,nZ))
		allocate(wsum2(nX,nZ))
		allocate(Wt(nX,nZ))
		allocate(x1(nX))
		allocate(x2(nY))
		allocate(x3(nZ))
		allocate(new(nZ))
	endif

	loc_n=(nZ-1)/proc
	bn=2+(my_rank*loc_n)
	en=bn+loc_n-1
	if(my_rank.eq.proc-1)en=nZ-1

	istart=bn
	iend=en
	if(my_rank.eq.0)istart=bn-1
	if(my_rank.eq.proc-1)iend=en+1

	print*,my_rank,istart,iend

	do i=1,nX
	read(12,*)
		do j=2,nY-1
		read(12,*)
			do k=istart,iend
				read(12,97)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Press(i,j,k),usum2(i,k),vsum2(i,k),wsum2(i,k),Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Wt(i,k),new(k)
			enddo
		enddo
	enddo
	close(12)
enddo

Ut(:,1,:,:)=Ut(:,nY-1,:,:)
Ut(:,nY,:,:)=Ut(:,2,:,:)

x2(1)=0.0d0
x2(nY)=4.0d0

open(unit=12,file='3D.tp',status='unknown')
write(12,03)'TITLE ="',time,nstep,'"'
write(12,*)'Variables = "x","y","z","U","V","W","P","Um","Vm","Wm","Uf","Vf","Wf","Wt","Nk"'
write(12,04)'ZONE k=',nX,',j=',nY,',i=',nZ,',DATAPACKING="POINT"'
do i=1,nX
	write(12,*)
	do j=1,nY
		write(12,*)
		do k=1,nZ
			WRITE(12,97)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Press(i,j,k),usum2(i,k),vsum2(i,k),wsum2(i,k),Ufluc(i,j,k,1),Ufluc(i,j,k,2),Ufluc(i,j,k,3),Wt(i,k),new(k)
		enddo
	enddo
enddo
close(12)
!*****************************************************************************
01 format(8X,F20.10,I9,1X)
02 format(7X,I4,3X,I4,3X,I6)
03 format(A,F20.10,I9,A)
04 format(A,I4,A,I4,A,I4,A)
97 FORMAT(3(1X,F10.6),11(1X,E22.15),1X,I4)

deallocate(Ut)
deallocate(x1)
deallocate(x2)
deallocate(x3)

end program restart
