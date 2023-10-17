program combineFile
!Author : Parvez Ahmad <pahmed333@gmail.com>
!It combines separate files into one
!Last modified on 02/10/2016

implicit none
DOUBLE PRECISION,dimension(:,:,:,:),allocatable::Ut
DOUBLE PRECISION,dimension(:),allocatable::x1,x2,x3
INTEGER,dimension(:),allocatable::b,e
DOUBLE PRECISION::time
INTEGER::id,nproc,nstep,it1,it2
INTEGER::i,j,k,nX,nY,nZ
CHARACTER(LEN=25)::filname

open(unit=12,file='Output/3D001.dat')
	read(12,101)time,nstep,nproc
	read(12,*)
	read(12,102)nX,nY,nZ
	nY=nY+2

	allocate(Ut(nX,nY,nZ,4))
	allocate(x1(nX),x2(nY),x3(nZ))
	ALLOCATE(b(0:nproc-1),e(0:nproc-1))
	
it1=nZ/nproc
it2=nproc-mod(nZ,nproc)

b(0)=1
do i=0,nproc-1
	if(i==it2) it1=it1+1
	e(i)=b(i)+it1-1
	if(i==nproc-1) exit
	b(i+1)=e(i)+1
enddo

101 FORMAT(8X,F20.10,I9,I3,1X)
102 format(7X,I4,3X,I4,3X,I4)

do id=0,nproc-1
	write(filname,'(a,i3.3,a)')'Output/3D',id+1,'.dat'
	open(unit=12,file=filname)
	print*,id,b(id),e(id)
	do i=1,nX
	read(12,*)
		do j=2,nY-1
		read(12,*)
			do k=b(id),e(id)
				read(12,104)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Ut(i,j,k,4)
			enddo
		enddo
	enddo
	close(12)
enddo

Ut(:,1,:,:)=Ut(:,nY-1,:,:)
Ut(:,nY,:,:)=Ut(:,2,:,:)

x2(1)=0.0d0
x2(nY)=4.0d0

open(unit=12,file='Input/3D.dat')
write(12,105)'TITLE ="',time,nstep,'"'
write(12,'(a)')'Variables = "x","y","z","U","V","W","P"'
write(12,106)'ZONE k=',nX,',j=',nY,',i=',nZ,',DATAPACKING=POINT'
do i=1,nX
	write(12,*)
	do j=1,nY
		write(12,*)
		do k=1,nZ
			write(12,104)x1(i),x2(j),x3(k),Ut(i,j,k,1),Ut(i,j,k,2),Ut(i,j,k,3),Ut(i,j,k,4)
		enddo
	enddo
enddo
close(12)
105 format(A,F20.10,I9,A)
106 format(A,I4,A,I4,A,I4,A)
104 FORMAT(3(1X,F10.6),5(1X,E22.15))

end program combineFile
