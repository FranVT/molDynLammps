
program readdumpspheres

implicit none

integer :: ij,ii,i,j,ntot,numbin,kmax
integer :: qmax,iq,jq,kq,sbin,ntotav,ball,b1,b2,b3,b4
real(8) :: ln1,ln2,lxn,lyn,lzn,xc,yc,zc,s,x2,y2,z2,timestp,sigma
real(8) :: dq0,bin0,sfsum1,sfsum2,qr(3),&
       rq,rq0,ki,smax,qmax0,times

character(60) :: dummy
character (60) :: filename
character (60) :: filedump, filedump0



real(8), parameter :: pi=acos(-1.0d00), pi2= 2.0d00* pi


real(8), dimension(:,:), allocatable :: r,sfhis
real(8), dimension(:), allocatable :: bc,kr


open(44,  file = "datareadsf",action = "read", status = "old" )
read(44,*) filedump0, filedump
read(44,*) bin0,qmax0
read(44,*) sigma
read(44,*) ball
!this is to turn on and off certain particle types
read(44,*) b1, b2, b3, b4
close(44)

filename=trim(filedump0)//trim("/")//trim(filedump)

open(45, file = filename,action = "read", status = "old" )

read(45,*) dummy
read(45,*) timestp
read(45,*) dummy
read(45,*) ntot
read(45,*) dummy
read(45,*) ln1,ln2
lxn=ln2-ln1
read(45,*) ln1,ln2
lyn=ln2-ln1
read(45,*) ln1,ln2
lzn=ln2-ln1
read(45,*) dummy
!Calculatin parameters
xc=lxn
yc=lyn
zc=lzn
s=sigma
x2=xc/2.0
y2=yc/2.0
z2=zc/2.0
!\Delta q values and maximum q vector.
dq0=pi2/xc
bin0=dq0*bin0
qmax=int(qmax0/dq0)
rq0=qmax*dq0
!Number of bins
numbin=int(qmax*dq0/bin0)+1
allocate(r(ntot,4),sfhis(numbin,2),bc(ntot),kr(ntot))	

r=0.0

!this is to multiply by 0 particles positions that are not needed for the calculation
	
if (ball==1) then
	bc=1
	ntotav=ntot
else
	bc=0
	ntotav=0
end if

do ij=1,ntot
	read(45,*) ii,i,j,r(ii,4),r(ii,1),r(ii,2),r(ii,3)
	if (ball==1) then
		bc(ii)=1
	else if (ball== 0 .and. r(ii,4)==1) then
		bc(ii)=b1
		ntotav=ntotav+b1
	else if (ball== 0 .and. r(ii,4)==2) then
		bc(ii)=b2		
		ntotav=ntotav+b2
	else if (ball== 0 .and. r(ii,4)==3) then
		bc(ii)=b3
		ntotav=ntotav+b3
	else if (ball== 0 .and. r(ii,4)==4) then
		bc(ii)=b4
		ntotav=ntotav+b4
	end if
end do
close(45)

do i=1, ntot
	if (r(i,1)==0 .and. r(i,2)==0 .and.r(i,3)==0 ) print*, "No data"
end do


! sf calculation
sfhis=0.0
do iq=-qmax,qmax
	qr(1)=iq*dq0
	do jq=-qmax,qmax
		qr(2)=jq*dq0
		do kq=0,qmax			
			if (iq==0 .and. jq==0 .and. kq==0) cycle
			times=2.0
			if (kq==0) times=1.0
			qr(3)=kq*dq0
                        rq=sqrt(dot_product(qr,qr))
                        if (rq> rq0 ) cycle
			sbin=int(rq/bin0)+1.0
			if (mod(rq,bin0)==0) sbin=sbin-1
			if (sbin> numbin) cycle
			kr=matmul(r(:,1:3),qr)
			sfsum1=sum(bc*dcos(kr/s))
			sfsum2=sum(bc*dsin(kr/s))
			sfhis(sbin,1)=sfhis(sbin,1)+(sfsum1*sfsum1+sfsum2*sfsum2)*times/(ntotav*1.0)
			sfhis(sbin,2)=sfhis(sbin,2)+times			
		end do
	end do
end do


sfhis(:,1)=sfhis(:,1)/sfhis(:,2)
kmax=int(maxloc(sfhis(:,1),dim=1))
smax=sfhis(kmax,1)

filename=trim(filedump0)//trim("/sf.dat")

open(22, file=filename, action = "write", status="unknown")
!q vector - S(q) - S(q)/S_max
!	write(22, *) "#q       S(q)            S(q)/Smax"
	do ij = 1,numbin-1
		ki=(ij*1.0)*bin0
    	write(22,133) ki," ",sfhis(ij,1)," ",sfhis(ij,1)/smax
	end do
close(22)

133 format(F6.3,A,F12.6,A,F12.6)


end program
