program main
implicit none
! Laplaces equation for heat conduction with non-homogenous boundary conditions
! Explicit Scheme
! Change
  integer, parameter :: PREC = kind( 1.0d0 ) ! double precision
  integer, parameter :: xnodes=1600 	 ! x grid points
  integer, parameter :: ynodes=xnodes 	 ! y grid points

  real(PREC) :: pi = 4 * atan( 1._PREC ) ! pi
  real(PREC) :: a=1.0                    ! x domain length
  real(PREC) :: b=1.0			 ! y domain length
  real(PREC) :: x			 ! x-location
  real(PREC) :: y			 ! y-location
  real(PREC) :: dx			 ! x grid spacing
  real(PREC) :: dy   			 ! y grid spacing
  real(PREC) :: R			 ! Cell Reynolds Number
  real(PREC) :: CFL = 0.3_PREC		 ! CLF Number  
  real :: start, finish
  real(PREC) :: Tn(xnodes,ynodes) = 0._PREC 
  real(PREC) :: Tnp1(xnodes,ynodes) = 0._PREC
  real(PREC) :: eps = epsilon(1._PREC)
  real(PREC) :: norm_np1=0.0
  real(PREC) :: norm_n=1.0
  real(PREC) :: resd
  
  integer :: un=11 
  integer :: ierror
  integer :: i 				 ! x spatial counter
  integer :: j				 ! y spatial counter
  integer ::u
  integer :: max_iter=5001		 ! maximum iteration count
  integer :: count=0
  open(unit=un,file="Tn_output_OpenMP.txt",status="replace",iostat=ierror)
  if (ierror==1) then
    write(*,*) "Error opening Tn output file."
  else
    write(un,100) xnodes,',',ynodes
  end if
100 format(I10,1A,I10)
  dx = a / (xnodes-1)
  dy = b / (ynodes-1) 
  
  call cpu_time(start)

! Main loop
resd=1
!do while (abs(norm_np1-norm_n)>1e-3.and.count<max_iter)
!do while (count<max_iter)
!$OMP PARALLEL
!$OMP DO SCHEDULE(DYNAMIC)
do u=1,max_iter
 Tnp1=Tn
! norm_np1=norm2(Tnp1)
 count = count+1
! write(*,*) count
 do i=1,xnodes
     do j=1,ynodes
      if(i==1) then
        Tn(i,j)=100.0
      else if (i==xnodes) then
        Tn(xnodes,j)=((Tn(xnodes-1,j)+Tn(xnodes-1,j))*dy*dy+(Tn(xnodes,j+1)+Tn(xnodes,j-1))*dx*dx)/(2*(dx*dx+dy*dy))
      else if (j==1) then
        Tn(i,1)=100.0*sin(pi*i*dx)
      else if  (j==ynodes) then
        Tn(:,ynodes)=0.0
      else 
        Tn(i,j)=((Tn(i+1,j)+Tn(i-1,j))*dy*dy+(Tn(i,j+1)+Tn(i,j-1))*dx*dx)/(2*(dx*dx+dy*dy))
      end if
    end do
 end do
 resd=0
 do i=2,xnodes
    do j=2,ynodes
      resd = resd + sqrt(abs(Tnp1(i,j)**2-Tn(i,j)**2))
    end do 
 end do
 resd=resd/sum(Tnp1)
! norm_n=norm2(Tn)
end do
!$OMP END DO
!$OMP END PARALLEL
 write(*,*) maxval(Tnp1(2:xnodes-1,2:ynodes-1))
! write(*,*) resd
  call cpu_time(finish)
  print '("Time = ",f9.3," seconds.")',finish-start
 do i=1,xnodes  
   do j=1,ynodes
      write(un,101) Tn(i,j)
   end do
 end do
101 format(f10.5,',')

close(un)
end program main
