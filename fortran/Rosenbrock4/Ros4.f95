program Ros4_DO
	implicit none
	integer, parameter :: eqNo = 2 !Number of equations in the system
	real*8, parameter :: dt = 1.0E-6, tfinal = 4.0, e = 0.1, w = 10.0 
	real*8, allocatable, dimension(:) :: x1,x2,t
	!---------Integration parameters------
	real*8, parameter :: c1 = 0.199293275701, c2 = 0.482645235674, c3 = 0.680614886256e-1, c4 = 0.25
	real*8, parameter :: chat1 = 0.346325833758, chat2 = 0.385693175712, chat3 = 0.367980990530
	real*8, parameter :: y = 0.395
	real*8, parameter :: a21 = 0.438, a31 = 0.796920457938, a32 = 0.730795420615e-1
	real*8, parameter :: y21 - -0.767672395484, y31 = -0.851675323742, y32 = 0.5229673, y41 = 0.28846311, y42 = 0.880214e-1, y43 = -0.33739  
	real*8, dimension(4,2) :: k
	real*8, external :: f1,f2
	integer*8 :: n,i
	n = int(tfinal/dt)
	
	
	allocate(x1(n),x2(n),t(n))
	
	
	print *,'Initializing time'
	do i=1,n
		t(i) = (i - 1)*dt
	end do
	
	print *,'Setting initial conditions'
	x1(1) = 0.02
	x2(1) = 0.0
	print *,'Starting integration'
	do i=2,n
		if (mod(i,10) == 0) then
			print *,"Time step", t(i)
		end if
		!Solving for the slopes
		k(1,1) = f1(x1(i-1),x2(i-1),e,w)
		k(1,2) = f1(x1(i-1),x2(i-1),e,w)
		k(2,1) = f1(x1(i-1) + k(1,1)*dt*0.5,x2(i-1) + k(1,2)*dt*0.5,e,w)
		k(2,2) = f2(x1(i-1) + k(1,1)*dt*0.5,x2(i-1) + k(1,2)*dt*0.5,e,w)
		k(3,1) = f1(x1(i-1) + k(2,1)*dt*0.5,x2(i-1) + k(2,2)*dt*0.5,e,w)
		k(3,2) = f2(x1(i-1) + k(2,1)*dt*0.5,x2(i-1) + k(2,2)*dt*0.5,e,w)
		k(4,1) = f1(x1(i-1) + k(3,1)*dt,x2(i-1) + k(3,2)*dt,e,w)
		k(4,2) = f2(x1(i-1) + k(3,1)*dt,x2(i-1) + k(3,2)*dt,e,w)
		
		x1(i) = x1(i-1) + (1.0/6.0)*dt*(k(1,1) + 2*k(2,1) + 2*k(3,1) + k(4,1))
		x2(i) = x2(i-1) + (1.0/6.0)*dt*(k(1,2) + 2*k(2,2) + 2*k(3,2) + k(4,2))
		
	end do
	
	print *,'Writing results'
	open(unit=1, file="data.dat", status="replace")
	
	10 format(f8.5,1X,f8.5,1X,f8.5)
	do i=1,n
		if (mod(i,10) == 0) then
			write(unit = 1, fmt = 10), t(i), x1(i), x2(i)
		end if
	end do
	
	close(unit=1)
	
	
	deallocate(x1,x2,t)
	
	
end program euler_DO

real*8 function f1(x1,x2,e,w)
	
	implicit none
	real*8, intent(in) :: x1,x2,e,w 
	f1 = x2

end function f1

real*8 function f2(x1,x2,e,w)
	implicit none
	real*8, intent(in) :: x1,x2,e,w
	
	f2 = -1.0*(2.0*e*w*x2 + x1*w**2)
	

end function f2
