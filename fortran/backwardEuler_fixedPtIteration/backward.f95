program backwardEul
	implicit none
	real*8, parameter :: h = 0.1, tfinal = 0.1
	real*8, allocatable, dimension(:) :: y,t
	integer*8 :: n,i
	real*8 :: tempSol, tol = 1E-8
	real*8, external :: fixedPtIteration
	
	n = int(tFinal/h) + 1
	print *,'N= ',n
	allocate(y(n),t(n))
	print *,'Initializing time'
	do i=1,n
		t(i) = (i-1)*h
	end do

	print *,'Setting initial condition'
	y(1) = 0.2
	print *,'Start integration...'

	do i=2,n
		tempSol = fixedPtIteration(y(i-1), tol, y(i-1), h)
		print *,'Updating solution ',i
		y(i) = tempSol
	end do
	
	print *,'Writing results'
	open(unit=1, file="data.dat", status="replace")
	
	10 format(f8.5,1X,f8.5)
	do i=1,n
			write(unit = 1, fmt = 10), t(i), y(i)
	end do
	close(unit=1)
	
	
	deallocate(y,t)
end program backwardEul

subroutine f(eval,x)
	implicit none
	real*8, intent(in) :: x
	real*8, intent(out) :: eval
	eval = 2.0*x*(1.0 - x)
end subroutine f

real*8 function fixedPtIteration(guess, tol, xn, h)
	implicit none
	real*8, intent(in) :: tol, xn, h
	real*8, intent(inout) :: guess
	real*8 :: feval
	real*8 :: finalSol = 0.0 !Works only when guess neq 0.0
	integer*8 :: i=0
	print *,'Starting FP iteration...'
	do 
		if (abs(finalSol - guess) .le. tol) then
			print *,'Convergence!'
			exit
		else
			if (i .ne. 0) then
				guess = finalSol
			end if
			call f(feval,guess)
			finalSol = xn + h*feval
			print *,'Final Sol: ', finalSol
			i = i+1
			print *,'FixedPt iteration: ', i
		end if
		if (i .ge. 100) then
			print *,'Max iterations exceeded!'
			exit
		end if
	end do
	fixedPtIteration = finalSol
end function fixedPtIteration
