program fptest
	implicit none
	real*8 :: x, xn
	real*8, external :: fixedPtIteration
	real*8 :: tol = 1E-4, h = 1.0
	x = 1.8
	xn = 0.0
	x = fixedPtIteration(x,tol,xn,h)


end program fptest


subroutine f(eval,x)
	implicit none
	real*8, intent(in) :: x
	real*8, intent(out) :: eval
	eval = (x + 10)**0.25
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
