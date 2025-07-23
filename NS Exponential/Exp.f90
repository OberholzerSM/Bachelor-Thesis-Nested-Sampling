program NS_Exp

	implicit none
	
	integer :: i, j, k, modus
	real (kind = 16) :: Z, oZ_H, oZ_X, oZ_B, oZ_100, Z_sum, Z2_sum
	
	call random_seed()
	
	modus = 2
	
	if(modus == 0) then
	
		call sampler(500, 1, Z, oZ_H, i)
		
		print *, Z
	
	else if (modus == 1) then
	
		open(8,file = 'Z_Exp.txt',status = 'old')
		open(9,file = 'i_Exp.txt',status = 'old')
		open(10,file = 'oZ_H_Exp.txt',status = 'old')
		open(11,file = 'oZ_X_Exp.txt',status = 'old')
		open(12,file = 'oZ_B_Exp.txt',status = 'old')
		open(13,file = 'oZ_100_Exp.txt',status = 'old')
	
		do k = 500,2000,10
	
			call sampler(k, 1, Z, i, oZ_H, oZ_X, oZ_B)
		
			write(8,*) Z
			write(9,*) i
			write(10,*) oZ_H
			write(11,*) oZ_X
			write(12,*) oZ_B
		
			Z_sum = real(0.0,16)
			Z2_sum = real(0.0,16)
		
			do j = 1,100
		
				call sampler(k, 1, Z, oZ_H, oZ_X, oZ_B, i)
			
				Z_sum = Z_sum + Z
				Z2_sum = Z2_sum + Z ** 2
		
			end do
		
			oZ_100 = real( sqrt( ( real(Z2_sum) - ( real(Z_sum **2) / 100.0 ) ) / 100.0 ) ,16)
			
			write(13,*) oZ_100
			
			print *, K
		
		end do
	
		close(8)
		close(9)
		close(10)
		close(11)
		close(12)
		close(13)

	else if (modus == 2) then
	
		open(8, file = 'Z_500_Exp.txt', status = 'old')
		open(9, file = 'oZ_500_Exp.txt', status = 'old')
		open(10, file = 'i_500_Exp.txt', status = 'old')
		
		open(11, file = 'Z_2000_Exp.txt', status = 'old')
		open(12, file = 'oZ_2000_Exp.txt', status = 'old')
		open(13, file = 'i_2000_Exp.txt', status = 'old')
	
		do j = 1,100
		
			call sampler(500, 1, Z, oZ_H, i)
			
			write(8,*) Z
			write(9,*) oZ_H
			write(10,*) i
			
			call sampler(2000, 1, Z, oZ_H, i)
			
			write(11,*) Z
			write(12,*) oZ_H
			write(13,*) i
			
			print*, j
		
		end do
		
	else if (modus == 3) then
	
		open(8,file = 'Z_Exp.txt',status = 'old')
		open(9,file = 'oZ_Exp.txt',status = 'new')
	
		do k = 500,2000,100
		
			Z_sum = real(0.0,16)
			Z2_sum = real(0.0,16)
		
			do j = 1,100
		
				call sampler(k, 1, Z, oZ_H, i)
			
				Z_sum = Z_sum + Z
				Z2_sum = Z2_sum + Z ** 2
		
			end do
			
			oZ_100 = real( sqrt( ( real(Z2_sum) - ( real(Z_sum **2) / 100.0 ) ) / 100.0 ) ,16)
			
			write(8,*) Z_sum / real(100.0,16)
			write(9,*) oZ_100
			
			print *, K
		
		end do
	
		close(8)
		close(9)

	end if

end program NS_Exp

subroutine prior(sample, d) !Creates a d-dimensional Sample from the Prior

	integer :: d
	real (kind = 16) :: u
	real (kind = 16), dimension(d) :: sample
	
	call random_number(u)
	sample(1) = u

end subroutine prior

real (kind = 16) function lnL(theta, d) !ln(likelihood)

	integer :: d
	real (kind = 16), dimension(d) :: theta
	lnL = theta(1)

end function lnL

subroutine boots(i, oZ_b)
	integer :: i, j, k, r
	real (kind = 16):: oZ_b, Z, Z_sum, Z2_sum, u
	real (kind = 16), dimension(i) :: w
	
	open(6,file = 'Weights_Exp.txt', status = 'old')
	
	do j = 1,i
	
		read(6,*) w(j)
	
	end do
	
	close(6)
	
	Z_sum = real(0.0,16)
	Z2_sum = real(0.0,16)
	
	do k = 1,100
	
		Z = real(0.0,16)
		do j = 1,i
		
			call random_number(u)
			
			r = 1 + floor( real(i)*u )
			
			Z = Z + w(r)
		
		end do
	
		Z_sum = Z_sum + Z
		Z2_sum = Z2_sum + ( Z ** 2 )
	
	end do
	
	oZ_b = real( sqrt( ( real(Z2_sum) - ( real(Z_sum **2) / 100.0 ) ) / 100.0 ) ,16)

end subroutine boots

subroutine X_Error(i, K, oZ_X)

	integer :: i, K, j, l
	real :: u
	real (kind = 16) :: oZ_X, t, w, Z, Z_sum, Z2_sum
	real (kind = 16), dimension(i) :: lnL, X
	
	open (7,file = 'lnL_Exp.txt', status = 'old')
	
	do j = 1,i
	
		read(7,*) lnL(j)
	
	end do
	
	close(7)
	
	Z_sum = real(0.0,16)
	Z2_sum = real(0.0,16)
	
	do l = 1,100
	
		Z = real(0.0,16)
		
		do j = 1,i
		
			call random_number(u)
		
			t = real( u ** ( 1 / real(i) ), 16)
			
			if(j > 2) then
				
				X(j) = t * X(j-1)
				w = real( ( exp(real( lnL(j-1) )) + exp(real( lnL(j) )) ) / 2.0 , 16 ) * ( X(j-2) - X(j) ) / real(2.0,16)
		
			else if(j == 2) then
		
				X(j) = t * X(j-1)
				w = real( ( exp(real( lnL(j-1) )) + exp(real( lnL(j) )) ) / 2.0 , 16 ) * ( real(1.0,16) - X(j) ) / real(2.0,16)
		
			else if(j == 1) then
		
				X(j) = t
				w = real( ( 0.0 + exp(real( lnL(j) )) ) / 2.0 , 16 ) * ( real(1.0,16) - X(j) ) / real(2.0,16)
		
			end if
			
			Z = Z + w
		
		end do
		
		Z_sum = Z_sum + Z
		Z2_sum = Z2_sum + ( Z ** 2 )
	
	end do
	
	oZ_X = real( sqrt( ( real(Z2_sum) - ( real(Z_sum **2) / 100.0 ) ) / 100.0 ) ,16)

end subroutine X_Error

!Inputs: K: Number of Live-Points. d: Dimension of one Live-Point. 
			!theta_max: Maximal values for theta theta_min: Minimal values for theta
!Outputs: Z: Evidence. i: Number of iterations.
subroutine sampler(K, d, Z, oZ_H, i)

	integer :: K, d 															!Dimensions
	integer :: i, i_new 														!Counters
	integer :: j, l 															!Do-Loops
	integer :: j_min, j_test 													!Positions
	integer :: halt, test1, test2												!Test
	
	real (kind = 16):: Z, w, H, oZ_H, oZ_X, oZ_B								!To calculate
	real (kind = 16):: lnL														!Functions
	real (kind = 16):: lnL_min, lnL_old, test_lnL								!Ln[L_min]
	real (kind = 16):: X_i, X_i2												!Prior-Volume
	real (kind = 16):: dZ, dlnZ, lnL_max										!Stopping-Criterion
	real (kind = 16):: sigma, u													!Other real variables
	real (kind = 16), dimension(d) :: sample, sigma_max, theta_max, theta_min	!Arrays
	real (kind = 16), dimension(K,d) :: livepoints								!Matrix
	
	open(1, file = 'Results_Exp.txt', status = 'old')
	open(2, file = 'Weights_Exp.txt', status = 'old')
	open(3, file = 'Deadpoints_Exp.txt', status = 'old')
	open(4, file = 'lnL_Exp.txt', status = 'old')
	open(5, file = 'PriorVolume_Exp.txt', status = 'old')
	
	
	!Generate K Live-Points from the Prior
	do j = 1,K
	
		call prior(sample,d)
		
		do l = 1,d 
		
			livepoints(j,l) = sample(l)
			
		end do
	
	end do
	
	i = 0
	Z = real(0.0,16)
	H = real(0.0,16)
	oZ_H = real(0.0,16)
	
	lnL_min = real(0.0,16)
	dlnZ = real(10000.0,16)
	test1 = 0
	
	! Begin main Loop
	halt = 0
	do while( halt == 0 )
	
		!i: Current iteration
		i = i + 1
	
		!Save the lnL_min from the previous Iteration
		lnL_old = lnL_min
	
		!Find the Live-Point with the smallest ln[L]
		lnL_min = real(100.0,16)
		do j = 1,K
		
			do l = 1,d 
			
				sample(l) = livepoints(j,l)
			
			end do
		
			test_lnL = lnL( sample, d )
			if( lnL_min > test_lnL ) then
			
				lnL_min = test_lnL
				j_min = j
			
			end if
		
		end do
		
		do l = 1,d
		
			write(3,*) livepoints(j_min,l)
		
		end do
		
		write(4,*) lnL_min
		
		!Estimate the Prior-Volume
		X_i = real( exp( - real(i) / real(K) ) , 16)
		X_i2 = real( exp( - real(i-2) / real(K) ) , 16)
		
		write(5,*) X_i
		
		!Determine the Weight
		
		if(i > 2) then
		
			w = real( ( exp(real(lnL_old)) + exp(real(lnL_min)) ) / 2.0 , 16 ) * ( X_i2 - X_i ) / real(2.0,16)
		
		else if(i == 2) then
		
			w = real( ( exp(real(lnL_old)) + exp(real(lnL_min)) ) / 2.0 , 16 ) * ( real(1.0,16) - X_i ) / real(2.0,16)
		
		else if(i == 1) then
		
			w = real( ( 0.0 + exp(real(lnL_min)) ) / 2.0 , 16 ) * ( real(1.0,16) - X_i ) / real(2.0,16)
		
		end if
		
		write(2,*) w
		
		!Determine Z, its error and H
		
		Z = Z + w
		
		H = H + ( w * lnL_min )
		
		!Replace the Live-Point with the smallest ln[L]
		
		!If directly sampling from the prior does not take too long
		if( test1 == 0 ) then
		
			i_new = 0
		
			test2 = 0
			do while( test2 == 0 )
		
				i_new = i_new + 1
			
				call prior(sample,d)
			
				if( lnL(sample, d) > lnL_min ) then
				
					do l = 1,d
					
						livepoints( j_min, l ) = sample(l)
						
					end do
				
					test2 = 1
				
				end if
				
				!If the loop takes too long, exist and use the other method
				if( i_new > 100000 ) then
				
					test1 = 1
					exit
				
				end if
		
			end do
		
		!If directly sampling from the prior does take too long
		else if (test1 == 1) then

			test2 = 0
			do while( test2 == 0 )
			
				!Choose a random Live-Point and determine sigma_max
				call random_number(u)
				j_test = 1 + floor( real(K) * real(u) )
				
				do l = 1,d
				
					sample(l) = livepoints(j_test,l)
					
					sigma_max(l) = real(0.5,16) * ( X_i ** ( real(1.0,16) / d ) )
				
				end do
			
				! Add j_max many random variations
				do j = 1,100
			
					do l = 1,d
				
						call random_number(u)
						sigma = real(2.0,16) * sigma_max(l) * u - sigma_max(l)
						
						sample(l) = sample(l) + sigma
				
					end do
					
				end do

				if( lnL(sample, d) >= lnL_min ) then
				
					do l = 1,d
					
						livepoints( j_min, l ) = sample(l)
						
					end do
				
					test2 = 1
				
				end if
			
			end do
		
		end if
		
		!Find the largest lnL
		lnL_max = 0
		do j = 1,K
		
			do l = 1,d 
			
				sample(l) = livepoints(j,l)
			
			end do
		
			test_lnL = lnL( sample, d )
			
			if( lnL_max < test_lnL ) then
			
				lnL_max = test_lnL
			
			end if
		
		end do
		
		!Evaluate the Stopping-Criterion
		
		dZ = exp(lnL_max) * X_i
		
		!Z may for the first iterations be numerically too small
		if ( Z > real(0.0,16) ) then
		
			dlnZ = log(Z + dZ) - log(Z)
			
			if( dlnZ < real(10.0 ** (-3),16) ) then
			
				halt = 1
			
			end if
		
		end if
	
	end do
	
	!Calculate the H-Error
	
	H = (H / Z) - real( log( real(Z) ) ,16)
	oZ_H = Z * real( sqrt( real(H) / real(K) ) ,16)
	
	!Print the Results
	
	write(1,*) Z
	write(1,*) oZ_H
	write(1,*) i
	
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	
end subroutine sampler