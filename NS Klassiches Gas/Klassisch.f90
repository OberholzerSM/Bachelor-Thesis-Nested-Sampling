program NS

	implicit none
	
	integer :: i, j, modus
	real (kind = 16) :: Z, oZ, T, v_max, v, E
	
	call random_seed()
	
	v_max = real(3500.0,16)
	
	modus = 2
	
	if(modus == 0) then
	
		open(6,file = 'Z_Klassisch.txt', status = 'old')
		open(7,file = 'oZ_Klassisch.txt', status = 'old')
		open(8,file = 'i_Klassisch.txt', status = 'old')
		open(9,file = 'v_Klassisch.txt',status = 'old')
		open(10,file = 'E_Klassisch.txt',status = 'old')
	
		do j = 1,40
		
			T = real( 10 * j, 16)
			call sampler(500, 1, T, v_max, Z, oZ, i, E, v)
			
			write(6,*) Z
			write(7,*) oZ
			write(8,*) i
			write(9,*) v
			write(10,*) E
			
			print *, floor(T) / 4, '%'
			
		end do
		
		close(6)
		close(7)
		close(8)
		close(9)
		close(10)
	
	else if(modus == 1) then
	
		open(6,file = 'Z_100K_Klassisch.txt', status = 'old')
		open(7,file = 'oZ_100K_Klassisch.txt', status = 'old')
		open(8,file = 'i_100K_Klassisch.txt', status = 'old')
		
		open(9,file = 'Z_200K_Klassisch.txt', status = 'old')
		open(10,file = 'oZ_200K_Klassisch.txt', status = 'old')
		open(11,file = 'i_200K_Klassisch.txt', status = 'old')
		
		open(12,file = 'Z_300K_Klassisch.txt', status = 'old')
		open(13,file = 'oZ_300K_Klassisch.txt', status = 'old')
		open(14,file = 'i_300K_Klassisch.txt', status = 'old')
		
		open(15,file = 'Z_400K_Klassisch.txt', status = 'old')
		open(16,file = 'oZ_400K_Klassisch.txt', status = 'old')
		open(17,file = 'i_400K_Klassisch.txt', status = 'old')
	
		do j = 500,10000,500
		
			T = real( 100, 16)
			call sampler(j, 1, T, v_max, Z, oZ, i, E, v)
			
			write(6,*) Z
			write(7,*) oZ
			write(8,*) i
			
			T = real( 200, 16)
			call sampler(j, 1, T, v_max, Z, oZ, i, E, v)
			
			write(9,*) Z
			write(10,*) oZ
			write(11,*) i
			
			T = real( 300, 16)
			call sampler(j, 1, T, v_max, Z, oZ, i, E, v)
			
			write(12,*) Z
			write(13,*) oZ
			write(14,*) i
			
			T = real( 400, 16)
			call sampler(j, 1, T, v_max, Z, oZ, i, E, v)
			
			write(15,*) Z
			write(16,*) oZ
			write(17,*) i
			
			print *, j
			
		end do
		
		close(6)
		close(7)
		close(8)
		close(9)
		close(10)
		close(11)
		close(12)
		close(13)
		close(14)
		close(15)
		close(16)
		close(17)
		
	else if(modus == 2) then
	
		open(6,file = 'Z_vmax_100K_Klassisch.txt', status = 'new')
		open(7,file = 'oZ_vmax_100K_Klassisch.txt', status = 'new')
		open(8,file = 'i_vmax_100K_Klassisch.txt', status = 'new')
		
		open(9,file = 'Z_vmax_200K_Klassisch.txt', status = 'new')
		open(10,file = 'oZ_vmax_200K_Klassisch.txt', status = 'new')
		open(11,file = 'i_vmax_200K_Klassisch.txt', status = 'new')
		
		open(12,file = 'Z_vmax_300K_Klassisch.txt', status = 'new')
		open(13,file = 'oZ_vmax_300K_Klassisch.txt', status = 'new')
		open(14,file = 'i_vmax_300K_Klassisch.txt', status = 'new')
		
		open(15,file = 'Z_vmax_400K_Klassisch.txt', status = 'new')
		open(16,file = 'oZ_vmax_400K_Klassisch.txt', status = 'new')
		open(17,file = 'i_vmax_400K_Klassisch.txt', status = 'new')
	
		do j = 1,10
		
			v_max = real( 500*j , 16 )
		
			T = real( 100, 16)
			call sampler(500, 1, T, v_max, Z, oZ, i, E, v)
			
			write(6,*) Z
			write(7,*) oZ
			write(8,*) i
			
			T = real( 200, 16)
			call sampler(500, 1, T, v_max, Z, oZ, i, E, v)
			
			write(9,*) Z
			write(10,*) oZ
			write(11,*) i
			
			T = real( 300, 16)
			call sampler(500, 1, T, v_max, Z, oZ, i, E, v)
			
			write(12,*) Z
			write(13,*) oZ
			write(14,*) i
			
			T = real( 400, 16)
			call sampler(500, 1, T, v_max, Z, oZ, i, E, v)
			
			write(15,*) Z
			write(16,*) oZ
			write(17,*) i
			
			print *, j
			
		end do
	
	end if
	
	
end program NS

subroutine prior(sample, d, v_max) !Creates a d-dimensional Sample from the Prior

	integer :: d
	real (kind = 16) :: u, v_max
	real (kind = 16), dimension(d) :: sample
	
	call random_number(u)
	sample(1) = v_max * u

end subroutine prior

real (kind = 16) function lnL(theta, d, T, v_max) !ln(likelihood)

	integer :: d
	real (kind = 16) :: T, v_max, pi, m, kB
	real (kind = 16), dimension(d) :: theta
	
	pi = real( 4.0 * atan(1.0), 16)
	kb = real( 1.380649 * 10.0 ** (-23) ,16)
	m = real( 4.002602 * 1.66053906661 * 10.0**(-27) ,16)
	
	lnL = real( log( real(v_max) * 4.0 * real(pi) * ( real(theta(1))**2 ) ) ,16) - ( (m * theta(1)**2) / ( real(2.0,16) * kB * T ) )

end function lnL

!Inputs: K: Number of Live-Points. d: Dimension of one Live-Point.
!Outputs: Z: Evidence. i: Number of iterations.
subroutine sampler(K, d, T, v_max, Z, oZ_H, i, E, v)

	integer :: K, d 											!Dimensions
	integer :: i, i_new 										!Counters
	integer :: j, l 											!Do-Loops
	integer :: j_min, j_test 									!Positions
	integer :: halt, test1, test2								!Test
	 
	real (kind = 16):: Z, w, H, oZ_H, E, v						!To calculate
	real (kind = 16):: lnL										!Functions
	real (kind = 16):: lnL_min, lnL_old, test_lnL				!Ln[L_min]
	real (kind = 16):: X_i, X_i2								!Prior-Volume
	real (kind = 16):: dZ, dlnZ, lnL_max						!Stopping-Criterion
	real (kind = 16):: sigma, u, T, v_max, m					!Other real variables
	real (kind = 16), dimension(d) :: sample, sigma_max			!Arrays
	real (kind = 16), dimension(K,d) :: livepoints				!Matrix
	
	open(1, file = 'Results_Klassisch.txt', status = 'old')
	open(2, file = 'Weights_Klassisch.txt', status = 'old')
	open(3, file = 'Deadpoints_Klassisch.txt', status = 'old')
	open(4, file = 'lnL_Klassisch.txt', status = 'old')
	open(5, file = 'PriorVolume_Klassisch.txt', status = 'old')
	
	m =  real( 4.002602 * 1.66053906661 * 10.0**(-27) ,16)
	
	!Generate K Live-Points from the Prior
	do j = 1,K
	
		call prior(sample,d, v_max)
		
		do l = 1,d 
		
			livepoints(j,l) = sample(l)
			
		end do
	
	end do
	
	i = 0
	Z = real(0.0,16)
	H = real(0.0,16)
	oZ_H = real(0.0,16)
	
	v = real(0.0,16)
	E = real(0.0,16)
	
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
		
			test_lnL = lnL( sample, d, T, v_max )
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
		
		!Determine the average Speed and Energy
		
		v = v + ( w * livepoints(j_min,1) ) 
		
		E = E + ( w * real(0.5,16) * m * ( livepoints(j_min,1)**2 ) )
		
		!Replace the Live-Point with the smallest ln[L]
		
		!If directly sampling from the prior does not take too long
		if( test1 == 0 ) then
		
			i_new = 0
		
			test2 = 0
			do while( test2 == 0 )
		
				i_new = i_new + 1
			
				call prior(sample,d, v_max)
			
				if( lnL(sample, d, T, v_max) > lnL_min ) then
				
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
					
					sigma_max(l) = v_max * X_i
				
				end do
			
				! Add j_max many random variations
				do j = 1,100
			
					do l = 1,d
				
						call random_number(u)
						sigma = real(2.0,16) * sigma_max(l) * u - sigma_max(l)
						
						sample(l) = sample(l) + sigma
				
					end do
					
				end do

				if( lnL(sample, d, T, v_max) >= lnL_min ) then
				
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
		
			test_lnL = lnL( sample, d, T, v_max )
			
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
	
	!Calculate E and v
	
	E = E / Z
	
	v = v / Z
	
	!Print the Results
	
	write(1,*) Z
	write(1,*) oZ_H
	write(1,*) i
	write(1,*) E
	write(1,*) v
	
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	
end subroutine sampler