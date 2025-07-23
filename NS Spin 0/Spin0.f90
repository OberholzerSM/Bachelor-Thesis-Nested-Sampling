program NS

	implicit none
	
	integer :: i, j, l, modus, d, N, K
	integer (kind = 8) :: k_max
	real (kind = 16) :: Z, oZ_H, T, p, E, pi, hbar, m, V, lnL, lnL_min, a, b, fact, v_max, kB, k2
	integer (kind = 8), dimension(30) :: sample
	
	call random_seed()
	
	modus = 1
	
	K = 500
	
	pi = real(4.0*atan(1.0),16)
	hbar = real( 6.62607015 * ( 10.0 ** -34 ) / (2.0 * real(pi) ) ,16)
	kb = real( 1.380649 * 10.0 ** (-23) ,16)
	m = real( 4.002602 * 1.66053906661 * 10.0**(-27) ,16)
	
	v_max = real(3500.0 ,16)
	
	if(modus == 0) then
		
		N = 10
		d = 3*N
		
		open(6, file = 'Z_N10_Spin0.txt',status = 'old')
		open(7, file = 'oZ_N10_Spin0.txt',status = 'old')
		open(8, file = 'E_N10_Spin0.txt',status = 'old')
		
		do j = 1,4
		
			T = real(100*j,16)
		
			do l = 1,4
			
				V = real(l,16) / real(4.0,16)
			
				k_max = floor( (m * ( V ** (1.0/3.0) ) * v_max ) / ( pi * hbar ) , 8 )
			
				call sampler(K, d, T, V, k_max, Z, oZ_H, i, E)
				
				write(6,*) Z
				write(7,*) oZ_H
				write(8,*) E
				
			end do
				
		end do
		
		close(6)
		close(7)
	
	else if(modus == 1) then
	
		N = 1
		d = 3*N
		
		T = 400
		V = 1.0
		
		k_max = floor( (m * ( V ** (1.0/3.0) ) * v_max ) / ( pi * hbar ) , 8 )
		
		call sampler(K, d, T, V, k_max, Z, oZ_H, i, E)
		
	else if(modus == 2) then
	
	
	end if
	
	
end program NS

subroutine prior(sample, d, k_max) !Creates a d-dimensional Sample from the Prior

	integer :: d, l
	integer (kind = 8) :: k_max
	real (kind = 16) :: u
	integer (kind = 8), dimension(d) :: sample
	
	do l = 1,d
	
		call random_number(u)
		sample(l) =  floor(1.0 + real(k_max) * u, 8 )
		
	end do

end subroutine prior

real (kind = 16) function lnL(theta, d, T, V, k_max) !ln(likelihood)

	integer :: d, j, l
	integer (kind = 8) :: k_max
	real :: fact
	real (kind = 16) :: T, V, pi, m, kB, hbar, k2, a
	integer (kind = 8), dimension(d) :: theta
	
	pi = real( 4.0 * atan(1.0), 16)
	hbar = real( 6.62607015 * ( 10.0 ** -34 ) / (2.0 * real(pi) ) ,16)
	kb = real( 1.380649 * 10.0 ** (-23) ,16)
	m = real( 4.002602 * 1.66053906661 * 10.0**(-27) ,16)
	
	a = ( (pi**2) * (hbar**2) ) / ( real(2.0,16) * m * kB * T * ( V ** (1.0/3.0) ) )
	
	k2 = real(0.0,16)
	fact = 1.0
	
	do j = 1,d/3
	
		fact = fact * real(j)
	
	end do
	
	do j = 1,d
		
		k2 = k2 + ( real(theta(j),16) ** 2 )
	
	end do
	
	lnL = real(d,16) * real( log( real(k_max) ) , 16) - real(log(fact),16) -  a * k2

end function lnL

subroutine lnL_new(d, k_max, a, b, sample)

	integer :: d, j, test1
	integer (kind = 8) :: k_max, k_test, k_max_j
	real (kind = 16) :: u, a, b, k2
	integer (kind = 8), dimension(d) :: sample
	
	k2 = real(0.0,16)
	j = 0
	
	k_max_j = k_max
	
	test1 = 0
	do while( test1 == 0)
	
		j = j + 1
	
		call random_number(u)
		k_test =  floor(1.0 + real(k_max_j) * u, 8 )
		
		sample(j) = k_test
		
		k2 = k2 + ( real(k_test,16) ** 2)
		
		if( k2 > b/a ) then
		
			k2 = k2 - ( real(k_test,16) ** 2)
			
			j = j - 1
			
			k_max_j = floor( ( (b/a) - k2 ) ** 0.5, 8)
		
		end if
		
		if(j == d) then
		
			test1 = 1
		
		end if
	
	end do

end subroutine lnL_new

!Inputs: K: Number of Live-Points. d: Dimension of one Live-Point.
!Outputs: Z: Evidence. i: Number of iterations.
subroutine sampler(K, d, T, V, k_max, Z, oZ_H, i, E)

	integer :: K, d													!Dimensions
	integer :: i, i_new 											!Counters
	integer :: j, l 												!Do-Loops
	integer :: j_min, j_test 										!Positions
	integer :: halt, test1, test2									!Test
	integer (kind = 8) :: k_max, sigma
	 
	real (kind = 16):: Z, w, H, oZ_H, E		 						!To calculate
	real (kind = 16):: lnL											!Functions
	real (kind = 16):: lnL_min, lnL_old, test_lnL					!Ln[L_min]
	real (kind = 16):: X_i, X_i1, X_i2								!Prior-Volume
	real (kind = 16):: dZ, dlnZ, lnL_max							!Stopping-Criterion
	real (kind = 16):: u, T, m, V, pi, hbar, kb, a, b, fact, k2		!Other real variables
	integer (kind = 8), dimension(d) :: sample, sigma_max			!Arrays
	integer (kind = 8), dimension(K,d) :: livepoints				!Matrix
	
	open(1, file = 'Results_Spin0_400K_1V.txt', status = 'new')
	open(2, file = 'Weights_Spin0_400K_1V.txt', status = 'new')
	open(3, file = 'Deadpoints_Spin0_400K_1V.txt', status = 'new')
	open(4, file = 'lnL_Spin0_400K_1V.txt', status = 'new')
	open(5, file = 'PriorVolume_Spin0_400K_1V.txt', status = 'new')
	
	pi = real( 4.0 * atan(1.0), 16)
	hbar = real( 6.62607015 * ( 10.0 ** -34 ) / (2.0 * real(pi) ) ,16)
	kb = real( 1.380649 * 10.0 ** (-23) ,16)
	m = real( 4.002602 * 1.66053906661 * 10.0**(-27) ,16)
	
	a = ( (pi**2) * (hbar**2) ) / ( real(2.0,16) * m * kB * T * ( V ** (2.0/3.0) ) )
	
	fact = 1.0
	
	do j = 1,d/3
	
		fact = fact * real(j)
	
	end do
	
	!Generate K Live-Points from the Prior
	do j = 1,K
	
		call prior(sample,d, k_max)
		
		do l = 1,d 
		
			livepoints(j,l) = sample(l)
			
		end do
	
	end do
	
	i = 0
	Z = real(0.0,16)
	H = real(0.0,16)
	oZ_H = real(0.0,16)

	E = real(0.0,16)
	
	lnL_min = real(0.0,16)
	dlnZ = real(10000.0,16)
	test1 = 1
	
	! Begin main Loop
	halt = 0
	do while( halt == 0 )
	
		!i: Current iteration
		i = i + 1
		
		print *, i, real(V), real(T), d/3
	
		!Save the lnL_min from the previous Iteration
		lnL_old = lnL_min
	
		!Find the Live-Point with the smallest ln[L]
		lnL_min = real(100000.0,16)
		do j = 1,K
		
			do l = 1,d 
			
				sample(l) = livepoints(j,l)
			
			end do
		
			test_lnL = lnL( sample, d, T, V, k_max )
			
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
		
		X_i2 = exp( - real(i-2,16) / real(K,16) ) 
		X_i = exp( - real(i,16) / real(K,16) )
		
		write(5,*) X_i
		
		!Determine the Weight
		
		if(i > 2) then
		
			w = real( ( exp(lnL_old) + exp(lnL_min) ) / 2.0 , 16 ) * ( X_i2 - X_i ) / real(2.0,16)
		
		else if(i == 2) then
		
			w = real( ( exp(lnL_old) + exp(lnL_min) ) / 2.0 , 16 ) * ( real(1.0,16) - X_i ) / real(2.0,16)
		
		else if(i == 1) then
		
			w = real( ( 0.0 + exp(lnL_min) ) / 2.0 , 16 ) * ( real(1.0,16) - X_i ) / real(2.0,16)
		
		end if
		
		write(2,*) w
		
		!Determine Z, its error and H
		
		Z = Z + w
		
		H = H + ( w * lnL_min )
		
		!Average Energy
		
		k2 = real(0.0,16)
		
		do l = 1,d
		
			k2 = k2 + ( real(livepoints(j_min,l),16) ** 2 )
	
		end do
		
		E = E + w * a * k2
		
		!Replace the Live-Point with the smallest ln[L]
		
		!If directly sampling from the prior does not take too long
		if( test1 == 0 ) then
		
			i_new = 0
		
			test2 = 0
			do while( test2 == 0 )
		
				i_new = i_new + 1
			
				call prior(sample,d, k_max)
			
				if( lnL(sample, d, T, V, k_max) > lnL_min ) then
				
					do l = 1,d
					
						livepoints( j_min, l ) = sample(l)
						
					end do
				
					test2 = 1
				
				end if
				
				!If the loop takes too long, exist and use the other method
				if( i_new > 10000 ) then
				
					test1 = 1
					exit
				
				end if
		
			end do
		
		!If directly sampling from the prior does take too long
		else if (test1 == 1) then

			test2 = 0
			do while ( test2 == 0)
			
				b = real(d,16) * real(log( real(k_max) ) ,16) - real( log(real(fact)) ,16) - lnL_min
			
				call lnL_new( d, k_max, a, b, sample)
			
				if( lnL(sample, d, T, V, k_max) > lnL_min ) then
				
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
		
			test_lnL = lnL( sample, d, T, V, k_max )
			
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
	
	H = (H / Z) - real( log( Z ) ,16)
	oZ_H = Z * sqrt( H / real(K,16) )
	
	
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