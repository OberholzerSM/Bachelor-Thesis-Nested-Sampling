program linearMH

	implicit none
	
	integer (kind = 8) :: samplesize
	integer :: i, j, datasize, u
	real (kind = 16) :: om, ob, ln_f, res, A, two, wj, wc
	real (kind = 16), dimension(100) :: x, y, dx, dy, time
	real (kind = 16), parameter :: pi = 4*atan(1.0)
	real (kind = 16), dimension(6) :: B
	real (kind = 16), dimension(6,6) :: C
	
	call random_seed()
	
	datasize = 100
	samplesize = 10**8
	u = 6
	
	!Lese die ersten 100 Daten ein. (Nutze nur 100 Daten, um das Programm zu beschleunigen.)
	open(1, file = 'CODE_X.txt', status = 'old')
	open(2, file = 'CODE_Y.txt', status = 'old')
	open(3, file = 'CODE_dX.txt', status = 'old')
	open(4, file = 'CODE_dY.txt', status = 'old')
	open(5, file = 'CODE_MJD.txt', status = 'old')
	do i = 1,datasize
		read(1,*) x(i)
		read(2,*) y(i)
		read(3,*) dx(i)
		read(4,*) dy(i)
		read(5,*) time(i)
		
		!Falls der Fehler als 0 angegeben ist, ersetzte diesen durch einen Fehler von 10^-6.
		if (dx(i) == 0 ) then
			dx(i) = real(10,16)**-6
		end if
		
		if (dy(i) == 0 ) then
			dy(i) = real(10,16)**-6
		end if
		
	end do
	close(1)
	close(2)
	close(3)
	close(4)
	close(5)
	
	!Bestimme die Bilinearfaktoren A, B, C
	
	two = real(2.0,16) !Zwei als real type 16
	wj = (two*pi) / real(365.25,16) !Winkelgeschwindigkeit des jährlichen Anteil j
	wc = (two*pi) / real(432.25,16) !Winkelgeschwindigkeit der Chandler-Periode c
	
	A = 0
	do i = 1,6
	
		B(i) = 0
		
		do j = 1,6
		
			C(i,j) = 0
			
		end do
	
	end do
	
	
	do i = 1,datasize
	
		A = A - real(Log( two*pi*dy(i)*dy(i) ),16) - ( x(i)**2 / (two*dx(i)**2) ) - ( y(i)**2 / (two*dy(i)**2) )
		
		B(1) = B(1) + (x(i) / (dx(i)**2) ) * real(cos(wj*time(i)),16) - (y(i) / (dy(i)**2) ) * real(cos(wj*time(i)),16)
		B(2) = B(2) - (x(i) / (dx(i)**2) ) * real(sin(wj*time(i)),16) - (y(i) / (dy(i)**2) ) * real(sin(wj*time(i)),16)
		B(3) = B(3) + (x(i) / (dx(i)**2) ) * real(cos(wc*time(i)),16) - (y(i) / (dy(i)**2) ) * real(cos(wc*time(i)),16)
		B(4) = B(4) - (x(i) / (dx(i)**2) ) * real(sin(wc*time(i)),16) - (y(i) / (dy(i)**2) ) * real(sin(wc*time(i)),16)
		B(5) = B(5) + (x(i) / (dx(i)**2) )
		B(6) = B(6) + (y(i) / (dy(i)**2) )
		
		C(1,1) = real(cos(wj*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (-1 / (two * dy(i)**2))  )
		C(1,2) = -two * real(cos(wj*time(i))**2,16) * real(sin(wj*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (1 / (two * dy(i)**2))  )
		C(1,3) = -two * real(cos(wj*time(i))**2,16) * real(cos(wc*time(i))**2,16) * ( (1 / (two * dx(i)**2)) + (1 / (two * dy(i)**2))  )
		C(1,4) = -two * real(cos(wj*time(i))**2,16) * real(sin(wc*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (1 / (two * dy(i)**2))  )
		C(1,5) = -two * real(cos(wj*time(i))**2,16) *  (1 / (two * dx(i)**2))
		C(1,6) = -two * real(cos(wj*time(i))**2,16) *  (-1 / (two * dy(i)**2))
		
		C(2,2) = real(sin(wj*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (-1 / (two * dy(i)**2))  )
		C(2,3) = -two * real(sin(wj*time(i))**2,16) * real(cos(wc*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (1 / (two * dy(i)**2))  )
		C(2,4) = -two * real(sin(wj*time(i))**2,16) * real(sin(wc*time(i))**2,16) * ( (1 / (two * dx(i)**2)) + (1 / (two * dy(i)**2))  )
		C(2,5) = -two * real(sin(wj*time(i))**2,16) *  (-1 / (two * dx(i)**2))
		C(2,6) = -two * real(sin(wj*time(i))**2,16) *  (-1 / (two * dy(i)**2))
		
		C(3,3) = real(cos(wc*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (-1 / (two * dy(i)**2))  )
		C(3,4) = -two * real(cos(wc*time(i))**2,16) * real(sin(wc*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (1 / (two * dy(i)**2))  )
		C(3,5) = -two * real(cos(wc*time(i))**2,16) *  (1 / (two * dx(i)**2))
		C(3,6) = -two * real(cos(wc*time(i))**2,16) *  (-1 / (two * dy(i)**2))
		
		C(4,4) = real(sin(wc*time(i))**2,16) * ( (-1 / (two * dx(i)**2)) + (-1 / (two * dy(i)**2))  )
		C(4,5) = -two * real(sin(wc*time(i))**2,16) *  (-1 / (two * dx(i)**2))
		C(4,6) = -two * real(sin(wc*time(i))**2,16) *  (-1 / (two * dy(i)**2))
		
		C(5,5) = ( (-1 / (two * dx(i)**2)) + (-1 / (two * dy(i)**2))  )
		C(6,6) = ( (-1 / (two * dx(i)**2)) + (-1 / (two * dy(i)**2))  )
	
	end do

	
	
	
	call MH_Sampler(datasize, samplesize, x, y, dx, dy, time, A, B, C)


end program linearMH


subroutine GaussSample(mu, sigma, x) !Generiert ein Sample aus einer Gauss-Verteilung. Inputs: mu, sigma. Output: x

	implicit none

	real (kind = 16) :: u1,u2,z, mu, sigma, x
	real (kind = 16), parameter :: pi = 4*atan(1.0)
	
	call random_number(u1)
	call random_number(u2)

	z = sqrt( -2*log(u1) )*cos(2*pi*u2)
	x = sigma*z + mu
	
end subroutine GaussSample

real (kind = 16) function ln_f(theta, A, B, C) !Berechnet den ln der Wahrscheinlichkeit, dass die Parameter bei gegeben Daten stimmen.

	integer :: i, j
	real (kind = 16) :: ln_L, A
	real (kind = 16), dimension(6) :: theta, B
	real (kind = 16), dimension(6,6) :: C
	
	
	ln_L = A
	
	do i = 1,6
	
		ln_L = ln_L + B(i)*theta(i)
		
		do j = 1, (7-i)
			
			ln_L = ln_L + C(i,j) * theta(i) * theta(j)

		end do
	
	end do
	
	!Füge dem log_L noch den ln vom Prior hinzu (pi = 1 / ( 1 - (-1) )^6 = 1/64)
	ln_f = ln_L - real(log(64.0),16)

end function ln_f

subroutine MH_Sampler(datasize, samplesize, x, y, dx, dy, time, A, B, C) !Erstellt ein .txt File, in den die Samplings aufgelistet sind.

	implicit none
	
	integer (kind = 8) :: i, samplesize
	integer :: j,  datasize, test, k
	integer , dimension (6) :: rate
	real (kind = 16) :: ln_f, r, A
	real (kind = 16), dimension(100) :: x, y, dx, dy, time
	real (kind = 16), dimension(6) :: theta1, theta2, o, B
	real (kind = 16), dimension(6,6) :: C
	
	do i = 1,6
		o(i) = real(0.00000000005,16) !Startwerte der Breiten
		rate(i) = 0
	end do
	!Startwerte von theta1 und theta2
	theta1(1) = real(0.034971,16)
	theta1(2) = real(0.0995418,16)
	theta1(3) = real(-0.0724511,16)
	theta1(4) = real(-0.0207436,16)
	theta1(5) = real(0.0825917,16)
	theta1(6) = real(0.34662,16)
	
	theta2(1) = real(0.034971,16)
	theta2(2) = real(0.0995418,16)
	theta2(3) = real(-0.0724511,16)
	theta2(4) = real(-0.0207436,16)
	theta2(5) = real(0.0825917,16)
	theta2(6) = real(0.34662,16)

	test = 0
	k = 0
	
	open(6, file = 'Pol_samplings.txt', status = 'old')
	
	do i = 1,samplesize
		
		do j = 1,6
		
			call random_number(r)
			call GaussSample(theta1(j), o(j), theta2(j))
			
			if ( ln_f(theta2, A, B, C) - ln_f(theta1, A, B, C)  > log(r) ) then
			
			theta1(j) = theta2(j)
			rate(j) = rate(j) + 1
			
			end if
		
		end do
		
		k = k + 1
		if( k == 100 ) then
		
			write(6,*) theta1(1), theta1(2), theta1(3), theta1(4), theta1(5) ,theta1(6)
			k = 0
			
		end if
			
			
		!Passe alle 1'000 Schritte die Breiten an.
		if (test == 0 .and. mod(i,1000) == 0 ) then
			
			do j = 1,6
			
				if (rate(j) < real(600,16) .and. rate(j) > 0) then
					o(j) = o(j) - real(0.000000000001,16)
				end if
		
				if (rate(j) > real(800,16) .or. o(j) <= 0) then
					o(j) = o(j) + real(0.000000000001,16)
				end if
			
			end do
		
			!Falls alle Raten zwischen 60% und 80% sind, passe die Breiten nicht mehr an.
			if( minval(rate) > 600 .and. maxval(rate) < 800 ) then
				
				test = 1
				print *, 'Daten koennen ab dem ', i, '. Schritt genutzt werden.'
			
			end if
			
			do j = 1,6
				rate(j) = 0
			end do
			
		end if
		!Gib den Fortschritt des Programmes an.
		if ( Mod(i, samplesize/100) == 0) then
			print *, (100*i/samplesize), '%'
		end if
		
	end do
	
	close(6)
	
	if (test == 0) then
		print *, 'Warnung: Breitenanpassung nicht abgeschlossen!'
	end if
	
	if (test == 1) then
		print *, "Raten:"
		do j = 1,6
			print *, real(rate(j)) / real(samplesize)
		end do
	end if

	print *, "Breiten:"
	do j = 1,6
			print *, o(j)
		end do

end subroutine MH_Sampler