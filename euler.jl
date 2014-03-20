# Step Calculator For The Selected Interval
# Basic stepCalculation, to be improved with a better way to find the step
function stepCalculator (first::FloatingPoint, last::FloatingPoint)

	# Step calculation
	step = abs(last - first) / 100000 
	
	return step

end


# Explicit Euler Method For First Order Ordinary Differential Equations
function eulerExplicit {T} (F::Function, tspan::AbstractVector, y0::AbstractVector{T})

	tfirst = tspan[1] # First Interval's Value
    tlast = tspan[end] # Last Interval's Value
	
	h = stepCalculator (tfirst, tlast) #Step Calculation
		
	return eulerExplicit (F, tspan, y0, h)
	
end


# Explicit Euler Method where the user defines his own step value
function eulerExplicit {T} (F::Function, tspan::AbstractVector, y0::AbstractVector{T}, h::FloatingPoint)

	tlast = tspan[end] # Last Interval's Value
	
	Tk = tspan[1] # First Interval's Value
	
	Yk = y0 # Starting Point	
	
	while Tk != tlast
	
		Yk1 = Yk + h*F(Tk, Yk) # Explicit Euler Calculation
		
		Tk = Tk + h # New Tk's Value
		
		Yk = Yk1 # New Yk's Value
	
	end
	
	tout = Tk - h # Last Time Value
	yout = Yk # Y Requested Value
	
	return (tout, yout)
	
end


# Implicit Euler Method For First Order Ordinary Differential Equations
function eulerImplicit {T} (F::Function, tspan::AbstractVector, y0::AbstractVector{T})

	tfirst = tspan[1] # First Interval's Value
    tlast = tspan[end] # Last Interval's Value
	
	h = stepCalculator (tfirst, tlast) #Step Calculation
		
	return eulerImplicit (F, tspan, y0, h)
	
end


# Implicit Euler Method where the user defines his own step value
function eulerImplicit {T} (F::Function, tspan::AbstractVector, y0::AbstractVector{T}, h::FloatingPoint)

	tlast = tspan[end] # Last Interval's Value
	
	Tk = tspan[1] + h # First Interval's Value
	
	Yk = y0 # Starting Point	
	
	while Tk != tlast
		
		# Solving G(Y*) = 0, To Find Y* = YK1	
		Yk1 = dicotomic (G, F, Tk, Yk, h) # Implicit Euler Calculation
		
		Tk = Tk + h # New Tk's Value 
		
		Yk = Yk1 # New Yk's Value
	
	end
	
	tout = Tk - h # Last Time Value
	yout = Yk # Y Requested Value
	
	return (tout, yout)
	
end


# G(Y*) = Y* -Y0 -h*F(Tk; Y*)
function G {T} (F::Function, t::FloatingPoint, y0::AbstractVector{T}, h::FloatingPoint) end

	return (Y - y0 -h*F(t, Y))

end


# Solves G(Y) = 0, where G = function (F, t, Y, h)
function dicotomic {T} (G::Function, F::Function, t::FloatingPoint, y0::AbstractVector{T}, h::FloatingPoint)
	
	
	precision = h/100
		
	m = y0 # First Root's Approssimation
	mold = y0 - precision - 1 # Fake Old Approssimation To Enter The Loop
	
	while (abs (m - mold) ) > precision
	
		a = y0 # First Extremum
		A = G (F, t, a, h) # a Image
	
		b = a + F(t+h, y0) # Second Extremum
		B = G (F, t, b, h) # b Image
		
		if A*B < 0 # Sign Check

			error ("Same Image Sign !!!")			
			
		end
			
		mold = m
	
		m = (a + b) / 2 # Middle Point
		M = G (F, t, m, h) # m Image
	
		if M == 0 # Root found, Y* = m
	
			break
		
		end
	
		if A*M < 0 # m Is Nearer To The Root Than b
	
			# Switchng m with b
			b = m
			B = M
		
		end
	
		if B*M < 0 # m Is Nearer To The Root Than a
		
			# Switchng m with a
			a = m 
			A = M
		
		end
	
	end
		
	return (m)
	
end	
