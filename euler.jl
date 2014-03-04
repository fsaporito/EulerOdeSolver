# Step Calculator For The Selected Interval
# Basic stepCalculation, to be improved with a better way to find the step
function stepCalculator (first::FloatingPoint, last::FloatingPoint)

	# Step calculation
	step = abs(last - first) / 1000000 
	
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
	
		Tk = Tk + h # New Tk's Value Calculation
	
		Yk = Yk1 # New Yk's Value Calculation
	
	end
	
	tout = Tk - h # Last Time Value
	yout = Yk # Y Requested Value
	
	return (tout, yout)
	
end
	
