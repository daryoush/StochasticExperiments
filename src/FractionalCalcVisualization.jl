
using Interact, Plots, DifferentialEquations, Cubature


### PLot of chars in http://elib.mi.sanu.ac.rs/files/journals/vm/57/vmn57p5-7.pdf

## ---
function plotFunAndIntegral(fun, from, to, range)
	@manipulate for from=slider(range, value=from), to=slider(range, value=to)

		fplot=plot(fun ,from, to, framestyle = :origin,  xlims=range) #,  fill = (to, from, :gray))

		areafrom(fun, x) = y-> hquadrature(fun, x, y ;
		                        reltol=1e-8, abstol=0, maxevals=0)

		plot!( y-> areafrom(fun, from)(y)[1], from, to)
		vbox(fplot)
	end
end

plotFunAndIntegral(x-> -x*(x-4), -0.8, 4.4, (-5:.1:5))

## ---
plotFunAndIntegral(x-> 1/x, 1, 5, (0.1:.1:5))

## ---
# Solve the 1/t differential eq with u0=0  observe that the integral of 1/x is the solution
# to the diff eq.

plot(x->1/x, 1, 5)

f(u,p,t) = 1/t
u0 = 0
tspan = (1,5)
prob = ODEProblem(f,u0,tspan)
sol=solve(prob)

plot!(sol)

## ---

f(x)=x
g(x)=exp(x)
q(t)=x->exp(t-x)
p(t)=x-> x*exp(t-x)

plot(f, 0,5, grid=true, gridlinewidth=1, ylims=(0,20), xlims=(0:5), labels="f")
plot!(g, 0,5, grid=true, gridlinewidth=1, ylims=(0,20), xlims=(0:5), labels="g")
plot!(q(3.2), 0,5, grid=true, gridlinewidth=1, ylims=(0,20), xlims=(0:5), labels="p")
plot!(p(3.2), 0,5, grid=true, gridlinewidth=1, ylims=(0,20), xlims=(0:5), labels="q")

# fun is a function of form:  f(t)(x), t would be the y
conv(fun, x) = y-> hquadrature(fun(y), x, y ;
						reltol=1e-8, abstol=0, maxevals=0)
@show conv(p, 0)(3.2)
plot!( y-> conv(p, 0)(y)[1], 0,5, grid=true, gridlinewidth=1, ylims=(0,20), xlims=(0:5), labels="convolution")



## --
using SpecialFunctions
@manipulate for a=slider(0:0.1:2, value=1.52)



f(x)=x
g(x)=.5x^2
ylims=(0,20)
xlims=(0:8)

# f(x)=sin(pi*x)
# g(x)=cos(pi*x)
# ylims=(-2,2)
# xlims=(0:8)

pp(t, h)= x-> f(x)*(t-x)^(h-1)
#pp(t, h)= x-> f(t-x)*((x)^(1-h))   ## this is what the document says, but I am not getting the correct result

fracint(fun, x) = (y,h)-> ((hquadrature(fun(y, h), x, y ;
						reltol=1e-8, abstol=0, maxevals=0))[1]/gamma(h))

fp=plot(f, 0,8, grid=true, gridlinewidth=1, ylims=ylims, xlims=xlims, labels="f")
plot!(g, 0,8, grid=true, gridlinewidth=1,ylims=ylims, xlims=xlims, labels="g")
plot!( y-> fracint(pp, 0)(y, a), 0,8, grid=true, gridlinewidth=1, ylims=ylims, xlims=xlims, labels="fractionalInt")
@show a
@show fracint(pp, 0)(3.2, a)
#do the integration by small increments
@show sum(pp(3.2,a)(t)*.0001/gamma(a) for t in 0:0.0001:3.2)

vbox(fp)
end

## ---
@manipulate for a=slider(-1:0.1:1, value=.868)



f(x)=x
fprime(x)=1
ylims=(-1,5)
xlims=(-1:10)

# f(x)=sin(pi*x)
# g(x)=cos(pi*x)
# ylims=(-2,2)
# xlims=(0:8)

pp(t, h)= x-> fprime(t-x)*(x^(-h))

fracdriv(fun, x) = (y,h)-> ((hquadrature(fun(y, h), x, y ;
						reltol=1e-8, abstol=0, maxevals=0))[1]/gamma(1-h))

fp=plot(f, 0,8, grid=true, gridlinewidth=1, ylims=ylims, xlims=xlims, labels="f")
plot!( y-> fracdriv(pp, 0)(y, a), 0,8, grid=true, gridlinewidth=1, ylims=ylims, xlims=xlims, labels="fractionaldriv")

@show fracdriv(pp, 0)(3.2, a)


vbox(fp)
end
