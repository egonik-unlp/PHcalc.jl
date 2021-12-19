using Optim, Plots
module pHcalc

export Neutral,Acid, System, α,pHsolve, minimise

struct Acid
	ka::Vector{Float64}
	conc::Float64
	charge::Vector{Int64}
	pKA::Vector{Float64}
	_ka::Vector{Float64}
	function Acid(ka, conc, charge::Int64)
			if typeof(ka) <: Number
				ka=[ka]
			end
			_ka=copy(ka)
			sort!(ka, rev=:true)
			prepend!(_ka,1.)
			charge=charge:-1:charge-(length(ka)) |> collect
			# charge=LinRange(charge,-length(ka), length(ka))
			new(ka, conc,charge , -log10.(ka),_ka )
		end
end
	
struct Neutral
	charge::Union{Vector{Int64},Int64}
	conc::Float64
end
	
struct System
	species
	function System(species...)
		new(species)
	end
end

function α(self::Acid,pH) # solo para un valor de pH, se podria pensar para un array de pHs.
		h3o=10.0^(-pH)
		power=length(self._ka)-1:-1:0
		h3o_power=(x->h3o^x).(power)
		Ka_prod=cumprod(self._ka)
		h3o_Ka=h3o_power.*(Ka_prod)
		h3o_Ka./sum(h3o_Ka)
	end

function α(self::Neutral,pH) # Devielve 1 porque no tiene disociacion, devuelvo en array para que se parezca a la otra implementación.
	[1]
end
	
function pHsolve(sys)
	function minimise(pH)
		h3o=10.0^(-pH)
     	oh = (10.0^(-14))/h3o
	 	x = (h3o - oh)
		for specie in sys.species
			x+=(specie.conc.*specie.charge.*α(specie,pH))|> sum
		end
	 	abs(x)
	end
	
	function phguess()
		phv=1:.1:14
		minimise.(phv) |> argmin |> x -> getindex(phv, x)
	end
	guess=phguess()
	optimize(x->minimise(first(x)), [guess], BFGS()).minimizer |> first
end


function minimise(sys::System,pH)
	h3o=10.0^(-pH)
     	oh = (10.0^(-14))/h3o
	 x = (h3o - oh)
	for specie in sys.species
		x+=(specie.conc.*specie.charge.*α(specie,pH))|> sum
	end
	 abs(x)
end

function plot_distribution(acid::Acid)
	interval=1:.01:14
	plot= (x-> α(acid, x)).(interval) |> x-> hcat(x...)' |> x-> Plots.plot(interval, x, legend=:false)
	xlabel!(plot, "pH")
	ylabel!(plot, "Fracción de concentración")

end




end
