using ITensors, ITensorMPS, ITensorCorrelators
using Dates, Printf
using HDF5
include("funcs.jl")

# number of spins
N = 150
sites = siteinds("S=1", N)

# parameters
J = 1.0

nw=now()
dirname_pa = pwd() * "/results/" * Dates.format(nw, "yyyymmdd_HHMMSS") * "_" * @sprintf("N%d_J%.1f", N, J)
mkpath(dirname_pa)


for theta in -pi:0.3:pi
    # directory
        # for results
        dirname = dirname_pa * "/" * @sprintf("theta_%.3f", theta)
        mkpath(dirname)
        # for MPS
        dirname_mps = pwd() * "/mps/" * @sprintf("theta%.3f_N%d_J%.1f", theta, N, J)
        mkpath(dirname_mps)

        # basic info
        for name in [dirname, dirname_mps]
            filename = name * "/info.txt"
            open(filename, "a") do fp
                println(fp, "*** basic info ***")
                @printf(fp, "theta = %.4f\n N = %d\n J = %.1f\n", theta, N, J)
            end
        end
    ##

    ## DMRG
    sites, psi, energy = run_dmrg(N, J, theta, dirname, dirname_mps)
    

    ###
    # Physical quantities
    i0 = Int(N/2)
    ##
        # energy per site
        filename = dirname_pa * "/energy_perSite.dat"
        open(filename, "a") do fp
            @printf(fp, "%f\t %f\n", theta, energy/N)
        end
        # expectation value
        for name in ["Sx", "Sy", "Sz"]
            Ss = expect(psi, name)
            filename = dirname * "/expval_" * name * ".dat"
            open(filename, "w") do fp
                [@printf(fp, "%d\t%f\n", i, s) for (i,s) in enumerate(Ss)]
            end
        end

        # correlation function
        for name in ["Sx", "Sy", "Sz"]
            cor = []
            for j in i0:N
                push!(cor, correlator(psi, (name,name), [(i0,j)]) )
            end
            filename = dirname * "/cor_" * name*name * ".dat"
            open(filename, "w") do fp
                for ele in cor
                    for (key, val) in ele
                        i,j = key[1], key[end]
                        @printf(fp, "%d\t%f\n", abs(j-i), val)
                    end
                end
            end
        end

        # string correlation
        function ITensors.op(::OpName"expSz", ::SiteType"S=1")
            o = zeros(ComplexF64, 3,3)
            phase = pi
            o[1,1] = exp(im*pi)
            o[2,2] = 1.0
            o[3,3] = exp(-im*pi)
            return o
        end
        function ITensors.op(::OpName"negativeSz", ::SiteType"S=1")
            o = zeros(3,3)
            o[1,1] = -1.0
            o[3,3] = 1.0
            return o
        end
        cor = []
        for j in i0+1:N
            # ops
            opseries = ["negativeSz"]
            [push!(opseries, "expSz") for b in i0+1:j-1]
            push!(opseries, "Sz")
            opseries = tuple(opseries...)
            # sites
            siteseries = [Tuple(i0:j)]

            push!(cor, correlator(psi, opseries, siteseries))
        end
        global stgop = 0.0+0.0*im
        filename = dirname * "/cor_stg.dat"
            open(filename, "w") do fp
            for ele in cor
                for (key, val) in ele
                    i,j = key[1], key[end]
                    @printf(fp, "%d\t%f\t%f\n", abs(j-i), real(val), imag(val))
                    if j == N 
                        stgop = real(val) 
                        @printf("theta = %f\t val = %f\t stgop = %f\n", theta, real(val), stgop) 
                    end
                end
            end
        end
        filename = dirname_pa * "/stgop.dat"
        open(filename, "a") do fp
            @printf(fp, "%f\t %f\t %f\n", theta, real(stgop), imag(stgop))
        end

        # entanglement entropy
        function myentropy(psi::MPS, b::Int)
            orthogonalize!(psi, i0)
            U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))
            SvN = 0.0
            for n in 1:dim(S, 1)
                p = S[n, n]^2
                if p > 1e-12
                    SvN -= p * log(p)
                end
            end
            return SvN
        end
        SvN = myentropy(psi, i0)
        for name in [dirname, dirname_mps]
            filename = name * "/info.txt"
            open(filename, "a") do fp
                println(fp, "*** SvN at the central site ***")
                println(fp, "SvN at i = $i0: $SvN")
            end
        end
        filename = dirname_pa * "/entropy.dat"
        open(filename, "a") do fp
            @printf(fp, "%f\t %f\n", theta, SvN)
        end

        # entanglement spectrum
        function entanglementspectrum(psi::MPS, b::Int)
            orthogonalize!(psi, i0)
            U, S, V = svd(psi[b], (linkind(psi, b-1), siteind(psi, b)))

            xi = Float64[]
            for n in 1:dim(S, 1)
                if S[n,n] > 1e-12
                    push!(xi, -2 * log(S[n,n]))
                end
            end
            sort!(xi)
            return xi
        end
        xi = entanglementspectrum(psi, i0)
        filename = dirname * "/entanglementspectrum.dat"
        open(filename, "w") do fp
            [@printf(fp, "%d\t%f\n", b,val) for (b,val) in enumerate(xi)]
        end
    ##


end
nothing
