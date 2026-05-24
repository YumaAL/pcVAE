using ITensors, ITensorMPS, ITensorCorrelators
using Dates, Printf



# number of spins
N = 300
sites = siteinds("S=1", N)

# parameters
theta = atan(1.0/3.0)
#theta = 0.0
J = 1.0

# directory
    # for results
    nw=now()
    dirname = pwd() * "/results/" * Dates.format(nw, "yyyymmdd_HHMMSS") * "_" * @sprintf("theta%.3f_N%d_J%.1f", theta, N, J)
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
    # MPO
        os = OpSum()
        for i=1:N-1
            # cos term
            os += J*cos(theta), "Sx",i, "Sx",i+1
            os += J*cos(theta), "Sy",i, "Sy",i+1
            os += J*cos(theta), "Sz",i, "Sz",i+1

            # sin term
            os += J*sin(theta), "Sx",i, "Sx",i, "Sx",i+1, "Sx",i+1
            os += J*sin(theta), "Sy",i, "Sy",i, "Sy",i+1, "Sy",i+1
            os += J*sin(theta), "Sz",i, "Sz",i, "Sz",i+1, "Sz",i+1
            os += J*sin(theta), "Sx",i, "Sy",i, "Sx",i+1, "Sy",i+1
            os += J*sin(theta), "Sy",i, "Sx",i, "Sy",i+1, "Sx",i+1
            os += J*sin(theta), "Sx",i, "Sz",i, "Sx",i+1, "Sz",i+1
            os += J*sin(theta), "Sz",i, "Sx",i, "Sz",i+1, "Sx",i+1
            os += J*sin(theta), "Sy",i, "Sz",i, "Sy",i+1, "Sz",i+1
            os += J*sin(theta), "Sz",i, "Sy",i, "Sz",i+1, "Sy",i+1
        end
        H = MPO(os, sites)
    ##

    #initial state
    psi0 = random_mps(sites)

    # DMRG
    nsweeps = 50
    maxdim = [1,5,10,20,30, 40, 50, 60, 70, 80, 100]
    cutoff = [1E-10]
    noise  = [0.0]
    energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff,noise)
    energy_perSite = energy/N

    # sweep schedule
    for name in [dirname, dirname_mps]
        filename = name * "/sweeps.txt"
        open(filename, "w") do fp
            println(fp, "nsweeps: $nsweeps")
            println(fp, "maxdim:  $maxdim")
            println(fp, "cutoff:  $cutoff")
            println(fp, "noise:   $noise")
        end
    end
    # basic info for DMRG
    for name in [dirname, dirname_mps]
        filename = name * "/info.txt"
        open(filename, "a") do fp
            println(fp, "*** DMRG calculation ***")
            println(fp, "ground energy: $energy")
            println(fp, "per site: $energy_perSite")
        end
    end
    println("per site: $energy_perSite")
    trueE0 = -1.0/3.0*2.0*cos(theta) # per site
    println("true energy per site: $trueE0")

##

###
# Physical quantities
i0 = Int(N/2)
##
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
    filename = dirname * "/cor_stg.dat"
        open(filename, "w") do fp
        for ele in cor
            for (key, val) in ele
                i,j = key[1], key[end]
                @printf(fp, "%d\t%f\t%f\n", abs(j-i), real(val), imag(val))
            end
        end
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



nothing
