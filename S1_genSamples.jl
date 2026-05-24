# file for generating samples from MPS


using ITensors, ITensorMPS, ITensorCorrelators
using Dates, Printf
using HDF5
include("funcs.jl")

# number of spins
N = 50

# parameters
J = 1.0

nw=now()
dirname_pa = pwd() * "/results/" * Dates.format(nw, "yyyymmdd_HHMMSS") * "_" * @sprintf("N%d_J%.1f", N, J)
mkpath(dirname_pa)


for theta in -pi:0.1:pi
    # directory
        # for results
        dirname = dirname_pa * "/" * @sprintf("theta_%.3f", theta)
        mkpath(dirname)
        # for MPS
        dirname_mps = pwd() * "/mps/" * @sprintf("theta%.3f_N%d_J%.1f", theta, N, J)
        mkpath(dirname_mps)
        filename_mps = dirname_mps * "/mps.h5"

        # basic info
        filename = dirname * "/info.txt"
        open(filename, "a") do fp
            println(fp, "*** basic info ***")
            @printf(fp, "theta = %.4f\n N = %d\n J = %.1f\n", theta, N, J)
        end
    ##

    sites, psi = run_dmrg(N, J, theta, dirname, dirname_mps)

    ### sampling
    num_samples = 500
    dirname_samples = dirname * "/samples_" * string(num_samples)
    mkpath(dirname_samples)
    
    # for S=1 siteinds, 1 --> Sz=+1, 2--> Sz=0. 3--> Sz=-1
    sz_map = Dict(1 => 1, 2 => 0, 3 => -1)

    samples = []
    normalize!(psi)
    for i in 1:num_samples
        vec = sample!(psi)
        push!(samples, [sz_map[s] for s in vec])
    end
    # output converted samples to file
    for i in 1:num_samples
        filename = dirname_samples * "/" * string(i) * ".dat"
        open(filename, "w") do fp
            for j in 1:N
                @printf(fp, "%d\n", samples[i][j])
            end
        end
    end

    ### benchmark for sampling
        function relative_error(vec_true::Vector, vec::Vector)
            return norm(vec_true - vec) / norm(vec_true) * 100
        end
        # actual string correlation
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
        function stringcor_true(psi::MPS)
            N = length(psi)
            i0 = 1
            cor_dict = []
            stgop = 0.0
            for j in i0+1:N
                # ops
                opseries = ["negativeSz"]
                [push!(opseries, "expSz") for b in i0+1:j-1]
                push!(opseries, "Sz")
                opseries = tuple(opseries...)
                # sites
                siteseries = [Tuple(i0:j)]

                push!(cor_dict, correlator(psi, opseries, siteseries))
            end
            cor = []
            for ele in cor_dict
                for (key, val) in ele
                push!(cor, real(val))
                end
            end
            return cor
        end
        # <Sz> from sampling
        function Sz_samples(samples::Vector)
            N = length(samples[1])
            num_samples = length(samples)
            
            expectSzs = zeros(N)
            for n in 1:N
                expectSzs[n] = sum(samples)[:][n]/num_samples
            end
            return expectSzs
        end
        # string cor from sampling
        function stringcor_samples(samples::Vector)
            N = length(samples[1])
            num_samples = length(samples)
            
            i0=1
            cor = zeros(N-1)
            for i in i0+1:N
                for n in 1:num_samples
                    cor[i-1] += real(-1.0*samples[n][i0]*exp(im*pi*sum(samples[n][i0+1:i-1]))*samples[n][i])
                end
                cor[i-1] /= num_samples
            end
            return cor
        end
        # true
        Sz_true = expect(psi, "Sz")
        stgcor_true = stringcor_true(psi)
        # from samples
        Sz_sample = Sz_samples(samples)
        stgcor_sample = stringcor_samples(samples)

        # relative error
        Sz_err = relative_error(Sz_true, Sz_sample)
        stgcor_err = relative_error(stgcor_true, stgcor_sample)

        # output 
        filename = dirname_samples * "/Sz.dat"
        open(filename, "w") do f
            for i in 1:length(Sz_true)
                @printf(f, "%f\t %f\n", Sz_true[i], Sz_sample[i])
            end
        end
        filename = dirname_samples * "/stgcor.dat"
        open(filename, "w") do f
            for i in 1:length(stgcor_true)
                @printf(f, "%f\t %f\n", stgcor_true[i], stgcor_sample[i])
            end
        end
        filename = dirname_samples * "/relative_error.txt"
        open(filename, "w") do f
            @printf(f, "relative error between true expectation value and value calculated with samples\n")
            @printf(f, "true expectation value calculated directly using operator and MPS\n")
            @printf(f, "value with samples calculated using samples generated from MPS\n\n")
            @printf(f, "<Sz> for each site: relative_error (percentage) = %f\n", Sz_err)
            @printf(f, "Cstg(r): relative_error (percentage) = %f\n", stgcor_err)
        end 
    ###

end
nothing
