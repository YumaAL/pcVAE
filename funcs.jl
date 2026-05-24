using ITensors, ITensorMPS
using Printf
using HDF5

function run_dmrg(N::Int, J::Float64, theta::Float64, dirname::String, dirname_mps::String)
    filename_mps = dirname_mps * "/mps.h5"

    if isfile(filename_mps) == false
        sites = siteinds("S=1", N)
        ## DMRG
            filename = dirname_mps * "/info.txt"
            open(filename, "a") do fp
                println(fp, "*** basic info ***")
                @printf(fp, "theta = %.4f\n N = %d\n J = %.1f\n", theta, N, J)
            end
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
            energy_tol = 1e-8
            nsweeps = 50
            maxdim = [1,5,10,20,30, 40, 50, 60, 70, 80, 100]
            cutoff = [1E-10]
            noise  = [0.0]
            observer = DMRGObserver(; energy_tol=energy_tol)
            energy,psi = dmrg(H,psi0;nsweeps,maxdim,cutoff,noise, observer=observer)
            energy_perSite = energy/N
            bonddim = maxlinkdim(psi)

            # sweep schedule
            for name in [dirname, dirname_mps]
                filename = name * "/sweeps.txt"
                open(filename, "w") do fp
                    println(fp, "nsweeps: $nsweeps")
                    println(fp, "maxdim:  $maxdim")
                    println(fp, "cutoff:  $cutoff")
                    println(fp, "noise:   $noise")
                    println(fp, "energy tol: $energy_tol")
                end
            end
            # basic info for DMRG
            for name in [dirname, dirname_mps]
                filename = name * "/info.txt"
                open(filename, "a") do fp
                    println(fp, "*** DMRG calculation ***")
                    println(fp, "ground energy: $energy")
                    println(fp, "per site: $energy_perSite")
                    println(fp, "max link dim of MPS: $bonddim")
                end
            end
            println("per site: $energy_perSite")

            # output of MPS
            h5open(filename_mps, "w") do file
                write(file, "psi", psi)
                write(file, "sites", sites)
            end
        ##
        return sites, psi, energy
    else
        println("MPS file exists")
        ### read MPS from file
        psi = h5open(filename_mps, "r") do fp
            read(fp, "psi", MPS)
        end

        sites = h5open(filename_mps, "r") do fp
            read(fp, "sites", Vector{Index})
        end

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
        energy = inner(psi', H, psi)
        return sites, psi, energy
    end
end