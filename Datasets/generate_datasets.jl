using XSim
import Random
Random.seed!(95616)

# Build genome and phenome
build_genome(n_chr=10, n_loci=100)
build_phenome([3, 8])

# Initialize founders
n_sires = 3
n_dams  = 20
sires   = Founders(n_sires)
dams    = Founders(n_dams)

# Define parameters
args     = Dict(# mating
                :nA               => 3,
                :nB_per_A         => 5,
                :n_per_mate       => 2,
                :ratio_malefemale => 1.0,
                # selection
                :h2               => [.8, .5],
                :weights          => [.6, .4],
                # breeding
                :n_gens           => 5,
                :n_select_A   => 3,
                :n_select_B => 20
                :select_all_gens  => true )

args_A  = Dict(# Mating
                :nA               => 3,
                :nB_per_A         => 5,
                :ratio_malefemale => 2,
                # Selection
                :is_random        => true,
                # Breeding
                :n_gens           => 5,
                :n_select_A   => 3,
                :select_all_gens=>true)
# Breeding program
males, females   = breed(sires, dams; args_A...)

ped = get_pedigree(males)

Genlist=[ females ]
 for i in Genlist 
 temped=XSim.get_pedigree(i)
 global ped=vcat(ped,temped)
 end

progeny = males + females
get_pedigree(progeny) 
