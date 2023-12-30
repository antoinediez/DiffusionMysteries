using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Symbolics, IfElse
using DifferentialEquations
using JLD2
using BenchmarkTools

println("Initialize everything...")

#-------- Add the scripts------------------------------------------------------#
include("reaction_taxis_growth.jl")
include("plot_save.jl")
include("ODE_func.jl")
#------------------------------------------------------------------------------#

simu_name = init_directory(simu_name="plot_RD_test")

#-------- Discretization parameters -------------------------------------------#
Nx = 100    # Number of cells in the x-direction (bulk and surface)
Ny = 100    # Number of cells in the y-direction (bulk)
Lx = 100.0     # Length of the domain in the x-direction (bulk and surface)
Ly = 100.0   # Length of the domain in the y-direction (bulk)
dx = Lx/Nx    # Cell length in the x-direction
dy = Ly/Ny    # Cell length in the y-direction

#------------------------------------------------------------------------------#


#-------- Interaction parameters ----------------------------------------------#


## Diffusion coefficients

# Du = 0.002
# Dv = 0.05
# Du = 1.0
# Dv = 20.0
# Dv = 1.0
# Dv = 100.0
# Du = 0.01
# Dv = 1.0
# Du = 0.1
# Dv = 0.1

#------------------------------------------------------------------------------#


#-------- Time parameters -----------------------------------------------------#

T = 380.0 # Simulation time 
# dt = 0.0001    # Discretisation step

#------------------------------------------------------------------------------#


#------------------- Model  ---------------------------------------------------#



#===================== REACTION ===========================#


# reac_name = "Linear_KS"
# # c_eq = 10.0
# c_eq = 12.0
# α_B = (a_S+b_S)^2/b_S*c_eq*(1-(a_S+b_S)/c_eq)
# β_B = (a_S+b_S)^2/b_S*c_eq
# uB_eq = (1 - α_B/β_B) * c_eq
# vB_eq = 1/β_B * c_eq
# param_reac = Dict(
#     :reac_name => "Linear_KS",
#     :α => α_B,
#     :β => β_B,
#     :γ => γ,
#     :c_eq => c_eq,
#     :sc => 1.0,
#     :taxis_decay_rate => 0.1
# )

# reac_name = "Custom Activator-Inhibitor"
# uB_eq = uS_eq
# vB_eq = vS_eq
# α_B = 8.0
# β_B = 0.1
# δ_B = 3.0
# c_eq = δ_B * vB_eq
# γ_B = α_B/(1 + β_B*vB_eq) * c_eq^4/(1 + c_eq^4)
# param_reac = Dict(
#     :α => α_B,
#     :β => β_B,
#     :γ => γ_B,
#     :δ => δ_B,
#     :c_eq => c_eq,
#     :sc => 1.0,
#     :taxis_decay_rate => 0.01
# )

# reac_name = "Activator-Inhibitor (u->c)"
# param_reac = Dict(
#     :act_synth => 0.2387,
#     :act_autocat => 0.08,
#     :inhib_act => 0.08,
#     :max_act_synth => 10000.0,
#     :act_deg => 0.03,
#     :inhib_synth => 0.15,
#     :inhib_autocat => 0.0,
#     :inhib_act_crea => 0.11,
#     :max_inhib_synth => 10000.0,
#     :inhib_deg => 0.06,
#     :sc => 0.2,
# )

reac_name = "Schnakenberg"
param_reac = (
    a = 0.2305,
    b = 0.7695,
    sc = 0.0
)

# ap = 1.97
# bp = 2.0

# param_reac = (
#     a = bp-ap,
#     b = ap+bp,
#     sc = 1.0
# )

u_eq = param_reac[:a] + param_reac[:b]
v_eq = param_reac[:b] / (param_reac[:a] + param_reac[:b])^2

# reac_name = "Fish"
# param_reac = Dict(
#     :a => 0.08,
#     :b => -0.08,
#     :c => 0.02,
#     :d => 0.03,
#     :e => 0.1,
#     :f => 0.06,
#     :g => -0.12,
#     :synUmax => 0.2,
#     :synVmax => 0.5,
#     :sc => 1.0
# )


#=======================================================#

#------------------- Initial condition  ---------------------------------------#

U_init = 0.5 .+ rand(Nx,Ny,2)
# The individual solutions can be recovered with the following index table
   
u_index = :,:,1 
v_index = :,:,2

# U_init[u_index...] = u_eq * (ones(Nx,Ny) + 0.1 * rand(Nx,Ny))
# U_init[v_index...] = v_eq * (ones(Nx,Ny) + 0.1 * rand(Nx,Ny))

U_init[u_index...] = u_eq * ones(Nx,Ny)
U_init[v_index...] = v_eq * ones(Nx,Ny)
G = 150
s2 = (Lx/20)^2
for j in 1:Ny
    for i in 1:Nx
        xi = dx/2 + (i-1)*dx
        yj = dy/2 + (j-1)*dy
        U_init[i,j,1] += G/(2*pi*s2) * exp(-((Lx/2 - xi)^2 + (Ly/2 - yj)^2)/(2*s2))
        U_init[i,j,2] += G/(2*pi*s2) * exp(-((Lx/2 - xi)^2 + (Ly/2 - yj)^2)/(2*s2))
    end
end

#------------------------------------------------------------------------------#


#--------------------- Create an ODE problem ----------------------------------#

sol = Vector{Any}(undef,4)
saveat = 0.25
algo = ROCK2()

for s in 1:4

    if s==1
        Du = 1.0
        Dv = 1.0
        param_reac = (
            a = 0.2305,
            b = 0.7695,
            sc = 0.0
        )
    elseif s==2
        Du = 1.0
        Dv = 20.0
        param_reac = (
            a = 0.2305,
            b = 0.7695,
            sc = 0.0
        )    
    elseif s==3
        Du = 1.0
        Dv = 1.0
        param_reac = (
            a = 0.2305,
            b = 0.7695,
            sc = 1.0
        )
    elseif s==4
        Du = 1.0
        Dv = 20.0
        param_reac = (
            a = 0.2305,
            b = 0.7695,
            sc = 1.0
        )
    end

    p = (Du,Dv,dx,dy,Nx,Ny,param_reac)    # Simulation parameters

    println("Initialize ODE problem...")

    dU_init = copy(U_init)
    jac_sparsity = Symbolics.jacobian_sparsity((du,u)->func!(du,u,p,0.0),dU_init,U_init)    # Automatic sparsity detection

    ode_func = ODEFunction(func!,jac_prototype=float.(jac_sparsity))
    prob = ODEProblem(ode_func,U_init,(0.,T),p)

    # prob = ODEProblem(func!,U_init,(0.,T),p)    # Without sparsity information

    println("Ok let's solve the ODE...")

    #------------------------------------------------------------------------------#


    #----------- Finallly run the simulation and save data ------------------------#
    # saveat = 3.0
    sol[s] = solve(prob,algo,alg_hints=[:stiff],saveat=saveat,dtmin=0.001,force_dtmin=true;progress=true,progress_steps=1)
    # sol = solve(prob,algo,saveat=10.0;dt=0.001,progress=true,progress_steps=1)    # For explicit Euler, specify dt
end

println("Solving is done, now plotting...")

umin = 0.5
umax = 1.5
vmin = 0.5
vmax = 1.5


plot_save_sol_4(
    sol,[U_init,U_init,U_init,U_init],u_index,v_index,dx,dy;
    colorrange_u=(umin,umax),
    colorrange_v=(vmin,vmax),
    dir=simu_name,video_name=simu_name,
    dt=saveat)

# println("Save data...")

# data = Dict(
#     "simu_name" => simu_name,
#     "reac_name" => reac_name,
#     "sol_final" => sol[end],
#     "T" => T,
#     "Du" => Du, 
#     "Dv" => Dv,
#     "dx" => dx,
#     "dy" => dy,
#     "Nx" => Nx,
#     "Ny" => Ny,
#     "param_reac" => param_reac
# )

# save(simu_name*"/"*simu_name*".jld2",data)

println("All done.")

#=======================================================#
#================ SECOND RUN ===========================#
#=======================================================#

# simu_name = init_directory(simu_name="simu")

# Du = 0.002
# Dv = 0.05
# U_init = sol[end]

# p = (Du,Dv,dx,dy,Nx,Ny,param_reac)    # Simulation parameters

# println("Initialize ODE problem...")

# dU_init = copy(U_init)
# jac_sparsity = Symbolics.jacobian_sparsity((du,u)->func!(du,u,p,0.0),dU_init,U_init)    # Automatic sparsity detection

# ode_func = ODEFunction(func!,jac_prototype=float.(jac_sparsity))
# prob = ODEProblem(ode_func,U_init,(0.,T),p)

# # prob = ODEProblem(func!,U_init,(0.,T),p)    # Without sparsity information

# println("Ok let's solve the ODE...")

# # algo = nothing    # Default choice
# # algo = Euler() 
# algo = ImplicitEuler()
# # algo = RadauIIA3()

# #------------------------------------------------------------------------------#


# #----------- Finallly run the simulation and save data ------------------------#

# sol = solve(prob,algo,alg_hints=[:stiff],saveat=10.0;progress=true,progress_steps=1)
# # sol = solve(prob,algo,saveat=10.0;dt=0.2,progress=true,progress_steps=1)    # For explicit Euler, specify dt

# println("Solving is done, now plotting...")

# plot_save_sol(
#     sol,U_init,u_index,v_index,dx,dy;
#     colorrange_u=(0.0,3.0),
#     colorrange_v=(0.0,3.0),
#     dir=simu_name,video_name=simu_name,
#     dt=10)

# println("Save data...")

# data = Dict(
#     "simu_name" => simu_name,
#     "reac_name" => reac_name,
#     "sol_final" => sol[end],
#     "T" => T,
#     "Du" => Du, 
#     "Dv" => Dv,
#     "dx" => dx,
#     "dy" => dy,
#     "Nx" => Nx,
#     "Ny" => Ny,
#     "param_reac" => param_reac
# )

# save(simu_name*"/"*simu_name*".jld2",data)

# println("All done.")

