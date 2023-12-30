using Logging: global_logger
using TerminalLoggers: TerminalLogger
global_logger(TerminalLogger())

using Symbolics, IfElse
using DifferentialEquations
using JLD2
using BenchmarkTools
using StaticArrays

println("Initialize everything...")

#-------- Add the scripts------------------------------------------------------#
include("scr/reaction_diffusion_solver.jl/reaction_taxis_growth.jl")
include("scr/reaction_diffusion_solver.jl/plot_save.jl")
include("scr/reaction_diffusion_solver.jl/ODE_func.jl")
#------------------------------------------------------------------------------#

simu_name = init_directory(simu_name="plot_RD")

#-------- Discretization parameters -------------------------------------------#
Nx = 512    # Number of cells in the x-direction (bulk and surface)
Ny = 512    # Number of cells in the y-direction (bulk)
Lx = 100.0     # Length of the domain in the x-direction (bulk and surface)
Ly = 100.0   # Length of the domain in the y-direction (bulk)
dx = Lx/Nx    # Cell length in the x-direction
dy = Ly/Ny    # Cell length in the y-direction

#------------------------------------------------------------------------------#


#-------- Time parameters -----------------------------------------------------#

T = 2500.0 # Simulation time 

#------------------------------------------------------------------------------#


#------------------- Model  ---------------------------------------------------#


#===================== REACTION ===========================#

# Choose from one of the models in `reaction_taxis_growth.jl`. 
# The parameters can be defined in a dictionary (slow) or in a named tuple (much faster)

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

# reac_name = "Activator-Inhibitor"
# param_reac = (
#     act_synth = 0.1,
#     act_autocat = 0.08,
#     inhib_act = 0.08,
#     max_act_synth = 0.2,
#     min_act_synth = 0.0,
#     act_deg = 0.03,
#     inhib_synth = -0.15,
#     inhib_autocat = 0.0,
#     inhib_act_crea = 0.11,
#     max_inhib_synth = 0.5,
#     inhib_deg = 0.06,
#     min_inhib_synth = 0.0,
#     sc = 1.0,
# )

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

reac_name = "Schnakenberg"

param_reac = (
    a = 0.01,
    b = 2.0,
    sc = 1.0
)

Du = 1.0
Dv = 130.0
u_eq = param_reac[:a] + param_reac[:b]
v_eq = param_reac[:b] / (param_reac[:a] + param_reac[:b])^2


#=======================================================#

#------------------- Initial condition  ---------------------------------------#

# Basic purely random initialization 
U_init = 0.0 .+ 3.0 .* rand(Nx,Ny,2)

# The individual solutions can be recovered with the following index table   
u_index = :,:,1 
v_index = :,:,2

# Small perturbation of the equilibrium
# U_init[u_index...] = u_eq * (ones(Nx,Ny) + 0.1 * rand(Nx,Ny))
# U_init[v_index...] = v_eq * (ones(Nx,Ny) + 0.1 * rand(Nx,Ny))

# Big spot in the middle of the domain
U_init[u_index...] = u_eq * ones(Nx,Ny)
U_init[v_index...] = v_eq * ones(Nx,Ny)
G = 1500
# G = 15.0    # Smaller amplitude
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

bc = 1    # Boundary condition: 0 for Neumann and 1 for periodic

p = (Du,Dv,dx,dy,Nx,Ny,param_reac,bc)    # Simulation parameters

println("Initialize ODE problem...")

dU_init = copy(U_init)
jac_sparsity = Symbolics.jacobian_sparsity((du,u)->func!(du,u,p,0.0),dU_init,U_init)    # Automatic sparsity detection

ode_func = ODEFunction(func!,jac_prototype=float.(jac_sparsity))
prob = ODEProblem(ode_func,U_init,(0.,T),p)

# prob = ODEProblem(func!,U_init,(0.,T),p)    # Without sparsity information (very slow)

println("Ok let's solve the ODE...")

# algo = nothing    # Default choice
# algo = Euler() 
# algo = ImplicitEuler()
# algo = RadauIIA3()
# algo = FBDF()
algo = ROCK2()   # Seems to be the best ans fastest one

#------------------------------------------------------------------------------#

# Events can be speficied as callbacks: here the diffusion is reduced at the times given in `dosetimes`. 

dosetimes = collect(200.0:100.0:1100.0)
append!(dosetimes,2000.0)
condition(u, t, integrator) = t ∈ dosetimes
# affect!(integrator) = integrator.p = (Du,Dv-1.0*(1.0+(integrator.t - 200.0)/100.0)^2,dx,dy,Nx,Ny,param_reac,bc)
affect!(integrator) = integrator.p = integrator.t<1999.0 ? (Du,Dv-1.0*(1.0+(integrator.t - 200.0)/100.0)^2,dx,dy,Nx,Ny,param_reac,bc) : (Du,20.0,dx,dy,Nx,Ny,param_reac,bc)
cb = DiscreteCallback(condition, affect!)

#----------- Finallly run the simulation and save data ------------------------#

saveat = 1.0
# sol = solve(prob,algo,alg_hints=[:stiff],saveat=saveat,dtmin=0.001,force_dtmin=true;progress=true,progress_steps=1)
sol = solve(prob,algo,alg_hints=[:stiff],saveat=saveat,dtmin=0.00001,force_dtmin=true;progress=true,progress_steps=1,callback=cb,tstops=dosetimes)

println("Solving is done, now plotting...")

view_param_u = (
    limits=(0,Lx,0,Ly,0,11),
    aspect=(1,1,0.2),
    elevation=0.2pi,
    xlabelvisible=false,xticklabelsvisible=false,
    ylabelvisible=false,yticklabelsvisible=false,
    zlabelvisible=false,zticklabelsvisible=false
    )

view_param_v = (
    limits=(0,Lx,0,Ly,0,11),
    aspect=(1,1,0.2),
    elevation=0.2pi,
    xlabelvisible=false,xticklabelsvisible=false,
    ylabelvisible=false,yticklabelsvisible=false,
    zlabelvisible=false,zticklabelsvisible=false
    )


plot_save_sol(
    sol,U_init,u_index,v_index,dx,dy;
    plot_3D=true,
    plot_v=true,
    plot_arrow=true,
    diffusion=Dv,
    diffusion_max=150.0,
    view_param_u=view_param_u,
    view_param_v=view_param_v,
    diffusion_time=dosetimes,diffusion_value=max.(Dv.-1.0.*(1.0.+(dosetimes .- 200.0)./100.0).^2,20.0),
    dir=simu_name,video_name=simu_name,
    dt=saveat)

# plot_save_sol(
#     sol,U_init,u_index,v_index,dx,dy;
#     colorrange_u=(umin,umax),
#     colorrange_v=(vmin,vmax),
#     dir=simu_name,video_name=simu_name,diffusion_value=0.0,
#     dt=saveat)

println("Save data...")

data = Dict(
    "simu_name" => simu_name,
    "reac_name" => reac_name,
    "sol_final" => sol[end],
    "T" => T,
    "Du" => Du, 
    "Dv" => Dv,
    "dx" => dx,
    "dy" => dy,
    "Nx" => Nx,
    "Ny" => Ny,
    "param_reac" => param_reac
)

save(simu_name*"/"*simu_name*".jld2",data)

println("All done.")
