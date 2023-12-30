function reaction(u,v,param)
    #==== Linear KS (with inhibitor) ====#
    # r_u = param[:sc]*(c - param[:α]*v - u)
    # r_v = param[:sc]*(c - param[:β]*v)
    #====================================#

    #==== Custom Activator-Inhibitor ====#
    # r_u = param[:sc]*(param[:α]*u/(1+param[:β]*v)*c^4/(1+c^4) - param[:γ]*u)
    # r_v = param[:sc]*(c - param[:δ]*v)
    #====================================#

    #== Activator-Inhibitor (no cell) ==#
    # act_synth = max(param[:min_act_synth],min(param[:act_synth]+param[:act_autocat]*u - param[:inhib_act]*v,param[:max_act_synth]))
    # inhib_synth = max(param[:min_inhib_synth],min(param[:inhib_synth]+param[:inhib_autocat]*v + param[:inhib_act_crea]*u,param[:max_inhib_synth]))
    # r_u = param[:sc]*(act_synth - param[:act_deg]*u)
    # r_v = param[:sc]*(inhib_synth - param[:inhib_deg]*v)
    #====================================#

    #== Activator-Inhibitor (u->c) ==#
    # act_synth = max(0.0,min(param[:act_synth]+param[:act_autocat] - param[:inhib_act]*v,param[:max_act_synth]))
    # inhib_synth = max(0.0,min(param[:inhib_synth]+param[:inhib_autocat]*v + param[:inhib_act_crea],param[:max_inhib_synth]))
    # r_u = param[:sc]*(act_synth - param[:act_deg]*u)
    # r_v = param[:sc]*(inhib_synth - param[:inhib_deg]*v)
    #====================================#

    #== Schnakenberg ====================#
    r_u = param[:sc]*(param[:a] + u^2*v - u)
    r_v = param[:sc]*(param[:b] - u^2*v)
    #====================================#

    #== Fish ============================#

    # act_synth = max(0.0, min(param[:a]*u + param[:b]*v + param[:c],param[:synUmax]))
    # inhib_synth = max(0.0, min(param[:e]*u + param[:g],param[:synVmax]))
    # r_u = act_synth - param[:d]*u
    # r_v = inhib_synth - param[:f]*v
    #====================================#




    #=========== No reaction ============#
    # r_u = 0.0
    # r_v = 0.0
    #====================================#
    return (r_u,r_v)
end

