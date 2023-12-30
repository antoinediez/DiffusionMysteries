# using CairoMakie
using GLMakie    #Faster than CairoMakie
using ProgressMeter

function init_plot(
    u,v;
    dx,dy,
    save_video=false,
    time=0.0,
    colorrange_u::Union{Tuple,Nothing}=nothing,
    colorrange_v::Union{Tuple,Nothing}=nothing,
    fig=nothing,
    bigfig=nothing,
)
    round_time = nice_float2string(time,2)
    Nx = size(u)[1]
    Ny = size(u)[2]
    Lx = dx * Nx
    Ly = dy * Ny    

    xx = dx .* collect(0:(Nx-1)) .+ dx/2
    yy = dy .* collect(0:(Ny-1)) .+ dy/2

    if isnothing(fig)
        fig = Figure(resolution=(1600,800))
    end

    if isnothing(bigfig)
        bigfig = fig
    end

    ax_u = Axis(fig[1, 1], title="u\ntime=$round_time")
    ax_u.aspect = AxisAspect(1)
    xlims!(ax_u,0,Lx)
    ylims!(ax_u,0,Ly)

    if isnothing(colorrange_u)
        @views hm_u = heatmap!(ax_u,xx,yy,u)
    else
        @views hm_u = heatmap!(ax_u,xx,yy,u,colorrange=colorrange_u)
    end
    Colorbar(fig[1,2],hm_u)

    ax_v = Axis(fig[1, 3], title="v\ntime=$round_time")
    ax_v.aspect = AxisAspect(1)
    xlims!(ax_v,0,Lx)
    ylims!(ax_v,0,Ly)
    if isnothing(colorrange_v)
        @views hm_v = heatmap!(ax_v,xx,yy,v)
    else
        @views hm_v = heatmap!(ax_v,xx,yy,v,colorrange=colorrange_v)
    end
    Colorbar(fig[1,4],hm_v)

    if save_video
        stream = VideoStream(bigfig,framerate=40)
        recordframe!(stream)
        return fig, ax_u, ax_v, hm_u, hm_v, stream
    else
        return fig, ax_u, ax_v, hm_u, hm_v, nothing
    end
end


function update_plot!(
    u,v;
    time,
    fig, ax_u, ax_v, hm_u, hm_v,
    stream=nothing,
)    
    round_time = nice_float2string(time,2)

    hm_u[3] = u
    ax_u.title = "u\ntime=$round_time"

    if !isnothing(hm_v)
        hm_v[3] = v
        ax_v.title = "v\ntime=$round_time"
    end

    # Specify the color range at the end for the simulation with decreasing diffusion (remove otherwise...).
    if time>1999.9
        hm_u.colorrange = (0.5203,3.49666)
        hm_v.colorrange = (0.31513,0.7015)
    end

    if !isnothing(stream)
        recordframe!(stream)
    end
end


function nice_float2string(x,K::Int)
    integer_part = trunc(Int,x)
    y = x - integer_part
    y10K = trunc(Int,y*10^K)
    if y10K >= 10
        decimal_part = rpad(y10K,K,"0")
    else
        decimal_part = lpad(y10K,K,"0")
    end
    return "$(integer_part).$decimal_part"
end


function plot_save_sol_4(
    sol,init,u_index,v_index,dx,dy;
    dt=0.1,
    colorrange_u=nothing,colorrange_v=nothing,colorrange_c=nothing,
    dir=pwd(),video_name="funny_video")
    
    println("Init plot...")

    bigfig = Figure(resolution=(1600,800))
    u = Vector{Any}(undef,4)
    v = Vector{Any}(undef,4)
    fig = Vector{Any}(undef,4)
    ax_u = Vector{Any}(undef,4)
    ax_v = Vector{Any}(undef,4)
    hm_u = Vector{Any}(undef,4)
    hm_v = Vector{Any}(undef,4)
    stream = nothing

    for s in 1:4

        if s==1
            i,j=1,1
            save_video=false
        elseif s==2
            i,j=1,2
            save_video=false
        elseif s==3
            i,j=2,1
            save_video=false
        elseif s==4
            i,j=2,2
            save_video=true
        end

        u[s] = @view init[s][u_index...]
        v[s] = @view init[s][v_index...]

        fig[s], ax_u[s], ax_v[s], hm_u[s], hm_v[s], stream = init_plot(
            u[s],v[s];
            dx,dy,
            save_video=save_video,
            time=0.0,
            colorrange_u=colorrange_u,
            colorrange_v=colorrange_v,
            fig=bigfig[i,j],bigfig=bigfig
        )
    end

    T = sol[4].t[end]
    K = floor(Int,T/dt)

    println("Plotting...")

    @showprogress for k in 1:K
        time = k*dt
        for s in 1:4
            if s==4
                stm = stream 
            else
                stm = nothing
            end
            U = sol[s](k*dt)
            u[s] = @view U[u_index...]
            v[s] = @view U[v_index...]
            update_plot!(
            u[s],v[s];
            time=time,
            fig=fig[s], ax_u=ax_u[s], ax_v=ax_v[s], hm_u=hm_u[s], hm_v=hm_v[s],
            stream=stm,
            )
        end
    end

    println("Save plot...")
    
    save(dir*"/"*video_name*".mp4", stream)
end


function plot_save_sol(
    sol,init,u_index,v_index,dx,dy;
    dt=0.1,
    colorrange_u=nothing,colorrange_v=nothing,plot_3D=false,
    view_param_u=nothing,view_param_v=nothing,plot_v=true,fig=nothing,bigfig=nothing,
    plot_arrow=true,diffusion=1.0,diffusion_max=10.0,
    diffusion_time=nothing,diffusion_value=nothing,
    dir=pwd(),video_name="funny_video")
    
    println("Init plot...")

    u = @view init[u_index...]
    v = @view init[v_index...]

    if plot_3D
        fig, ax_u, ax_v, hm_u, hm_v, ar_fig, stream = init_plot_3D(
            u,v;
            dx,dy,
            save_video=true,
            time=0.0,
            view_param_u=view_param_u,
            view_param_v=view_param_v,
            plot_v=plot_v,
            plot_arrow=plot_arrow,
            diffusion=diffusion,diffusion_max=diffusion_max,
            colorrange_u=colorrange_u,
            colorrange_v=colorrange_v,
            fig=fig,
            bigfig=bigfig,
        )
    else
        fig, ax_u, ax_v, hm_u, hm_v, stream = init_plot(
            u,v;
            dx,dy,
            save_video=true,
            time=0.0,
            colorrange_u=colorrange_u,
            colorrange_v=colorrange_v,
        )
    end

    T = sol.t[end]
    K = floor(Int,T/dt)
    diff_plot = 1

    println("Plotting...")

    @showprogress for k in 1:K

        if plot_3D
            if !isnothing(ar_fig)
                if diff_plot<=length(diffusion_time)
                    if k*dt>=diffusion_time[diff_plot]
                        x0 = ar_fig[4][1].val[1][1]
                        ar_fig[4][1].val = [[x0,diffusion_value[diff_plot]]]
                        notify(ar_fig[4][1])
                        diff_plot += 1
                    end
                end
            end
        end

        U = sol(k*dt)
        time = k*dt
        u = @view U[u_index...]
        v = @view U[v_index...]
        update_plot!(
        u,v;
        time,
        fig, ax_u, ax_v, hm_u, hm_v,
        stream=stream,
        )

    end

    println("Save plot...")
    
    save(dir*"/"*video_name*".mp4", stream)
end

function init_directory(;simu_name="simu")
    if ispath(simu_name)
        k=1
        while ispath(simu_name*"_$k")
            k+=1
        end
        mkdir(simu_name*"_$k")
        return simu_name*"_$k"
    else
        mkdir(simu_name)
        return simu_name
    end
end

function init_plot_3D(
    u,v;
    dx,dy,
    save_video=false,
    time=0.0,
    view_param_u=nothing,
    view_param_v=nothing,
    plot_v=false,
    plot_arrow=true,
    diffusion=1.0,
    diffusion_max=100.0,
    colorrange_u::Union{Tuple,Nothing}=nothing,
    colorrange_v::Union{Tuple,Nothing}=nothing,
    fig=nothing,
    bigfig=nothing,
)

    round_time = nice_float2string(time,2)
    Nx = size(u)[1]
    Ny = size(u)[2]
    Lx = dx * Nx
    Ly = dy * Ny    

    xx = dx .* collect(0:(Nx-1)) .+ dx/2
    yy = dy .* collect(0:(Ny-1)) .+ dy/2

    if isnothing(fig)
        if plot_v
            if plot_arrow
                fig = Figure(resolution=(1800,600))
            else
                fig = Figure(resolution=(1600,800))
            end
        else
            if plot_arrow
                fig = Figure(resolution=(1000,800))
            else
                fig = Figure(resolution=(800,800))
            end
        end
    end

    if isnothing(bigfig)
        bigfig = fig
    end

    if plot_arrow
        ar_fig = plot_arrow_point(
            diffusion;ymax=diffusion_max,ymin=0,
            txt="diffusion",
            markersize=54,
            fig=fig[1,1],bigfig=bigfig)
        colsize!(fig.layout, 1, Auto(0.25))
    else
        ar_fig = nothing
    end

    if isnothing(view_param_u)
        ax_u = Axis3(fig[1, 2], title="u\ntime=$round_time")
    else
        ax_u = Axis3(fig[1, 2]; title="u\ntime=$round_time", view_param_u...)
    end

    if isnothing(colorrange_u)
        sf_u = surface!(ax_u,xx,yy,u)
    else
        sf_u = surface!(ax_u,xx,yy,u,colorrange=colorrange_u)
    end

    if plot_v
        if isnothing(view_param_v)
            ax_v = Axis3(fig[1, 3], title="v\ntime=$round_time")
        else
            ax_v = Axis3(fig[1, 3]; title="v\ntime=$round_time", view_param_v...)
        end
        if isnothing(colorrange_v)
            sf_v = surface!(ax_v,xx,yy,v)
        else
            sf_v = surface!(ax_v,xx,yy,v,colorrange=colorrange_v)
        end
    else
        sf_v = nothing
        ax_v = nothing
    end

    if save_video
        stream = VideoStream(bigfig,framerate=40)
        recordframe!(stream)
        return fig, ax_u, ax_v, sf_u, sf_v, ar_fig, stream
    else
        return fig, ax_u, ax_v, sf_u, sf_v, ar_fig, nothing
    end

end

function plot_arrow_point(
    y;ymax=nothing,ymin=0,
    txt=nothing,
    markersize=54,
    fig=nothing,bigfig=nothing)

    if isnothing(fig)
        fig = Figure(resolution=(200,800))
    end

    if isnothing(bigfig)
        bigfig = fig
    end

    if isnothing(ymax)
        ymax = 2*y
    end

    ax = Axis(
        fig[1,1],
        limits=(0,ymax-ymin,ymin,1.1*ymax),
        xlabelvisible=false,xticksvisible=false,xticklabelsvisible=false,xgridvisible=false,
        ylabelvisible=false,yticksvisible=false,yticklabelsvisible=false,ygridvisible=false,
        leftspinevisible=false,rightspinevisible=false,bottomspinevisible=false,topspinevisible=false
        )
    x0 = (ymax-ymin)/2.0
    ar = arrows!(ax,[x0],[ymin],[0.0],[ymax-ymin],linewidth=8,arrowsize=36)
    if !isnothing(txt)
        tx = text!(ax,x0,1.07*ymax,text=txt,align=(:center,:center),fontsize=42)
    else
        tx = nothing
    end

    pt = scatter!(ax,[x0],[y],markersize=markersize)
    return fig, ax, ar, pt, tx
end