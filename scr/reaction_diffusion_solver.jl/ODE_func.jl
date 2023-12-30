
#--------------------------------------------------#
#                     ODE function                 #
#--------------------------------------------------#

function func!(dU,U,p,t)

    Du,Dv,dx,dy,Nx,Ny,param_reac,bc = p 



    ############################ BULK ###############################

    @inbounds for j in 1:Ny
        @simd for i in 1:Nx

            if bc==0    #Neumann
                ip1 = min(i+1,Nx)
                im1 = max(i-1,1)
                jp1 = min(j+1,Ny)
                jm1 = max(j-1,1)
            elseif bc==1    #Periodic
                ip1 = i<Nx ? i+1 : 1
                im1 = i>1 ? i-1 : Nx
                jp1 = j<Ny ? j+1 : 1
                jm1 = j>1 ? j-1 : Nx
            end
            
            ########## Laplace #############################
            Δu = (U[ip1,j,1] + U[im1,j,1] - 2*U[i,j,1])/dx^2 + (U[i,jp1,1] + U[i,jm1,1] - 2*U[i,j,1])/dy^2
            Δv = (U[ip1,j,2] + U[im1,j,2] - 2*U[i,j,2])/dx^2 + (U[i,jp1,2] + U[i,jm1,2] - 2*U[i,j,2])/dy^2
            ################################################

            r_u,r_v = reaction(U[i,j,1],U[i,j,2],param_reac)
            dU[i,j,1] = Du * Δu + r_u
            dU[i,j,2] = Dv * Δv + r_v
        end
    end

end

