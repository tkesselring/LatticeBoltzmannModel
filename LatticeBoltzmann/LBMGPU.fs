/// Simplistic 2D9 Lattice Boltzman Model implementation for GPU.
module LattieBoltzmann.GPU
////////////////////////////////////////////////////////////////////////////////
//
// Inspired by the "Crude 2D Lattice Boltzmann Demo program" by Graham Pullan, 
// see: http://www.many-core.group.cam.ac.uk/projects/LBdemo.shtml
//
//      f6  f2   f5
//        \  |  /
//         \ | /
//          \|/
//      f3---|--- f1
//          /|\
//         / | \       and f0 for the rest (zero) velocity
//        /  |  \
//      f7  f4   f8
//
///////////////////////////////////////////////////////////////////////////////

open Alea.CUDA
open Alea.CUDA.Utilities

[<Literal>]
let TILE_I = 16
[<Literal>]
let TILE_J = 8

[<ReflectedDefinition>]
let I2D ni i j = ni*j + i

[<AOTCompile(AOTOnly = true, SpecificArchs = "sm20;sm30;sm35")>]
type SimulatorModule(target, ni, nj, roout_, vxin_, tau_) =
    inherit GPUModule(target)
    let size = ni*nj
    let vxIn = vxin_
    let roOut = roout_
    let tau = tau_

    let faceq1 = 4.f / 9.f
    let faceq2 = 1.f / 9.f
    let faceq3 = 1.f / 36.f


    let F0 = Array.init size (fun _ -> faceq1*roOut*(1.f                             - 1.5f*vxIn*vxIn) )
    let F1 = Array.init size (fun _ -> faceq2*roOut*(1.f + 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
    let F2 = Array.init size (fun _ -> faceq2*roOut*(1.f                             - 1.5f*vxIn*vxIn) )
    let F3 = Array.init size (fun _ -> faceq2*roOut*(1.f - 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
    let F4 = Array.init size (fun _ -> faceq2*roOut*(1.f                             - 1.5f*vxIn*vxIn) )
    let F5 = Array.init size (fun _ -> faceq3*roOut*(1.f + 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
    let F6 = Array.init size (fun _ -> faceq3*roOut*(1.f - 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
    let F7 = Array.init size (fun _ -> faceq3*roOut*(1.f - 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
    let F8 = Array.init size (fun _ -> faceq3*roOut*(1.f + 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )


    let velocity = Array.zeroCreate size

    /// we do not have access to registers at the moment. For the first attempt read and write form global memory.
    // todo: only works for all in oune block by now.
    [<Kernel;ReflectedDefinition>]
    member this.StreamKernel (f1:deviceptr<float32>) (f2:deviceptr<float32>) (f3:deviceptr<float32>) (f4:deviceptr<float32>) (f5:deviceptr<float32>) (f6:deviceptr<float32>) (f7:deviceptr<float32>) (f8:deviceptr<float32>) =
        let i = blockIdx.x*TILE_I + threadIdx.x
        let j = blockIdx.y*TILE_J + threadIdx.y

        let i2d = I2D ni i j

        let sharedPos = __shared__.ExternArray<float32>()

        let i0 = i2D ni i j

        let jm1 = if j = 0 then 0
                  else j - 1
        let jp1 = if j = nj - 1 then nj - 1
                  else j + 1
    
        let im1 = if i = 0 then 0
                  else i - 1
        let ip1 = if i = ni - 1 then ni - 1
                  else i + 1

        sharedPos.[i0] <- f1.[i2D ni im1 j]

        __syncthreads()

        f1.[i0] <- sharedPos.[i0]

    member this.Stream() =
        let numBlocks = dim3(ni/TILE_I, nj/TILE_J)
        let blockSize = dim3(TILE_I, TILE_J)
        let sharedMemSize = ni*nj*__sizeof<float>()
        let lp = LaunchParam(numBlocks, blockSize, sharedMemSize)
        this.GPULaunch <@ this.StreamKernel @> lp f0 f1 f2 f3 f4 f5 f6 f7 f8


    [<Kernel;ReflectedDefinition>]
    member this.CollideKernel pitch tau faceq1 faceq2 faceq3 (f0:deviceptr<float32>) (f1:deviceptr<float32>) (f2:deviceptr<float32>) (f3:deviceptr<float32>) (f4:deviceptr<float32>) (f5:deviceptr<float32>) (f6:deviceptr<float32>) (f7:deviceptr<float32>) (f8:deviceptr<float32>) (plotData:deviceptr<float32>) =
        let i = blockIdx.x*TILE_I + threadIdx.x
        let j = blockIdx.y*TILE_J + threadIdx.y

        let i2d = i + j*pitch/__sizeof<float32>()

        let rtau = 1.f/tau;
        let rtau1 = 1.f - rtau;    

        // Read all f's and store in registers
        let f0now = f0.[i2d]
        let f1now = f1.[i2d]
        let f2now = f2.[i2d]
        let f3now = f3.[i2d]
        let f4now = f4.[i2d]
        let f5now = f5.[i2d]
        let f6now = f6.[i2d]
        let f7now = f7.[i2d]
        let f8now = f8.[i2d]

        // Macroscopic flow properties
        let ro =  f0now + f1now + f2now + f3now + f4now + f5now + f6now + f7now + f8now
        let vx = (f1now - f3now + f5now - f6now - f7now + f8now)/ro
        let vy = (f2now - f4now + f5now + f6now - f7now - f8now)/ro

        // Set plotting variable to velocity magnitude
        plotData.[i2d] <- sqrt (vx*vx + vy*vy)
    
        // Calculate equilibrium functions
        let v_sq_term = 1.5f*(vx*vx + vy*vy)
        let f0eq = ro*faceq1*(1.f - v_sq_term)
        let f1eq = ro*faceq2*(1.f + 3.f*vx + 4.5f*vx*vx - v_sq_term)
        let f2eq = ro*faceq2*(1.f + 3.f*vy + 4.5f*vy*vy - v_sq_term)
        let f3eq = ro*faceq2*(1.f - 3.f*vx + 4.5f*vx*vx - v_sq_term)
        let f4eq = ro*faceq2*(1.f - 3.f*vy + 4.5f*vy*vy - v_sq_term)
        let f5eq = ro*faceq3*(1.f + 3.f*( vx + vy) + 4.5f*( vx + vy)*( vx + vy) - v_sq_term)
        let f6eq = ro*faceq3*(1.f + 3.f*(-vx + vy) + 4.5f*(-vx + vy)*(-vx + vy) - v_sq_term)
        let f7eq = ro*faceq3*(1.f + 3.f*(-vx - vy) + 4.5f*(-vx - vy)*(-vx - vy) - v_sq_term)
        let f8eq = ro*faceq3*(1.f + 3.f*( vx - vy) + 4.5f*( vx - vy)*( vx - vy) - v_sq_term)

        // collide
        f0.[i2d] <- rtau1*f0now + rtau*f0eq
        f1.[i2d] <- rtau1*f1now + rtau*f1eq
        f2.[i2d] <- rtau1*f2now + rtau*f2eq
        f3.[i2d] <- rtau1*f3now + rtau*f3eq
        f4.[i2d] <- rtau1*f4now + rtau*f4eq
        f5.[i2d] <- rtau1*f5now + rtau*f5eq
        f6.[i2d] <- rtau1*f6now + rtau*f6eq
        f7.[i2d] <- rtau1*f7now + rtau*f7eq
        f8.[i2d] <- rtau1*f8now + rtau*f8eq

    member this.Collide() =
        let blockSize = dim3(TILE_I, TILE_J)
        let numBlocks = dim3(ni/TILE_I, nj/TILE_J)
        let lp = LaunchParam(numBlocks, blockSize)
        this.GPULaunch <@ this.CollideKernel @> lp tau faceq1 faceq2 faceq3 f0 f1 f2 f3 f4 f5 f6 f7 f8 plotData

// todo nj i snot used, mistake?
    [<Kernel;ReflectedDefinition>]
    member this.ApplyBCskernel ni (nj:int) vxin roout faceq2 faceq3 (f0:deviceptr<float32>) (f1:deviceptr<float32>) (f2:deviceptr<float32>) (f3:deviceptr<float32>) (f4:deviceptr<float32>) (f5:deviceptr<float32>,  f6:deviceptr<float32>) (f7:deviceptr<float32>) (f8:deviceptr<float32>) (solid_data:deviceptr<int>) =
        let i = blockIdx.x*TILE_I + threadIdx.x
        let j = blockIdx.y*TILE_J + threadIdx.y
    
        let i2d = i + j*pitch/__sizeof<float>()
    
        // Solid BC: "bounce-back"
        if solid_data.[i2d] = 0 then
            let f1old = f1.[i2d]
            let f2old = f2.[i2d]
            let f3old = f3.[i2d]
            let f4old = f4.[i2d]
            let f5old = f5.[i2d]
            let f6old = f6.[i2d]
            let f7old = f7.[i2d]
            let f8old = f8.[i2d]
        
            f1.[i2d] <- f3old
            f2.[i2d] <- f4old
            f3.[i2d] <- f1old
            f4.[i2d] <- f2old
            f5.[i2d] <- f7old
            f6.[i2d] <- f8old
            f7.[i2d] <- f5old
            f8.[i2d] <- f6old
    
        // Inlet BC
        if i = 0 then
            let v_sq_term = 1.5f*(vxin * vxin)
        
            f1.[i2d] <- roout*faceq2*(1.f + 3.f*vxin + 3.f*v_sq_term)
            f5.[i2d] <- roout*faceq3*(1.f + 3.f*vxin + 3.f*v_sq_term)
            f8.[i2d] <- roout*faceq3*(1.f + 3.f*vxin + 3.f*v_sq_term)

        // Exit BC
        if i = (ni-1) then
            let i2d2 = i2d - 1
            f3.[i2d] <- f3.[i2d2]
            f6.[i2d] <- f6.[i2d2]
            f7.[i2d] <- f7.[i2d2]

    member this.ApplyBCs() =
        let numBlocks = dim3(ni/TILE_I, nj/TILE_J)
        let blockSize = dim3(TILE_I, TILE_J)
        let lp = LaunchParam(numBlocks, blockSize)
        this.GPULaunch <@ this.ApplyBCskernel @> lp ni nj vxin roout faceq2 faceq3 f0 f1 f2 f3 f4 f5 f6 f7 f8 solid_data

    [<Kernel;ReflectedDefinition>]
    member this.ApplyPeriodicBCkernel ni nj (f2:deviceptr<float32>) (f4:deviceptr<float32>) (f5:deviceptr<float32>) (f6:deviceptr<float32>) (f7:deviceptr<float32>) (f8:deviceptr<float32>) =
        let i = blockIdx.x*TILE_I + threadIdx.x
        let j = blockIdx.y*TILE_J + threadIdx.y

        let i2d = i + j/__sizeof<float32>()

        if j = 0 then
            let i2d2 = i + (nj-1)/__sizeof<float32>()
            f2.[i2d] <- f2.[i2d2]
            f5.[i2d] <- f5.[i2d2]
            f6.[i2d] <- f6.[i2d2]
        if j = (nj-1) then
            let i2d2 = i
            f4.[i2d] <- f4.[i2d2]
            f7.[i2d] <- f7.[i2d2]
            f8.[i2d] <- f8.[i2d2]

    member this.ApplyPeriodicBC() =
        let numBlocks = dim3(ni/TILE_I, nj/TILE_J)
        let blockSize = dim3(TILE_I, TILE_J)
        let lp = LaunchParam(numBlocks, blockSize)
        this.GPULaunch <@ this.ApplyPeriodicBCkernel @> lp ni nj f2 f4 f5 f6 f7 f8

    [<Kernel;ReflectedDefinition>]
    member this.GetRgbaKernel ncol minvar maxvar plotData (plot_rgba_data:deviceptr<float32>) (cmap_rgba_data:deviceptr<float32>) (solid_data:deviceptr<float32>) =
        let i = blockIdx.x*TILE_I + threadIdx.x
        let j = blockIdx.y*TILE_J + threadIdx.y

        let i2d = i + j/__sizeof<float>()

        let frac = (plotData[i2d]-minvar)/(maxvar-minvar)
        let icol = (int)(frac * (float)ncol)
        plot_rgba_data.[i2d] <- solid_data.[i2d] * cmap_rgba_data.[icol]

    member this.GetRgba() =
        let numBlocks = dim3(ni/TILE_I, nj/TILE_J)
        let blockSize = dim3(TILE_I, TILE_J)
        let lp = LaunchParam(numBlocks, blockSize)
        this.GPULaunch <@ this.GetRgbaKernel @> lp ncol minvar maxvar plotData plot_rgba_data cmap_rgba_data solid_data