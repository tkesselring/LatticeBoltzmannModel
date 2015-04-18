/// Simplistic 2D9 Lattice Boltzman Model implementation for CPU.
module LattieBoltzmann
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
let i2D ni i j = ni*j + i

type LBMLattice = 
    {
        Ni : int
        Nj : int
        F0 : float32[]
        F1 : float32[]
        F2 : float32[]
        F3 : float32[]
        F4 : float32[]
        F5 : float32[]
        F6 : float32[]
        F7 : float32[]
        F8 : float32[]
    }

type LatticeBoltzmannModel(ni, nj, roout_, vxin_, tau_) =
    let size = ni*nj
    let vxIn = vxin_
    let roOut = roout_
    let tau = tau_

    let faceq1 = 4.f / 9.f
    let faceq2 = 1.f / 9.f
    let faceq3 = 1.f / 36.f

    let lattice = 
        {
            Ni = ni
            Nj = nj
            F0 = Array.init size (fun _ -> faceq1*roOut*(1.f                             - 1.5f*vxIn*vxIn) )
            F1 = Array.init size (fun _ -> faceq2*roOut*(1.f + 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
            F2 = Array.init size (fun _ -> faceq2*roOut*(1.f                             - 1.5f*vxIn*vxIn) )
            F3 = Array.init size (fun _ -> faceq2*roOut*(1.f - 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
            F4 = Array.init size (fun _ -> faceq2*roOut*(1.f                             - 1.5f*vxIn*vxIn) )
            F5 = Array.init size (fun _ -> faceq3*roOut*(1.f + 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
            F6 = Array.init size (fun _ -> faceq3*roOut*(1.f - 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
            F7 = Array.init size (fun _ -> faceq3*roOut*(1.f - 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
            F8 = Array.init size (fun _ -> faceq3*roOut*(1.f + 3.f*vxIn + 4.5f*vxIn*vxIn - 1.5f*vxIn*vxIn) )
        }
    let velocity = Array.zeroCreate size

    /// Let the particel densities flow.
    member this.stream () =
        let mutable jm1 = 0
        let mutable jp1 = 0
        let mutable im1 = 0
        let mutable ip1 = 0
        let ni, nj = lattice.Ni, lattice.Nj
        let size = ni*nj

        // create separate copies
        let (tmpf1 : float32[]) = Array.zeroCreate size
        let (tmpf2 : float32[]) = Array.zeroCreate size
        let (tmpf3 : float32[]) = Array.zeroCreate size
        let (tmpf4 : float32[]) = Array.zeroCreate size
        let (tmpf5 : float32[]) = Array.zeroCreate size
        let (tmpf6 : float32[]) = Array.zeroCreate size
        let (tmpf7 : float32[]) = Array.zeroCreate size
        let (tmpf8 : float32[]) = Array.zeroCreate size

        for j = 0 to nj - 1 do
            jm1 <- j - 1
            jp1 <- j + 1
            if j = 0 then jm1 <- 0
            if j = nj - 1 then  jp1 <- nj - 1
            for i = 1 to ni - 1 do
                let i0  = i2D ni i j 
                im1 <- i - 1
                ip1 <- i + 1
                if i = 0 then im1 <- 0
                if i = ni - 1 then ip1 <- ni - 1
                tmpf1.[i0] <- lattice.F1.[i2D ni im1 j]
                tmpf2.[i0] <- lattice.F2.[i2D ni i   jm1]
                tmpf3.[i0] <- lattice.F3.[i2D ni ip1 j]
                tmpf4.[i0] <- lattice.F4.[i2D ni i   jp1]
                tmpf5.[i0] <- lattice.F5.[i2D ni im1 jm1]
                tmpf6.[i0] <- lattice.F6.[i2D ni ip1 jm1]
                tmpf7.[i0] <- lattice.F7.[i2D ni ip1 jp1]
                tmpf8.[i0] <- lattice.F8.[i2D ni im1 jp1]
        
        // copy temporary arrays to lattice.
        for j = 0 to nj - 1 do
            for i = 1 to ni - 1 do
                let i0 = i2D ni i j
                lattice.F1.[i0] <- tmpf1.[i0]
                lattice.F2.[i0] <- tmpf2.[i0]
                lattice.F3.[i0] <- tmpf3.[i0]
                lattice.F4.[i0] <- tmpf4.[i0]
                lattice.F5.[i0] <- tmpf5.[i0]
                lattice.F6.[i0] <- tmpf6.[i0]
                lattice.F7.[i0] <- tmpf7.[i0]
                lattice.F8.[i0] <- tmpf8.[i0]

    /// Lattice Boltzmann collisions for 2D9 Lattice.
    member this.collide () =
        for j = 0 to lattice.Nj - 1 do
            for i = 0 to lattice.Ni - 1 do

                let i0 = i2D lattice.Ni i j
            
                let ro = lattice.F0.[i0] + lattice.F1.[i0] + lattice.F2.[i0] + lattice.F3.[i0] + 
                         lattice.F4.[i0] + lattice.F5.[i0] + lattice.F6.[i0] + lattice.F7.[i0] + lattice.F8.[i0]
                let rovx = lattice.F1.[i0] - lattice.F3.[i0] + lattice.F5.[i0] - lattice.F6.[i0] - lattice.F7.[i0] + lattice.F8.[i0]
                let rovy = lattice.F2.[i0] - lattice.F4.[i0] + lattice.F5.[i0] + lattice.F6.[i0] - lattice.F7.[i0] - lattice.F8.[i0]
                let vx = rovx / ro
                let vy = rovy / ro
            
                velocity.[i0] <- sqrt (vx*vx + vy*vy)

                let v_sq_term = 1.5f * (vx*vx + vy*vy)
            
                // Calculate equilibrium functions:
                let f0eq = ro*faceq1*(1.f - v_sq_term)
                let f1eq = ro*faceq2*(1.f + 3.f*vx + 4.5f*vx*vx - v_sq_term)
                let f2eq = ro*faceq2*(1.f + 3.f*vy + 4.5f*vy*vy - v_sq_term)
                let f3eq = ro*faceq2*(1.f - 3.f*vx + 4.5f*vx*vx - v_sq_term)
                let f4eq = ro*faceq2*(1.f - 3.f*vy + 4.5f*vy*vy - v_sq_term)
                let f5eq = ro*faceq3*(1.f + 3.f*( vx + vy) + 4.5f*( vx + vy)*( vx + vy) - v_sq_term)
                let f6eq = ro*faceq3*(1.f + 3.f*(-vx + vy) + 4.5f*(-vx + vy)*(-vx + vy) - v_sq_term)
                let f7eq = ro*faceq3*(1.f + 3.f*(-vx - vy) + 4.5f*(-vx - vy)*(-vx - vy) - v_sq_term)
                let f8eq = ro*faceq3*(1.f + 3.f*( vx - vy) + 4.5f*( vx - vy)*( vx - vy) - v_sq_term)

                let rtau = 1.0f / tau
                let rtau1 = 1.0f - rtau
                // Simulate collisions by "relaxing" toward the local equilibrium
                lattice.F0.[i0] <- rtau1*lattice.F0.[i0] + rtau*f0eq
                lattice.F1.[i0] <- rtau1*lattice.F1.[i0] + rtau*f1eq
                lattice.F2.[i0] <- rtau1*lattice.F2.[i0] + rtau*f2eq
                lattice.F3.[i0] <- rtau1*lattice.F3.[i0] + rtau*f3eq
                lattice.F4.[i0] <- rtau1*lattice.F4.[i0] + rtau*f4eq
                lattice.F5.[i0] <- rtau1*lattice.F5.[i0] + rtau*f5eq
                lattice.F6.[i0] <- rtau1*lattice.F6.[i0] + rtau*f6eq
                lattice.F7.[i0] <- rtau1*lattice.F7.[i0] + rtau*f7eq
                lattice.F8.[i0] <- rtau1*lattice.F8.[i0] + rtau*f8eq

    /// Apply bounce back boundry conditions if the variable `solid`is 0.0
    member this.solid_BC (solid:float32[]) =
        let ni, nj = lattice.Ni, lattice.Nj
        for j = 0 to nj - 1 do
            for i = 0 to ni - 1 do
                let i0 = i2D ni i j
                if solid.[i0] = 0.0f then
                    let f1old = lattice.F1.[i0]
                    let f2old = lattice.F2.[i0]
                    let f3old = lattice.F3.[i0]
                    let f4old = lattice.F4.[i0]
                    let f5old = lattice.F5.[i0]
                    let f6old = lattice.F6.[i0]
                    let f7old = lattice.F7.[i0]
                    let f8old = lattice.F8.[i0]
                
                    lattice.F1.[i0] <- f3old
                    lattice.F2.[i0] <- f4old
                    lattice.F3.[i0] <- f1old
                    lattice.F4.[i0] <- f2old
                    lattice.F5.[i0] <- f7old
                    lattice.F6.[i0] <- f8old
                    lattice.F7.[i0] <- f5old
                    lattice.F8.[i0] <- f6old

    /// Add periodic boundry conditions in the flow direction (x-axis/i-direction).
    member this.per_BC () =
        let ni, nj = lattice.Ni, lattice.Nj
        for i = 0 to ni - 1 do
            let i0 = i2D ni i 0
            let i1 = i2D ni i (nj - 1)
            lattice.F2.[i0] <- lattice.F2.[i1]
            lattice.F5.[i0] <- lattice.F5.[i1]
            lattice.F6.[i0] <- lattice.F6.[i1]
            lattice.F4.[i1] <- lattice.F4.[i0]
            lattice.F7.[i1] <- lattice.F7.[i0]
            lattice.F8.[i1] <- lattice.F8.[i0]

    /// Add inflow boundry conditions (drives the system).
    member this.in_BC () =
        let vx_term = 1.0f + 3.0f*vxIn + 3.0f*vxIn*vxIn
        let f1new = roOut * faceq2 * vx_term
        let f5new = roOut * faceq3 * vx_term
        let f8new = f5new

        for j = 0 to lattice.Nj - 1 do
          let i0 = i2D  lattice.Ni 0 j
          lattice.F1.[i0] <- f1new
          lattice.F5.[i0] <- f5new
          lattice.F8.[i0] <- f8new

/// Simplistic outflow boundary conditions. All the f values pointing
/// into the domain at the exit (ni-1) are set equal to those one node into
/// the domain (ni-2)
    member this.ex_BC_crude () =
        let ni, nj = lattice.Ni, lattice.Nj
        for j = 0 to nj - 1 do
            let i0 = i2D ni (ni - 1) j
            let i1 = i0 - 1
            lattice.F3.[i0] <- lattice.F3.[i1]
            lattice.F6.[i0] <- lattice.F6.[i1]
            lattice.F7.[i0] <- lattice.F7.[i1]

/// Combine all boundary conditions. 
    member this.apply_BCs solid =
        this.per_BC ()
        this.solid_BC solid
        this.in_BC ()
        this.ex_BC_crude ()

/// Applies a full Lattice Boltzmann Step.
    member this.applyLbmStep solid =
        this.stream ()
        this.apply_BCs solid
        this.collide ()

    member this.getVelocity () =
        velocity