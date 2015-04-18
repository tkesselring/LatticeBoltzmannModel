module Tests

open NUnit.Framework
open FsUnit
open LattieBoltzmann


[<Test>]
let ``poiseuille flow`` () =
    let ni = 20      // system size in i-direction
    let nj = 10      // system size in j-direction
    let vxIn = 0.04f // influx of fluid
    let roOut = 1.0f // outflow of fluid
    let tau = 0.6f   // collision tau
    let lbm = LatticeBoltzmannModel(ni, nj, roOut, vxIn, tau)

    let solid = Array.init (ni*nj) (fun i -> if i/ni = 0 || i/ni = (nj-1) then 0.0f
                                             else 1.0f)
    for i = 0 to 500 do
        lbm.applyLbmStep solid

    let velocities = lbm.getVelocity ()
    let velocities_ = Array.init nj (fun j -> velocities.[ni/2 + ni*j])
    let expected = [|0.0003876934; 0.009466321; 0.0188583955; 0.0253776275; 0.02863464; 0.02863464; 0.0253776275; 0.0188583955; 0.009466321; 0.0003876934|]
    velocities_ |> should (equalWithin 1e-5) expected