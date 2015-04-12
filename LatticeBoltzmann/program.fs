module program

open System
open LattieBoltzmann

[<EntryPoint>]
let main argv =
    LBM.Display.runSim()
    0