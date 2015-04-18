module program

open System
open LattieBoltzmann

[<EntryPoint>]
let main argv =
//    LBM.Display.runSim()
    Tests.``poiseuille flow``()
    0