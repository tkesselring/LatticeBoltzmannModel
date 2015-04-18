#r @"packages/FAKE/tools/FakeLib.dll"

open Fake

let resultsDir = "release"

Target "Clean" (fun _ ->
    DeleteDirs [resultsDir]
)

Target "Build" (fun _ ->
    !! "/**/*.*proj"
    |> MSBuildRelease resultsDir "Build"
    |> Log "Build-Output: "
)

Target "Tests" (fun _ ->
    !! "/**/*.exe"
//    ++ "/**/Test.*.dll"
    |> SetBaseDir resultsDir
    |> NUnitParallel (fun defaults -> { defaults with Framework = "net-4.5"})
)

Target "Default" DoNothing

"Clean"
    ==> "Build"
    =?> ("Tests", not <| hasBuildParam "NoTests")
    ==> "Default"

RunTargetOrDefault "Default"
