﻿// Colormap for displaying velocities.
module LBM.PlottingColours

    let nCol = 236
    let cmap_rgb = 
        [|         
            0.0000, 0.0000, 1.0000
            0.0000, 0.0157, 1.0000
            0.0000, 0.0314, 1.0000
            0.0000, 0.0471, 1.0000
            0.0000, 0.0627, 1.0000
            0.0000, 0.0784, 1.0000
            0.0000, 0.0941, 1.0000
            0.0000, 0.1098, 1.0000
            0.0000, 0.1255, 1.0000
            0.0000, 0.1412, 1.0000
            0.0000, 0.1569, 1.0000
            0.0000, 0.1725, 1.0000
            0.0000, 0.1882, 1.0000
            0.0000, 0.2039, 1.0000
            0.0000, 0.2196, 1.0000
            0.0000, 0.2353, 1.0000
            0.0000, 0.2510, 1.0000
            0.0000, 0.2667, 1.0000
            0.0000, 0.2824, 1.0000
            0.0000, 0.2980, 1.0000
            0.0000, 0.3137, 1.0000
            0.0000, 0.3294, 1.0000
            0.0000, 0.3451, 1.0000
            0.0000, 0.3608, 1.0000
            0.0000, 0.3765, 1.0000
            0.0000, 0.3922, 1.0000
            0.0000, 0.4078, 1.0000
            0.0000, 0.4235, 1.0000
            0.0000, 0.4392, 1.0000
            0.0000, 0.4549, 1.0000
            0.0000, 0.4706, 1.0000
            0.0000, 0.4863, 1.0000
            0.0000, 0.5020, 1.0000
            0.0000, 0.5176, 1.0000
            0.0000, 0.5333, 1.0000
            0.0000, 0.5490, 1.0000
            0.0000, 0.5647, 1.0000
            0.0000, 0.5804, 1.0000
            0.0000, 0.5961, 1.0000
            0.0000, 0.6118, 1.0000
            0.0000, 0.6275, 1.0000
            0.0000, 0.6431, 1.0000
            0.0000, 0.6588, 1.0000
            0.0000, 0.6745, 1.0000
            0.0000, 0.6902, 1.0000
            0.0000, 0.7059, 1.0000
            0.0000, 0.7216, 1.0000
            0.0000, 0.7373, 1.0000
            0.0000, 0.7529, 1.0000
            0.0000, 0.7686, 1.0000
            0.0000, 0.7843, 1.0000
            0.0000, 0.8000, 1.0000
            0.0000, 0.8157, 1.0000
            0.0000, 0.8314, 1.0000
            0.0000, 0.8471, 1.0000
            0.0000, 0.8627, 1.0000
            0.0000, 0.8784, 1.0000
            0.0000, 0.8941, 1.0000
            0.0000, 0.9098, 1.0000
            0.0000, 0.9255, 1.0000
            0.0000, 0.9412, 1.0000
            0.0000, 0.9569, 1.0000
            0.0000, 0.9725, 1.0000
            0.0000, 0.9882, 1.0000
            0.0000, 1.0000, 0.9961
            0.0000, 1.0000, 0.9804
            0.0000, 1.0000, 0.9647
            0.0000, 1.0000, 0.9490
            0.0000, 1.0000, 0.9333
            0.0000, 1.0000, 0.9176
            0.0000, 1.0000, 0.9020
            0.0000, 1.0000, 0.8863
            0.0000, 1.0000, 0.8706
            0.0000, 1.0000, 0.8549
            0.0000, 1.0000, 0.8392
            0.0000, 1.0000, 0.8235
            0.0000, 1.0000, 0.8078
            0.0000, 1.0000, 0.7922
            0.0000, 1.0000, 0.7765
            0.0000, 1.0000, 0.7608
            0.0000, 1.0000, 0.7451
            0.0000, 1.0000, 0.7294
            0.0000, 1.0000, 0.7137
            0.0000, 1.0000, 0.6980
            0.0000, 1.0000, 0.6824
            0.0000, 1.0000, 0.6667
            0.0000, 1.0000, 0.6510
            0.0000, 1.0000, 0.6353
            0.0000, 1.0000, 0.6196
            0.0000, 1.0000, 0.6039
            0.0000, 1.0000, 0.5882
            0.0000, 1.0000, 0.5725
            0.0000, 1.0000, 0.5569
            0.0000, 1.0000, 0.5412
            0.0000, 1.0000, 0.5255
            0.0000, 1.0000, 0.5098
            0.0000, 1.0000, 0.4941
            0.0000, 1.0000, 0.4784
            0.0000, 1.0000, 0.4627
            0.0000, 1.0000, 0.4471
            0.0000, 1.0000, 0.4314
            0.0000, 1.0000, 0.4157
            0.0000, 1.0000, 0.4000
            0.0000, 1.0000, 0.3843
            0.0000, 1.0000, 0.3686
            0.0000, 1.0000, 0.3529
            0.0000, 1.0000, 0.3373
            0.0000, 1.0000, 0.3216
            0.0000, 1.0000, 0.3059
            0.0000, 1.0000, 0.2745
            0.0000, 1.0000, 0.2431
            0.0000, 1.0000, 0.2118
            0.0000, 1.0000, 0.1804
            0.0000, 1.0000, 0.1490
            0.0000, 1.0000, 0.1176
            0.0000, 1.0000, 0.0863
            0.0000, 1.0000, 0.0549
            0.0000, 1.0000, 0.0235
            0.0078, 1.0000, 0.0000
            0.0392, 1.0000, 0.0000
            0.0706, 1.0000, 0.0000
            0.1020, 1.0000, 0.0000
            0.1333, 1.0000, 0.0000
            0.1647, 1.0000, 0.0000
            0.1961, 1.0000, 0.0000
            0.2275, 1.0000, 0.0000
            0.2588, 1.0000, 0.0000
            0.2902, 1.0000, 0.0000
            0.3216, 1.0000, 0.0000
            0.3373, 1.0000, 0.0000
            0.3529, 1.0000, 0.0000
            0.3686, 1.0000, 0.0000
            0.3843, 1.0000, 0.0000
            0.4000, 1.0000, 0.0000
            0.4157, 1.0000, 0.0000
            0.4314, 1.0000, 0.0000
            0.4471, 1.0000, 0.0000
            0.4627, 1.0000, 0.0000
            0.4784, 1.0000, 0.0000
            0.4941, 1.0000, 0.0000
            0.5098, 1.0000, 0.0000
            0.5255, 1.0000, 0.0000
            0.5412, 1.0000, 0.0000
            0.5569, 1.0000, 0.0000
            0.5725, 1.0000, 0.0000
            0.5882, 1.0000, 0.0000
            0.6039, 1.0000, 0.0000
            0.6196, 1.0000, 0.0000
            0.6353, 1.0000, 0.0000
            0.6510, 1.0000, 0.0000
            0.6667, 1.0000, 0.0000
            0.6824, 1.0000, 0.0000
            0.6980, 1.0000, 0.0000
            0.7137, 1.0000, 0.0000
            0.7294, 1.0000, 0.0000
            0.7451, 1.0000, 0.0000
            0.7608, 1.0000, 0.0000
            0.7765, 1.0000, 0.0000
            0.7922, 1.0000, 0.0000
            0.8078, 1.0000, 0.0000
            0.8235, 1.0000, 0.0000
            0.8392, 1.0000, 0.0000
            0.8549, 1.0000, 0.0000
            0.8706, 1.0000, 0.0000
            0.8863, 1.0000, 0.0000
            0.9020, 1.0000, 0.0000
            0.9176, 1.0000, 0.0000
            0.9333, 1.0000, 0.0000
            0.9490, 1.0000, 0.0000
            0.9647, 1.0000, 0.0000
            0.9804, 1.0000, 0.0000
            0.9961, 1.0000, 0.0000
            1.0000, 0.9882, 0.0000
            1.0000, 0.9725, 0.0000
            1.0000, 0.9569, 0.0000
            1.0000, 0.9412, 0.0000
            1.0000, 0.9255, 0.0000
            1.0000, 0.9098, 0.0000
            1.0000, 0.8941, 0.0000
            1.0000, 0.8784, 0.0000
            1.0000, 0.8627, 0.0000
            1.0000, 0.8471, 0.0000
            1.0000, 0.8314, 0.0000
            1.0000, 0.8157, 0.0000
            1.0000, 0.8000, 0.0000
            1.0000, 0.7843, 0.0000
            1.0000, 0.7686, 0.0000
            1.0000, 0.7529, 0.0000
            1.0000, 0.7373, 0.0000
            1.0000, 0.7216, 0.0000
            1.0000, 0.7059, 0.0000
            1.0000, 0.6902, 0.0000
            1.0000, 0.6745, 0.0000
            1.0000, 0.6588, 0.0000
            1.0000, 0.6431, 0.0000
            1.0000, 0.6275, 0.0000
            1.0000, 0.6118, 0.0000
            1.0000, 0.5961, 0.0000
            1.0000, 0.5804, 0.0000
            1.0000, 0.5647, 0.0000
            1.0000, 0.5490, 0.0000
            1.0000, 0.5333, 0.0000
            1.0000, 0.5176, 0.0000
            1.0000, 0.5020, 0.0000
            1.0000, 0.4863, 0.0000
            1.0000, 0.4706, 0.0000
            1.0000, 0.4549, 0.0000
            1.0000, 0.4392, 0.0000
            1.0000, 0.4235, 0.0000
            1.0000, 0.4078, 0.0000
            1.0000, 0.3922, 0.0000
            1.0000, 0.3765, 0.0000
            1.0000, 0.3608, 0.0000
            1.0000, 0.3451, 0.0000
            1.0000, 0.3294, 0.0000
            1.0000, 0.3137, 0.0000
            1.0000, 0.2980, 0.0000
            1.0000, 0.2824, 0.0000
            1.0000, 0.2667, 0.0000
            1.0000, 0.2510, 0.0000
            1.0000, 0.2353, 0.0000
            1.0000, 0.2196, 0.0000
            1.0000, 0.2039, 0.0000
            1.0000, 0.1882, 0.0000
            1.0000, 0.1725, 0.0000
            1.0000, 0.1569, 0.0000
            1.0000, 0.1412, 0.0000
            1.0000, 0.1255, 0.0000
            1.0000, 0.1098, 0.0000
            1.0000, 0.0941, 0.0000
            1.0000, 0.0784, 0.0000
            1.0000, 0.0627, 0.0000
            1.0000, 0.0471, 0.0000
            1.0000, 0.0314, 0.0000
            1.0000, 0.0157, 0.0000
            1.0000, 0.0000, 0.0000
        |] |> Array.map (fun (rcol,gcol,bcol) -> ((int)(255.0f) <<< 24) ||| ((int)(bcol * 255.0) <<< 16) ||| ((int)(gcol * 255.0) <<<  8) ||| ((int)(rcol * 255.0) <<<  0) )

    let nColGray = 64

    let cmap_gray = 
        [|         
            0.0000
            0.0157
            0.0314
            0.0471
            0.0627
            0.0784
            0.0941
            0.1098
            0.1255
            0.1412
            0.1569
            0.1725
            0.1882
            0.2039
            0.2196
            0.2353
            0.2510
            0.2667
            0.2824
            0.2980
            0.3137
            0.3294
            0.3451
            0.3608
            0.3765
            0.3922
            0.4078
            0.4235
            0.4392
            0.4549
            0.4706
            0.4863
            0.5020
            0.5176
            0.5333
            0.5490
            0.5647
            0.5804
            0.5961
            0.6118
            0.6275
            0.6431
            0.6588
            0.6745
            0.6902
            0.7059
            0.7216
            0.7373
            0.7529
            0.7686
            0.7843
            0.8000
            0.8157
            0.8314
            0.8471
            0.8627
            0.8784
            0.8941
            0.9098
            0.9255
            0.9412
            0.9569
            0.9725
            0.9882
        |] |> Array.map (fun gray -> (uint32)(gray * 255.0) )