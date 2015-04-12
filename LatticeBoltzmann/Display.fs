/// running Lattice Boltzmann model and displaying it using openTK/openGL.
module LBM.Display

open System

open OpenTK
open OpenTK.Graphics
open OpenTK.Graphics.OpenGL
open OpenTK.Input

open LattieBoltzmann

let h = 112
let l = 320

#nowarn "9"
#nowarn "51"

type LBSimWindow() =
    inherit GameWindow(l, h, GraphicsMode.Default, "Lattice Bolzmann simulation")
    let ni = l          // system size in i-direction
    let nj = h          // system size in j-direction
    let vxIn = 0.04f    // influx of fluid
    let roOut = 1.0f    // outflow of fluid
    let tau = 0.51f     // collision tau

    let lbm = LatticeBoltzmannModel(ni, nj, roOut, vxIn, tau)

    let mutable glTex = 0
    let mutable glPBO = 0

    let mutable xPosOld, yPosOld = 0,0

    let solid = Array.init (ni*nj) (fun i -> if i / ni = 0 || i / ni = (nj-1) then 0.0f
                                             else 1.0f)

    do        
        GL.ClearColor(0.0f, 0.0f, 0.0f, 0.0f)
        GL.MatrixMode(MatrixMode.Projection)
        GL.LoadIdentity()
        GL.Ortho(0.0, float ni, 0.0, float nj, -200.0, 200.0)

        GL.Enable(EnableCap.Texture2D)
        GL.GenTextures(1, &&glTex)
        GL.BindTexture(TextureTarget.Texture2D, glTex)
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapS, (int) TextureWrapMode.Clamp)      
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureWrapT, (int) TextureWrapMode.Clamp)      
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMagFilter, (int) TextureMagFilter.Linear)
        GL.TexParameter(TextureTarget.Texture2D, TextureParameterName.TextureMinFilter, (int) TextureMagFilter.Linear)

        GL.TexImage2D(TextureTarget.Texture2D, 0, PixelInternalFormat.Rgba8, ni, nj, 0, PixelFormat.Rgba, PixelType.UnsignedByte, (null : _ []))

        GL.GenBuffers(1, &&glPBO)
        GL.BindBuffer(BufferTarget.PixelUnpackBuffer, glPBO)

    override this.Dispose(disposing) =
        base.Dispose(disposing)

    override this.OnLoad e =
        base.OnLoad(e)
        this.VSync <- VSyncMode.Off

    override this.OnRenderFrame e =
        base.OnRenderFrame(e)

        let mutable isol = 0.0f

        /// Do one Lattice Boltzmann step (stream, BC, collide):
        lbm.applyLbmStep solid

        let velocities = lbm.getVelocity ()

        let minvar = 0.0f
        let maxvar = 0.2f
        
        // Create colourful array for plotting using the color dataset in data.fd
        let ncol = LBM.PlottingColours.nCol - 1
        let plot_rgba = Array.init (ni*nj) (fun i -> 
                let frac = (velocities.[i] - minvar) / (maxvar - minvar)
                let icol = (int) (frac*float32 ncol)
                let isol = (int) solid.[i]
                isol*LBM.PlottingColours.cmap_rgb.[max (min icol ncol) 0]
            )

        // Fill the pixel buffer with the plot_rgba array.
        GL.BindBuffer(BufferTarget.ArrayBuffer, 0 )
        GL.BufferData(BufferTarget.PixelUnpackBuffer, IntPtr (sizeof<uint32>*ni*nj), plot_rgba, BufferUsageHint.StreamCopy)

        // Copy the pixel buffer to the texture, ready to display.
        GL.TexSubImage2D(TextureTarget.Texture2D, 0, 0, 0, ni, nj, PixelFormat.Rgba, PixelType.UnsignedByte, IntPtr 0)

        // Use a quad to display texture:
        // i.e. showing the volocities in color code from data.fs.
        GL.Clear(ClearBufferMask.ColorBufferBit)
        GL.Begin(PrimitiveType.Quads)
        GL.TexCoord2(0.0, 0.0)
        GL.Vertex3(0.0, 0.0, 0.0)
        GL.TexCoord2(1.0, 0.0)
        GL.Vertex3(float ni, 0.0, 0.0)
        GL.TexCoord2(1.0, 1.0)
        GL.Vertex3(float ni,float nj, 0.0)
        GL.TexCoord2(0.0, 1.0)
        GL.Vertex3(0.0, float nj, 0.0)
        GL.End()

        GL.PopMatrix()
        this.SwapBuffers()

    override this.OnMouseDown e =
        base.OnMouseDown(e)
        if base.Mouse.Item MouseButton.Left then
            let x = ni*e.Position.X/this.ClientRectangle.Width - 1
            let y = nj - nj*e.Position.Y/this.ClientRectangle.Height - 1
            if x >= 0 && x < ni && y >= 0 && y < nj then
                solid.[ni*y + x] <- 0.0f

        if base.Mouse.Item MouseButton.Right then
            let x = ni*e.Position.X/this.ClientRectangle.Width - 1
            let y = nj - nj*e.Position.Y/this.ClientRectangle.Height - 1
            if x >= 0 && x < ni && y >= 0 && y < nj then
                solid.[ni*y + x] <- 1.0f

    override this.OnResize e =
        base.OnResize(e)
        GL.Viewport(0, 0, this.ClientRectangle.Width, this.ClientRectangle.Height)
        GL.MatrixMode(MatrixMode.Projection)
        GL.LoadIdentity()
        GL.Ortho(0.0, float ni, 0.0, float nj, -200.0 ,200.0)
        GL.MatrixMode(MatrixMode.Modelview)
        GL.LoadIdentity()

(**
Defaults to use `SimWindow`.
*)
let runSim() =
    use window = new LBSimWindow()
    window.Run()