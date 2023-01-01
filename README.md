# SkyImages.jl

Load astronomical images of the sky and process them with convenient, general, and composable functions.

`SkyImages.jl` focuses on the functionality to **load images and transform them** to the shape most suitable for further analysis --- whatever this shape is in your particular case.

Loaded images are directly **represented as keyed arrays** (from `AxisKeys.jl`) with axis keys being special objects defined here in `SkyImages.jl`. These axis key objects know how to translate between pixels and sky coordinates on the basis of FITS WCS or other projection (see e.g. a healpix example below).

It's most convenient to work with regular rectangular array, and `SkyImages` gives multiple way to obtain image in such form. It can **project onto any reasonable coordinate grid**, or **extract the original data array** as stored in FITS with its coordinates.

_Extra features such as plotting are explicitly out of scope for `SkyImages`. There are multiple commonly used plotting libraries, and we don't impose any specific choice. Examples in this notebook use `matplotlib` through `PyPlot.jl`._

See the [Pluto notebook](https://aplavin.github.io/SkyImages.jl/test/notebook.jl.html) for more details and examples.
