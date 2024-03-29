### A Pluto.jl notebook ###
# v0.19.11

using Markdown
using InteractiveUtils

# ╔═╡ e3c72be9-8488-4c01-843c-16553057a0b5
begin
	using Revise
	import Pkg
	eval(:(Pkg.develop(path="..")))
	Pkg.resolve()
	using SkyImages
end

# ╔═╡ 37bd1aba-72ba-4b42-b8aa-789318691328
using Unitful, UnitfulAngles

# ╔═╡ 8e437f0f-20fd-40b7-bfac-e660d8dfa7b8
using AxisKeys

# ╔═╡ bae108cc-f9d5-4667-aa51-6fc3481b02c5
using RectiGrids: grid

# ╔═╡ f27b825e-b90f-4b8b-a04c-bc7df379383c
using SkyCoords

# ╔═╡ eb8eff01-cdaf-4536-8bd9-b6b99665a5ac
using PyPlotUtils; pyplot_style!()

# ╔═╡ c1f7876c-9262-445f-b199-6d7f46434294
# ╠═╡ show_logs = false
using DataPipes

# ╔═╡ 36ccd512-b91d-47a7-a86d-a0baeb7faf32
# ╠═╡ show_logs = false
using DisplayAs: Text as AsText, DisplayAs

# ╔═╡ 9d343fd3-0258-41e4-9e06-938282491ecf
# ╠═╡ show_logs = false
using StatsBase

# ╔═╡ 8dc78491-a634-4290-9ecc-c13e9be46c2e
# ╠═╡ show_logs = false
using SkipNan

# ╔═╡ f7bb31e5-9a69-4a9d-bd62-02c3b7d0d332
# ╠═╡ show_logs = false
using AccessorsExtra

# ╔═╡ 1c17fd5e-f047-4bd4-8b00-56fa077c7bb0
# ╠═╡ show_logs = false
using PyFormattedStrings

# ╔═╡ 5ebbcfb6-f534-490f-947f-b7c039156908
using PlutoUI

# ╔═╡ d27f13ff-72de-4136-b145-d8a973f23047
md"""
!!! info "SkyImages.jl"
	Load astronomical images of the sky and process them with convenient, general, and composable functions.
"""

# ╔═╡ a0a180c0-15c5-408d-82a2-3e8d5f1ff484
md"""
`SkyImages.jl` focuses on the functionality to **load images and transform them** to the shape most suitable for further analysis --- whatever this shape is in your particular case.

Loaded images are directly **represented as keyed arrays** (from `AxisKeys.jl`) with axis keys being special objects defined here in `SkyImages.jl`. These axis key objects know how to translate between pixels and sky coordinates on the basis of FITS WCS or other projection (see e.g. a healpix example below).

It's most convenient to work with regular rectangular array, and `SkyImages` gives multiple way to obtain image in such form. It can **project onto any reasonable coordinate grid**, or **extract the original data array** as stored in FITS with its coordinates.

_Extra features such as plotting are explicitly out of scope for `SkyImages`. There are multiple commonly used plotting libraries, and we don't impose any specific choice. Examples in this notebook use `matplotlib` through `PyPlot.jl`._
"""

# ╔═╡ aa6da229-640c-4eae-a0fe-9bd2742fa849
md"""
Notable differences from `AstroImages.jl`:
- `AstroImages` makes image a custom type and implements plotting functions specific to this type. \
  `SkyImages` represents images as regular keyed arrays, so general plotting packages can plot them just fine. One doesn't easily get the non-uniform grid lines plotted though, unlike `AstroImages`.
- `SkyImages` directly treats images as having coordinates (`SkyCoords.jl`), not just pixel indices. For example, one can:
  - Directly extract values at coordinates of interest
  - Project an image onto an arbitrary coordinate grid, see various examples below
  - Extract the image bounding box and easily generate a coordinate grid on it (`RectiGrids.jl` interface)
  Using `AstroImages` involves explicitly working with WCS objects more often.
- `SkyImages` support healpix FITS images, in addition to FITS + WCS images
- For now, we only support celestial WCS. Other axes such as the polarization or frequency aren't implemented yet, but this is not a fundamental limitation.
"""

# ╔═╡ b02f575f-bbd5-42c0-aab3-37415b66bd76
md"""
# Detailed example
"""

# ╔═╡ 72a56fff-84ab-4fdb-882f-e0ee6efada45
md"""
Let's load one of the images. This is the map of the whole sky in gamma rays, as see by Fermi LAT.
"""

# ╔═╡ 3d581238-56b7-44b9-a9c3-9af9530ae136
# uncomment to download:
# download("https://fermi.gsfc.nasa.gov/ssc/data/access/lat/12yr_catalog/intens_scaled_ait_144m_gt1000_psf3_gal_0p1.fits.gz", "data/intens_scaled_ait_144m_gt1000_psf3_gal_0p1.fits")

# ╔═╡ 7de9fc6e-f652-48bc-babe-af0e56a8c427
simg = SkyImages.load("data/intens_scaled_ait_144m_gt1000_psf3_gal_0p1.fits");

# ╔═╡ 41f3ff07-8ebe-4224-ad1b-9342ee4e7e07
simg |> AsText

# ╔═╡ b2649127-99b3-4ee4-aa54-d11f46f2c04e
md"""
Images are represented in `SkyImages` as keyed arrays (from `AxisKeys.jl`), with one "axis" corresponding to celestial coordinates. Other axes _(not present in this example)_ could represent frequencies, polarizations, or other parameters.

Displaying keyed arrays prints the axis keys (coordinates) as the first column, and array values as the second. This looks cleaner in the Julia REPL, while Pluto doesn't show colors.
"""

# ╔═╡ f26a1360-ba03-4df3-8ba0-ffd707a34bfb
md"""
All shown values are zero, but this is only true for array edges. `extrema()` confirm that:
"""

# ╔═╡ e42ac814-ca50-4464-8a19-59695b1f9ef0
extrema(simg)

# ╔═╡ ec6d0165-18e5-41f6-91e0-382840cde0f0
md"""
Coordinate values are shown here, but they are never materialized. Instead, the image axis computes them on demand when needed. \
For example, in this case coordinates are defined by the World Coordinate System read from the FITS file. The WCS object can directly be accessed, if necessary:
"""

# ╔═╡ b0c2a790-231a-4ebc-a755-627cca858325
axiskeys(simg, :coords).wcs

# ╔═╡ 170f4c04-0961-44d5-876e-9e9f51834e32
md"""
We can extract image values closest to any given coordinates using the `Near` selector (defined in `AxisKeys.jl`). \
`Near` supports both the native image coordinates -- galactic in this case:
"""

# ╔═╡ 4285e435-549e-487f-9013-95f8d2958cbe
simg[Near(GalCoords(0, 0))]

# ╔═╡ de270142-3c0c-45b9-968f-0ef68c5074b0
md"""
... and any other kind of coordinates as well:
"""

# ╔═╡ 62fd3d3c-d8a6-48e4-8d03-c27d2ef581f0
simg[Near(ICRSCoords(0, 0))]

# ╔═╡ 6ed9fc21-08b6-4f39-b9c3-29abbcd3cee1
md"""
Alternatively, instead of extracting the nearest value, we can use interpolation of different orders. \
Zeroth order -- constant (same as the `Near` selector), first order -- linear, and so on:
"""

# ╔═╡ 17671ee5-bda8-4ffd-b1ef-dbb98eb7339a
simg(SkyImages.Interp(ICRSCoords(0, 0), order=0))

# ╔═╡ 16bc419e-1c2a-4617-ba06-fcaa943066ca
simg(SkyImages.Interp(ICRSCoords(0, 0), order=1))

# ╔═╡ 6cb2f138-5bba-406a-a505-4c60645d1dcb
md"""
Arrays of `Near` or `Interp` are selectors themselves, which is the base `AxisKeys` semantic. Such vectorized selections are implemented efficiently in `SkyImages`, and can directly be used to obtain a regular rectangular image in any coordinates.

First, let's create a _(lazy)_ 500×500 coordinate grid in equatorial coordinates:
"""

# ╔═╡ 6c7b4e2d-7f67-428a-8d22-c2bfb0a96a42
coords = grid(ICRSCoords, range(π, -π; length=500), range(-π/2, π/2; length=500));

# ╔═╡ e052b45d-ef9e-41de-aeee-c0fce166f4f6
coords |> DisplayAs.withcontext(:displaysize => (20, 200))

# ╔═╡ a6d02948-0e29-4fa6-a2b6-a8f7d605fb97
md"""
This is a regular keyed array, with keys being RA and Dec, and values being the actual coordinates.
"""

# ╔═╡ 9905274c-40b2-4f8c-bf87-9f0f939699bf
md"""
Select the closest image value to each coordinate from this grid:
"""

# ╔═╡ 7ac37418-514a-4b41-beae-2bec303c9737
eq_img = simg(Near.(coords)) |> copy;

# ╔═╡ 83946094-3f59-4d98-92cc-b79558f6a504
eq_img |> AsText

# ╔═╡ 8efcbb46-9685-479f-96b9-816aacd2a851
md"""
And transform both axes keys to degrees for convenience:
"""

# ╔═╡ 59c20159-4614-47d0-98ab-bf6f4991caed
eq_deg_img = @modify(xs -> xs .* rad2deg(1), axiskeys(eq_img) |> Elements());

# ╔═╡ e16e262d-8366-41cc-ba5c-085188dbe3c7
md"""
That's the visualization of what we got:
"""

# ╔═╡ e9b9da8e-80c3-47dc-b1e4-db6d490e7d12
begin
	plt.figure()
	imshow_ax(eq_deg_img)
	xylabels("RA (°)", "Dec (°)")
	plt.gcf()
end

# ╔═╡ d177d489-072c-4dcd-a72f-0eeab7fedc2a
md"""
# Images of compact regions
"""

# ╔═╡ 3b5fbd3c-cfff-4f21-9736-589e9c01e656
md"""
## Extract bounding box
"""

# ╔═╡ a5ed3228-f068-4784-9d7a-61530cbd40e6
md"""
The first example demonstrated an image that covers the whole sky, and we manually generated the target coordinate grid. \
It was easy for the whole sky, but what to do when the image only covers a small part -- some specific object?

First, let's load such an image:
"""

# ╔═╡ 6eaf40d7-5513-4a03-a77d-a1f38aa47af2
# uncomment to download:
# download("http://www.astro.uvic.ca/~wthompson/astroimages/fits/656nmos.fits", "data/656nmos.fits")

# ╔═╡ 18d26c89-e953-46bf-9a96-512d55c75a7d
small_img = SkyImages.load("data/656nmos.fits");

# ╔═╡ 0c7b5add-fc89-4633-8c1b-9f472cb5707e
md"""
We want to extract its bounding box, and generate a grid over that box. The bbox can be calculated from the corresponding WCS axiskeys:
"""

# ╔═╡ 8f6a723c-bc5c-46e9-aaeb-af13a25ba4d8
axk = axiskeys(small_img, :coords)

# ╔═╡ 5ab12071-3a83-452c-b85e-7b475b5a86e2
md"""
and here's the bbox itself:
"""

# ╔═╡ dfb02cfb-5451-44da-bb74-716aa0fb9cb8
bbox = boundingbox(axk)

# ╔═╡ cb61bd15-0c0c-43f8-abc5-4585b7544c16
md"""
It's always a rectangle in some coordinates. The default choice for the coordinate type is taken from the WCS itself: in this case, the image is defined in equatorial FK5 coordinates.

Now, create a _(again, lazy)_ 1000×1000 grid of coordinates over this bbox:
"""

# ╔═╡ 50c0fd27-9f9d-4fd5-a053-8a9cbdc06f81
small_coords = grid(bbox; lengths=1000);

# ╔═╡ 385d8b97-acb4-4cc0-b258-a4237bbf296e
small_coords |> DisplayAs.withcontext(:displaysize => (20, 200))

# ╔═╡ c7a78987-6c06-41e9-bcc3-067468efe961
md"""
... and compute the image values on this grid, same as we did in the first example:
"""

# ╔═╡ d0de65a3-82e1-44f7-be84-5c09511adbda
eagle_plt = begin
	plt.figure()
	@p begin
		small_img(SkyImages.Interp.(small_coords, order=1))
		@modify(xs -> xs .* rad2deg(1), __ |> axiskeys(_) |> Elements())
		imshow_ax(__; norm=matplotlib.colors.Normalize(vmin=0, vmax=quantile(skipnan(__), 1-3e-4)))
	end
	xylabels("RA (°)", "Dec (°)")
	plt.xticks(rotation=15)
	plt.gcf()
end

# ╔═╡ e2e336f8-26da-456f-8feb-546086b97861
md"""
## Projected coordinates
"""

# ╔═╡ 7180a98a-1639-4647-a973-5cb4a654776c
md"""
It's common to work with astronomical images in so-called projected coordinates. They have their origin at a certain point on the celestial sphere, and define local axes that typically go along the corresponding latitude and longitude. The main difference from regular coordinates is that vertical and horizontal scales are the same: one arcsecond of local longitude is always the same separation as one degree of local latitude. This greatly simplifies working with images of compact sky regions, making them effectively flat.
"""

# ╔═╡ 7b1bc5cd-15a2-4b6e-bc23-983201ea6f88
md"""
`SkyImages` let you compute the image bounding box in projected coordinates as well. Here we use projected Galactic coordinates:
"""

# ╔═╡ 8b2bc998-e4d6-42eb-88b7-abec6d3c5b1b
bbox_proj = boundingbox(ProjectedCoords{GalCoords}, axk)

# ╔═╡ f4c9635b-86b9-4ea7-8433-6d42f347695f
md"""
Everything else stays the same as before:
"""

# ╔═╡ 918e65d6-2e7f-4db5-9937-ebcd719b42ff
begin
	plt.figure()
	@p begin
		grid(bbox_proj; lengths=1000)
		small_img(SkyImages.Interp.(__, order=2))
		@modify(xs -> xs .* rad2deg(1), axiskeys(__) |> Elements());
		imshow_ax(__; norm=matplotlib.colors.Normalize(vmin=0, vmax=quantile(skipnan(__), 1-3e-4)))
	end
	xylabels("Rel Gal lon (°)", "Rel Gal lat (°)")
	plt.xticks(rotation=15)
	plt.gcf()
end

# ╔═╡ e7952445-dbfc-4e06-a00e-d23038f6435b
md"""
# All-sky images
"""

# ╔═╡ 6555eb04-d630-4d03-800b-5d023c8e0522
md"""
We've shown an all-sky image in the first example. However, there are a few other notable details.
"""

# ╔═╡ 2baea1f0-2a27-4aca-991b-6f3e7bb13c64
md"""
## Spherical projections
"""

# ╔═╡ 5e631aff-d7a4-487c-b21a-03c47be08d8e
md"""
Celestial images are often shown in a common spherical projection instead of a simple rectangle. \
Here we show how to do this with `matplotlib` (`PyPlot.jl`):
"""

# ╔═╡ e98ccebd-b7fd-4de3-8e16-4056512ed1a9
fermi_galplt = let
	plt.figure(figsize=(10, 6))
	plt.subplot(111, projection="mollweide")
	@p begin
		SkyImages.load("data/intens_scaled_ait_144m_gt1000_psf3_gal_0p1.fits")
		__(Near.(grid(GalCoords, range(-π, π; length=1500), range(-π/2, π/2; length=1500))))
		@modify(-, __ |> axiskeys(_, :l) |> Elements())
		pcolormesh_ax(__; cmap=:gnuplot2)
	end
	plt.gcf()
end

# ╔═╡ a97f89d2-b8fc-4b7b-b0ab-c4170017db27
md"""
Note that base `matplotlib` cannot invert the horizontal axis of a spherical projection itself. Here the inversion is done manually before actual plotting, by negating the corresponding axis keys.
"""

# ╔═╡ 42be67e7-bc56-49d4-a391-fededc1bd1b6
md"""
## Healpix

Healpix is also a common format for astronomical images covering the whole sky. \
`SkyImages.jl` transparently supports them as well, with the same interface as FITS + WCS images shown above:
"""

# ╔═╡ 123564f2-8a13-47cb-9f9c-5af7fd19b6b8
# uncomment to download:
# download("https://faun.rc.fas.harvard.edu/dfink/skymaps/halpha/data/v1_1/healpix/Halpha_fwhm06_1024.fits", "data/Halpha_fwhm06_1024.fits")

# ╔═╡ 999c382b-5de9-4324-bea4-f6ddb1f37735
let 
	plt.figure(figsize=(10, 6))
	plt.subplot(111, projection="mollweide")
	@p begin
		SkyImages.load("data/Halpha_fwhm06_1024.fits")
		__(Near.(grid(ICRSCoords, range(-π, π; length=1500), range(-π/2, π/2; length=1500))))
		@modify(-, __ |> axiskeys(_, :ra) |> Elements())
		pcolormesh_ax(__; norm=matplotlib.colors.SymLogNorm(vmin=0, linthresh=1), cmap=:gnuplot2)
	end
	plt.gcf()
end

# ╔═╡ c03ddad8-d708-4d4f-8e08-ae7ce5a05bb9
md"""
# Advanced
"""

# ╔═╡ 3d868917-38a8-41cf-acb5-d5a1a2615448
md"""
## Multiple images of M87
"""

# ╔═╡ 487f817d-d67f-4152-bae6-7ec593f253fe
md"""
Here, we show a neat example of visualizing multiple images of the same astronomical object --- the M87 active galaxy.

Images span many orders of magnitude in extent and resolution, from kiloparsec to astronomical units. Together they provide a beautiful view of the M87 jet and the immediate surroundings of the central black hole.
"""

# ╔═╡ 31c8cf4b-cc86-4323-8587-16082443d3ae
# uncomment to download:
# download("https://github.com/yuanjea/eht-demo/blob/master/eht.fits", "data/m87_eht.fits")
# download("https://www.cv.nrao.edu/2cmVLBA/data/1228+126/2021_05_01/1228+126.u.2021_05_01.icn.fits.gz", "data/1228+126.u.2021_05_01.icn.fits.gz") # and unpack
# download("https://osf.io/download/te6jb/?view_only=9e8124dcdc7643cda6d980c142de23a3", "data/1-VLASS__2.2.ql.T14t19.J123049.42+122328.0.fits")
# download("https://osf.io/download/3jdpz/?view_only=9e8124dcdc7643cda6d980c142de23a3", "data/NVSS_1995-02-27_J123049.42+122328.0.fits")
# download("https://osf.io/download/v9gzy/?view_only=9e8124dcdc7643cda6d980c142de23a3", "data/hst_05122_01_wfpc2_f547m_wf_sci.fits")

# ╔═╡ 594c91e1-565a-49ad-9b77-20aabe7e71f4
imspecs = @p [
	(file="NVSS_1995-02-27_J123049.42+122328.0.fits", size=300u"as", dynrange=10, label="VLA NVSS"),
	(file="1-VLASS__2.2.ql.T14t19.J123049.42+122328.0.fits", size=45u"as", dynrange=7, label="VLASS"),
	(file="1228+126.u.2021_05_01.icn.fits", size=30u"mas", dynrange=300, label="VLBI MOJAVE"),
	(file="hst_05122_01_wfpc2_f547m_wf_sci.fits", size=25u"as", dynrange=100, label="HST"),
	(file="hst_05122_01_wfpc2_f547m_wf_sci.fits", size=2u"as", dynrange=10, label="HST zoom"),
	(file="m87_eht.fits", size=1u"mas", dynrange=1, label="EHT"),
	(file="intens_scaled_ait_144m_gt1000_psf3_gal_0p1.fits", size=2u"°", dynrange=1, label="Fermi"),
] |> @set(__ |> Elements() |> _.file |> dirname = "data")

# ╔═╡ 09a7db9c-0afb-41a1-a66b-924469ea003c
m87_coord = ICRSCoords(deg2rad.((187.70593075416664, 12.391123236111111))...)

# ╔═╡ 1c4de053-f8d5-4363-8224-af5a71594978
function center_at_peak(timg)
	peakix = @p timg |> skipnan |> argmax
	timg = @modify(x -> x .- x[peakix[1]], axiskeys(timg)[1])
	timg = @modify(x -> x .- x[peakix[2]], axiskeys(timg)[2])
end

# ╔═╡ 95b0e957-1deb-4176-b819-9d716a09e5b1
m87_images = map(imspecs) do p
	timg = @p begin
		SkyImages.load(p.file)
		# if crval is zero - no actual position in FITS, replace crval with m87 coords
		simg = @modify(__ |> axiskeys(_, :coords).wcs.crval |> If(iszero)) do _
			[SkyCoords.lon(m87_coord), SkyCoords.lat(m87_coord)] .|> rad2deg
		end
		# extract the bbox
		axiskeys(__, :coords)
		boundingbox(ProjectedCoords{ICRSCoords})
		# the original projection center can be arbitrary - we want (0, 0) to be M87
		SkyImages.project(m87_coord)
		# cut the bbox to twice the specified image size
		@modify(__ |> Properties() |> _.xy |> Elements() |> abs) do x
			min(x, ustrip(u"rad", 2*p.size))
		end
		# compute the image on a grid over this bbox
		grid(lengths=500)
		simg(SkyImages.Interp.(__, 1))
	end
	if p.size < 5u"as"
		# for small sub-as images, center them on the peak: stored image coordinates can be noticeably off at these scale
		timg = center_at_peak(timg)
	end
	# finally, cut the image to exactly the specified size
	timg = timg(0±p.size, 0±p.size)
	@insert p.image = timg
end;

# ╔═╡ 2bd8cb76-4c90-48f6-80d5-6c790045f314
let
	fig, ax = plt.subplots(1, length(m87_images); figsize=(length(m87_images) * 3, 3.5))
	@p begin
		m87_images
		# should go smallest to largest
		sort(by=__ -> __.image |> axiskeys(__, :ra) |> extrema |> abs(__[2] - __[1]))
		zip(__, ax)
		map() do (p, ax)
			plt.sca(ax)
			plt.xticks([])
			plt.yticks([])
			plt.title(p.label)
			@p p.image |>
				imshow_ax(__;
					norm=matplotlib.colors.SymLogNorm(vmin=max(0, __ |> skipnan |> minimum), linthresh=maximum(skipnan(__))/p.dynrange, vmax=quantile(skipnan(__), 1-3e-4)),
					interpolation=:gaussian
				)
			plt.gca().add_artist(ScalebarArtist(
				[
					(u"rad"/u"as", x -> x < 1e-4 ? f"{10^6*x:d} μas" : x < 0.1 ? f"{1000x:d} mas" : x < 500 ? f"{x}\"" : f"{x/60:d}'"),
					(80*u"rad"/u"as", x -> x < 1e-2 ? f"{round(x * 206265, sigdigits=1):d} AU" : x < 1e3 ? f"{x} pc" : f"{x / 1000:d} kpc")
				];
				loc="lower right", color=:w,
			))
			plt.gca().invert_xaxis()
		end
	end
	for (a, b) in zip(ax[1:end-1], ax[2:end])
		add_zoom_patch(a, b, :horizontal; color=:"0.5")
	end
	plt.tight_layout()
	plt.gcf()
end

# ╔═╡ a9cba3ab-275a-4250-a207-88506d796945
md"""
## Extracting original array
"""

# ╔═╡ 16ac4fe7-009b-4be5-8ce0-1c59385bec0a
md"""
In all the above example, we created an arbitrary coordinate grid, and projected the loaded image onto that grid.
This is often the right way to go, as image analysis tends to be most convenient in a specific fixed coordinate system --- typically, in projected coordinates for compact images.

However, sometimes it's useful or even necessary to work with the original image exactly as it's stored in the FITS file. `SkyImages.native_rect_image` function extracts the original rectangular image, putting it into a keyed array with axis keys derived on a best-effort basis. The image data is not copied, so this operation is efficient.
"""

# ╔═╡ 4d7e1f0f-327a-4316-8595-22d149a13939
md"""
Here is the same image as in the first example (Fermi LAT sky map), shown directly in pixel coordinates:
"""

# ╔═╡ cb93aba9-0fd6-4716-b9a4-dde69c365570
begin
	plt.figure()
	plt.imshow(SkyImages.native_rect_image(simg) |> permutedims)
	plt.gcf()
end

# ╔═╡ 0fe458e1-c021-4115-a941-8a1eeb67f63b
md"""
We can display the same array taking axis keys into account:
"""

# ╔═╡ 167eec1b-7f36-4ff0-9333-b086a09fb6b2
begin
	plt.figure()
	pcolormesh_ax(SkyImages.native_rect_image(simg))
	plt.gca().set_aspect(:equal)
	plt.gcf()
end

# ╔═╡ ec6d6033-0fd9-4b03-a82d-16f4258dd60b
md"""
Of course, this works with compact images as well:
"""

# ╔═╡ 8dcbf3c5-0de4-4235-8330-45886c562535
begin
	plt.figure()
	@p SkyImages.native_rect_image(small_img) |>
		@modify(xs -> xs .* rad2deg(1), axiskeys(__) |> Elements()) |>
		pcolormesh_ax(__; norm=matplotlib.colors.Normalize(vmin=0, vmax=quantile(skipnan(__), 1-3e-4)))
	plt.gca().set_aspect(:equal)
	plt.gcf()
end

# ╔═╡ 6bcc1479-30f6-44e1-9745-04217b9961d9
md"""
And with projected coordinates:
"""

# ╔═╡ 2086df08-3915-46d7-aafa-2700145265c4
begin
	plt.figure()
	@p SkyImages.native_rect_image(ProjectedCoords, small_img) |>
		@modify(xs -> xs .* rad2deg(1), axiskeys(__) |> Elements()) |>
		pcolormesh_ax(__; norm=matplotlib.colors.Normalize(vmin=0, vmax=quantile(skipnan(__), 1-3e-4)))
	plt.gca().set_aspect(:equal)
	plt.gcf()
end

# ╔═╡ cc5439ff-17b7-4a8f-98fb-3bf82bc1544f
md"""
Compare this plot to the same image in Galactic and equatorial coordinates from previous examples.
"""

# ╔═╡ c0ef8625-3d6a-457b-a07c-cc03ecb31a5d
md"""
## `AstroImages.jl` plot comparison
"""

# ╔═╡ ed5d3e9d-5b11-4c9e-b701-c02ece625d01


# ╔═╡ 2e9afe63-0048-4e91-b9b3-778957d3df00


# ╔═╡ 2f1996bb-514b-4666-ae47-c84fdca266b6


# ╔═╡ a0a45805-ec74-4df4-8c6f-d0a7d1534103


# ╔═╡ 8d66003b-3123-42a1-8f3e-4ca58eab0c32


# ╔═╡ 1265dd87-7395-4637-8b07-8dd5df20e43b
import AstroImages, Plots

# ╔═╡ 0962b3ae-4183-4e5c-8e5d-1c0817e548b7
[
	fermi_galplt;
	AstroImages.implot(AstroImages.load("data/intens_scaled_ait_144m_gt1000_psf3_gal_0p1.fits"))
]

# ╔═╡ 9fbb2cba-005c-4770-97d9-6a4e279cb2a7
[
	eagle_plt;
	AstroImages.implot(AstroImages.load("data/656nmos.fits"))
]

# ╔═╡ e28ee6fb-dbab-4884-867f-c60e9bf55e65
[
	let
		plt.figure()
		@p begin
			SkyImages.load("data/1228+126.u.2021_05_01.icn.fits")
			SkyImages.native_rect_image()
			@modify(xs -> (xs .- mean(xs)) .* (rad2deg(1) * 3600e3), axiskeys(__) |> Elements())
			imshow_ax(__; norm=matplotlib.colors.SymLogNorm(vmin=0, linthresh=maximum(skipnan(__))/1e3, vmax=maximum(skipnan(__))))
		end
		set_xylims((0 ± 30)^2; inv=:x)
		plt.colorbar()
		xylabels("Rel RA (mas)", "Rel Dec (mas)")
		plt.gcf()
	end,
	AstroImages.implot(AstroImages.load("data/1228+126.u.2021_05_01.icn.fits")[Z=1, AstroImages.Dim{:4}(1)])
]

# ╔═╡ 9515317f-3fdb-4b8c-8686-3d0b1a42ffe3
TableOfContents()

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
AccessorsExtra = "33016aad-b69d-45be-9359-82a41f556fd4"
AstroImages = "fe3fc30c-9b16-11e9-1c73-17dabf39f4ad"
AxisKeys = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
DataPipes = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
DisplayAs = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
Pkg = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
PyFormattedStrings = "5f89f4a4-a228-4886-b223-c468a82ed5b9"
PyPlotUtils = "5384e752-6c47-47b3-86ac-9d091b110b31"
RectiGrids = "8ac6971d-971d-971d-971d-971d5ab1a71a"
Revise = "295af30f-e4ad-537b-8983-00126c2a3abe"
SkipNan = "aed68c70-c8b0-4309-8cd1-d392a74f991a"
SkyCoords = "fc659fc5-75a3-5475-a2ea-3da92c065361"
SkyImages = "2d546a2e-713c-402d-bee5-ba90cc43b84b"
StatsBase = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"
UnitfulAngles = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"

[compat]
AccessorsExtra = "~0.1.9"
AstroImages = "~0.3.1"
AxisKeys = "~0.2.7"
DataPipes = "~0.2.17"
DisplayAs = "~0.1.6"
Plots = "~1.31.7"
PlutoUI = "~0.7.39"
PyFormattedStrings = "~0.1.10"
PyPlotUtils = "~0.1.17"
RectiGrids = "~0.1.13"
Revise = "~3.4.0"
SkipNan = "~0.2.0"
SkyCoords = "~1.0.1"
SkyImages = "~0.1.0"
StatsBase = "~0.33.21"
Unitful = "~1.11.0"
UnitfulAngles = "~0.6.2"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.8.0"
manifest_format = "2.0"
project_hash = "805b274e7f4e397bdcf18adf6c5d06a93f01029c"

[[deps.AbstractFFTs]]
deps = ["ChainRulesCore", "LinearAlgebra"]
git-tree-sha1 = "69f7020bd72f069c219b5e8c236c1fa90d2cb409"
uuid = "621f4979-c628-5d54-868e-fcf4e3e8185c"
version = "1.2.1"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "8eaf9f1b4921132a4cff3f36a1d9ba923b14a481"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.1.4"

[[deps.Accessors]]
deps = ["Compat", "CompositionsBase", "ConstructionBase", "Dates", "InverseFunctions", "LinearAlgebra", "MacroTools", "Requires", "Test"]
git-tree-sha1 = "8557017cfc7b58baea05a43ed35538857e6d35b4"
uuid = "7d9f7c33-5ae7-4f3b-8dc6-eff91059b697"
version = "0.1.19"

[[deps.AccessorsExtra]]
deps = ["Accessors", "ConstructionBase", "InverseFunctions", "Reexport", "Requires"]
git-tree-sha1 = "ef45a3c71f3a7e98a107ec66222e04250185c7bb"
uuid = "33016aad-b69d-45be-9359-82a41f556fd4"
version = "0.1.9"

[[deps.Adapt]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "195c5505521008abea5aee4f96930717958eac6f"
uuid = "79e6a3ab-5dfb-504d-930d-738a2a938a0e"
version = "3.4.0"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.ArrayInterface]]
deps = ["ArrayInterfaceCore", "Compat", "IfElse", "LinearAlgebra", "Static"]
git-tree-sha1 = "0582b5976fc76523f77056e888e454f0f7732596"
uuid = "4fba245c-0d91-5ea0-9b3e-6abc04ee57a9"
version = "6.0.22"

[[deps.ArrayInterfaceCore]]
deps = ["LinearAlgebra", "SparseArrays", "SuiteSparse"]
git-tree-sha1 = "40debc9f72d0511e12d817c7ca06a721b6423ba3"
uuid = "30b0a656-2188-435a-8636-2ec0e6a096e2"
version = "0.1.17"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.AstroAngles]]
git-tree-sha1 = "41621fa5ed5f7614b75eea8e0b3cfd967b284c87"
uuid = "5c4adb95-c1fc-4c53-b4ea-2a94080c53d2"
version = "0.1.3"

[[deps.AstroImages]]
deps = ["AbstractFFTs", "AstroAngles", "ColorSchemes", "DimensionalData", "FITSIO", "FileIO", "ImageAxes", "ImageBase", "ImageIO", "ImageShow", "MappedArrays", "PlotUtils", "Printf", "RecipesBase", "Statistics", "Tables", "UUIDs", "WCS"]
git-tree-sha1 = "796b7afebb3577d70c9c91caa530b32e42b353fc"
uuid = "fe3fc30c-9b16-11e9-1c73-17dabf39f4ad"
version = "0.3.1"

[[deps.AxisAlgorithms]]
deps = ["LinearAlgebra", "Random", "SparseArrays", "WoodburyMatrices"]
git-tree-sha1 = "66771c8d21c8ff5e3a93379480a2307ac36863f7"
uuid = "13072b0f-2c55-5437-9ae7-d433b7a33950"
version = "1.0.1"

[[deps.AxisArrays]]
deps = ["Dates", "IntervalSets", "IterTools", "RangeArrays"]
git-tree-sha1 = "1dd4d9f5beebac0c03446918741b1a03dc5e5788"
uuid = "39de3d68-74b9-583c-8d2d-e117c070f3a9"
version = "0.4.6"

[[deps.AxisKeys]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "IntervalSets", "InvertedIndices", "LazyStack", "LinearAlgebra", "NamedDims", "OffsetArrays", "Statistics", "StatsBase", "Tables"]
git-tree-sha1 = "88cc6419032d0e3ea69bc65d012aa82302774ab8"
uuid = "94b1ba4f-4ee9-5380-92f1-94cde586c3c5"
version = "0.2.7"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "19a35467a82e236ff51bc17a3a44b69ef35185a2"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+0"

[[deps.CEnum]]
git-tree-sha1 = "eb4cb44a499229b3b8426dcfb5dd85333951ff90"
uuid = "fa961155-64e5-5f13-b03f-caf6b980ea82"
version = "0.4.2"

[[deps.CFITSIO]]
deps = ["CFITSIO_jll"]
git-tree-sha1 = "8425c47db102577eefb93cb37b4480e750116b0d"
uuid = "3b1b4be9-1499-4b22-8d78-7db3344d1961"
version = "1.4.1"

[[deps.CFITSIO_jll]]
deps = ["Artifacts", "JLLWrappers", "LibCURL_jll", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "9c91a9358de42043c3101e3a29e60883345b0b39"
uuid = "b3e40c51-02ae-5482-8a39-3ace5868dcf4"
version = "4.0.0+0"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "4b859a208b2397a7a623a03449e4636bdb17bcf2"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.16.1+1"

[[deps.ChainRulesCore]]
deps = ["Compat", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "80ca332f6dcb2508adba68f22f551adb2d00a624"
uuid = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
version = "1.15.3"

[[deps.ChangesOfVariables]]
deps = ["ChainRulesCore", "LinearAlgebra", "Test"]
git-tree-sha1 = "38f7a08f19d8810338d4f5085211c7dfa5d5bdd8"
uuid = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
version = "0.1.4"

[[deps.CodeTracking]]
deps = ["InteractiveUtils", "UUIDs"]
git-tree-sha1 = "1833bda4a027f4b2a1c984baddcf755d77266818"
uuid = "da1fd8a2-8d9e-5ec2-8556-3022fb5608a2"
version = "1.1.0"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "ded953804d019afa9a3f98981d99b33e3db7b6da"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.0"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Random"]
git-tree-sha1 = "1fd869cc3875b57347f7027521f561cf46d1fcd8"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.19.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "eb7f0f8307f71fac7c606984ea5fb2817275d6e4"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.4"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "SpecialFunctions", "Statistics", "TensorCore"]
git-tree-sha1 = "d08c20eef1f2cbc6e60fd3612ac4340b89fea322"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.9.9"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "417b0ed7b8b838aa6ca0a87aadf1bb9eb111ce40"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.8"

[[deps.Compat]]
deps = ["Dates", "LinearAlgebra", "UUIDs"]
git-tree-sha1 = "5856d3031cdb1f3b2b6340dfdc66b6d9a149a374"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.2.0"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "0.5.2+0"

[[deps.CompositeTypes]]
git-tree-sha1 = "d5b014b216dc891e81fea299638e4c10c657b582"
uuid = "b152e2b5-7a66-4b01-a709-34e65c35f657"
version = "0.1.2"

[[deps.CompositionsBase]]
git-tree-sha1 = "455419f7e328a1a2493cabc6428d79e951349769"
uuid = "a33af91c-f02d-484b-be07-31d278c5ca2b"
version = "0.1.1"

[[deps.Conda]]
deps = ["Downloads", "JSON", "VersionParsing"]
git-tree-sha1 = "6e47d11ea2776bc5627421d59cdcc1296c058071"
uuid = "8f4d0f93-b110-5947-807f-2305c1781a2d"
version = "1.7.0"

[[deps.ConstructionBase]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "fb21ddd70a051d882a1686a5a550990bbe371a95"
uuid = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
version = "1.4.1"

[[deps.Contour]]
git-tree-sha1 = "d05d9e7b7aedff4e5b51a029dced05cfb6125781"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.2"

[[deps.CovarianceEstimation]]
deps = ["LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "3c8de95b4e932d76ec8960e12d681eba580e9674"
uuid = "587fd27a-f159-11e8-2dae-1979310e6154"
version = "0.2.8"

[[deps.DataAPI]]
git-tree-sha1 = "fb5f5316dd3fd4c5e7c30a24d50643b73e37cd40"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.10.0"

[[deps.DataPipes]]
deps = ["Accessors", "SplitApplyCombine"]
git-tree-sha1 = "ab6b5bf476e9111b0166cc3f8373638204d7fafd"
uuid = "02685ad9-2d12-40c3-9f73-c6aeda6a7ff5"
version = "0.2.17"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "d1fff3a548102f48987a52a2e0d114fa97d730f0"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.13"

[[deps.DataValueInterfaces]]
git-tree-sha1 = "bfc1187b79289637fa0ef6d4436ebdfe6905cbd6"
uuid = "e2d170a0-9d28-54be-80f0-106bbe20a464"
version = "1.0.0"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"

[[deps.Dictionaries]]
deps = ["Indexing", "Random", "Serialization"]
git-tree-sha1 = "96dc5c5c8994be519ee3420953c931c55657a3f2"
uuid = "85a47980-9c8c-11e8-2b9f-f7ca1fa99fb4"
version = "0.3.24"

[[deps.DimensionalData]]
deps = ["Adapt", "ArrayInterface", "ConstructionBase", "Dates", "Extents", "IntervalSets", "LinearAlgebra", "Random", "RecipesBase", "SparseArrays", "Statistics", "Tables"]
git-tree-sha1 = "ee34aa915d13091b0357e4746fe825aac616a30a"
uuid = "0703355e-b756-11e9-17c0-8b28908087d0"
version = "0.20.11"

[[deps.DirectionalStatistics]]
deps = ["AccessorsExtra", "IntervalSets", "InverseFunctions", "LinearAlgebra", "Statistics", "StatsBase"]
git-tree-sha1 = "156365de4369a6cf587d0d59ce52fe688f2b5f92"
uuid = "e814f24e-44b0-11e9-2fd5-aba2b6113d95"
version = "0.1.19"

[[deps.DisplayAs]]
git-tree-sha1 = "43c017d5dd3a48d56486055973f443f8a39bb6d9"
uuid = "0b91fe84-8a4c-11e9-3e1d-67c38462b6d6"
version = "0.1.6"

[[deps.Distributed]]
deps = ["Random", "Serialization", "Sockets"]
uuid = "8ba89e20-285c-5b6f-9357-94700520ee1b"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "5158c2b41018c5f7eb1470d558127ac274eca0c9"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.1"

[[deps.DomainSets]]
deps = ["CompositeTypes", "IntervalSets", "LinearAlgebra", "Random", "StaticArrays", "Statistics"]
git-tree-sha1 = "dc45fbbe91d6d17a8e187abad39fb45963d97388"
uuid = "5b8099bc-c8ec-5219-889f-1d9e522a28bf"
version = "0.5.13"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EarCut_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3f3a2501fa7236e9b911e0f7a588c657e822bb6d"
uuid = "5ae413db-bbd1-5e63-b57d-d24a61df00f5"
version = "2.2.3+0"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bad72f730e9e91c08d9427d5e8db95478a3c323d"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.4.8+0"

[[deps.Extents]]
git-tree-sha1 = "5e1e4c53fa39afe63a7d356e30452249365fba99"
uuid = "411431e0-e8b7-467b-b5e0-f676ba4f2910"
version = "0.1.1"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "Pkg", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "ccd479984c7838684b3ac204b716c89955c76623"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.2+0"

[[deps.FITSIO]]
deps = ["CFITSIO", "Printf", "Reexport", "Tables"]
git-tree-sha1 = "e6033823834ec0070125120d4d4a1234f1826a47"
uuid = "525bcba6-941b-5504-bd06-fd0dc1a4d2eb"
version = "0.16.12"

[[deps.FileIO]]
deps = ["Pkg", "Requires", "UUIDs"]
git-tree-sha1 = "94f5101b96d2d968ace56f7f2db19d0a5f592e28"
uuid = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
version = "1.15.0"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Formatting]]
deps = ["Printf"]
git-tree-sha1 = "8339d61043228fdd3eb658d86c926cb282ae72a8"
uuid = "59287772-0a20-5a39-b81b-1366585eb4c0"
version = "0.4.2"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "87eb71354d8ec1a96d4a7636bd57a7347dde3ef9"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.10.4+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Pkg", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "d972031d28c8c8d9d7b41a536ad7bb0c2579caca"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.8+0"

[[deps.GR]]
deps = ["Base64", "DelimitedFiles", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Printf", "Random", "RelocatableFolders", "Serialization", "Sockets", "Test", "UUIDs"]
git-tree-sha1 = "cf0a9940f250dc3cb6cc6c6821b4bf8a4286cf9c"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.66.2"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Pkg", "Qt5Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "2d908286d120c584abbe7621756c341707096ba4"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.66.2+0"

[[deps.GeoInterface]]
deps = ["Extents"]
git-tree-sha1 = "fb28b5dc239d0174d7297310ef7b84a11804dfab"
uuid = "cf35fbd7-0cd7-5166-be24-54bfbe79505f"
version = "1.0.1"

[[deps.GeometryBasics]]
deps = ["EarCut_jll", "GeoInterface", "IterTools", "LinearAlgebra", "StaticArrays", "StructArrays", "Tables"]
git-tree-sha1 = "a7a97895780dab1085a97769316aa348830dc991"
uuid = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
version = "0.4.3"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "a32d672ac2c967f3deb8a81d828afc739c838a06"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.68.3+2"

[[deps.Graphics]]
deps = ["Colors", "LinearAlgebra", "NaNMath"]
git-tree-sha1 = "d61890399bc535850c4bf08e4e0d3a7ad0f21cbd"
uuid = "a2bd30eb-e257-5431-a919-1863eab51364"
version = "1.1.2"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "Dates", "IniFile", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "59ba44e0aa49b87a8c7a8920ec76f8afe87ed502"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.3.3"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Healpix]]
deps = ["CFITSIO", "LazyArtifacts", "Libsharp", "LinearAlgebra", "Pkg", "Printf", "RecipesBase", "StaticArrays", "Test"]
git-tree-sha1 = "e4583be5e3fc6436d3eb00f163c799cb46712aa6"
uuid = "9f4e344d-96bc-545a-84a3-ae6b9e1b672b"
version = "4.0.1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "8d511d5b81240fc8e6802386302675bdf47737b9"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.4"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "c47c5fa4c5308f27ccaac35504858d8914e102f9"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.4"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "f7be53659ab06ddc986428d3a9dcc95f6fa6705a"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.2"

[[deps.IfElse]]
git-tree-sha1 = "debdd00ffef04665ccbb3e150747a77560e8fad1"
uuid = "615f187c-cbe4-4ef1-ba3b-2fcf58d6d173"
version = "0.1.1"

[[deps.ImageAxes]]
deps = ["AxisArrays", "ImageBase", "ImageCore", "Reexport", "SimpleTraits"]
git-tree-sha1 = "c54b581a83008dc7f292e205f4c409ab5caa0f04"
uuid = "2803e5a7-5153-5ecf-9a86-9b4c37f5f5ac"
version = "0.6.10"

[[deps.ImageBase]]
deps = ["ImageCore", "Reexport"]
git-tree-sha1 = "b51bb8cae22c66d0f6357e3bcb6363145ef20835"
uuid = "c817782e-172a-44cc-b673-b171935fbb9e"
version = "0.1.5"

[[deps.ImageCore]]
deps = ["AbstractFFTs", "ColorVectorSpace", "Colors", "FixedPointNumbers", "Graphics", "MappedArrays", "MosaicViews", "OffsetArrays", "PaddedViews", "Reexport"]
git-tree-sha1 = "acf614720ef026d38400b3817614c45882d75500"
uuid = "a09fc81d-aa75-5fe9-8630-4744c3626534"
version = "0.9.4"

[[deps.ImageIO]]
deps = ["FileIO", "IndirectArrays", "JpegTurbo", "LazyModules", "Netpbm", "OpenEXR", "PNGFiles", "QOI", "Sixel", "TiffImages", "UUIDs"]
git-tree-sha1 = "342f789fd041a55166764c351da1710db97ce0e0"
uuid = "82e4d734-157c-48bb-816b-45c225c6df19"
version = "0.6.6"

[[deps.ImageShow]]
deps = ["Base64", "FileIO", "ImageBase", "ImageCore", "OffsetArrays", "StackViews"]
git-tree-sha1 = "b563cf9ae75a635592fc73d3eb78b86220e55bd8"
uuid = "4e3cecfd-b093-5904-9786-8bbb286a6a31"
version = "0.3.6"

[[deps.Imath_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "87f7662e03a649cffa2e05bf19c303e168732d3e"
uuid = "905a6f67-0a94-5f89-b386-d35d92009cd1"
version = "3.1.2+0"

[[deps.Indexing]]
git-tree-sha1 = "ce1566720fd6b19ff3411404d4b977acd4814f9f"
uuid = "313cdc1a-70c2-5d6a-ae34-0150d3930a38"
version = "1.1.1"

[[deps.IndirectArrays]]
git-tree-sha1 = "012e604e1c7458645cb8b436f8fba789a51b257f"
uuid = "9b13fd28-a010-5f03-acff-a1bbcff69959"
version = "1.0.0"

[[deps.Inflate]]
git-tree-sha1 = "5cd07aab533df5170988219191dfad0519391428"
uuid = "d25df0c9-e2be-5dd7-82c8-3ad0b3e990b9"
version = "0.1.3"

[[deps.IniFile]]
git-tree-sha1 = "f550e6e32074c939295eb5ea6de31849ac2c9625"
uuid = "83e8ac13-25f8-5344-8a64-a9f2b223428f"
version = "0.5.1"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.Interpolations]]
deps = ["Adapt", "AxisAlgorithms", "ChainRulesCore", "LinearAlgebra", "OffsetArrays", "Random", "Ratios", "Requires", "SharedArrays", "SparseArrays", "StaticArrays", "WoodburyMatrices"]
git-tree-sha1 = "64f138f9453a018c8f3562e7bae54edc059af249"
uuid = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
version = "0.14.4"

[[deps.IntervalSets]]
deps = ["Dates", "Random", "Statistics"]
git-tree-sha1 = "076bb0da51a8c8d1229936a1af7bdfacd65037e1"
uuid = "8197267c-284f-5f27-9208-e0e47529a953"
version = "0.7.2"

[[deps.InverseFunctions]]
deps = ["Test"]
git-tree-sha1 = "b3364212fb5d870f724876ffcd34dd8ec6d98918"
uuid = "3587e190-3f89-42d0-90ee-14403ec27112"
version = "0.1.7"

[[deps.InvertedIndices]]
git-tree-sha1 = "bee5f1ef5bf65df56bdd2e40447590b272a5471f"
uuid = "41ab1584-1d38-5bbf-9106-f11c6c58b48f"
version = "1.1.0"

[[deps.IrrationalConstants]]
git-tree-sha1 = "7fd44fd4ff43fc60815f8e764c0f352b83c49151"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.1.1"

[[deps.IterTools]]
git-tree-sha1 = "fa6287a4469f5e048d763df38279ee729fbd44e5"
uuid = "c8e1da08-722c-5040-9ed9-7db0dc04731e"
version = "1.4.0"

[[deps.IteratorInterfaceExtensions]]
git-tree-sha1 = "a3f24677c21f5bbe9d2a714f95dcd58337fb2856"
uuid = "82899510-4779-5014-852e-03e436cf321d"
version = "1.0.0"

[[deps.JLLWrappers]]
deps = ["Preferences"]
git-tree-sha1 = "abc9885a7ca2052a736a600f7fa66209f96506e1"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.4.1"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "3c837543ddb02250ef42f4738347454f95079d4e"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.3"

[[deps.JpegTurbo]]
deps = ["CEnum", "FileIO", "ImageCore", "JpegTurbo_jll", "TOML"]
git-tree-sha1 = "a77b273f1ddec645d1b7c4fd5fb98c8f90ad10a5"
uuid = "b835a17e-a41a-41e7-81f0-2f016b05efe0"
version = "0.1.1"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b53380851c6e6664204efb2e62cd24fa5c47e4ba"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "2.1.2+0"

[[deps.JuliaInterpreter]]
deps = ["CodeTracking", "InteractiveUtils", "Random", "UUIDs"]
git-tree-sha1 = "0f960b1404abb0b244c1ece579a0ec78d056a5d1"
uuid = "aa1ae85d-cabe-5617-a682-6adf51b2e16a"
version = "0.9.15"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "f2355693d6778a178ade15952b7ac47a4ff97996"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.0"

[[deps.Latexify]]
deps = ["Formatting", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "Printf", "Requires"]
git-tree-sha1 = "1a43be956d433b5d0321197150c2f94e16c0aaa0"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.15.16"

[[deps.LazyArtifacts]]
deps = ["Artifacts", "Pkg"]
uuid = "4af54fe1-eca0-43a8-85a7-787d91b784e3"

[[deps.LazyModules]]
git-tree-sha1 = "a560dd966b386ac9ae60bdd3a3d3a326062d3c3e"
uuid = "8cdb02fc-e678-4876-92c5-9defec4f444e"
version = "0.3.1"

[[deps.LazyStack]]
deps = ["ChainRulesCore", "LinearAlgebra", "NamedDims", "OffsetArrays"]
git-tree-sha1 = "2eb4a5bf2eb0519ebf40c797ba5637d327863637"
uuid = "1fad7336-0346-5a1a-a56f-a06ba010965b"
version = "0.0.8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.3"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "7.84.0+0"

[[deps.LibGit2]]
deps = ["Base64", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.10.2+0"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "7739f837d6447403596a75d19ed01fd08d6f56bf"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.3.0+3"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "42b62845d70a619f063a7da093d995ec8e15e778"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.16.1+1"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9c30530bf0effd46e15e0fdcf2b8636e78cbbd73"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.35.0+0"

[[deps.Libsharp]]
deps = ["Libdl", "libsharp2_jll"]
git-tree-sha1 = "e09051a3f95b83091fc9b7a26e6c585c6691b5bc"
uuid = "ac8d63fe-4615-43ae-9860-9cd4a3820532"
version = "0.2.0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "Pkg", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "3eb79b0ca5764d4799c06699573fd8f533259713"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.4.0+0"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "7f3efec06033682db852f8b3bc3c1d2b0a0ab066"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.36.0+0"

[[deps.LinearAlgebra]]
deps = ["Libdl", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["ChainRulesCore", "ChangesOfVariables", "DocStringExtensions", "InverseFunctions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "94d9c52ca447e23eac0c0f074effbcd38830deb5"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.18"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "5d4d2d9904227b8bd66386c1138cf4d5ffa826bf"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "0.4.9"

[[deps.LoweredCodeUtils]]
deps = ["JuliaInterpreter"]
git-tree-sha1 = "dedbebe234e06e1ddad435f5c6f4b85cd8ce55f7"
uuid = "6f1432cf-f94c-5a45-995e-cdbf5db27b0b"
version = "2.2.2"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "3d3e902b31198a27340d0bf00d6ac452866021cf"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.9"

[[deps.MappedArrays]]
git-tree-sha1 = "e8b359ef06ec72e8c030463fe02efe5527ee5142"
uuid = "dbb5928d-eab1-5f90-85c2-b9b0edb7c900"
version = "0.4.1"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "Random", "Sockets"]
git-tree-sha1 = "2f0be365951a88dfb084f754005177e6dfb00ed0"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.4"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.0+0"

[[deps.Measures]]
git-tree-sha1 = "e498ddeee6f9fdb4551ce855a46f54dbd900245f"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.1"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "bf210ce90b6c9eed32d25dbcae1ebc565df2687f"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.0.2"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MosaicViews]]
deps = ["MappedArrays", "OffsetArrays", "PaddedViews", "StackViews"]
git-tree-sha1 = "b34e3bc3ca7c94914418637cb10cc4d1d80d877d"
uuid = "e94cdb99-869f-56ef-bcf0-1ae2bcbe0389"
version = "0.3.3"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2022.2.1"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "a7c3d1da1189a1c2fe843a3bfa04d18d20eb3211"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.1"

[[deps.NamedDims]]
deps = ["AbstractFFTs", "ChainRulesCore", "CovarianceEstimation", "LinearAlgebra", "Pkg", "Requires", "Statistics"]
git-tree-sha1 = "f39537cbe1cf4f407e65bdf7aca6b04f5877fbb1"
uuid = "356022a1-0364-5f58-8944-0da4b18d706f"
version = "1.1.0"

[[deps.Netpbm]]
deps = ["FileIO", "ImageCore"]
git-tree-sha1 = "18efc06f6ec36a8b801b23f076e3c6ac7c3bf153"
uuid = "f09324ee-3d7c-5217-9330-fc30815ba969"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.NonNegLeastSquares]]
deps = ["Distributed", "LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "1271344271ffae97e2855b0287356e6ea5c221cc"
uuid = "b7351bd1-99d9-5c5d-8786-f205a815c4d7"
version = "0.4.0"

[[deps.OffsetArrays]]
deps = ["Adapt"]
git-tree-sha1 = "1ea784113a6aa054c5ebd95945fa5e52c2f378e7"
uuid = "6fe1bfb0-de20-5000-8ca7-80f57d26f881"
version = "1.12.7"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.20+0"

[[deps.OpenEXR]]
deps = ["Colors", "FileIO", "OpenEXR_jll"]
git-tree-sha1 = "327f53360fdb54df7ecd01e96ef1983536d1e633"
uuid = "52e1d378-f018-4a11-a4be-720524705ac7"
version = "0.3.2"

[[deps.OpenEXR_jll]]
deps = ["Artifacts", "Imath_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "923319661e9a22712f24596ce81c54fc0366f304"
uuid = "18a262bb-aa17-5467-a713-aee519bc75cb"
version = "3.1.1+0"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+0"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e60321e3f2616584ff98f0a4f18d98ae6f89bbb3"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "1.1.17+0"

[[deps.OpenSpecFun_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "13652491f6856acfd2db29360e1bbcd4565d04f1"
uuid = "efe28fd5-8261-553b-a9e1-b2916fc3738e"
version = "0.5.5+0"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "85f8e6578bf1f9ee0d11e7bb1b1456435479d47c"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.4.1"

[[deps.PCRE_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b2a7af664e098055a7529ad1a900ded962bca488"
uuid = "2f80f16e-611a-54ab-bc61-aa92de5b98fc"
version = "8.44.0+0"

[[deps.PNGFiles]]
deps = ["Base64", "CEnum", "ImageCore", "IndirectArrays", "OffsetArrays", "libpng_jll"]
git-tree-sha1 = "e925a64b8585aa9f4e3047b8d2cdc3f0e79fd4e4"
uuid = "f57f5aa1-a3ce-4bc8-8ab9-96f992907883"
version = "0.3.16"

[[deps.PaddedViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "03a7a85b76381a3d04c7a1656039197e70eda03d"
uuid = "5432bcbf-9aad-5242-b902-cca2824c8663"
version = "0.5.11"

[[deps.Parsers]]
deps = ["Dates"]
git-tree-sha1 = "3d5bf43e3e8b412656404ed9466f1dcbf7c50269"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.4.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "b4f5d02549a10e20780a24fce72bea96b6329e29"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.40.1+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.8.0"

[[deps.PkgVersion]]
deps = ["Pkg"]
git-tree-sha1 = "f6cf8e7944e50901594838951729a1861e668cb8"
uuid = "eebad327-c553-4316-9ea0-9fa01ccd7688"
version = "0.3.2"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "8162b2f8547bc23876edd0c5181b27702ae58dce"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.0.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "9888e59493658e476d3073f1ce24348bdc086660"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.3.0"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "GeometryBasics", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "Unzip"]
git-tree-sha1 = "a19652399f43938413340b2068e11e55caa46b65"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.31.7"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "Markdown", "Random", "Reexport", "UUIDs"]
git-tree-sha1 = "8d1f54886b9037091edf146b517989fc4a09efec"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.39"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "47e5f437cc0e7ef2ce8406ce1e7e24d44915f88d"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.3.0"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.ProgressMeter]]
deps = ["Distributed", "Printf"]
git-tree-sha1 = "d7a7aef8f8f2d537104f170139553b14dfe39fe9"
uuid = "92933f4c-e287-5a05-a399-4b506db050ca"
version = "1.7.2"

[[deps.PyCall]]
deps = ["Conda", "Dates", "Libdl", "LinearAlgebra", "MacroTools", "Serialization", "VersionParsing"]
git-tree-sha1 = "53b8b07b721b77144a0fbbbc2675222ebf40a02d"
uuid = "438e738f-606a-5dbb-bf0a-cddfbfd45ab0"
version = "1.94.1"

[[deps.PyFormattedStrings]]
deps = ["Printf"]
git-tree-sha1 = "427aadbcb003bd2b39c87caa5177dded65f81ddf"
uuid = "5f89f4a4-a228-4886-b223-c468a82ed5b9"
version = "0.1.10"

[[deps.PyPlot]]
deps = ["Colors", "LaTeXStrings", "PyCall", "Sockets", "Test", "VersionParsing"]
git-tree-sha1 = "f9d953684d4d21e947cb6d642db18853d43cb027"
uuid = "d330b81b-6aea-500a-939a-2ce795aea3ee"
version = "2.11.0"

[[deps.PyPlotUtils]]
deps = ["Accessors", "AxisKeys", "Colors", "DataPipes", "DirectionalStatistics", "DomainSets", "IntervalSets", "LinearAlgebra", "NonNegLeastSquares", "PyCall", "PyPlot", "StatsBase", "Unitful"]
git-tree-sha1 = "37b811018ec1ebde8e4698c118b72f90ed2ec275"
uuid = "5384e752-6c47-47b3-86ac-9d091b110b31"
version = "0.1.17"

[[deps.QOI]]
deps = ["ColorTypes", "FileIO", "FixedPointNumbers"]
git-tree-sha1 = "18e8f4d1426e965c7b532ddd260599e1510d26ce"
uuid = "4b34888f-f399-49d4-9bb3-47ed5cae4e65"
version = "1.0.0"

[[deps.Qt5Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Pkg", "Xorg_libXext_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "xkbcommon_jll"]
git-tree-sha1 = "c6c0f690d0cc7caddb74cef7aa847b824a16b256"
uuid = "ea2cea3b-5b76-57ae-a6ef-0a8af62496e1"
version = "5.15.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA", "Serialization"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RangeArrays]]
git-tree-sha1 = "b9039e93773ddcfc828f12aadf7115b4b4d225f5"
uuid = "b3c3ace0-ae52-54e7-9d0b-2c1406fd6b9d"
version = "0.3.2"

[[deps.Ratios]]
deps = ["Requires"]
git-tree-sha1 = "dc84268fe0e3335a62e315a3a7cf2afa7178a734"
uuid = "c84ed2f1-dad5-54f0-aa8e-dbefe2724439"
version = "0.4.3"

[[deps.RecipesBase]]
git-tree-sha1 = "6bf3f380ff52ce0832ddd3a2a7b9538ed1bcca7d"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.2.1"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "RecipesBase"]
git-tree-sha1 = "e7eac76a958f8664f2718508435d058168c7953d"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.3"

[[deps.RectiGrids]]
deps = ["AxisKeys", "ConstructionBase", "Random", "StaticArraysCore"]
git-tree-sha1 = "940a23a7472b7352ff0ebbb661da8bbb5d7f932f"
uuid = "8ac6971d-971d-971d-971d-971d5ab1a71a"
version = "0.1.13"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "22c5201127d7b243b9ee1de3b43c408879dff60f"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "0.3.0"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.Revise]]
deps = ["CodeTracking", "Distributed", "FileWatching", "JuliaInterpreter", "LibGit2", "LoweredCodeUtils", "OrderedCollections", "Pkg", "REPL", "Requires", "UUIDs", "Unicode"]
git-tree-sha1 = "dad726963ecea2d8a81e26286f625aee09a91b7c"
uuid = "295af30f-e4ad-537b-8983-00126c2a3abe"
version = "3.4.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "f94f779c94e58bf9ea243e77a37e16d9de9126bd"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.1.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.SharedArrays]]
deps = ["Distributed", "Mmap", "Random", "Serialization"]
uuid = "1a1011a3-84de-559e-8e89-a11a2f7dc383"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.SimpleTraits]]
deps = ["InteractiveUtils", "MacroTools"]
git-tree-sha1 = "5d7e3f4e11935503d3ecaf7186eac40602e7d231"
uuid = "699a6c99-e7fa-54fc-8d76-47d257e15c1d"
version = "0.9.4"

[[deps.Sixel]]
deps = ["Dates", "FileIO", "ImageCore", "IndirectArrays", "OffsetArrays", "REPL", "libsixel_jll"]
git-tree-sha1 = "8fb59825be681d451c246a795117f317ecbcaa28"
uuid = "45858cf5-a6b0-47a3-bbea-62219f50df47"
version = "0.1.2"

[[deps.SkipNan]]
git-tree-sha1 = "b07be17ad1c4dd3e2d11aff5aa06157838ee6a6a"
uuid = "aed68c70-c8b0-4309-8cd1-d392a74f991a"
version = "0.2.0"

[[deps.SkyCoords]]
deps = ["AstroAngles", "ConstructionBase", "StaticArrays"]
git-tree-sha1 = "3302abfbde42db0c029e86d6155775b474a536d3"
uuid = "fc659fc5-75a3-5475-a2ea-3da92c065361"
version = "1.0.1"

[[deps.SkyImages]]
deps = ["AccessorsExtra", "AxisKeys", "ConstructionBase", "DataPipes", "DirectionalStatistics", "FITSIO", "Healpix", "Interpolations", "IntervalSets", "RectiGrids", "SkyCoords", "WCS"]
path = "../../home/aplavin/.julia/dev/SkyImages"
uuid = "2d546a2e-713c-402d-bee5-ba90cc43b84b"
version = "0.1.4"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "b3363d7460f7d098ca0912c69b082f75625d7508"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.0.1"

[[deps.SparseArrays]]
deps = ["LinearAlgebra", "Random"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"

[[deps.SpecialFunctions]]
deps = ["ChainRulesCore", "IrrationalConstants", "LogExpFunctions", "OpenLibm_jll", "OpenSpecFun_jll"]
git-tree-sha1 = "d75bda01f8c31ebb72df80a46c88b25d1c79c56d"
uuid = "276daf66-3868-5448-9aa4-cd146d93841b"
version = "2.1.7"

[[deps.SplitApplyCombine]]
deps = ["Dictionaries", "Indexing"]
git-tree-sha1 = "48f393b0231516850e39f6c756970e7ca8b77045"
uuid = "03a91e81-4c3e-53e1-a0a4-9c0c8f19dd66"
version = "1.2.2"

[[deps.StackViews]]
deps = ["OffsetArrays"]
git-tree-sha1 = "46e589465204cd0c08b4bd97385e4fa79a0c770c"
uuid = "cae243ae-269e-4f55-b966-ac2d0dc13c15"
version = "0.1.1"

[[deps.Static]]
deps = ["IfElse"]
git-tree-sha1 = "f94f9d627ba3f91e41a815b9f9f977d729e2e06f"
uuid = "aedffcd0-7271-4cad-89d0-dc628f76c6d3"
version = "0.7.6"

[[deps.StaticArrays]]
deps = ["LinearAlgebra", "Random", "StaticArraysCore", "Statistics"]
git-tree-sha1 = "dfec37b90740e3b9aa5dc2613892a3fc155c3b42"
uuid = "90137ffa-7385-5640-81b9-e52037218182"
version = "1.5.6"

[[deps.StaticArraysCore]]
git-tree-sha1 = "ec2bd695e905a3c755b33026954b119ea17f2d22"
uuid = "1e83bf80-4336-4d27-bf5d-d5a4f845583c"
version = "1.3.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "f9af7f195fb13589dd2e2d57fdb401717d2eb1f6"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.5.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "d1bf48bfcc554a3761a133fe3a9bb01488e06916"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.33.21"

[[deps.StructArrays]]
deps = ["Adapt", "DataAPI", "StaticArraysCore", "Tables"]
git-tree-sha1 = "8c6ac65ec9ab781af05b08ff305ddc727c25f680"
uuid = "09ab397b-f2b6-538f-b94a-2f83cf4a842a"
version = "0.6.12"

[[deps.SuiteSparse]]
deps = ["Libdl", "LinearAlgebra", "Serialization", "SparseArrays"]
uuid = "4607b0f0-06f3-5cda-b6b1-a6196a1729e9"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.0"

[[deps.TableTraits]]
deps = ["IteratorInterfaceExtensions"]
git-tree-sha1 = "c06b2f539df1c6efa794486abfb6ed2022561a39"
uuid = "3783bdb8-4a98-5b6b-af9a-565f29a5fe9c"
version = "1.0.1"

[[deps.Tables]]
deps = ["DataAPI", "DataValueInterfaces", "IteratorInterfaceExtensions", "LinearAlgebra", "OrderedCollections", "TableTraits", "Test"]
git-tree-sha1 = "5ce79ce186cc678bbb5c5681ca3379d1ddae11a1"
uuid = "bd369af6-aec1-5ad0-b16a-f7cc5008161c"
version = "1.7.0"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TiffImages]]
deps = ["ColorTypes", "DataStructures", "DocStringExtensions", "FileIO", "FixedPointNumbers", "IndirectArrays", "Inflate", "Mmap", "OffsetArrays", "PkgVersion", "ProgressMeter", "UUIDs"]
git-tree-sha1 = "70e6d2da9210371c927176cb7a56d41ef1260db7"
uuid = "731e570b-9d59-4bfa-96dc-6df516fadf69"
version = "0.6.1"

[[deps.TranscodingStreams]]
deps = ["Random", "Test"]
git-tree-sha1 = "ed5d390c7addb70e90fd1eb783dcb9897922cbfa"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.9.8"

[[deps.Tricks]]
git-tree-sha1 = "6bac775f2d42a611cdfcd1fb217ee719630c4175"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.6"

[[deps.URIs]]
git-tree-sha1 = "e59ecc5a41b000fa94423a578d29290c7266fc10"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.4.0"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["ConstructionBase", "Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "b649200e887a487468b71821e2644382699f1b0f"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.11.0"

[[deps.UnitfulAngles]]
deps = ["Dates", "Unitful"]
git-tree-sha1 = "d6cfdb6ddeb388af1aea38d2b9905fa014d92d98"
uuid = "6fb2a4bd-7999-5318-a3b2-8ad61056cd98"
version = "0.6.2"

[[deps.Unzip]]
git-tree-sha1 = "34db80951901073501137bdbc3d5a8e7bbd06670"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.1.2"

[[deps.VersionParsing]]
git-tree-sha1 = "58d6e80b4ee071f5efd07fda82cb9fbe17200868"
uuid = "81def892-9a0e-5fdd-b105-ffc91e053289"
version = "1.3.0"

[[deps.WCS]]
deps = ["WCS_jll"]
git-tree-sha1 = "905da53f09c677385c786135b0774de48227e219"
uuid = "15f3aee2-9e10-537f-b834-a6fb8bdb944d"
version = "0.6.1"

[[deps.WCS_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "947bfa11fcd65dac9e9b2e963504fba6b4971d31"
uuid = "550c8279-ae0e-5d1b-948f-937f2608a23e"
version = "7.7.0+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "3e61f0b86f90dacb0bc0e73a0c5a83f6a8636e23"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.19.0+0"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4528479aa01ee1b3b4cd0e6faef0e04cf16466da"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.25.0+0"

[[deps.WoodburyMatrices]]
deps = ["LinearAlgebra", "SparseArrays"]
git-tree-sha1 = "de67fa59e33ad156a590055375a30b23c40299d3"
uuid = "efce3f68-66dc-5838-9240-27a6d6f5f9b6"
version = "0.5.5"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "58443b63fb7e465a8a7210828c91c08b92132dff"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.9.14+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "5be649d550f3f4b95308bf0183b82e2582876527"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.6.9+4"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4e490d5c960c314f33885790ed410ff3a94ce67e"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.9+4"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fe47bd2247248125c428978740e18a681372dd4"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.3+4"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "6783737e45d3c59a4a4c4091f5f88cdcf0908cbb"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.0+3"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "daf17f441228e7a3833846cd048892861cff16d6"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.13.0+3"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "926af861744212db0eb001d9e40b5d16292080b2"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.0+4"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "4bcbf660f6c2e714f87e960a171b119d06ee163b"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.2+4"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "5c8424f8a67c3f2209646d4425f3d415fee5931d"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.27.0+4"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "79c31e7844f6ecf779705fbc12146eb190b7d845"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.4.0+3"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.12+3"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e45044cd873ded54b6a5bac0eb5c971392cf1927"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.2+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl", "OpenBLAS_jll"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.1.1+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "94d180a6d2b5e55e447e2d27a29ed04fe79eb30c"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.38+0"

[[deps.libsharp2_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9a7b8663b021ba776c89d9de7cdb2069cc27c00a"
uuid = "180207a7-b08e-5162-af94-7d62a04fe081"
version = "1.0.2+1"

[[deps.libsixel_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "78736dab31ae7a53540a6b752efc61f77b304c5b"
uuid = "075b6546-f08a-558a-be8f-8157d0f608a5"
version = "1.8.6+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.48.0+0"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+0"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9ebfc140cc56e8c2156a15ceac2f0302e327ac0a"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+0"
"""

# ╔═╡ Cell order:
# ╟─d27f13ff-72de-4136-b145-d8a973f23047
# ╟─a0a180c0-15c5-408d-82a2-3e8d5f1ff484
# ╟─aa6da229-640c-4eae-a0fe-9bd2742fa849
# ╟─b02f575f-bbd5-42c0-aab3-37415b66bd76
# ╟─72a56fff-84ab-4fdb-882f-e0ee6efada45
# ╠═3d581238-56b7-44b9-a9c3-9af9530ae136
# ╟─41f3ff07-8ebe-4224-ad1b-9342ee4e7e07
# ╠═7de9fc6e-f652-48bc-babe-af0e56a8c427
# ╟─b2649127-99b3-4ee4-aa54-d11f46f2c04e
# ╟─f26a1360-ba03-4df3-8ba0-ffd707a34bfb
# ╠═e42ac814-ca50-4464-8a19-59695b1f9ef0
# ╟─ec6d0165-18e5-41f6-91e0-382840cde0f0
# ╠═b0c2a790-231a-4ebc-a755-627cca858325
# ╟─170f4c04-0961-44d5-876e-9e9f51834e32
# ╠═4285e435-549e-487f-9013-95f8d2958cbe
# ╟─de270142-3c0c-45b9-968f-0ef68c5074b0
# ╠═62fd3d3c-d8a6-48e4-8d03-c27d2ef581f0
# ╟─6ed9fc21-08b6-4f39-b9c3-29abbcd3cee1
# ╠═17671ee5-bda8-4ffd-b1ef-dbb98eb7339a
# ╠═16bc419e-1c2a-4617-ba06-fcaa943066ca
# ╟─6cb2f138-5bba-406a-a505-4c60645d1dcb
# ╟─e052b45d-ef9e-41de-aeee-c0fce166f4f6
# ╠═6c7b4e2d-7f67-428a-8d22-c2bfb0a96a42
# ╟─a6d02948-0e29-4fa6-a2b6-a8f7d605fb97
# ╟─9905274c-40b2-4f8c-bf87-9f0f939699bf
# ╟─83946094-3f59-4d98-92cc-b79558f6a504
# ╠═7ac37418-514a-4b41-beae-2bec303c9737
# ╟─8efcbb46-9685-479f-96b9-816aacd2a851
# ╠═59c20159-4614-47d0-98ab-bf6f4991caed
# ╟─e16e262d-8366-41cc-ba5c-085188dbe3c7
# ╠═e9b9da8e-80c3-47dc-b1e4-db6d490e7d12
# ╟─d177d489-072c-4dcd-a72f-0eeab7fedc2a
# ╟─3b5fbd3c-cfff-4f21-9736-589e9c01e656
# ╟─a5ed3228-f068-4784-9d7a-61530cbd40e6
# ╠═6eaf40d7-5513-4a03-a77d-a1f38aa47af2
# ╠═18d26c89-e953-46bf-9a96-512d55c75a7d
# ╟─0c7b5add-fc89-4633-8c1b-9f472cb5707e
# ╠═8f6a723c-bc5c-46e9-aaeb-af13a25ba4d8
# ╟─5ab12071-3a83-452c-b85e-7b475b5a86e2
# ╠═dfb02cfb-5451-44da-bb74-716aa0fb9cb8
# ╟─cb61bd15-0c0c-43f8-abc5-4585b7544c16
# ╟─385d8b97-acb4-4cc0-b258-a4237bbf296e
# ╠═50c0fd27-9f9d-4fd5-a053-8a9cbdc06f81
# ╟─c7a78987-6c06-41e9-bcc3-067468efe961
# ╠═d0de65a3-82e1-44f7-be84-5c09511adbda
# ╟─e2e336f8-26da-456f-8feb-546086b97861
# ╟─7180a98a-1639-4647-a973-5cb4a654776c
# ╟─7b1bc5cd-15a2-4b6e-bc23-983201ea6f88
# ╠═8b2bc998-e4d6-42eb-88b7-abec6d3c5b1b
# ╟─f4c9635b-86b9-4ea7-8433-6d42f347695f
# ╠═918e65d6-2e7f-4db5-9937-ebcd719b42ff
# ╟─e7952445-dbfc-4e06-a00e-d23038f6435b
# ╟─6555eb04-d630-4d03-800b-5d023c8e0522
# ╟─2baea1f0-2a27-4aca-991b-6f3e7bb13c64
# ╟─5e631aff-d7a4-487c-b21a-03c47be08d8e
# ╠═e98ccebd-b7fd-4de3-8e16-4056512ed1a9
# ╟─a97f89d2-b8fc-4b7b-b0ab-c4170017db27
# ╟─42be67e7-bc56-49d4-a391-fededc1bd1b6
# ╠═123564f2-8a13-47cb-9f9c-5af7fd19b6b8
# ╠═999c382b-5de9-4324-bea4-f6ddb1f37735
# ╟─c03ddad8-d708-4d4f-8e08-ae7ce5a05bb9
# ╟─3d868917-38a8-41cf-acb5-d5a1a2615448
# ╟─487f817d-d67f-4152-bae6-7ec593f253fe
# ╠═37bd1aba-72ba-4b42-b8aa-789318691328
# ╠═31c8cf4b-cc86-4323-8587-16082443d3ae
# ╠═594c91e1-565a-49ad-9b77-20aabe7e71f4
# ╠═09a7db9c-0afb-41a1-a66b-924469ea003c
# ╠═95b0e957-1deb-4176-b819-9d716a09e5b1
# ╠═1c4de053-f8d5-4363-8224-af5a71594978
# ╠═2bd8cb76-4c90-48f6-80d5-6c790045f314
# ╟─a9cba3ab-275a-4250-a207-88506d796945
# ╟─16ac4fe7-009b-4be5-8ce0-1c59385bec0a
# ╟─4d7e1f0f-327a-4316-8595-22d149a13939
# ╠═cb93aba9-0fd6-4716-b9a4-dde69c365570
# ╟─0fe458e1-c021-4115-a941-8a1eeb67f63b
# ╠═167eec1b-7f36-4ff0-9333-b086a09fb6b2
# ╟─ec6d6033-0fd9-4b03-a82d-16f4258dd60b
# ╠═8dcbf3c5-0de4-4235-8330-45886c562535
# ╟─6bcc1479-30f6-44e1-9745-04217b9961d9
# ╠═2086df08-3915-46d7-aafa-2700145265c4
# ╟─cc5439ff-17b7-4a8f-98fb-3bf82bc1544f
# ╟─c0ef8625-3d6a-457b-a07c-cc03ecb31a5d
# ╠═0962b3ae-4183-4e5c-8e5d-1c0817e548b7
# ╠═9fbb2cba-005c-4770-97d9-6a4e279cb2a7
# ╠═e28ee6fb-dbab-4884-867f-c60e9bf55e65
# ╠═ed5d3e9d-5b11-4c9e-b701-c02ece625d01
# ╠═2e9afe63-0048-4e91-b9b3-778957d3df00
# ╠═2f1996bb-514b-4666-ae47-c84fdca266b6
# ╠═a0a45805-ec74-4df4-8c6f-d0a7d1534103
# ╠═8d66003b-3123-42a1-8f3e-4ca58eab0c32
# ╠═8e437f0f-20fd-40b7-bfac-e660d8dfa7b8
# ╠═1265dd87-7395-4637-8b07-8dd5df20e43b
# ╠═e3c72be9-8488-4c01-843c-16553057a0b5
# ╠═bae108cc-f9d5-4667-aa51-6fc3481b02c5
# ╠═f27b825e-b90f-4b8b-a04c-bc7df379383c
# ╠═eb8eff01-cdaf-4536-8bd9-b6b99665a5ac
# ╠═c1f7876c-9262-445f-b199-6d7f46434294
# ╠═36ccd512-b91d-47a7-a86d-a0baeb7faf32
# ╠═9d343fd3-0258-41e4-9e06-938282491ecf
# ╠═8dc78491-a634-4290-9ecc-c13e9be46c2e
# ╠═f7bb31e5-9a69-4a9d-bd62-02c3b7d0d332
# ╠═1c17fd5e-f047-4bd4-8b00-56fa077c7bb0
# ╠═5ebbcfb6-f534-490f-947f-b7c039156908
# ╠═9515317f-3fdb-4b8c-8686-3d0b1a42ffe3
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
