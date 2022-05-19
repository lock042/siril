#ifndef SRC_CORE_COMMAND_DEF_H_
#define SRC_CORE_COMMAND_DEF_H_

#define STR_NONE ""

#define STR_ADDMAX N_("Computes a new image combining the image in memory with the image <b>filename</b>. At each pixel location, the new value is determined as the max of value in current image and in <b>filename</b>")
#define STR_AUTOSTRETCH N_("Auto-stretches the currently loaded image, with different parameters for each channel (unlinked) unless <b>-linked</b> is passed. Arguments are optional, <b>shadowclip</b> is the shadows clipping point, measured in sigma units from the main histogram peak (default is -2.8), <b>targetbg</b> is the target background value, giving a final brightness to the image, range [0, 1], default is 0.25")
#define STR_ASINH N_("Stretches the image to show faint objects using an hyperbolic arcsin transformation. The mandatory argument <b>stretch</b>, typically between 1 and 1000, will give the strength of the stretch. The black point can be offset by providing an optional <b>offset</b> argument in the normalized pixel value of [0, 1]. Finally the option <b>-human</b> enables using human eye luminous efficiency weights to compute the luminance used to compute the stretch value for each pixel, instead of the simple mean of the channels pixel values")

#define STR_BG N_("Returns the background level of the image loaded in memory")
#define STR_BGNOISE N_("Returns the background noise level of the image loaded in memory")
#define STR_BOXSELECT N_("Make a selection area in the currently loaded image with the arguments <b>x</b>, <b>y</b>, <b>width</b> and <b>height</b>, with <b>x</b> and <b>y</b> being the coordinates of the top left corner starting at (0, 0), and <b>width</b> and <b>height</b>, the size of the selection. The <b>-clear</b> argument deletes any selection area. If no argument is passed, the current selection is printed.")

#define STR_CD N_("Sets the new current working directory.\n\nThe argument <b>directory</b> can contain the ~ token, expanded as the home directory, directories with spaces in the name can be protected using single or double quotes")
#define STR_CDG N_("Returns the coordinates of the center of gravity of the image")
#define STR_CLAHE N_("Equalizes the histogram of an image using Contrast Limited Adaptive Histogram Equalization.\n\n<b>cliplimit</b> sets the threshold for contrast limiting.\n<b>tilesize</b> sets the size of grid for histogram equalization. Input image will be divided into equally sized rectangular tiles")
#define STR_CLEAR N_("Clears the graphical output logs")
#define STR_CLEARSTAR N_("Clear all the stars saved in memory and displayed on the screen")
#define STR_CLOSE N_("Properly closes the opened image and the opened sequence, if any")
#define STR_CONVERT N_("Converts all images in a known format into Siril's FITS images.\n\nThe argument <b>basename</b> is the basename of the new sequence. For FITS images, Siril will try to make a symbolic link. If not possible, files will be copied.\nThe flags <b>-fitseq</b> and <b>-ser</b> can be used to specify an alternative output format, other than the default FITS.\nThe option <b>-debayer</b> applies demosaicing to images. In this case no symbolic link are done.\n<b>-start=index</b> sets the starting index number and the <b>-out=</b> option converts files into the directory <b>out</b>")
#define STR_CONVERTRAW N_("Converts DSLR RAW files into Siril's FITS images.\n\nThe argument <b>basename</b> is the basename of the new sequence. For FITS images, Siril will try to make a symbolic link. If not possible, files will be copied.\nThe flags <b>-fitseq</b> and <b>-ser</b> can be used to specify an alternative output format, other than the default FITS.\nThe option <b>-debayer</b> applies demosaicing to images. In this case no symbolic link are done.\n<b>-start=index</b> sets the starting index number and the <b>-out=</b> option converts files into the directory <b>out</b>")
#define STR_COSME N_("Applies the local mean to a set of pixels on the in-memory image (cosmetic correction). The coordinates of these pixels are in an ASCII file [.lst file]. COSME is adapted to correct residual hot and cold pixels after preprocessing")
#define STR_COSME_CFA N_("Same function as COSME but applying to RAW CFA images")
#define STR_CROP N_("Crops a selection of the loaded image.\n\nIn the GUI, if a selection is active, no further arguments are required. Otherwise, or in scripts, arguments have to be given, with <b>x</b> and <b>y</b> being the coordinates of the top left corner, and <b>width</b> and <b>height</b> the size of the selection")

#define STR_DDP N_("Performs a DDP (digital development processing) as described first by Kunihiko Okano. This implementation is the one described in IRIS.\n\nIt combines a linear distribution on low levels (below <b>level</b>) and a non-linear one on high levels.\nIt uses a Gaussian filter of standard deviation <b>sigma</b> and multiplies the resulting image by <b>coef</b>. Typical values for <b>sigma</b> are within 0.7 and 2")
#define STR_DUMPHEADER N_("Dumps the FITS header")

#define STR_ENTROPY N_("Computes the entropy of the opened image on the displayed layer, only in the selected area if one has been selected or in the whole image. The entropy is one way of measuring the noise or the details in an image")
#define STR_EXIT N_("Quits the application")
#define STR_EXTRACT N_("Extracts <b>NbPlans</b> planes of wavelet domain")
#define STR_EXTRACTHA N_("Extracts Ha signal from a CFA image. The output file name starts with the prefix \"Ha_\"")
#define STR_EXTRACTHAOIII N_("Extracts Ha and OIII signals from a CFA image. The output files names start with the prefix \"Ha_\" and \"OIII_\"")
#define STR_EXTRACTGREEN N_("Extracts green signal from a CFA image. The output file name starts with the prefix \"Green_\"")

#define STR_FDIV N_("Divides the image in memory by the image given in argument. The resulting image is multiplied by the value of the <b>scalar</b> argument. See also IDIV")
#define STR_FFTD N_("Applies a Fast Fourier Transform to the image loaded in memory. <b>Modulus</b> and <b>phase</b> given in argument are saved in FITS files")
#define STR_FFTI N_("Retrieves corrected image applying an inverse transformation. The <b>modulus</b> and <b>phase</b> used are the files given in argument")
#define STR_FILL N_("Fills the whole current image (or selection) with pixels having the <b>value</b> intensity expressed in ADU")
#define STR_FILL2 N_("Same command as FILL but this is a symmetric fill of a region defined by the mouse. Used to process an image in the Fourier (FFT) domain")
#define STR_FIND_COSME N_("Applies an automatic detection of cold and hot pixels following the thresholds written in arguments")
#define STR_FIND_COSME_CFA N_("Same command as FIND_COSME but for monochromatic CFA images")
#define STR_FIND_HOT N_("Saves a list file <b>filename</b> (text format) in the working directory which contains the coordinates of the pixels which have an intensity <b>hot_sigma</b> times higher and <b>cold_sigma</b> lower than standard deviation. We generally use this command on a master-dark file")
#define STR_FINDSTAR N_("Detects stars having a level greater than a threshold computed by Siril.\n\nThe algorithm is based on the publication of Mighell, K. J. 1999, in ASP Conf. Ser., Vol. 172, Astronomical Data Analysis Software and Systems VIII, eds. D. M. Mehringer, R. L. Plante, and D. A. Roberts (San Francisco: ASP), 317.\nAfter that, a PSF is applied and Siril rejects all detected structures that don't fulfill a set of prescribed detection criteria. Finally, a circle is drawn around detected stars.\n\nSee also the command CLEARSTAR")
#define STR_FIXBANDING N_("Tries to remove the canon banding.\n\nArgument <b>amount</b> define the amount of correction.\n<b>Sigma</b> defines a protection level of the algorithm, higher sigma gives higher protection")
#define STR_FIXXTRANS N_("Fixes the Fujifilm X-Trans Auto Focus pixels.\n\nIndeed, because of the phase detection auto focus system, the photosites used for auto focus get a little less light than the surrounding photosites. The camera compensates for this and increases the values from these specific photosites giving a visible square in the middle of the dark/bias frames")
#define STR_FMEDIAN N_("Performs a median filter of size <b>ksize</b> x <b>ksize</b> (<b>ksize</b> MUST be odd) to the original image with a modulation parameter <b>modulation</b>.\n\nThe output pixel is computed as : out=mod x m + (1 − mod) x in, where m is the median-filtered pixel value. A modulation's value of 1 will apply no modulation")
#define STR_FMUL N_("Multiplies the loaded image by the <b>scalar</b> given in argument")

#define STR_GAUSS N_("Performs a Gaussian filter with the given <b>sigma</b>")
#define STR_GREY_FLAT N_("Equalizes the mean intensity of RGB layers in a CFA image")

#define STR_HELP N_("Lists the available commands or help for one command")
#define STR_HISTO N_("Calculates the histogram of the image channel in memory and produces file histo_[channel name].dat in the working directory")

#define STR_IADD N_("Adds the image in memory to the image <b>filename</b> given in argument")
#define STR_IDIV N_("Divides the image in memory by the image <b>filename</b> given in argument. See also FDIV")
#define STR_IMUL N_("Multiplies the image in memory by the image <b>filename</b> given in argument")
#define STR_ISUB N_("Subtracts the image in memory by the image <b>filename</b> given in argument")

#define STR_LINK N_("Links all FITS images in the working directory with the <b>basename</b> given in argument.\n\nIf no symbolic links could be created, files are copied. It is possible to convert files in another directory with the <b>-out=</b> option")
#define STR_LMATCH N_("Computes a linear function between a <b>reference</b> image and the image in memory.\n\nThe function is then applied to the current image to match it to the reference one. The algorithm will ignore all reference pixels whose values are outside of the [<b>low</b>, <b>high</b>] range")
#define STR_LOAD N_("Loads the image <b>filename</b>\n\nIt first attempts to load <b>filename</b>, then <b>filename</b>.fit, finally <b>filename</b>.fits and finally all supported formats, aborting if none of these are found.\nThis scheme is applicable to every Siril command that involves reading files.\nFits headers MIPS-HI and MIPS-LO are read and their values given to the current viewing levels.\nWriting a known extension <b>.ext</b> at the end of <b>filename</b> will load specifically the image <b>filename.ext</b>: this is used when numerous files have the same name but not the same extension")
#define STR_LOG N_("Computes and applies a logarithmic scale to the current image")
#define STR_LS N_("Lists files and directories in the working directory")

#define STR_MERGE N_("Merges several sequences into one")
#define STR_MIRRORX N_("Flips the image about the vertical axis")
#define STR_MIRRORY N_("Flips the image about the horizontal axis")
#define STR_MTF N_("Applies midtones transfer function to the current loaded image.\n\nThree parameters are needed, <b>low</b>, <b>midtones</b> and <b>high</b> where midtones balance parameter defines a nonlinear histogram stretch in the [0,1] range")

#define STR_NEG N_("Shows the negative view of the current image")
#define STR_NEW N_("Creates a new image filled with zeros with a size of <b>width</b> x <b>height</b>.\n\nThe image is in 32-bit format, and it contains <b>nb_channel</b> channels, <b>nb_channel</b> being 1 or 3. It is not saved, but displayed and can be saved afterwards")
#define STR_NOZERO N_("Replaces null values by <b>level</b> values. Useful before an idiv or fdiv operation")

#define STR_OFFSET N_("Adds the constant <b>value</b> (specified in ADU) to the current image. This constant can take a negative value.\n\nIn 16-bit mode, values of pixels that fall outside of [0, 65535] are clipped. In 32-bit mode, no clipping occurs")

#define STR_PREPROCESS N_("Preprocesses the sequence <b>sequencename</b> using bias, dark and flat given in argument.\n\nFor bias, a uniform level can be specified instead of an image, by entering a quoted expression starting with an = sign, such as -bias=\"=256\" or -bias=\"=64*$OFFSET\".\n\nBy default, cosmetic correction is not activated. If you wish to apply some, you will need to specify it with <b>-cc=</b> option.\nYou can use <b>-cc=dark</b> to detect hot and cold pixels from the masterdark (a masterdark must be given with the <b>-dark=</b> option), optionally followed by <b>siglo</b> and <b>sighi</b> for cold and hot pixels respectively. A value of 0 deactivates the correction. If sigmas are not provided, only hot pixels detection with a sigma of 3 will be applied.\nAlternatively, you can use <b>-cc=bpm</b> followed by the path to your Bad Pixel Map to specify which pixels must be corrected. An example file can be obtained with a <i>find_hot</i> command on a masterdark.\n\nIt is possible to specify if images are CFA for cosmetic correction purposes with the option <b>-cfa</b> and also to demosaic images at the end of the process with <b>-debayer</b>.\nThe <b>-fix_xtrans</b> option is dedicated to X-Trans files by applying a correction on darks and biases to remove an ugly square pattern.\nThe <b>-equalize_cfa</b> option equalizes the mean intensity of RGB layers of the CFA flat master.\nIt is also possible to optimize the dark subtraction with <b>-opt</b>.\nThe output sequence name starts with the prefix \"pp_\" unless otherwise specified with option <b>-prefix=</b>.\nIf <b>-fitseq</b> is provided, the output sequence will be a FITS sequence (single file).\n")
#define STR_PREPROCESS_SINGLE N_("Preprocesses the image <b>imagename</b> using bias, dark and flat given in argument.\n\nFor bias, a uniform level can be specified instead of an image, by entering a quoted expression starting with an = sign, such as -bias=\"=256\" or -bias=\"=64*$OFFSET\".\n\nBy default, cosmetic correction is not activated. If you wish to apply some, you will need to specify it with <b>-cc=</b> option.\nYou can use <b>-cc=dark</b> to detect hot and cold pixels from the masterdark (a masterdark must be given with the <b>-dark=</b> option), optionally followed by <b>siglo</b> and <b>sighi</b> for cold and hot pixels respectively. A value of 0 deactivates the correction. If sigmas are not provided, only hot pixels detection with a sigma of 3 will be applied.\nAlternatively, you can use <b>-cc=bpm</b> followed by the path to your Bad Pixel Map to specify which pixels must be corrected. An example file can be obtained with a <i>find_hot</i> command on a masterdark.\n\nIt is possible to specify if images are CFA for cosmetic correction purposes with the option <b>-cfa</b> and also to demosaic images at the end of the process with <b>-debayer</b>.\nThe <b>-fix_xtrans</b> option is dedicated to X-Trans files by applying a correction on darks and biases to remove an ugly square pattern.\nThe <b>-equalize_cfa</b> option equalizes the mean intensity of RGB layers of the CFA flat master.\nIt is also possible to optimize the dark subtraction with <b>-opt</b>.\nThe output filename starts with the prefix \"pp_\" unless otherwise specified with option <b>-prefix=</b>.\n")
#define STR_PSF N_("Performs a PSF (Point Spread Function) on the selected star and display the results. If provided, the <b>channel</b> argument selects the image channel on which the star will be analyzed. It can be omitted for monochrome images or when run from the GUI with one of the channels active in the view.")

#define STR_REGISTER N_("Performs geometric transforms on images of the sequence given in argument so that they may be superimposed on the reference image. Using stars for registration, this algorithm only works with deepsky images.\n\nThe output sequence name starts with the prefix <b>\"r_\"</b> unless otherwise specified with <b>-prefix=</b> option.\nThe option <b>-drizzle</b> activates the sub-pixel stacking, either by up-scaling by 2 the images created in the rotated sequence or by setting a flag that will proceed to the up-scaling during stacking if <b>-noout</b> is passed.\nThe option <b>-transf=</b> specifies the use of either <b>\"shift\"</b>, <b>\"similarity\"</b>, <b>\"affine\"</b> or <b>\"homography\"</b> transformations respectively, homography being the default.\nThe option <b>-minpairs=</b> will specify the minimum number of star pairs a frame must have with the reference frame, otherwise the frame will be dropped.\nThe option <b>-maxstars=</b> will specify the maximum number of star to find within each frame (must be between 100 and 2000). A larger value will enable to find more stars and perform a more accurate registration but will take more time to run.\nThe registration is done on the green layer for RGB images unless specified by <b>-layer=</b> option (0, 1 or 2 for R, G and B respectively).\n")
#define STR_RELOADSCRIPTS N_("Rescans the scripts folders and updates scripts menu")
#define STR_REQUIRES N_("Returns an error if the version of Siril is older than the one passed in argument")
#define STR_RESAMPLE N_("Resamples image with a factor <b>factor</b>. This is generally used to resize images, a factor of 0.5 divides size by 2.\nIn the graphical user interface, we can see that several interpolation algorithms are proposed. Here, ''Pixel Area Relation'' is used and cannot be changed.")
#define STR_RGBCOMP N_("Create an RGB composition using three independent images, or an LRGB composition using the optional luminance image and three monochrome images or a color image. Result image is called composed_rgb.fit or composed_lrgb.fit unless another name is provided in the optional argument")
#define STR_RGRADIENT N_("Creates two images, with a radial shift (<b>dR</b> in pixels) and a rotational shift (<b>dalpha</b> in degrees) with respect to the point (<b>xc</b>, <b>yc</b>).\n\nBetween these two images, the shifts have the same amplitude, but an opposite sign. The two images are then added to create the final image. This process is also called Larson Sekanina filter")
#define STR_RL N_("Restores an image using the Richardson-Lucy method.\n\n<b>Sigma</b> is the size of the kernel to be applied, while <b>corner_radius_boost</b> is a value which is added to Gaussian sigma for the tiles in the corners of an image.\n<b>Iterations</b> is the number of iterations to be performed")
#define STR_RMGREEN N_("Applies a chromatic noise reduction filter. It removes green tint in the current image. This filter is based on PixInsight's SCNR Average Neutral algorithm and it is the same filter used by HLVG plugin in Photoshop.\nWith the command, lightness is always preserved. For image processing without L* preservation, use the graphical tool and uncheck the corresponding box.\n\n<b>Type</b> can take values 0 for Average Neutral Protection or 1 for Maximum Neutral Protection, defaulting to 0")
#define STR_ROTATE N_("Rotates the image by an angle of <b>degree</b> value. The option <b>-nocrop</b> can be added to avoid the cropping")
#define STR_ROTATEPI N_("Rotates the image of an angle of 180° around its center. This is equivalent to the command \"ROTATE 180\" or \"ROTATE -180\"")

#define STR_SATU N_("Enhances the color saturation of the loaded image. Try iteratively to obtain best results.\n<b>amount</b> can be a positive number to increase color saturation, negative to decrease it, 0 would do nothing, 1 would increase it by 100%\n<b>background_factor</b> is a factor to (median + sigma) used to set a threshold for which only pixels above it would be modified. This allows background noise to not be color saturated, if chosen carefully. Defaults to 1. Setting 0 disables the threshold.\n<b>hue_range_index</b> can be [0, 6], meaning: 0 for pink to orange, 1 for orange to yellow, 2 for yellow to cyan, 3 for cyan, 4 for cyan to magenta, 5 for magenta to pink, 6 for all (default)")
#define STR_SAVE N_("Saves current image to <b>filename</b>.fit (or .fits, depending on your preferences, see SETEXT). Fits headers MIPS-HI and MIPS-LO are added with values corresponding to the current viewing levels")
#define STR_SAVEBMP N_("Saves current image under the form of a bitmap file with 8-bit per channel: <b>filename</b>.bmp (BMP 24-bit)")
#define STR_SAVEJPG N_("Saves current image into a JPG file: <b>filename</b>.jpg.\n\nYou have the possibility to adjust the quality of the compression. A value 100 for <b>quality</b> parameter offers best fidelity while a low value increases the compression ratio. If no value is specified, a default value of 100 is applied")
#define STR_SAVEPNG N_("Saves current image into a PNG file: <b>filename</b>.png, with 16 bits per channel if the loaded image is 16 or 32 bits, and 8 bits per channel if the loaded image is 8 bits")
#define STR_SAVEPNM N_("Saves current image under the form of a Netpbm file format with 16-bit per channel.\n\nThe extension of the output will be <b>filename</b>.ppm for RGB image and <b>filename</b>.pgm for gray-level image")
#define STR_SAVETIF N_("Saves current image under the form of a uncompressed TIFF file with 16-bit per channel: <b>filename</b>.tif")
#define STR_SAVETIF32 N_("Same command as SAVETIF but the output file is saved in 32-bit per channel: <b>filename</b>.tif")
#define STR_SAVETIF8 N_("Same command as SAVETIF but the output file is saved in 8-bit per channel: <b>filename</b>.tif")
#define STR_SELECT N_("This command allows easy mass selection of images in the loaded sequence (from <b>from</b> to <b>to</b> included)")
#define STR_SEQCLEAN N_("This command clears registration and/or statistics data stored in <b>sequencename</b>.\n\nYou can specify to clear only registration or statistics with <b>-reg</b> and <b>-stat</b> options respectively. Both are cleared if no option is passed")
#define STR_SEQCOSME N_("Same command as COSME but for the the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"cosme_\" unless otherwise specified with option <b>-prefix=</b>")
#define STR_SEQCOSME_CFA N_("Same command as COSME_CFA but for the the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"cosme_\" unless otherwise specified with option <b>-prefix=</b>")
#define STR_SEQCROP N_("Crops the sequence given in argument <b>sequencename</b>.\n\nThe crop selection is specified by the upper left corner position <b>x</b> and <b>y</b> and the selection <b>width</b> and <b>height</b>.\nThe output sequence name starts with the prefix \"cropped_\" unless otherwise specified with <b>-prefix=</b> option")
#define STR_SEQEXTRACTHA N_("Same command as EXTRACT_HA but for the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"Ha_\" unless otherwise specified with option <b>-prefix=</b>")
#define STR_SEQEXTRACTGREEN N_("Same command as EXTRACT_GREEN but for the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"Green_\" unless otherwise specified with option <b>-prefix=</b>")
#define STR_SEQEXTRACTHAOIII N_("Same command as EXTRACT_HAOIII but for the sequence <b>sequencename</b>.\n\nThe output sequences names start with the prefixes \"Ha_\" and \"OIII_\"")
#define STR_SEQFIND_COSME N_("Same command as FIND_COSME but for the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"cc_\" unless otherwise specified with <b>-prefix=</b> option")
#define STR_SEQFIND_COSME_CFA N_("Same command as FIND_COSME_CFA but for the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"cc_\" unless otherwise specified with <b>-prefix=</b> option")
#define STR_SEQMTF N_("Same command as MTF but for the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"mtf_\" unless otherwise specified with <b>-prefix=</b> option")
#define STR_SEQPSF N_("Same command as PSF but for a sequence.\n\nResults are displayed in the plots tab if used from the GUI, else it is printed in the console in a form that can be used to produce brightness variation curves. For headless operation, arguments are mandatory and the center of the search box in pixels can be provided with the <b>-at=</b> argument")
#define STR_SEQSPLIT_CFA N_("Same command as SPLIT_CFA but for the sequence <b>sequencename</b>.\n\nThe output sequences names start with the prefix \"CFA_\" and a number unless otherwise specified with <b>-prefix=</b> option")
#define STR_SEQSTAT N_("Same command as STAT for sequence <b>sequencename</b>.\n\nThe <b>output</b> is saved as a csv file given in second argument.\nThe optional parameter defines the number of statistical values computed: <b>basic</b>, <b>main</b> or <b>full</b> (more detailed but longer to compute)")
#define STR_SEQSUBSKY N_("Same command as SUBSKY but for the sequence <b>sequencename</b>.\n\nThe output sequence name starts with the prefix \"bkg_\" unless otherwise specified with <b>-prefix=</b> option")
#define STR_SEQTILT N_("Same command as TILT but for the loaded sequence or the sequence <b>sequencename</b>.\n\nIt generally gives better result")
#define STR_SET16 N_("Disables images to be saved with 32 bits per channel on processing. It uses 16 bits instead")
#define STR_SET32 N_("Allows images to be saved with 32 bits per channel on processing")
#define STR_SETCOMPRESS N_("Defines if images are compressed or not.\n\n<b>0</b> means no compression while <b>1</b> enables compression.\nIf compression is enabled, the type must be explicitly written in the option <b>-type=</b> (\"rice\", \"gzip1\", \"gzip2\").\nAssociated to the compression, the quantization value must be within [0, 256] range. For example, \"setcompress 1 -type=rice 16\" sets the rice compression with a quantization of 16")
#define STR_SETCPU N_("Defines the number of processing threads used for calculation.\n\nCan be as high as the number of virtual threads existing on the system, which is the number of CPU cores or twice this number if hyperthreading (Intel HT) is available")
#define STR_SETEXT N_("Sets the extension used and recognized by sequences.\n\nThe argument <b>extension</b> can be \"fit\", \"fts\" or \"fits\"")
#define STR_SETFINDSTAR N_("Defines stars detection parameters for FINDSTAR and REGISTER commands.\n\n<b>-radius=</b> defines the radius of the initial search box and must be between 3 and 50.\n<b>-sigma=</b> defines the threshold above noise and must be greater or equal to 0.05.\n<b>-roundness=</b> defines minimum star roundness and must between 0 and 0.95\n<b>-auto=</b> is used if set to \"on\", in which case the <b>radius</b> option is corrected to account for the actual sampling defined by <b>focal</b> and <b>pixelsize</b> options (unless already defined in the image header). This option can be deactivated by setting <b>-auto=</b>off.\n<b>-focal=</b> defines the focal length of the telescope.\n<b>-pixelsize=</b> defines the pixel size of the sensor.")
#define STR_SETMAG N_("Calibrates the magnitude by selecting a star and giving the known apparent magnitude.\n\nAll PSF computations will return the calibrated apparent magnitude afterwards, instead of an apparent magnitude relative to ADU values.\nTo reset the magnitude constant see UNSETMAG")
#define STR_SETMAGSEQ N_("Same as SETMAG command but for the loaded sequence. \n\nThis command is only valid after having run SEQPSF or its graphical counterpart (select the area around a star and launch the PSF analysis for the sequence, it will appear in the graphs).\nThis command has the same goal as SETMAG but recomputes the reference magnitude for each image of the sequence where the reference star has been found.\nWhen running the command, the last star that has been analysed will be considered as the reference star. Displaying the magnitude plot before typing the command makes it easy to understand.\nTo reset the reference star and magnitude offset, see UNSETMAGSEQ")
#define STR_SETMEM N_("Sets a new ratio of free memory on memory used for stacking.\n\n<b>Ratio</b> value should be between 0.05 and 2, depending on other activities of the machine. A higher ratio should allow siril to stack faster, but setting the ratio of memory used for stacking above 1 will require the use of on-disk memory, which is very slow and unrecommended")
#define STR_SETPHOT N_("Gets or sets photometry settings, mostly used by SEQPSF. If arguments are provided, they will update the settings. None are mandatory, any can be provided, default values are shown in the command's syntax. At the end of the command, the active configuration will be printed. Aperture is dynamic unless forced, the <b>aperture</b> value from settings is not used if dynamic, FWHM is used instead")
#define STR_SETREF N_("Sets the reference image of the sequence given in first argument")
#define STR_SPLIT N_("Splits the color image into three distinct files (one for each color) and save them in <b>fileR</b>.fit, <b>fileG</b>.fit and <b>fileB</b>.fit files")
#define STR_SPLIT_CFA N_("Splits the CFA image into four distinct files (one for each channel) and save them in files")
#define STR_STACK N_("Stacks the <b>sequencename</b> sequence, using options.\n\nRejection type:\nThe allowed types are: \"sum\", \"max\", \"min\", \"med\" (or \"median\") and \"rej\" (or \"mean\"). If no argument other than the sequence name is provided, sum stacking is assumed.\n\nStacking with rejection:\nTypes rej or mean require the use of additional arguments for rejection type and values. The rejection type is one of { n[one] | p[ercentile] | s[igma] | m[edian] | w[insorized] | l[inear] | g[eneralized] | [m]a[d] } for Percentile, Sigma, Median, Winsorized, Linear-Fit, Generalized Extreme Studentized Deviate Test or k-MAD clipping. If omitted, the default (Winsorized) is used.\nThe <b>sigma low</b> and <b>sigma high</b> parameters of rejection are mandatory unless <b>none</b> is selected.\n\nNormalization of input images:\nFor med|median and rej|mean stacking types, different types of normalization are allowed: <b>-norm=add</b> for additive, <b>-norm=mul</b> for multiplicative. Options <b>-norm=addscale</b> and <b>-norm=mulscale</b> apply same normalization but with scale operations. <b>-nonorm</b> is the option to disable normalization. Otherwise addtive with scale method is applied by default.\n<b>-fastnorm</b> option specifies to use faster estimators for location and scale than the default IKSS.\n\nOther options for rejection stacking:\n<b>-weight_from_noise</b> is an option to add larger weights to frames with lower background noise.\n<b>-weight_from_nbstack</b> weights input images based on how many images were used to create them, useful for live stacking.\n<b>-output_norm</b> applies a normalization at the end of the stacking to rescale result in the [0, 1] range.\n<b>-rgb_equal</b> equalizes the RGB channels of the stacked image (color-only).\n\nOutputs:\nResult image name can be set with the <b>-out=</b> option. Otherwise, it will be named as <b>sequencename</b>_stacked.fit\n\nFiltering out images:\nImages to be stacked can be selected based on some filters, like manual selection or best FWHM, with some of the <b>-filter-*</b> options.\nSee the command reference for the complete documentation on this command")
#define STR_STACKALL N_("Opens all sequences in the current directory and stacks them with the optionally specified stacking type and filtering or with sum stacking. See STACK command for options description")
#define STR_STAT N_("Returns global statistics of the current image. If a selection is made, the command returns statistics within the selection")
#define STR_SUBSKY N_("Computes a synthetic background gradient using either the polynomial function model of <b>degree</b> degrees or the RBF model (if <b>-rbf</b> is provided instead) and subtracts it from the image. The number of samples per horizontal line and the tolerance to exclude brighter areas can be adjusted with the optional arguments. Tolerance is in mad units: median + tolerance * mad. For RBF, the additional smoothing parameter is also available")

#define STR_THRESHLO N_("Replaces values below <b>level</b> with <b>level</b>")
#define STR_THRESHHI N_("Replaces values above <b>level</b> with <b>level</b>")
#define STR_THRESH N_("Replaces values below <b>lo</b> with <b>lo</b> and values above <b>hi</b> with <b>hi</b>")
#define STR_TILT N_("Computes the sensor tilt as the fwhm difference between the best and worst corner truncated mean values. The <b>clear</b> option allows to clear the drawing")

#define STR_UNSELECT N_("Allows easy mass unselection of images in the loaded sequence (from <b>from</b> to <b>to</b> included). See SELECT")
#define STR_UNSETMAG N_("Reset the magnitude calibration to 0. See SETMAG")
#define STR_UNSETMAGSEQ N_("Resets the magnitude calibration and reference star for the sequence. See SETMAGSEQ")
#define STR_UNSHARP N_("Applies to the working image an unsharp mask with sigma <b>sigma</b> and coefficient <b>multi</b>")

#define STR_VISU N_("Displays an image with <b>low</b> and <b>high</b> as the low and high threshold")

#define STR_WAVELET N_("Computes the wavelet transform on <b>nbr_plan</b> plans using linear (<b>type</b>=1) or bspline (<b>type</b>=2) version of the 'à trous' algorithm. The result is stored in a file as a structure containing the planes, ready for weighted reconstruction with WRECONS")
#define STR_WRECONS N_("Reconstructs to current image from the planes previously computed with wavelets and weighted with coefficients <b>c1</b>, <b>c2</b>, ..., <b>cn</b> according to the number of planes used for wavelet transform")

#define STR_START_LS N_("Initialize a livestacking session, using the optional calibration files and wait for input files to be provided by the LIVESTACK command until STOP_LS is called")
#define STR_LIVESTACK N_("Process the provided image for live stacking. Only possible after START_LS")
#define STR_STOP_LS N_("Stop the live stacking session. Only possible after START_LS")


#endif /* SRC_CORE_COMMAND_DEF_H_ */
