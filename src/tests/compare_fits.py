from astropy.io import fits
import numpy as np

def compare_fits_files(file1, file2):
    """
    Compare the data and headers of two FITS files.
    """
    retval = 1
    with fits.open(file1) as hdul1, fits.open(file2) as hdul2:
        # Check if the number of HDUs is the same
        if len(hdul1) != len(hdul2):
            print(f"Number of HDUs differs: {len(hdul1)} vs {len(hdul2)}")
            return

        # Compare each HDU
        for i, (hdu1, hdu2) in enumerate(zip(hdul1, hdul2)):
            print(f"--- HDU {i} ---")

            # Compare data
            data1 = hdu1.data
            data2 = hdu2.data

            if data1 is None and data2 is None:
                print("Both HDUs have no data.")
            elif data1 is None or data2 is None:
                print("One HDU has data, the other does not.")
            else:
                if data1.shape != data2.shape:
                    print("Data arrays have different shapes.")
                else:
                    diff_pixels = np.sum(data1 != data2)
                if diff_pixels == 0:
                    print("Data arrays are identical.")
                    retval = 0
                else:
                    print(f"Data arrays differ in {diff_pixels} pixels.")
                    if diff_pixels < 3:
                        retval = 0

                # Compare headers
                header1 = hdu1.header
                header2 = hdu2.header

                keys1 = set(header1.keys())
                keys2 = set(header2.keys())

                common_keys = keys1 & keys2
                unique_to_1 = keys1 - keys2
                unique_to_2 = keys2 - keys1

                if unique_to_1:
                    print("Keys only in file1:", unique_to_1)
                if unique_to_2:
                    print("Keys only in file2:", unique_to_2)

                differing_values = []
                for key in common_keys:
                    if header1[key] != header2[key]:
                        differing_values.append((key, header1[key], header2[key]))

                if differing_values:
                    print("Keys with different values:")
                    for key, val1, val2 in differing_values:
                        print(f"  {key}: {val1} (file1) vs {val2} (file2)")
                else:
                    print("All common header keys have identical values.")
    return retval

if __name__ == "__main__":
    import sys
    if len(sys.argv) != 3:
        print("Usage: python compare_fits.py file1.fits file2.fits")
        sys.exit(1)
    file1 = sys.argv[1]
    file2 = sys.argv[2]
    sys.exit(compare_fits_files(file1, file2))

