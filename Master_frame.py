# -*- coding: utf-8 -*-
"""
Created on Fri Jan 31 13:01:20 2025

@author: sofie
"""

import numpy as np
import matplotlib.pyplot as plt
from astropy.nddata import CCDData
from astropy.stats import mad_std
from astropy.io import fits
import ccdproc as ccdp
from pathlib import Path
from astropy import units as u
from matplotlib.gridspec import GridSpec


from skimage.measure import block_reduce

plt.rc("font", family=["Arial"]) 
plt.rc("axes", labelsize=16)   
plt.rc("xtick", labelsize=16, top=True, direction="in")  
plt.rc("ytick", labelsize=16, right=True, direction="in")
plt.rc("axes", titlesize=18)
plt.rc("legend", fontsize=16)


# path = Path('C:/Users/sofie/OneDrive - Aarhus universitet/Speciale')

path = Path("D:/")
date = path/"2025-02-17"

raw_path = date /"Raw"
raw_images = ccdp.ImageFileCollection(raw_path)

master_path = date / "Masters"
master_path.mkdir(parents=True, exist_ok=True)

cal_path = date/"Calibrated/No_bin"
cal_path.mkdir(parents=True, exist_ok=True)

log_path = path / "log"
log_path.mkdir(parents=True, exist_ok=True)


def Huber_weight(x, delta=1):
    return np.where(np.abs(x) < delta, 1, delta / np.abs(x))
    
    
def robust_mean(image, delta=1):
    # image = image.data.ravel()
    # image = image[~np.isnan(image)]
    
    median = np.median(image)
    residuals = image - median
    weight = Huber_weight(residuals, delta)
    mean = np.average(image, weights=weight)
    
    MAD = np.median(np.abs(image - mean))
    std = 1.4826 * MAD
    
    return mean, std, MAD


def get_date_from_header(filename):
    with fits.open(filename) as hdul:
        date_obs = hdul[0].header.get("Date-obs", None)
    if date_obs:
        return date_obs.split("T")[0]
    return None


def histogram(masterframe, types):
    fig, ax = plt.subplots(figsize=(8,6))
    hist, bin_edges = np.histogram(masterframe, bins=500)
    ax.bar(bin_edges[:-1], hist, width=np.diff(bin_edges), align="edge")
    ax.set_yscale('log')
    ax.set_xlabel("Pixel value [adu]")
    ax.set_ylabel("Number of pixels")
    fig.tight_layout()
    fig.savefig(master_path/f"Histogram_{types}.png")


def master_frame_png(masterframe, filename, date, types):
    mean = robust_mean(masterframe)[0]
    std = robust_mean(masterframe)[1]
    
    xlen = masterframe.shape[1]
    ylen = masterframe.shape[0]
    x_cross_section = xlen//2  
    y_cross_section = ylen//2
    
    x_profile = masterframe[x_cross_section, :]
    x_profile = np.array(x_profile)
    y_profile = masterframe[:, y_cross_section]
    y_profile = np.array(y_profile)
    

    # Create figure and gridspec layout
    fig = plt.figure(figsize=(8, 8*ylen/xlen))
    gs = GridSpec(2, 3, width_ratios=[4, 1, 0.2], height_ratios=[4, 1], wspace=0.03, hspace=0.01)
    
    
    # Main CCD image
    ax_main = plt.subplot(gs[0, 0])
    im = ax_main.imshow(masterframe, cmap='gray', vmin=mean-2*std, vmax=mean+2*std, origin="lower")
    ax_main.axhline(y_cross_section, color='r', linestyle='-')
    ax_main.axvline(x_cross_section, color='b', linestyle='-')
    ax_main.set_ylabel("Y pixels")
    ax_main.text(xlen//15,ylen-ylen//10, f"Median: {mean:.1f}", color="purple")
    ax_main.text(xlen//15,ylen-1.5*ylen//10, f"St.dev: {std:.1f}", color="purple")
    
    # Bottom plot (X-axis cut)
    ax_bottom = plt.subplot(gs[1, 0], sharex=ax_main)
    ax_bottom.plot(np.arange(len(x_profile)),x_profile, color='r')
    ax_bottom.set_ylabel("Pixel value")
    ax_bottom.set_xlabel("X pixels")
    
    # Right plot (Y-axis cut)
    ax_right = plt.subplot(gs[0, 1], sharey=ax_main)
    ax_right.plot(y_profile, np.arange(len(y_profile)), color='b')
    ax_right.set_xlabel("Pixel value")
    
    # Colorbar
    cbar_ax = plt.subplot(gs[0, 2])
    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("Pixel intensity")
        
    # Making the plot easier to read
    ax_bottom.tick_params(labeltop=False)
    ax_right.tick_params(labelleft=False)
    ax_main.tick_params(labelbottom=False)
    
    if types == "flat":
        ax_right.set_xlim(mean-2*std, mean+2*std)
        ax_bottom.set_ylim(mean-2*std, mean+2*std)
        
    
    fig.savefig(master_path/f"masterframe_{filename}_{date}.png")
    
    
def logg(file_path, data_type="raw"):
    file = ccdp.CCDData.read(file_path, unit="adu")
    header = file.header 
    
    JD = header.get("JD", "JD not found")
    Filter = header.get("Filter", "Filter not found")
    Imagetyp = header.get("imagetyp", "Imagetype not found")
    Object = header.get("object", "No object")
    Exposure = header.get("exposure", "Exposure not found")
    mean, std, MAD = robust_mean(file)
    filename = Path(file_path).name
    
    
    columns = [
        ("Filename", 50, "<"),
        ("Julian Date", 20, "<"),
        ("Imagetype", 15, "<"),  
        ("Object", 20, "<"),  
        ("Filter", 10, "<"),  
        ("Exposure", 15, "<"),  
        ("Mean", 15, "<"),  
        ("Std", 15, "<"),  
    ]
    
    header_line = "".join(f"{col[0]:{col[2]}{col[1]}}" for col in columns) + "\n"
    
    log_file = log_path / f"log_{data_type}_{date.name}.txt"
    if not log_file.exists():
        with open(log_file, "w") as f:
            f.write(header_line)
    
    with open(log_file, "a") as f:
        f.write(
            f"{filename:<50}{JD:<20}{Imagetyp:<15}{Object:<20}{Filter:<10}"
            f"{Exposure:<15}{mean:<15.2f}{std:<15.2f}\n"
        )
        return 0


def master_bias_frame(file_collection):
    bias_files = file_collection.files_filtered(imagetyp='bias', include_path=True)
    
    if len(bias_files) == 0:
        bias_path = path / "Bias"
        bias_images = ccdp.ImageFileCollection(bias_path)
        bias_files = bias_images.files_filtered(imagetyp="bias", include_path=True)
    
    bias_obss = [get_date_from_header(f) for f in bias_files]
    bias_obss = sorted(set(filter(None, bias_obss)))
    
    if bias_obss:
        bias_obs = bias_obss[0] 
    else:
        bias_obs = date.name
    
    master_bias = ccdp.combine(bias_files,
                                 method='average',
                                 sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                                 sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                                 mem_limit=350e6, unit="adu"
                                )
    
    if not isinstance(master_bias, CCDData):
        master_bias = CCDData(master_bias, unit="adu")
        
    master_bias.meta['combined'] = True
    master_bias_path = master_path / f'Master_bias_{bias_obs}.fits'
    master_bias.write(master_bias_path, overwrite=True)
    
    logg(master_bias_path, data_type="cal")
    histogram(master_bias, f"bias_{bias_obs}")
    master_frame_png(master_bias, "bias", f"{bias_obs}", "bias")
    return master_bias



def master_flat_frame(filter_name, file_collection):
    bias_master_path = master_path / f'Master_bias_{date.name}.fits'
    
    if bias_master_path.exists():
        bias_master = CCDData.read(bias_master_path, unit="adu")
    else:
        bias_master = master_bias_frame(file_collection)
  
    flat_files = file_collection.files_filtered(imagetyp="flat", filter=filter_name, include_path=True)
    
    if len(flat_files) == 0:
        flat_path = path / "Flats"
        flat_images = ccdp.ImageFileCollection(flat_path)
        flat_files = flat_images.files_filtered(imagetyp="flat", filter=filter_name, include_path=True)
    
    flat_obss = [get_date_from_header(f) for f in flat_files]
    flat_obss = sorted(set(filter(None, flat_obss)))
    
    if flat_obss:
        flat_obs = flat_obss[0] #Use the oldest flat date
    else:
        flat_obs = date.name
    
    
    master_flat = ccdp.combine((ccdp.subtract_bias(CCDData.read(f, unit="adu"), bias_master) for f in flat_files),
                           method='average',
                           sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                           sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                           mem_limit=1e9, unit="adu"
                          )
    master_flat.data = master_flat.data.astype(np.float32)

    if not isinstance(master_flat, CCDData):
        master_flat = CCDData(master_flat, unit="adu")

    master_flat.data /= robust_mean(master_flat.data)[0]
    
    master_flat.meta['combined'] = True
    master_flat_path = master_path / f"Master_flat_{filter_name}_{flat_obs}.fits"
    master_flat.write(master_flat_path, overwrite=True)

    logg(master_flat_path, data_type="cal")
    histogram(master_flat, f"flat_{filter_name}_{flat_obs}")
    master_frame_png(master_flat, f"flat_{filter_name}", f"{flat_obs}", "flat")
    return master_flat



def bad_pixel_frame(file_collection):
    dark_files = file_collection.files_filtered(imagetyp='dark', include_path=True)
    
    if len(dark_files) == 0:
        dark_path = path / "Darks"
        dark_images = ccdp.ImageFileCollection(dark_path)
        dark_files = dark_images.files_filtered(imagetyp="dark", include_path=True)
    
    dark_obss = [get_date_from_header(f) for f in dark_files]
    dark_obss = sorted(set(filter(None, dark_obss)))
    
    if dark_obss:
        dark_obs = dark_obss[0] #Use the oldest flat date
    else:
        dark_obs = date.name
    
    dark_images = []
    
    for dark_file in dark_files:
        dark = CCDData.read(dark_file, unit="adu", format="fits")
        exposure = dark.header.get("exposure", "Exposure not found")
        
        if exposure is None or exposure <= 0:
            return ValueError(f"Invalid exposure time in {dark_file}")
        
        dark_corrected = CCDData(dark.data / exposure, unit="adu", meta=dark.header)
        dark_images.append(dark_corrected)
        
    
    combined_dark = ccdp.combine(dark_images,
                               method='average',
                               sigma_clip=True, sigma_clip_low_thresh=5, sigma_clip_high_thresh=5,
                               sigma_clip_func=np.ma.median, sigma_clip_dev_func=mad_std,
                               mem_limit=350e6, unit="adu"
                               )

    mean, std, MAD = robust_mean(combined_dark, delta=5)
    
    bad_pixels = np.logical_or(combined_dark > mean+5*MAD, combined_dark < mean-5*MAD)
    bad_pixels = bad_pixels.astype(np.uint8)
    
    bad_pixel_sum = bad_pixels.sum()
    percent = bad_pixel_sum/combined_dark.data.sum()
    if not isinstance(bad_pixels, CCDData):
        bad_pixels = CCDData(bad_pixels, unit=u.dimensionless_unscaled, meta=combined_dark.meta)
        
    bad_pixels.meta['combined'] = True
    bad_pixels_path = master_path / f'Bad_pixels_{dark_obs}.fits'
    bad_pixels.write(bad_pixels_path, overwrite=True)
    logg(bad_pixels_path, data_type="cal")

    return bad_pixels, percent, mean, MAD



def calibration(file_path, filter_name, bin_factor=1):
    file_collection = ccdp.ImageFileCollection(raw_path)
    
    raw_image = CCDData.read(file_path, unit="adu")    
    
    bias_master_files = list(master_path.glob("Master_bias_*.fits"))
    if bias_master_files:
        bias_master_path = bias_master_files[0]
        print("check master bias")
        bias_master = CCDData.read(bias_master_path, unit="adu")
    else:
        bias_master = master_bias_frame(file_collection)
        bias_master_path =  master_path / "Master_bias_*.fits"
        print("check master bias")
        
    master_flat_files = list(master_path.glob(f"Master_flat_{filter_name}_*.fits"))
    if master_flat_files:
        flat_master_path = master_flat_files[0]
        print("check flat bias")
        flat_master = CCDData.read(flat_master_path, unit="adu")
    else:
        flat_master = master_flat_frame(filter_name, file_collection)
        flat_master_path = master_path / f"Master_flat_{filter_name}_*.fits"
        print("check flat bias")
        
    bad_pixel_files = list(master_path.glob('Bad_pixels_*.fits'))
    if bad_pixel_files:
        bad_pixel_path = bad_pixel_files[0]
        print("check bad pixels")
        bad_pixels = CCDData.read(bad_pixel_path, unit="adu")
    else:
        bad_pixels, percent, mean, MAD = bad_pixel_frame(file_collection)
        bad_pixel_path = master_path / 'Bad_pixels_*.fits'
        print("check bad pixels")
      
        
    reduced = ccdp.subtract_bias(raw_image, bias_master)
    reduced = ccdp.flat_correct(reduced, flat_master)
    
    
    row_indices, col_indices = np.where(bad_pixels)    
    for i in range (0,len(row_indices)):
        row = row_indices[i]
        col = col_indices[i]
        
        value = []
        
        for r in range(row-1,row+2):
            for c in range(col-1,col+2):
                if (r == row and c == col) or r<0 or c<0 or r>= reduced.shape[0] or c>=reduced.shape[1]:
                    continue
                else:
                    value.append(reduced.data[r, c])

        #print(corrected_value)
        if value:
            reduced.data[row, col] = sum(value) / len(value)
    
    if bin_factor > 1:
        reduced.data = block_reduce(reduced.data, (bin_factor, bin_factor), np.mean)
        if reduced.mask is not None:
            reduced.mask = block_reduce(reduced.mask, (bin_factor, bin_factor), np.min)  

        reduced.header['NAXIS1'] = reduced.data.shape[1]
        reduced.header['NAXIS2'] = reduced.data.shape[0]

    
    
    reduced.data = reduced.data.astype(np.float32)
    
    reduced.meta['BITPIX'] = -32
    reduced.meta['Calibrated'] = True
    reduced.meta['M_bias'] =  Path(bias_master_path).name
    reduced.meta['M_flat'] =  Path(flat_master_path).name
    reduced.meta['Bad_pixel'] =  Path(bad_pixel_path).name
    reduced.meta['Binning'] = bin_factor
    reduced_path = cal_path / f'{Path(file_path).stem}_cal.fits'
    reduced.write(reduced_path, overwrite=True)
    
    print("Check")
    logg(reduced_path, data_type="cal")
    return reduced




raw_bias_files = raw_images.files_filtered(imagetyp="bias", include_path=True)
raw_dark_files = raw_images.files_filtered(imagetyp="dark", include_path=True)
raw_image_files = raw_images.files_filtered(imagetyp="object", include_path=True)


master_flat_frame("V", raw_images)
master_flat_frame("I", raw_images)
master_flat_frame("B", raw_images)


for raw_image in raw_bias_files:
    single_raw_image = raw_image 
    logg(single_raw_image, "raw")
    
for raw_image in raw_dark_files:
    single_raw_image = raw_image 
    logg(single_raw_image, "raw")

for i in range(0, len(raw_image_files)):
    single_raw_image = raw_image_files[i]
    logg(single_raw_image, "raw")
    
    
    raw_image = CCDData.read(single_raw_image, unit="adu") 
    meta = raw_image.header
    filter_name = meta["Filter"]
    calibrated_image = calibration(single_raw_image, f"{filter_name}", bin_factor=1)


    
# fig, ax = plt.subplots(figsize=(8,6))
# im = ax.imshow(calibrated_image.data, cmap='gray', origin="lower", vmin=100, vmax=400)
# fig.colorbar(im, ax=ax)

