import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser
from tqdm import tqdm
import os
import pandas as pd

def calculate_ppm_range(mz, ppm):
    """Calculate valid m/z range for given ppm tolerance."""
    tolerance = mz * ppm / 1e6
    mz_min = mz - tolerance
    mz_max = mz + tolerance
    return mz_min, mz_max

def find_spectra_in_range(parser, target_mz, ppm, show_progress=True):
    """Find all spectra within ppm error range."""
    mz_min, mz_max = calculate_ppm_range(target_mz, ppm)
    matching_coordinates = []
    matching_intensities = []
    iterator = tqdm(range(len(parser.coordinates)), desc="Processing pixels") if show_progress else range(len(parser.coordinates))
    for coord in iterator:
        mz, intensity = parser.getspectrum(coord)
        mask = (mz >= mz_min) & (mz <= mz_max)
        if np.any(mask):
            matching_coordinates.append(parser.coordinates[coord])
            matching_intensities.append(np.sum(intensity[mask]))
    return matching_coordinates, matching_intensities

def plot_imaging_map(coordinates, intensities, title="m/z Imaging", all_coords=None, vmin=None, vmax=None, save_path=None):
    """Plot MS imaging map, fill zero for missing regions. Save image if save_path is specified."""
    coordinates = np.array(coordinates)
    intensities = np.array(intensities)
    if all_coords is not None:
        all_coords = np.array(all_coords)
        x_min = int(np.min(all_coords[:, 0]))
        x_max = int(np.max(all_coords[:, 0]))
        y_min = int(np.min(all_coords[:, 1]))
        y_max = int(np.max(all_coords[:, 1]))
    else:
        x_min = int(np.min(coordinates[:, 0]))
        x_max = int(np.max(coordinates[:, 0]))
        y_min = int(np.min(coordinates[:, 1]))
        y_max = int(np.max(coordinates[:, 1]))
    image = np.zeros((y_max - y_min + 1, x_max - x_min + 1))
    for coord, intensity in zip(coordinates, intensities):
        x = int(coord[0]) - x_min
        y = int(coord[1]) - y_min
        image[y, x] = intensity
    plt.figure(figsize=(8, 6))
    plt.imshow(image, cmap='hot', origin='lower', vmin=vmin, vmax=vmax)
    plt.colorbar(label="Intensity")
    plt.title(title)
    plt.xlabel("X Coordinate")
    plt.ylabel("Y Coordinate")
    if save_path:
        plt.savefig(save_path, dpi=300, bbox_inches='tight')
        plt.close()
    else:
        plt.show()
    return image

def search_mz_and_plot(imzml_file, target_mz, ppm, vmin=None, vmax=None, verbose=True, show_progress=True, save_dir=None):
    """Search spectra by m/z and plot imaging map."""
    parser = ImzMLParser(imzml_file)
    all_coords = parser.coordinates
    matching_coordinates, matching_intensities = find_spectra_in_range(parser, target_mz, ppm, show_progress)
    if matching_coordinates:
        if save_dir:
            os.makedirs(save_dir, exist_ok=True)
            base_name = f"mz_{target_mz:.4f}_ppm_{ppm:.1f}"
            base_name = base_name.replace('.', '_')
            csv_filename = f"{base_name}.csv"
            png_filename = f"{base_name}.png"
            csv_save_path = os.path.join(save_dir, csv_filename)
            png_save_path = os.path.join(save_dir, png_filename)
        else:
            csv_save_path = None
            png_save_path = None
        image = plot_imaging_map(
            matching_coordinates, matching_intensities, 
            title=f"m/z = {target_mz} Imaging (ppm = {ppm})", 
            all_coords=all_coords,
            vmin=vmin, vmax=vmax,
            save_path=png_save_path
        )
        if csv_save_path:
            np.savetxt(csv_save_path, image, delimiter=',', fmt='%.4f')
            if verbose:
                print(f"Matrix data saved to: {csv_save_path}")
        if png_save_path and verbose:
            print(f"Image saved to: {png_save_path}")
    elif verbose:
        print(f"No spectra found within the m/z = {target_mz} ({ppm} ppm error) range.")

def read_mz_ppm_from_csv(csv_file):
    """Read m/z and ppm values from CSV file."""
    try:
        df = pd.read_csv(csv_file)
        if 'm/z' not in df.columns or 'ppm' not in df.columns:
            print("CSV file must contain 'm/z' and 'ppm' columns")
            return []
        mz_ppm_list = []
        for _, row in df.iterrows():
            mz_ppm_list.append((float(row['m/z']), float(row['ppm'])))
        return mz_ppm_list
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return []

if __name__ == "__main__":
    imzml_file = input("Enter imzML file path: ")
    save_dir = input("Enter output folder path (leave blank to not save): ")
    mode = input("Select mode: 1 - Manual m/z and ppm input, 2 - Batch from CSV: ")
    use_custom_scale = input("Set custom intensity range? (y/n): ").lower() == 'y'
    vmin, vmax = None, None
    if use_custom_scale:
        try:
            vmin = float(input("Enter minimum intensity value (vmin): "))
            vmax = float(input("Enter maximum intensity value (vmax): "))
            print(f"Intensity range set: {vmin} - {vmax}")
        except ValueError:
            print("Invalid input, using automatic intensity range")
            vmin, vmax = None, None
    if mode == "1":
        while True:
            target_mz_input = input("Please enter the m/z value (type 'q' to quit): ")
            if target_mz_input.lower() == 'q':
                break
            try:
                target_mz = float(target_mz_input)
            except ValueError:
                print("Invalid m/z value, please enter a number.")
                continue
            ppm_input = input("Please enter the ppm value: ")
            try:
                ppm = float(ppm_input)
            except ValueError:
                print("Invalid ppm value, please enter a number.")
                continue
            mz_min, mz_max = calculate_ppm_range(target_mz, ppm)
            print(f"Selected m/z range: {mz_min} to {mz_max}")
            search_mz_and_plot(imzml_file, target_mz, ppm, vmin=vmin, vmax=vmax, save_dir=save_dir)
    elif mode == "2":
        csv_file = input("Enter CSV file path: ")
        verbose = input("Show detailed processing info? (y/n): ").lower() == 'y'
        mz_ppm_list = read_mz_ppm_from_csv(csv_file)
        if mz_ppm_list:
            print(f"Read {len(mz_ppm_list)} m/z-ppm pairs, processing...")
            for target_mz, ppm in tqdm(mz_ppm_list, desc="Processing m/z items", unit="item"):
                if verbose:
                    print(f"Processing: m/z = {target_mz}, ppm = {ppm}")
                    mz_min, mz_max = calculate_ppm_range(target_mz, ppm)
                    print(f"Selected m/z range: {mz_min} to {mz_max}")
                search_mz_and_plot(imzml_file, target_mz, ppm, 
                                vmin=vmin, vmax=vmax, 
                                verbose=verbose, 
                                show_progress=False, save_dir=save_dir)
            print("All m/z-ppm pairs processed")
        else:
            print("No valid m/z-ppm pairs read from CSV file")
    else:
        print("Invalid selection, exiting")