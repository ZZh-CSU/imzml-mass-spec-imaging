import numpy as np
import matplotlib.pyplot as plt
from pyimzml.ImzMLParser import ImzMLParser
from tqdm import tqdm
import pandas as pd
import os
import glob
import re

def calculate_tolerance_range(target, tolerance):
    """
    Calculate m/z range based on absolute tolerance:
      target - tolerance ~ target + tolerance
    """
    return target - tolerance, target + tolerance

def get_spectrum_metadata(parser, idx):
    """
    Extract spectrum metadata from ImzMLParser (compatible with different pyimzML versions).
    Returns a dictionary containing standard CV parameters.
    """
    try:
        # Method 1: Try accessing _spectra attribute directly
        if hasattr(parser, '_spectra'):
            spectrum = parser._spectra[idx]
            if hasattr(spectrum, 'cvParams'):
                return {'cvParams': spectrum.cvParams}
            elif isinstance(spectrum, dict) and 'cvParams' in spectrum:
                return {'cvParams': spectrum['cvParams']}
        
        # Method 2: Try using getMetadata
        if hasattr(parser, 'getMetadata'):
            meta = parser.getMetadata(idx)
            if 'cvParams' in meta:
                return meta
        
        # Method 3: Fallback
        return {'cvParams': []}
    except Exception:
        return {'cvParams': []}

def extract_intensity(parser, ms_level, target=None, target_fragment=None, tolerance=0.01):
    """
    Enhanced signal extraction function for MS1 or MS2 imaging:
    
    - MS1 mode: Extract signal based on target m/z.
    - MS2 mode: Since imzML files are filtered to contain only the target parent ion during conversion,
                extract signal directly for the target fragment m/z.
    
    Parameters:
      - ms_level: 1 for MS1, 2 for MS2
      - target: Target m/z for MS1 mode
      - target_fragment: Target fragment m/z for MS2 mode
      - tolerance: Absolute m/z tolerance
    """
    results = {}
    
    for idx in tqdm(range(len(parser.coordinates)), desc="Processing pixels"):
        try:
            mz_array, intensity_array = parser.getspectrum(idx)
        except Exception as e:
            print(f"Error reading spectrum {idx}: {str(e)}")
            continue
        
        if ms_level == 2 and target_fragment is not None:
            # For MS2 data, extract signal for target fragment m/z
            frag_low, frag_high = calculate_tolerance_range(target_fragment, tolerance)
            frag_mask = (mz_array >= frag_low) & (mz_array <= frag_high)
            intensity_val = np.max(intensity_array[frag_mask]) if np.any(frag_mask) else 0.0
        else:
            # MS1 mode: Extract signal for target m/z
            mz_low, mz_high = calculate_tolerance_range(target, tolerance)
            mz_mask = (mz_array >= mz_low) & (mz_array <= mz_high)
            intensity_val = np.max(intensity_array[mz_mask]) if np.any(mz_mask) else 0.0

        # Record results (keep only first two dimensions of coordinates)
        coords = tuple(parser.coordinates[idx][:2])
        results[coords] = intensity_val

    return results

def read_ms_data_from_csv(csv_file, ms_level):
    """Read MS1 or MS2 data from CSV file."""
    try:
        df = pd.read_csv(csv_file)
        data_list = []
        
        if ms_level == 1:
            if 'm/z' not in df.columns or 'tolerance' not in df.columns:
                print("CSV file must contain 'm/z' and 'tolerance' columns")
                return []
            for _, row in df.iterrows():
                data_list.append((float(row['m/z']), float(row['tolerance'])))
        elif ms_level == 2:
            if 'fragment' not in df.columns or 'tolerance' not in df.columns:
                print("CSV file must contain 'fragment' and 'tolerance' columns")
                return []
            for _, row in df.iterrows():
                data_list.append((float(row['fragment']), float(row['tolerance'])))
        
        return data_list
    except Exception as e:
        print(f"Error reading CSV file: {e}")
        return []

def plot_imaging(results, title="Imaging", cmap="hot", save_path=None, show=False, 
                 tick_fontsize=10, label_fontsize=12, title_fontsize=14, colorbar_tick_fontsize=30,
                 vmin=None, vmax=None):
    """
    tick_fontsize: Font size for axis ticks
    label_fontsize: Font size for axis labels (X, Y)
    title_fontsize: Font size for image title
    colorbar_tick_fontsize: Font size for colorbar ticks
    vmin/vmax: Optional display range
    """
    import numpy as _np
    import matplotlib.pyplot as _plt
    from matplotlib.ticker import MaxNLocator, FuncFormatter
    import math

    coords = _np.array([tuple(key[:2]) for key in results.keys()]) if results else _np.array([])
    if coords.size == 0:
        print("No data to plot")
        return None
    x_coords = coords[:, 0]
    y_coords = coords[:, 1]

    x_min, x_max = int(_np.min(x_coords)), int(_np.max(x_coords))
    y_min, y_max = int(_np.min(y_coords)), int(_np.max(y_coords))

    img = _np.zeros((y_max - y_min + 1, x_max - x_min + 1))
    for key, val in results.items():
        x, y = key[:2]
        img[int(y) - y_min, int(x) - x_min] = val

    _plt.figure(figsize=(8, 6))
    im = _plt.imshow(img, cmap=cmap, origin='lower', aspect='equal', vmin=vmin, vmax=vmax)

    # Hide axes, ticks, and labels
    _plt.gca().axis('off')

    vmin_used, vmax_used = im.get_clim()

    # Force colorbar to start from 0
    vmin_used = 0.0
    im.set_clim(vmin_used, vmax_used)

    if vmax_used is None:
        vmax_used = 0.0
    if vmax_used > 0:
        # Calculate scale factor n
        n = int(math.floor(math.log10(vmax_used)))
        factor = 10 ** n

        # Use MaxNLocator to generate ticks
        locator = MaxNLocator(nbins=6, prune='both')
        ticks = locator.tick_values(vmin_used, vmax_used)

        if len(ticks) < 2:
            ticks = _np.linspace(max(0.0, vmin_used), vmax_used, num=5)

        ticks = _np.unique(_np.concatenate(([vmin_used], ticks)))
        ticks = _np.sort(ticks)

        cb = _plt.colorbar(im, pad=0.01, shrink=0.95)
        cb.set_ticks(ticks)

        vals = (ticks / factor) if factor != 0 else ticks

        def needed_decimals(x, max_check=3, tol=1e-8):
            for d in range(0, max_check + 1):
                if abs(x - round(x, d)) < tol:
                    return d
            return max_check

        decs = [needed_decimals(float(v)) for v in vals]
        max_d = min(3, max(decs)) if decs else 0

        fmt = FuncFormatter(lambda x, pos, d=max_d, f=factor: f"{(x / f):.{d}f}" if f != 0 else f"{x:.{d}f}")
        cb.ax.yaxis.set_major_formatter(fmt)
        cb.ax.tick_params(labelsize=colorbar_tick_fontsize)
        
        for lbl in cb.ax.get_yticklabels():
            lbl.set_fontweight('bold')

        if n != 0:
            sup_map = {'-':'⁻','0':'⁰','1':'¹','2':'²','3':'³','4':'⁴','5':'⁵','6':'⁶','7':'⁷','8':'⁸','9':'⁹'}
            s = str(n)
            sup = ''.join(sup_map[ch] if ch in sup_map else ch for ch in s)
            label = f"×10{sup}"
            cb.ax.text(0, 1.02, label, transform=cb.ax.transAxes,
                       va='bottom', ha='left', fontsize=colorbar_tick_fontsize, fontweight='bold')
    else:
        cb = _plt.colorbar(im, pad=0.01, shrink=0.95)
        if float(vmax_used) == 0.0:
            cb.set_ticks([0.0])
            try:
                from matplotlib.ticker import FuncFormatter
                cb.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, pos: f"{x:.0f}"))
            except Exception:
                cb.ax.set_yticklabels(['0'])
            cb.ax.tick_params(labelsize=colorbar_tick_fontsize)
            for lbl in cb.ax.get_yticklabels():
                lbl.set_fontweight('bold')
        else:
            cb.ax.tick_params(labelsize=colorbar_tick_fontsize)
            for lbl in cb.ax.get_yticklabels():
                lbl.set_fontweight('bold')

    _plt.tight_layout()

    if save_path:
        _plt.savefig(save_path, dpi=300, bbox_inches='tight')

    if show:
        _plt.show()
    else:
        _plt.close()

    return img

def extract_parent_mz_from_filename(filename):
    """Extract parent ion m/z from imzML filename."""
    base_name = os.path.basename(filename)
    match = re.search(r'(\d+\.\d+)', base_name)
    if match:
        return float(match.group(1))
    return None

def get_fragment_csv_path(parent_mz, fragment_dir):
    """Construct corresponding fragment CSV path based on parent m/z."""
    return os.path.join(fragment_dir, f"{parent_mz:.4f}.csv")

def batch_process_imzml_files(imzml_dir, fragment_dir, output_base_dir):
    """Batch process multiple imzML files and their corresponding fragment lists."""
    imzml_files = glob.glob(os.path.join(imzml_dir, "*.imzML"))
    
    if not imzml_files:
        print(f"No imzML files found in {imzml_dir}")
        return
    
    print(f"Found {len(imzml_files)} imzML files, starting processing...")
    
    for imzml_file in imzml_files:
        parent_mz = extract_parent_mz_from_filename(imzml_file)
        if not parent_mz:
            print(f"Could not extract parent m/z from filename {os.path.basename(imzml_file)}, skipping")
            continue
            
        fragment_csv = get_fragment_csv_path(parent_mz, fragment_dir)
        if not os.path.exists(fragment_csv):
            print(f"Corresponding fragment CSV file not found: {fragment_csv}, skipping")
            continue
            
        parent_output_dir = os.path.join(output_base_dir, f"parent_{parent_mz:.4f}")
        os.makedirs(parent_output_dir, exist_ok=True)
        
        print(f"\nProcessing parent ion {parent_mz:.4f} imzML file: {os.path.basename(imzml_file)}")
        print(f"Using fragment CSV file: {fragment_csv}")
        print(f"Output directory: {parent_output_dir}")
        
        try:
            parser = ImzMLParser(imzml_file)
            
            data_list = read_ms_data_from_csv(fragment_csv, ms_level=2)
            if not data_list:
                print(f"Failed to read fragment data from CSV file {fragment_csv}, skipping")
                continue
                
            print(f"Read {len(data_list)} fragments, starting processing...")
            
            for i, (fragment, tolerance) in enumerate(data_list, start=1):
                print(f"Processing {i}/{len(data_list)}: fragment = {fragment}, tolerance = {tolerance}")
                results = extract_intensity(parser, ms_level=2, target_fragment=fragment, tolerance=tolerance)
                title = f"MS/MS Imaging\nParent: {parent_mz:.4f}  Fragment: {fragment:.4f}\nTolerance: ±{tolerance}"
                
                base_name = f"MS2_parent_{parent_mz:.4f}_fragment_{fragment:.4f}_tolerance_{tolerance:.4f}"
                base_name = base_name.replace('.', '_')
                
                csv_filename = f"{base_name}.csv"
                csv_save_path = os.path.join(parent_output_dir, csv_filename)
                
                img_filename = f"{base_name}.png"
                img_save_path = os.path.join(parent_output_dir, img_filename)
                
                img = plot_imaging(results, title, save_path=img_save_path)
                if img is not None:
                    np.savetxt(csv_save_path, img, delimiter=',', fmt='%.4f')
                    
            print(f"All fragments for parent ion {parent_mz:.4f} processed")
            
        except Exception as e:
            print(f"Error processing file {imzml_file}: {str(e)}")
    
    print("\nAll imzML files processed")

if __name__ == "__main__":
    print("Please select operation mode:")
    print("1 - Process single imzML file")
    print("2 - Batch process multiple imzML files")
    mode_choice = input("Enter choice (1/2): ")
    
    if mode_choice == "1":
        imzml_file = input("Enter imzML file path: ")
        if not os.path.exists(imzml_file):
            print(f"File does not exist: {imzml_file}")
            exit(1)
            
        parser = ImzMLParser(imzml_file)
        
        try:
            ms_level = int(input("Enter MS level (1 for MS1; 2 for MS2): "))
        except ValueError:
            print("Invalid number, exiting")
            exit(1)
        
        mode = input("Select mode: 1 - Manual input, 2 - Batch from CSV: ")
        
        save_dir = input("Enter save directory (leave empty to skip): ")
        if save_dir and not os.path.exists(save_dir):
            os.makedirs(save_dir, exist_ok=True)
        
        if mode == "1":
            if ms_level == 2:
                target_fragment = float(input("Enter target fragment m/z: "))
                tolerance = float(input("Enter m/z tolerance (absolute, e.g., 0.01): "))
                results = extract_intensity(parser, ms_level, target_fragment=target_fragment, tolerance=tolerance)
                title = f"MS/MS Imaging\nParent: 610.5781  Fragment: {target_fragment:.4f}\nTolerance: ±{tolerance}"
                
                if save_dir:
                    base_name = f"MS2_fragment_{target_fragment:.4f}_tolerance_{tolerance:.4f}"
                    base_name = base_name.replace('.', '_')
                    csv_filename = f"{base_name}.csv"
                    csv_save_path = os.path.join(save_dir, csv_filename)
                    img_filename = f"{base_name}.png"
                    img_save_path = os.path.join(save_dir, img_filename)
                    img = plot_imaging(results, title, save_path=img_save_path)
                    if img is not None:
                        np.savetxt(csv_save_path, img, delimiter=',', fmt='%.4f')
                        print(f"Matrix data saved to: {csv_save_path}")
            else:
                target_mz = float(input("Enter target m/z: "))
                tolerance = float(input("Enter m/z tolerance (absolute, e.g., 0.01): "))
                results = extract_intensity(parser, ms_level, target=target_mz, tolerance=tolerance)
                title = f"MS Imaging\nTarget: {target_mz:.4f}\nTolerance: ±{tolerance}"
                
                if save_dir:
                    base_name = f"MS1_target_{target_mz:.4f}_tolerance_{tolerance:.4f}"
                    base_name = base_name.replace('.', '_')
                    csv_filename = f"{base_name}.csv"
                    csv_save_path = os.path.join(save_dir, csv_filename)
                    img_filename = f"{base_name}.png"
                    img_save_path = os.path.join(save_dir, img_filename)
                    img = plot_imaging(results, title, save_path=img_save_path)
                    if img is not None:
                        np.savetxt(csv_save_path, img, delimiter=',', fmt='%.4f')
                        print(f"Matrix data saved to: {csv_save_path}")
        
        elif mode == "2":
            csv_file = input("Enter CSV file path: ")
            data_list = read_ms_data_from_csv(csv_file, ms_level)
            
            if data_list:
                print(f"Read {len(data_list)} items, starting processing...")
                for i, (value, tolerance) in enumerate(data_list, start=1):
                    if ms_level == 2:
                        print(f"Processing {i}/{len(data_list)}: fragment = {value}, tolerance = {tolerance}")
                        results = extract_intensity(parser, ms_level, target_fragment=value, tolerance=tolerance)
                        title = f"MS/MS Imaging\nParent: 800.6609  Fragment: {value:.4f}\nTolerance: ±{tolerance}"
                        
                        if save_dir:
                            base_name = f"MS2_fragment_{value:.4f}_tolerance_{tolerance:.4f}"
                            base_name = base_name.replace('.', '_')
                            csv_filename = f"{base_name}.csv"
                            csv_save_path = os.path.join(save_dir, csv_filename)
                            img_filename = f"{base_name}.png"
                            img_save_path = os.path.join(save_dir, img_filename)
                            img = plot_imaging(results, title, save_path=img_save_path)
                            if img is not None:
                                np.savetxt(csv_save_path, img, delimiter=',', fmt='%.4f')
                    else:
                        print(f"Processing {i}/{len(data_list)}: m/z = {value}, tolerance = {tolerance}")
                        results = extract_intensity(parser, ms_level, target=value, tolerance=tolerance)
                        title = f"MS Imaging\nTarget: {value:.4f}\nTolerance: ±{tolerance}"
                        
                        if save_dir:
                            base_name = f"MS1_target_{value:.4f}_tolerance_{tolerance:.4f}"
                            base_name = base_name.replace('.', '_')
                            csv_filename = f"{base_name}.csv"
                            csv_save_path = os.path.join(save_dir, csv_filename)
                            img_filename = f"{base_name}.png"
                            img_save_path = os.path.join(save_dir, img_filename)
                            img = plot_imaging(results, title, save_path=img_save_path)
                            if img is not None:
                                np.savetxt(csv_save_path, img, delimiter=',', fmt='%.4f')
                print("All items processed")
            else:
                print("Failed to read valid data from CSV")
        
    elif mode_choice == "2":
        imzml_dir = input("Enter directory containing imzML files: ")
        fragment_dir = input("Enter directory containing fragment CSV files: ")
        output_dir = input("Enter main output directory: ")
        
        if not os.path.exists(imzml_dir):
            print(f"Directory does not exist: {imzml_dir}")
            exit(1)
            
        if not os.path.exists(fragment_dir):
            print(f"Directory does not exist: {fragment_dir}")
            exit(1)
            
        if output_dir and not os.path.exists(output_dir):
            os.makedirs(output_dir, exist_ok=True)
            
        batch_process_imzml_files(imzml_dir, fragment_dir, output_dir)
    
    else:
        print("Invalid choice, exiting")
