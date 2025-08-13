import base64
import struct
import zlib
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
from tqdm import tqdm
from pyimzml.ImzMLWriter import ImzMLWriter

def decode_binary_data(binary_data: str, compression: str, dtype: str) -> np.ndarray:
    """Decode Base64-encoded binary mass spectrometry data with validation."""
    if not binary_data:
        return np.array([])
    try:
        decoded_data = base64.b64decode(binary_data)
        if compression == 'zlib':
            decoded_data = zlib.decompress(decoded_data)
        elif compression not in ['none', None]:
            raise ValueError(f"Unsupported compression format: {compression}")
        dtype_map = {
            '32-bit float': ('f', 4),
            '64-bit float': ('d', 8)
        }
        if dtype not in dtype_map:
            raise ValueError(f"Unsupported data type: {dtype}")
        format_char, num_bytes = dtype_map[dtype]
        num_values = len(decoded_data) // num_bytes
        if len(decoded_data) != num_bytes * num_values:
            raise ValueError("Binary data length mismatch")
        return np.array(struct.unpack('<' + format_char * num_values, decoded_data))
    except Exception as e:
        print(f"Decoding error: {str(e)}")
        return np.array([])

def parse_mzml(file_path: str, ms_level: int = 1) -> list:
    """Parse mzML file and filter by MS level."""
    namespaces = {'ns': 'http://psi.hupo.org/ms/mzml'}
    spectra = []
    valid_rt_count = 0
    total_spectra = 0
    try:
        tree = ET.parse(file_path)
        root = tree.getroot()
        for spectrum in tqdm(root.iterfind('.//ns:spectrum', namespaces),
                            desc="Parsing spectra",
                            unit="spectrum"):
            ms_level_param = spectrum.find('.//ns:cvParam[@accession="MS:1000511"]', namespaces)
            if ms_level_param is None or int(ms_level_param.get('value')) != ms_level:
                continue
            total_spectra += 1
            rt_sec = None
            scan_list = spectrum.find('ns:scanList', namespaces)
            if scan_list is not None:
                rt_param = scan_list.find('.//ns:cvParam[@accession="MS:1000016"]', namespaces)
                if rt_param is not None:
                    rt_value = float(rt_param.get('value'))
                    unit_accession = rt_param.get('unitAccession', 'UO:0000031')
                    if unit_accession == 'UO:0000031':
                        rt_sec = rt_value * 60
                    elif unit_accession == 'UO:0000010':
                        rt_sec = rt_value
                    else:
                        print(f"Unknown time unit: {unit_accession}")
                    valid_rt_count += 1
            mz_values, inten_values = None, None
            for array in spectrum.findall('ns:binaryDataArrayList/ns:binaryDataArray', namespaces):
                array_type = None
                compression = 'none'
                dtype = '32-bit float'
                for cv_param in array.findall('ns:cvParam', namespaces):
                    accession = cv_param.get('accession')
                    if accession == 'MS:1000514':
                        array_type = 'mz'
                    elif accession == 'MS:1000515':
                        array_type = 'intensity'
                    elif accession == 'MS:1000574':
                        compression = 'zlib'
                    elif accession == 'MS:1000521':
                        dtype = '32-bit float'
                    elif accession == 'MS:1000523':
                        dtype = '64-bit float'
                if (binary_elem := array.find('ns:binary', namespaces)) is not None:
                    decoded = decode_binary_data(binary_elem.text, compression, dtype)
                    if array_type == 'mz':
                        mz_values = decoded
                    elif array_type == 'intensity':
                        inten_values = decoded
            if mz_values is not None and inten_values is not None:
                valid_mask = inten_values > 0
                spectra.append({
                    'rt_sec': rt_sec,
                    'mz': mz_values[valid_mask],
                    'intensity': inten_values[valid_mask]
                })
        print(f"\nParsing summary:")
        print(f"Total spectra (ms{ms_level}): {total_spectra}")
        print(f"Spectra with valid retention time: {valid_rt_count} ({valid_rt_count/total_spectra:.1%})")
    except Exception as e:
        print(f"Error parsing file: {str(e)}")
        raise
    return spectra

def create_imzml(mzml_file: str,
                location_file: str,
                output_file: str,
                time_window: float = None,
                start_rt: float = None,
                ms_level: int = 1):
    """Main conversion function with MS level filtering."""
    if None in (time_window, start_rt):
        try:
            time_window = float(input("Enter time window per pixel (seconds): "))
            start_rt = float(input("Enter start retention time (seconds): "))
            ms_level = int(input("Enter MS level to process (1/2): "))
        except ValueError:
            print("Invalid input, please enter numbers")
            return
    print(f"\n[1/3] Parsing mzML file (MS{ms_level})...")
    spectra_data = parse_mzml(mzml_file, ms_level=ms_level)
    print("\n[2/3] Loading coordinate file...")
    location_df = pd.read_csv(location_file)
    if not {'X', 'Y'}.issubset(location_df.columns):
        raise ValueError("Coordinate file must contain X and Y columns")
    print("\n[3/3] Preprocessing data...")
    sorted_spectra = sorted(
        [s for s in spectra_data if s['rt_sec'] is not None],
        key=lambda x: x['rt_sec']
    )
    if not sorted_spectra:
        raise ValueError("No valid time series data found")
    rt_array = np.array([s['rt_sec'] for s in sorted_spectra])
    max_rt = rt_array[-1]
    num_pixels = len(location_df)
    total_required_time = start_rt + num_pixels * time_window
    if total_required_time > max_rt + 1:
        print(f"\nWarning: Required total time {total_required_time:.1f}s exceeds max data time {max_rt:.1f}s")
        if input("Auto adjust time window? (y/n) ").lower() == 'y':
            time_window = (max_rt - start_rt) / num_pixels
            print(f"Auto-adjusted time window: {time_window:.2f} seconds")
        else:
            if input("Continue may result in many empty data, continue? (y/n) ").lower() != 'y':
                return
    print("\nGenerating imzML file...")
    has_valid_data = False
    with ImzMLWriter(output_file) as writer:
        progress_bar = tqdm(location_df.iterrows(), total=num_pixels, desc="Processing pixels")
        for idx, (_, row) in enumerate(progress_bar):
            x, y = int(row['X']), int(row['Y'])
            window_start = start_rt + idx * time_window
            window_end = window_start + time_window
            start_idx = np.searchsorted(rt_array, window_start, side='left')
            end_idx = np.searchsorted(rt_array, window_end, side='right')
            window_spectra = sorted_spectra[start_idx:end_idx]
            if window_spectra:
                all_mz = np.concatenate([s['mz'] for s in window_spectra])
                all_inten = np.concatenate([s['intensity'] for s in window_spectra])
                unique_mz, inverse_idx, counts = np.unique(all_mz, 
                                                          return_inverse=True,
                                                          return_counts=True)
                summed_inten = np.bincount(inverse_idx, weights=all_inten)
                avg_inten = summed_inten / counts
                writer.addSpectrum(unique_mz, avg_inten, (x, y))
                has_valid_data = True
            else:
                writer.addSpectrum(np.array([0.0]), np.array([0.0]), (x, y))
                progress_bar.write(f"Warning: Pixel({x}, {y}) uses empty data placeholder")
    if not has_valid_data:
        raise ValueError(
            "No valid data included in output file, please check:\n"
            f"1. Start time setting (current: {start_rt}s)\n"
            f"2. Time window setting (current: {time_window}s)\n"
            f"3. Data time range ({rt_array[0]:.1f}s - {max_rt:.1f}s)"
        )
    print(f"\nConversion complete! File saved to: {output_file}")

if __name__ == "__main__":
    input_mzml = input("Enter mzML file path: ")
    coordinate_csv = input("Enter coordinate file path: ")
    output_imzml = input("Enter output imzML file path: ")
    time_window = float(input("Enter time window per pixel (seconds): "))
    start_rt = float(input("Enter start retention time (seconds): "))
    ms_level = int(input("Enter MS level to process (1/2): "))
    create_imzml(
        mzml_file=input_mzml,
        location_file=coordinate_csv,
        output_file=output_imzml,
        time_window=time_window,
        start_rt=start_rt,
        ms_level=ms_level
    )