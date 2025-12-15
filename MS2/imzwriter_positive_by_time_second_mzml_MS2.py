import base64
import struct
import zlib
import numpy as np
import xml.etree.ElementTree as ET
import pandas as pd
from tqdm import tqdm
from pyimzml.ImzMLWriter import ImzMLWriter
import os
import time

def decode_binary_data(binary_data: str, compression: str, dtype: str) -> np.ndarray:
    """Decode Base64-encoded binary mass spectrometry data (with data validation)"""
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
    """
    Parse mzML file, support filtering by MS level and extract precursor information.
    For MS2 data extraction, precursor m/z is in cvParam (accession="MS:1000744")
    """
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
            precursor_mz = None
            precursor_charge = None
            
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

            if ms_level == 2:
                precursor_list = spectrum.find('ns:precursorList', namespaces)
                if precursor_list is not None:
                    precursor = precursor_list.find('ns:precursor', namespaces)
                    if precursor is not None:
                        precursor_mz_param = precursor.find('.//ns:cvParam[@accession="MS:1000744"]', namespaces)
                        if precursor_mz_param is not None:
                            precursor_mz = float(precursor_mz_param.get('value'))
                        precursor_charge_param = precursor.find('.//ns:cvParam[@accession="MS:1000041"]', namespaces)
                        if precursor_charge_param is not None:
                            precursor_charge = int(precursor_charge_param.get('value'))

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
                    'intensity': inten_values[valid_mask],
                    'precursor_mz': precursor_mz,
                    'precursor_charge': precursor_charge
                })

        print(f"\nParsing completed statistics:")
        print(f"Total spectra (MS{ms_level}): {total_spectra}")
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
                 ms_level: int = 1,
                 target_parent_mz: float = None,
                 parent_tolerance: float = 0.01
                ):
    """
    Main conversion function (using standard CV parameters for MS2 metadata)

    For MS2 mode, if target_parent_mz is specified, only spectra with precursor within the target range are retained during merging,
    ensuring that subsequent plotting of target precursor fragments does not mix signals from other precursors.
    """
    if None in (time_window, start_rt):
        try:
            time_window = float(input("Enter time window per pixel (seconds): "))
            start_rt = float(input("Enter start retention time (seconds): "))
            ms_level = int(input("Enter MS level to process (1/2): "))
            if ms_level == 2:
                inp = input("Filter by precursor (retain only specified precursor data)? (y/n): ").lower()
                if inp == 'y':
                    target_parent_mz = float(input("Enter target precursor m/z value: "))
                    parent_tolerance = float(input("Enter precursor m/z tolerance (e.g., 0.01): "))
        except ValueError:
            print("Invalid input, please ensure numeric values")
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
        raise ValueError("No valid time-series data found")
    
    rt_array = np.array([s['rt_sec'] for s in sorted_spectra])
    max_rt = rt_array[-1]
    num_pixels = len(location_df)
    
    total_required_time = start_rt + num_pixels * time_window
    if total_required_time > max_rt + 1:
        print(f"\nWarning: Required total time {total_required_time:.1f}s exceeds max data time {max_rt:.1f}s")
        if input("Auto-adjust time window? (y/n) ").lower() == 'y':
            time_window = (max_rt - start_rt) / num_pixels
            print(f"Auto-adjusted time window to: {time_window:.2f} seconds")
        else:
            if input("Continue may result in many empty data points, proceed? (y/n) ").lower() != 'y':
                return

    print("\nGenerating imzML file...")
    has_valid_data = False
    
    with ImzMLWriter(output_file, mode="processed") as writer:
        progress_bar = tqdm(location_df.iterrows(), total=num_pixels, desc="Processing pixels")
        
        for idx, (_, row) in enumerate(progress_bar):
            x, y = int(row['X']), int(row['Y'])
            window_start = start_rt + idx * time_window
            window_end = window_start + time_window
            
            start_idx = np.searchsorted(rt_array, window_start, side='left')
            end_idx = np.searchsorted(rt_array, window_end, side='right')
            window_spectra = sorted_spectra[start_idx:end_idx]

            if ms_level == 2 and target_parent_mz is not None:
                filtered_spectra = []
                for s in window_spectra:
                    pmz = s.get('precursor_mz')
                    if pmz is not None and abs(pmz - target_parent_mz) <= parent_tolerance:
                        filtered_spectra.append(s)
                window_spectra = filtered_spectra

            if window_spectra:
                all_mz = np.concatenate([s['mz'] for s in window_spectra])
                all_inten = np.concatenate([s['intensity'] for s in window_spectra])

                if all_mz.size == 0 or all_inten.size == 0:
                    writer.addSpectrum(
                        np.array([0.0]),
                        np.array([0.0]),
                        coords=(x, y),
                        userParams=[]
                    )
                    progress_bar.write(f"Warning: Pixel ({x}, {y}) has empty merged data")
                    continue

                unique_mz, inverse_idx, counts = np.unique(
                    all_mz, 
                    return_inverse=True,
                    return_counts=True
                )
                summed_inten = np.bincount(inverse_idx, weights=all_inten)
                avg_inten = summed_inten / counts
                
                avg_inten = np.nan_to_num(avg_inten, nan=0.0, posinf=0.0, neginf=0.0)
                
                user_params = []
                if ms_level == 2:
                    precursor_mz = window_spectra[0]['precursor_mz']
                    precursor_charge = window_spectra[0]['precursor_charge']
                    if precursor_mz is not None:
                        user_params.extend([
                            {
                                'accession': 'MS:1000744',
                                'value': str(precursor_mz),
                                'cvRef': 'MS',
                                'name': 'selected ion m/z'
                            },
                            {
                                'accession': 'MS:1000041',
                                'value': str(precursor_charge) if precursor_charge else '0',
                                'cvRef': 'MS',
                                'name': 'charge state'
                            }
                        ])
                
                writer.addSpectrum(
                    mzs=unique_mz,
                    intensities=avg_inten,
                    coords=(x, y),
                    userParams=user_params
                )
                has_valid_data = True
            else:
                writer.addSpectrum(
                    np.array([0.0]),
                    np.array([0.0]),
                    coords=(x, y),
                    userParams=[]
                )
                progress_bar.write(f"Warning: Pixel ({x}, {y}) using empty data placeholder")

    if not has_valid_data:
        raise ValueError("Generated file contains no valid data")
    print(f"\nConversion completed! File saved to: {output_file}")

def batch_process_parents(mzml_file: str,
                         location_file: str,
                         output_dir: str,
                         parent_mz_list: list,
                         time_window: float,
                         start_rt: float,
                         parent_tolerance: float = 0.02):
    """
    Batch process multiple precursor m/z values, generating separate imzML files for each precursor
    
    Parameters:
        mzml_file: mzML file path
        location_file: coordinate file path
        output_dir: output directory
        parent_mz_list: list of precursor m/z values
        time_window: time window per pixel (seconds)
        start_rt: start retention time (seconds)
        parent_tolerance: precursor m/z tolerance
    """
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    total_start_time = time.time()
    
    print(f"Starting batch processing of {len(parent_mz_list)} precursors...")
    
    for i, parent_mz in enumerate(parent_mz_list, 1):
        print(f"\n[{i}/{len(parent_mz_list)}] Processing precursor m/z = {parent_mz:.4f}")
        
        base_name = os.path.basename(mzml_file).split('.')[0]
        output_file = os.path.join(output_dir, f"{base_name}-{parent_mz:.4f}.imzML")
        
        try:
            print(f"Generating imzML file: {os.path.basename(output_file)}")
            create_imzml(
                mzml_file=mzml_file,
                location_file=location_file,
                output_file=output_file,
                time_window=time_window,
                start_rt=start_rt,
                ms_level=2,
                target_parent_mz=parent_mz,
                parent_tolerance=parent_tolerance
            )
            print(f"Precursor m/z = {parent_mz:.4f} processing completed")
        except Exception as e:
            print(f"Error processing precursor m/z = {parent_mz:.4f}: {str(e)}")
    
    total_time = time.time() - total_start_time
    print(f"\nBatch processing completed! Total time: {total_time:.1f} seconds")

def read_parent_mz_from_csv(csv_file: str) -> list:
    """Read precursor m/z values list from CSV file"""
    try:
        df = pd.read_csv(csv_file)
        
        if 'parent_mz' in df.columns:
            parent_column = 'parent_mz'
        elif 'm/z' in df.columns:
            parent_column = 'm/z'
        else:
            print("CSV file must contain 'parent_mz' or 'm/z' column")
            return []
        
        parent_mz_list = df[parent_column].astype(float).tolist()
        return parent_mz_list
    
    except Exception as e:
        print(f"Error reading CSV file: {str(e)}")
        return []

if __name__ == "__main__":
    # Select operation mode
    print("Select operation mode:")
    print("1 - Process single precursor")
    print("2 - Batch process multiple precursors")
    mode = input("Enter choice (1/2): ")
    
    # Get basic parameters
    input_mzml = input("Enter mzML file path: ")
    coordinate_csv = input("Enter coordinate file path: ")
    time_window = float(input("Enter time window per pixel (seconds): "))
    start_rt = float(input("Enter start retention time (seconds): "))
    
    if mode == "1":
        # Single precursor processing mode
        output_imzml = input("Enter imzML output file path: ")
        target_parent_mz = float(input("Enter target precursor m/z value: "))
        parent_tolerance = float(input("Enter precursor m/z tolerance (e.g., 0.02): "))
        
        create_imzml(
            mzml_file=input_mzml,
            location_file=coordinate_csv,
            output_file=output_imzml,
            time_window=time_window,
            start_rt=start_rt,
            ms_level=2,
            target_parent_mz=target_parent_mz,
            parent_tolerance=parent_tolerance
        )
    
    elif mode == "2":
        # Batch process multiple precursors
        output_dir = input("Enter output directory path: ")
        parent_tolerance = float(input("Enter precursor m/z tolerance (e.g., 0.02): "))
        
        print("\nSelect input method for precursor list:")
        print("1 - Read from CSV file")
        print("2 - Manual input of multiple precursor m/z values")
        input_method = input("Enter choice (1/2): ")
        
        parent_mz_list = []
        
        if input_method == "1":
            # Read from CSV
            csv_file = input("Enter CSV file path containing precursor m/z values: ")
            parent_mz_list = read_parent_mz_from_csv(csv_file)
            if not parent_mz_list:
                print("Failed to read valid precursor m/z values from CSV, exiting")
                exit(1)
            print(f"Read {len(parent_mz_list)} precursor m/z values from CSV")
        
        elif input_method == "2":
            # Manual input
            try:
                num_parents = int(input("Enter number of precursors to process: "))
                print("Enter each precursor m/z value:")
                for i in range(1, num_parents + 1):
                    parent_mz = float(input(f"Precursor {i}: "))
                    parent_mz_list.append(parent_mz)
            except ValueError:
                print("Invalid input, exiting")
                exit(1)
        else:
            print("Invalid choice, exiting")
            exit(1)
        
        # Execute batch processing
        if parent_mz_list:
            batch_process_parents(
                mzml_file=input_mzml,
                location_file=coordinate_csv,
                output_dir=output_dir,
                parent_mz_list=parent_mz_list,
                time_window=time_window,
                start_rt=start_rt,
                parent_tolerance=parent_tolerance
            )
        else:
            print("No precursor m/z values to process, exiting")
    
    else:
        print("Invalid choice, exiting")
