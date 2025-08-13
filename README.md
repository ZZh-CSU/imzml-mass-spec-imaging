# imzml-mass-spec-imaging

This repository provides two Python scripts for mass spectrometry imaging (MSI) data processing and visualization. The workflow consists of converting mzML files to imzML format and then visualizing the imaging data based on m/z and ppm values.

## Features

- **mzML to imzML conversion**: Convert time-series mass spectrometry data (mzML) into spatial imaging data (imzML) using pixel coordinates and time windows.
- **MSI visualization**: Read imzML files and plot intensity maps for specified m/z values and ppm tolerances. Supports batch processing from CSV.

## Requirements

- Python 3.7+
- numpy
- pandas
- tqdm
- matplotlib
- pyimzml

Install dependencies with:

```bash
pip install numpy pandas tqdm matplotlib pyimzml
```

## Workflow

1. **Convert mzML to imzML**

   Use `imzwriter_by_time_second_from_mzml.py` to convert your mzML file to imzML format.

   ```bash
   python imzwriter_by_time_second_from_mzml.py
   ```

   You will be prompted to enter:

   - mzML file path
   - coordinate CSV file path (must contain columns `X` and `Y`)
   - output imzML file path
   - time window per pixel (seconds)
   - start retention time (seconds)
   - MS level (1 or 2)

   The script will process the mzML file and generate an imzML file for imaging.

2. **Visualize imzML data**

   Use `read_imzml_and_draw_picture.py` to plot MSI images for specific m/z and ppm values.

   ```bash
   python read_imzml_and_draw_picture.py
   ```

   You can choose:

   - Manual mode: input m/z and ppm values one by one.
   - Batch mode: provide a CSV file with columns `m/z` and `ppm` for batch processing.

   The script will generate intensity maps and optionally save them as PNG and CSV files.

## Example Usage

### Step 1: Convert mzML to imzML

```bash
python imzwriter_by_time_second_from_mzml.py
```

### Step 2: Visualize imzML

```bash
python read_imzml_and_draw_picture.py
```

## File Descriptions

- `imzwriter_by_time_second_from_mzml.py`: Converts mzML files to imzML format using pixel coordinates and time windows.
- `read_imzml_and_draw_picture.py`: Reads imzML files and plots MSI images for specified m/z and ppm values.

## Notes

- The coordinate CSV file must contain columns named `X` and `Y`.
- For batch visualization, the CSV file must contain columns `m/z` and `ppm`.
- Make sure your input files are accessible and paths are correct.

## License

MIT License
