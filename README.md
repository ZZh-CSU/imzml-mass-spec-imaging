# imzml-mass-spec-imaging

This repository provides Python scripts for mass spectrometry imaging (MSI) data processing and visualization. The workflow consists of converting mzML files to imzML format and then visualizing the imaging data based on m/z and ppm values. It supports both MS1 and MS2 (Tandem MS) data workflows.

## Features

- **mzML to imzML conversion**: Convert time-series mass spectrometry data (mzML) into spatial imaging data (imzML) using pixel coordinates and time windows.
- **MSI visualization**: Read imzML files and plot intensity maps for specified m/z values and ppm tolerances. Supports batch processing from CSV.
- **MS2 Support**: Dedicated tools for processing MS2 data, allowing filtering by precursor ions and visualizing specific fragment ions.

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

### MS1 Workflow (Standard)

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

2. **Visualize imzML data**

   Use `read_imzml_and_draw_picture.py` to plot MSI images for specific m/z and ppm values.

   ```bash
   python read_imzml_and_draw_picture.py
   ```

   You can choose:

   - Manual mode: input m/z and ppm values one by one.
   - Batch mode: provide a CSV file with columns `m/z` and `ppm` for batch processing.

### MS2 Workflow (Tandem MS)

For MS2 analysis, use the scripts located in the `MS2` folder.

1. **Convert mzML to imzML (MS2)**

   Use `MS2/imzwriter_by_time_second_mzml_MS2.py` to convert mzML files. This script allows you to filter spectra based on precursor m/z values.

   ```bash
   python MS2/imzwriter_by_time_second_mzml_MS2.py
   ```

   **Modes:**

   - **Single Precursor**: Extract data for one specific precursor ion.
   - **Batch Processing**: Process a list of precursor ions (manually or from CSV) to generate multiple imzML files automatically.

2. **Visualize MS2 imzML data**

   Use `MS2/read_imzml_and_draw_picture_MS2.py` to visualize fragment ions.

   ```bash
   python MS2/read_imzml_and_draw_picture_MS2.py
   ```

   **Modes:**

   - **Single File**: Process a single imzML file for specific fragments (manual input or CSV).
   - **Batch Processing**: Process a directory of imzML files against a directory of fragment CSV files. This is useful when you have generated multiple precursor-specific imzML files in step 1.

## File Descriptions

- `imzwriter_by_time_second_from_mzml.py`: Converts mzML files to imzML format (Standard/MS1).
- `read_imzml_and_draw_picture.py`: Visualizes MSI images for specified m/z and ppm values (Standard/MS1).
- `MS2/imzwriter_by_time_second_mzml_MS2.py`: Converts mzML to imzML with support for MS2 precursor filtering and batch processing.
- `MS2/read_imzml_and_draw_picture_MS2.py`: Visualizes MS2 fragment ions, supporting batch processing of multiple precursor files.

## Notes

- The coordinate CSV file must contain columns named `X` and `Y`.
- For batch visualization (MS1), the CSV file must contain columns `m/z` and `ppm`.
- For MS2 batch visualization, fragment CSV files should contain `fragment` and `tolerance` columns.
- Make sure your input files are accessible and paths are correct.

## License

MIT License
