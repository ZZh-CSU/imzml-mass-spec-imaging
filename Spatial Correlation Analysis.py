import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import pearsonr
import seaborn as sns
import os

def calculate_moran_I(image):
    """Calculate Global Moran's I spatial autocorrelation coefficient"""
    image = np.array(image, dtype=float)
    rows, cols = image.shape
    N = rows * cols
    
    mean = np.nanmean(image)
    z = image - mean
    sum_sq_dev = np.nansum(z**2)
    
    if sum_sq_dev == 0: return 0.0 
    
    numerator = 0.0
    w_sum = 0.0 
    shifts = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    
    for dy, dx in shifts:
        shifted_z = np.roll(z, shift=(dy, dx), axis=(0, 1))
        
        mask = np.ones_like(z, dtype=bool)
        if dy == 1: mask[0, :] = False
        if dy == -1: mask[-1, :] = False
        if dx == 1: mask[:, 0] = False
        if dx == -1: mask[:, -1] = False
        
        term = z * shifted_z
        valid = mask & ~np.isnan(z) & ~np.isnan(shifted_z)
        
        numerator += np.sum(term[valid])
        w_sum += np.sum(valid) 
        
    if w_sum == 0: return 0.0
    moran_I = (N / w_sum) * (numerator / sum_sq_dev)
    return moran_I

class SpatialCorrelationAnalysis:
    def __init__(self, img1, img2):
        self.img1 = np.array(img1, dtype=float)
        self.img2 = np.array(img2, dtype=float)
        
    def diagnose_autocorrelation(self):
        print("--- 1. Diagnose Spatial Autocorrelation (Moran's I) ---")
        mi = calculate_moran_I(self.img1)
        ml = calculate_moran_I(self.img2)
        print(f"Global Moran's I (Image 1): {mi:.4f}")
        print(f"Global Moran's I (Image 2): {ml:.4f}")
        return mi, ml

    def _get_blocks(self, img, block_size):
        h, w = img.shape
        pad_h = (block_size - h % block_size) % block_size
        pad_w = (block_size - w % block_size) % block_size
        
        padded_img = np.pad(img, ((0, pad_h), (0, pad_w)), mode='constant', constant_values=np.nan)
        
        n_rows = padded_img.shape[0] // block_size
        n_cols = padded_img.shape[1] // block_size
        
        blocks = padded_img.reshape(n_rows, block_size, n_cols, block_size).transpose(0, 2, 1, 3)
        blocks_flat = blocks.reshape(-1, block_size, block_size)
        return blocks_flat, (h, w)

    def _reconstruct_from_blocks(self, blocks, original_shape, block_size):
        h, w = original_shape
        pad_h = (block_size - h % block_size) % block_size
        pad_w = (block_size - w % block_size) % block_size
        full_h = h + pad_h
        full_w = w + pad_w
        
        n_cols = full_w // block_size
        n_rows = full_h // block_size
        
        reshaped = blocks.reshape(n_rows, n_cols, block_size, block_size)
        transposed = reshaped.transpose(0, 2, 1, 3)
        reconstructed = transposed.reshape(full_h, full_w)
        return reconstructed[:h, :w]

    def run_permutation_test(self, block_size=3, n_perms=999):
        print(f"\n--- 2. Start Spatial Block Permutation Test (Block Size={block_size}, N={n_perms}) ---")
        mask_obs = ~np.isnan(self.img1) & ~np.isnan(self.img2)
        v1 = self.img1[mask_obs]
        v2 = self.img2[mask_obs]
        
        if len(v1) < 2:
            print("Error: Not enough valid data to calculate correlation")
            return None
            
        r_obs, _ = pearsonr(v1, v2)
        print(f"Observed Pearson r: {r_obs:.4f}")
        
        blocks1, shape = self._get_blocks(self.img1, block_size)
        n_blocks = blocks1.shape[0]
        
        r_perms = []
        for i in range(n_perms):
            perm_indices = np.random.permutation(n_blocks)
            shuffled_blocks = blocks1[perm_indices]
            img1_perm = self._reconstruct_from_blocks(shuffled_blocks, shape, block_size)
            
            mask_perm = ~np.isnan(img1_perm) & ~np.isnan(self.img2)
            v1_p = img1_perm[mask_perm]
            v2_p = self.img2[mask_perm]
            
            if len(v1_p) > 2:
                try:
                    r_p, _ = pearsonr(v1_p, v2_p)
                    if np.isnan(r_p): r_p = 0
                    r_perms.append(r_p)
                except:
                    r_perms.append(0)
            else:
                 r_perms.append(0)

        r_perms = np.array(r_perms)
        p_value = (np.sum(np.abs(r_perms) >= np.abs(r_obs)) + 1) / (n_perms + 1)
        print(f"Corrected P-value: {p_value:.5f}")
        
        self.result = {'r_obs': r_obs, 'p_value': p_value, 'r_perms': r_perms}
        return self.result

    def plot_results(self):
        if not hasattr(self, 'result'): return
        r_obs = self.result['r_obs']
        r_perms = self.result['r_perms']
        p_val = self.result['p_value']
        
        fig, axes = plt.subplots(1, 2, figsize=(12, 5))
        mask = ~np.isnan(self.img1) & ~np.isnan(self.img2)
        
        sns.scatterplot(x=self.img1[mask], y=self.img2[mask], ax=axes[0], alpha=0.6, color='#1f77b4')
        axes[0].set_title(f"Intensity Correlation\nObserved r = {r_obs:.3f}")
        axes[0].set_xlabel("Ion 1 (e.g., LPC-PC)")
        axes[0].set_ylabel("Ion 2 (e.g., PLA2)")
        
        axes[1].hist(r_perms, bins=30, color='gray', density=True, alpha=0.7)
        axes[1].axvline(r_obs, color='red', linestyle='--', linewidth=2, label=f'Observed r')
        axes[1].set_title(f"Null Distribution\nCorrected P-value = {p_val:.5f}")
        axes[1].legend()
        plt.tight_layout()
        plt.show()

script_dir = os.path.dirname(os.path.abspath(__file__))
print(f"Searching in directory: {script_dir}")

file1_path = os.path.join(script_dir, '1.xlsx')
file2_path = os.path.join(script_dir, '2.xlsx')

try:
    if not os.path.exists(file1_path) or not os.path.exists(file2_path):
        print("\n================ Error: File Not Found ================")
        print(f"Please check directory: {script_dir}")
        print("Please ensure files are renamed to '1.xlsx' and '2.xlsx'")
        print("=======================================================")
        raise FileNotFoundError

    print("Reading Excel files, please wait...")
    
    df1 = pd.read_excel(file1_path, header=None, engine='openpyxl')
    df2 = pd.read_excel(file2_path, header=None, engine='openpyxl')

    print(f"Data loaded successfully! Dimensions: {df1.shape} & {df2.shape}")

    img1 = df1.apply(pd.to_numeric, errors='coerce').values
    img2 = df2.apply(pd.to_numeric, errors='coerce').values
    
    analysis = SpatialCorrelationAnalysis(img1, img2)
    analysis.diagnose_autocorrelation()

    analysis.run_permutation_test(block_size=3, n_perms=999)
    
    analysis.plot_results()
    
except Exception as e:
    print(f"\nProgram terminated: {e}")