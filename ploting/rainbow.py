import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.cm as cm

# 1. Load your data
# Change 'your_data.csv' to your actual file name
df = pd.read_csv('your_data.csv') 

# Adjust based on your Prism structure
wavelength = df.iloc[:, 0]
intensity_data = df.iloc[:, 1:] 
num_timepoints = intensity_data.shape[1]

# 2. Setup the Plot with explicit Axes
fig, ax = plt.subplots(figsize=(10, 6))

# Use 'nipy_spectral' or 'rainbow_r' for Red -> Violet
cmap = plt.get_cmap('rainbow_r') 

# 3. Plot each column
for i, column_name in enumerate(intensity_data.columns):
    color = cmap(i / (num_timepoints - 1))
    ax.plot(wavelength, intensity_data[column_name], color=color)

# 4. Formatting
ax.set_xlabel('Wavelength (nm)')
ax.set_ylabel('Intensity')
ax.set_title('Spectral Data Over Time (Red=0s, Violet=End)')
ax.grid(True, linestyle='--', alpha=0.6)

# 5. Fix the Colorbar error
# We explicitly tell the colorbar to steal space from 'ax'
sm = cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=0, vmax=num_timepoints))
fig.colorbar(sm, ax=ax, label='Timepoint Index')

plt.tight_layout()
plt.show()
