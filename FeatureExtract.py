import numpy as np
import pandas as pd
import sys
from skimage import measure
from sklearn.decomposition import PCA
import imageio


# Voxel size (scale) in micrometers (z, y, x)
voxel_size = (0.554, 0.554, 0.554)

toreadfrom = sys.argv[1] 
toprintto = sys.argv[2] 

label_image_stack = imageio.volread(toreadfrom)

'''
stardist_labels_layer = skio.imread(toreadfrom, plugin='tifffile')

#access stardist labels layer
#stardist_labels_layer = toreadfrom

# Extract data from the layer
label_image_stack = stardist_labels_layer.data
'''

# Function to perform PCA and get orientation axes and lengths
def get_pca_results(coords):
    pca = PCA(n_components=3)
    pca.fit(coords)
    # Eigenvalues give us measures for the axes lengths (variance along each axis)
    axes_lengths = np.sqrt(pca.explained_variance_)
    orientation_axes = pca.components_
    return orientation_axes, axes_lengths

# Prepare a list to store region properties
data = []

# Compute region properties for the 3D labeled image
regions = measure.regionprops(label_image_stack)

# Gather region properties
for i, region in enumerate(regions):
    # Calculate volume
    volume = region.area * np.prod(voxel_size)  # Convert voxel count to micrometers^3

    # Calculate orientation and axes lengths using PCA
    coordsu = region.coords * voxel_size  # Scale the coordinates
    orientation_axes, axes_lengths = get_pca_results(coordsu)

    # Get centroid
    centroid = np.array(region.centroid) * voxel_size
    
    coords_shift = coordsu - centroid
    xvals = coords_shift[:, 0]
    yvals = coords_shift[:, 1]
    zvals = coords_shift[:, 2]
    rvals = np.sqrt(xvals * xvals + yvals * yvals + zvals * zvals)
    xmax = xvals[np.argmax(rvals)]
    ymax = yvals[np.argmax(rvals)]
    zmax = zvals[np.argmax(rvals)]
    
    #test if orientations axes is in same direction or not, by dot product
    dotprod = orientation_axes[0][0] * xmax + orientation_axes[0][1] * ymax + orientation_axes[0][2] * zmax
    if (dotprod < 0):
        orientation_axes[0] = - orientation_axes[0]
        
    dotprod = orientation_axes[1][0] * xmax + orientation_axes[1][1] * ymax + orientation_axes[1][2] * zmax
    if (dotprod < 0):
        orientation_axes[1] = - orientation_axes[1]
        
    dotprod = orientation_axes[2][0] * xmax + orientation_axes[2][1] * ymax + orientation_axes[2][2] * zmax
    if (dotprod < 0):
        orientation_axes[2] = - orientation_axes[2]


    # Append data
    data.append({
        'Object': i + 1,
        'Volume (micrometers^3)': volume,
        'Orientation Axes0': orientation_axes[0].tolist(),  # Convert ndarray to list for CSV
        'Orientation Axes1': orientation_axes[1].tolist(),  # Convert ndarray to list for CSV
        'Orientation Axes2': orientation_axes[2].tolist(),  # Convert ndarray to list for CSV
        'Axes Lengths (micrometers)': axes_lengths.tolist(),  # Record the axes lengths
        'Centroid (micrometers)': centroid.tolist()     # Convert ndarray to list for CSV
    })

# Convert the list to a pandas DataFrame
df = pd.DataFrame(data)

# Write the DataFrame to a CSV file
csv_filename = toprintto
df.to_csv(csv_filename, index=False)

print(f"Region properties with axes lengths (without area) have been saved to {csv_filename}")
