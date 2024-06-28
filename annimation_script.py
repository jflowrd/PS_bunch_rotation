import os
import matplotlib.pyplot as plt
import matplotlib.animation as animation
from PIL import Image

# Directory containing the images
image_dir = 'fig'

# Create a list of image filenames
image_files = [f'long_distr_{i}.png' for i in range(5, 146, 5)]

# Ensure all the images exist
image_files = [os.path.join(image_dir, f) for f in image_files if os.path.isfile(os.path.join(image_dir, f))]

# Create a figure and axis
fig, ax = plt.subplots()

# Remove axes
ax.axis('off')

# Load the first image to get the dimensions
first_image = Image.open(image_files[0])
im = ax.imshow(first_image, aspect='auto')

def update(frame):
    """Update the image for the animation."""
    image = Image.open(frame)
    im.set_array(image)
    return [im]

# Create the animation
ani = animation.FuncAnimation(fig, update, frames=image_files, blit=True, repeat=True)

# Save the animation as a high-quality GIF
ani.save('animation.gif', writer='pillow', fps=3, dpi=300)

plt.show()



