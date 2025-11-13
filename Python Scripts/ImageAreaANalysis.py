# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 17:41:33 2025

@author: tik105
"""

import numpy as np
import cv2
import matplotlib.pyplot as plt
from aicspylibczi import CziFile

# Load the .czi file
czi_path = "C://Users//tik105//Desktop//Litter 1 pup 1 slide 1-1.czi"
czi = CziFile(czi_path)

# Extract image array
image_array = czi.asarray()

# Check dimensions
print(image_array.shape)  # Typically (T, Z, C, Y, X) or (C, Y, X)

# Extract specific channels (assuming C dimension exists)
red_channel = image_array[2, :, :]  # Adjust based on channel order
green_channel = image_array[1, :, :]
purple_channel = image_array[3, :, :]

# Display extracted channels
plt.figure(figsize=(12, 4))
plt.subplot(1, 3, 1)
plt.imshow(red_channel, cmap="Reds")
plt.title("Amylase (Red)")

plt.subplot(1, 3, 2)
plt.imshow(green_channel, cmap="Greens")
plt.title("Insulin (Green)")

plt.subplot(1, 3, 3)
plt.imshow(purple_channel, cmap="Purples")
plt.title("Glucagon (Purple)")

plt.show()
