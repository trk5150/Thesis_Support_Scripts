import os
import re
import cv2
import numpy as np
import tifffile as tiff
import pandas as pd

def process_images_in_subdir(subdir_path):
    """Processes images in a given subdirectory and returns a DataFrame with results."""
    
    # Regex pattern to match filenames and extract the channel number
    pattern = re.compile(r'.*-Image Export-\d+_c(\d+)')

    # Define independent threshold values for each channel
    channel_thresholds = {
        1: 30,  # Example threshold for channel 1
        2: 25,  # Example threshold for channel 2
        3: 50   # Example threshold for channel 3
    }

    results = []

    # Create a folder for processed images
    processed_dir = os.path.join(subdir_path, "Processed_Masks")
    os.makedirs(processed_dir, exist_ok=True)

    # Iterate through files in the given subdirectory
    for filename in os.listdir(subdir_path):
        match = pattern.match(filename)
        if match:
            channel = int(match.group(1))
            if channel in channel_thresholds:
                file_path = os.path.join(subdir_path, filename)

                # Load the image
                image = tiff.imread(file_path)

                # Convert to grayscale (if necessary)
                if len(image.shape) > 2:
                    image = cv2.cvtColor(image, cv2.COLOR_RGB2GRAY)

                # Apply channel-specific thresholding
                threshold_value = channel_thresholds[channel]
                _, binary_mask = cv2.threshold(image, threshold_value, 255, cv2.THRESH_BINARY)

                # Save the binary mask as a new TIFF file
                mask_filename = f"Mask_{filename}"
                mask_path = os.path.join(processed_dir, mask_filename)
                tiff.imwrite(mask_path, binary_mask.astype(np.uint8))  # Save in grayscale (uint8)

                # Calculate stained area (count nonzero pixels)
                stained_area = np.count_nonzero(binary_mask)

                results.append({
                    'Subdirectory': os.path.basename(subdir_path),  # Include subdirectory name
                    'Filename': filename,
                    'Channel': channel,
                    'Threshold': threshold_value, 
                    'Stained Area': stained_area,
                    'Mask File': mask_filename
                })

    # Convert results to DataFrame
    df = pd.DataFrame(results)
    return df

def process_all_subdirectories(main_directory):
    """Cycles through all subdirectories in a given directory and runs analysis for each."""
    
    # Initialize an empty DataFrame to store results from all subdirectories
    all_results = pd.DataFrame()

    # Iterate through all subdirectories in the main directory
    for subdir in os.listdir(main_directory):
        subdir_path = os.path.join(main_directory, subdir)

        # Ensure it's a directory before processing
        if os.path.isdir(subdir_path):
            print(f"Processing: {subdir_path}")

            # Process the subdirectory and get results
            df_results = process_images_in_subdir(subdir_path)

            # Append to cumulative DataFrame
            all_results = pd.concat([all_results, df_results], ignore_index=True)

    # Save final DataFrame to an Excel file
    excel_path = os.path.join(main_directory, "Final_Analysis_Results.xlsx")
    all_results.to_excel(excel_path, index=False)
    
    print(f"Final results saved to: {excel_path}")

    return all_results

# Example usage
main_directory = "D://amylas, ins and gcg costains"
df_final_results = process_all_subdirectories(main_directory)
