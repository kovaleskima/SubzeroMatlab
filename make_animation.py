import cv2
import os
from glob import glob

# === CONFIGURATION ===
input_folder = './figs/'   # <- change this to your folder
output_file = 'animation.mp4'
frame_rate = 10  # frames per second

# === LOAD PNG FILES ===
image_files = sorted(glob(os.path.join(input_folder, '*.jpg')))

# Check if we found any files
if not image_files:
    raise FileNotFoundError(f"No JPG files found in {input_folder}")

# Read the first image to get the dimensions
frame = cv2.imread(image_files[0])
height, width, layers = frame.shape

# === VIDEO WRITER ===
fourcc = cv2.VideoWriter_fourcc(*'mp4v')  # Use 'avc1' or 'H264' if needed
video = cv2.VideoWriter(output_file, fourcc, frame_rate, (width, height))

for image_file in image_files:
    frame = cv2.imread(image_file)
    if frame is None:
        print(f"Warning: Could not read {image_file}, skipping.")
        continue
    video.write(frame)

video.release()
print(f"Video saved to {output_file}")
