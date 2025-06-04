import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import io
from scipy.signal import butter, filtfilt
import cv2
import os
from tqdm import tqdm

# Configure matplotlib for better performance on M1
matplotlib.use('agg')  # Non-interactive backend, better for video generation
plt.rcParams['figure.dpi'] = 100

def create_apple_eeg_video(mat_file, out_file='eeg_video.mp4', window_size=6, fps=30):
    """
    Create an Apple Health-style EEG video from a .mat file
    
    Parameters:
    -----------
    mat_file : str
        Path to the .mat file containing EEG data
    out_file : str
        Path to save the output video file (default: 'eeg_video.mp4')
    window_size : float
        Visible time window in seconds (default: 6)
    fps : int
        Frames per second (default: 30)
    """
    print(f"Loading data from {mat_file}...")
    
    # Load MAT file
    try:
        data = io.loadmat(mat_file)
        if 'raw_data' in data:
            # Get raw_data and extract the first row (channel)
            raw_signal = data['raw_data']
            if raw_signal.ndim > 1:
                raw_signal = raw_signal[0, :]  # Extract first row/channel
            # Ensure signal is a 1D array
            raw_signal = np.squeeze(raw_signal)
        else:
            raise ValueError(f"Variable 'raw_data' not found in {mat_file}")
    except Exception as e:
        print(f"Error loading MAT file: {e}")
        return
    
    # Constants
    fs = 1000  # Given sampling frequency (Hz)
    
    # Signal processing
    print("Processing signal...")
    # Apply bandpass filter (1-100 Hz)
    nyquist = 0.5 * fs
    b, a = butter(4, [1/nyquist, 100/nyquist], btype='band')
    filtered_signal = filtfilt(b, a, raw_signal)
    
    # Normalize signal to center it in the middle of the grid
    # Use amplitude scaling to keep it visible but not too large
    signal_min = np.min(filtered_signal)
    signal_max = np.max(filtered_signal)
    signal_range = signal_max - signal_min
    
    # Center around 0.5 with amplitude scaled to use 80% of the vertical space
    vertical_scale = 0.40  # Use 40% above and below center
    signal = 0.5 + ((filtered_signal - (signal_min + signal_range/2)) / signal_range) * vertical_scale
    
    # Create scale bar image with same scale as video
    scale_bar_file = out_file.replace('.mp4', '_scalebar.png')
    create_scale_bar_image(signal_min, signal_max, vertical_scale, 
                          width=200, height=height, out_file=scale_bar_file)
    
    # Calculate duration and sample counts
    n_samples = len(signal)
    duration = (n_samples - 1) / fs
    print(f"Signal duration: {duration:.2f} seconds ({n_samples} samples)")
    
    # Video settings
    # Base width for 10.2:1 aspect ratio (compatible with 1260x720 format)
    base_width = 1260
    base_height = int(base_width / 10.2)
    # Ensure dimensions are even (required by some codecs)
    width = base_width - (base_width % 2)
    height = base_height - (base_height % 2)
    print(f"Video dimensions: {width}x{height} (10.2:1 aspect ratio)")
    
    # Grid parameters - matching createAppleGrid.m exactly
    major = 0.2   # 200ms
    minor = 0.04  # 40ms
    second = 1.0  # 1 second
    
    # Grid colors - matching createAppleGrid.m exactly
    major_color = np.array([173, 172, 176]) / 255
    minor_color = np.array([236, 236, 236]) / 255
    second_color = np.array([140, 140, 140]) / 255
    
    # Line colors
    line_color = np.array([232, 60, 91]) / 255  # Apple Health pink
    line_width = 2.7  # Reduced by 33% from 4
    
    # Create time vector
    time = np.arange(n_samples) / fs
    
    # Calculate frames
    samples_per_frame = int(fs / fps)
    total_frames = int(np.ceil(n_samples / samples_per_frame))
    
    # Initialize video writer
    fourcc = cv2.VideoWriter_fourcc(*'mp4v')
    video = cv2.VideoWriter(out_file, fourcc, fps, (width, height))
    
    print(f"Generating {total_frames} frames...")
    
    # Generate frames
    for frame_idx in tqdm(range(total_frames)):
        # Calculate sample indices for current frame
        end_sample = min(n_samples, (frame_idx + 1) * samples_per_frame)
        current_time = time[end_sample - 1]
        
        # Calculate visible window
        if current_time > window_size:
            visible_start = current_time - window_size
            visible_end = current_time
        else:
            visible_start = 0
            visible_end = window_size
        
        # Create figure
        fig, ax = plt.figure(figsize=(width/100, height/100), dpi=100, 
                           facecolor='white'), plt.axes([0, 0, 1, 1])
        ax.set_xlim(visible_start, visible_end)
        ax.set_ylim(0, 1)
        ax.set_aspect('auto')
        ax.axis('off')
        
        # Draw grid exactly like createAppleGrid.m
        # 1. Minor grid lines first (bottom layer)
        for x in np.arange(np.floor(visible_start/minor)*minor, visible_end+minor, minor):
            if abs(x % major) >= 1e-12:  # Not a major gridline
                ax.axvline(x, color=minor_color, linewidth=0.6, zorder=1)
                
        for y in np.arange(0, 1+minor, minor):
            if abs(y % major) >= 1e-12:  # Not a major gridline
                ax.axhline(y, color=minor_color, linewidth=0.6, zorder=1)
        
        # 2. Major grid lines (middle layer)
        for x in np.arange(np.floor(visible_start/major)*major, visible_end+major, major):
            if abs(x % second) >= 1e-12:  # Not a 1-second gridline
                ax.axvline(x, color=major_color, linewidth=1.2, zorder=2)
                
        for y in np.arange(0, 1+major, major):
            ax.axhline(y, color=major_color, linewidth=1.2, zorder=2)
        
        # 3. Second grid lines (top layer, thickest)
        for x in np.arange(np.floor(visible_start/second)*second, visible_end+second, second):
            ax.axvline(x, color=second_color, linewidth=1.8, zorder=3)
        
        # Find visible signal segment
        visible_start_idx = np.searchsorted(time, visible_start)
        visible_signal = time[visible_start_idx:end_sample]
        visible_values = signal[visible_start_idx:end_sample]
        
        # Plot signal
        ax.plot(visible_signal, visible_values, color=line_color, linewidth=line_width, zorder=4)
        
        # Convert to image - using the new buffer_rgba method instead of deprecated tostring_rgb
        fig.canvas.draw()
        w, h = fig.canvas.get_width_height()
        buf = fig.canvas.buffer_rgba()
        img = np.frombuffer(buf, dtype=np.uint8).reshape(h, w, 4)
        
        # Convert RGBA to RGB by removing the alpha channel
        img = img[:, :, :3]
        
        # Convert RGB to BGR (for OpenCV)
        img = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        
        # Write to video
        video.write(img)
        
        # Close figure to free memory
        plt.close(fig)
    
    # Finalize video
    video.release()
    print(f"Video saved to {out_file}")

def create_scale_bar_image(signal_min, signal_max, vertical_scale, width=200, height=600, out_file=None):
    """
    Create a simplified scale bar image showing only the essential amplitude references
    
    Parameters:
    -----------
    signal_min : float
        Minimum value of the filtered signal
    signal_max : float
        Maximum value of the filtered signal
    vertical_scale : float
        Vertical scaling factor (e.g., 0.40)
    width : int
        Width of the scale bar image (default: 200)
    height : int
        Height of the scale bar image (default: 600)
    out_file : str
        Output file path for the scale bar image
    """
    # Calculate the actual voltage range that corresponds to the normalized display
    signal_range = signal_max - signal_min
    signal_center = signal_min + signal_range/2
    
    # The normalized range is 0.5 ± vertical_scale
    # So 0.1 to 0.9 if vertical_scale = 0.4
    norm_min = 0.5 - vertical_scale
    norm_max = 0.5 + vertical_scale
    norm_center = 0.5  # Zero baseline
    
    # Create figure for scale bar
    fig, ax = plt.subplots(figsize=(width/100, height/100), dpi=100, facecolor='black')
    ax.set_xlim(0, 1)
    ax.set_ylim(norm_min, norm_max)
    ax.set_facecolor('black')
    
    # Calculate the three key voltage values
    # Top value (maximum)
    top_value = signal_center + (vertical_scale * signal_range)
    # Zero baseline (center)
    zero_value = signal_center  
    # Bottom value (minimum)
    bottom_value = signal_center - (vertical_scale * signal_range)
    
    # Draw only three essential reference lines and labels
    reference_points = [
        (norm_max, top_value, "TOP"),
        (norm_center, zero_value, "ZERO"),
        (norm_min, bottom_value, "BOTTOM")
    ]
    
    for norm_pos, voltage_val, label_type in reference_points:
        # Draw horizontal reference line across the width
        ax.plot([0.05, 0.4], [norm_pos, norm_pos], color='white', linewidth=2)
        
        # Add voltage value label
        ax.text(0.45, norm_pos, f'{voltage_val:.1f}', color='white', 
               fontsize=14, va='center', fontweight='bold')
        
        # Add descriptive label for zero baseline
        if label_type == "ZERO":
            ax.text(0.75, norm_pos, '0', color='yellow', 
                   fontsize=16, va='center', fontweight='bold')
    
    # Draw the main scale line
    ax.plot([0.2, 0.2], [norm_min, norm_max], color='white', linewidth=3)
    
    # Add title
    ax.text(0.5, norm_max + 0.05, 'Amplitude Scale', color='white', 
           fontsize=14, ha='center', fontweight='bold')
    
    # Add unit label (assuming microvolts, adjust as needed)
    ax.text(0.5, norm_min - 0.05, 'μV', color='white', 
           fontsize=12, ha='center', fontweight='bold')
    
    ax.axis('off')
    
    # Convert to image
    fig.canvas.draw()
    w, h = fig.canvas.get_width_height()
    buf = fig.canvas.buffer_rgba()
    img = np.frombuffer(buf, dtype=np.uint8).reshape(h, w, 4)
    
    # Convert RGBA to RGB
    img = img[:, :, :3]
    
    # Save as image if output file specified
    if out_file:
        # Convert RGB to BGR for OpenCV
        img_bgr = cv2.cvtColor(img, cv2.COLOR_RGB2BGR)
        cv2.imwrite(out_file, img_bgr)
        print(f"Simplified scale bar image saved to {out_file}")
    
    plt.close(fig)
    return img

def main():
    """Command line interface"""
    import argparse
    parser = argparse.ArgumentParser(description='Create Apple Health-style EEG video')
    parser.add_argument('mat_file', help='Path to .mat file containing EEG data')
    parser.add_argument('--output', '-o', default='eeg_video.mp4', help='Output video file')
    parser.add_argument('--window', '-w', type=float, default=6, help='Visible window size in seconds')
    parser.add_argument('--fps', '-f', type=int, default=30, help='Frames per second')
    args = parser.parse_args()
    
    create_apple_eeg_video(args.mat_file, args.output, args.window, args.fps)

if __name__ == '__main__':
    main()