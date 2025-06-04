# This .py was implemented in the 'SI: STIM condition' video, in order to harmonically unify the voltage (uV) scale bar (y-axis).
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy import io
from scipy.signal import butter, filtfilt
import cv2
import os
from tqdm import tqdm
import glob
import h5py

# Configure matplotlib for better performance on M1
matplotlib.use('agg')  # Non-interactive backend, better for video generation
plt.rcParams['figure.dpi'] = 1600

def calculate_scale_positions(signal_min, signal_max, target_bottom, target_top, adaptive_scale=True):
    """
    Calculate scale bar positions and values for both adaptive and fixed scaling modes.
    
    Args:
        signal_min, signal_max: Actual signal range in µV
        target_bottom, target_top: Y-coordinate range for display (0.2 to 0.8)
        adaptive_scale: If True, adapt to signal range. If False, use fixed -5000 to +5000 µV range
    
    Returns:
        scale_positions: Y-coordinates for scale bar marks
        scale_values: µV values for scale bar labels
    """
    # Fixed harmonic scale values
    fixed_harmonic_scale_values = np.array([-5000, -2500, 0, 2500, 5000])
    
    scale_positions = []
    scale_values = []
    
    if adaptive_scale:
        # ADAPTIVE MODE: Only show scale values that are within the visible signal range
        signal_range = signal_max - signal_min
        
        for uv_value in fixed_harmonic_scale_values:
            if signal_range > 0:
                # Calculate where this µV value would appear in the scaled coordinate system
                scale_factor = (target_top - target_bottom) / signal_range
                y_pos = target_bottom + (uv_value - signal_min) * scale_factor
            else:
                y_pos = 0.5  # Center if signal is flat
            
            # Only include values that are within the visible range (with small margin)
            if target_bottom - 0.1 <= y_pos <= target_top + 0.1:
                scale_positions.append(y_pos)
                scale_values.append(uv_value)
    else:
        # FIXED MODE: Always show all 5 reference points at fixed positions
        # Map -5000 µV to target_bottom, +5000 µV to target_top
        for uv_value in fixed_harmonic_scale_values:
            # Linear mapping: -5000 → target_bottom, +5000 → target_top
            y_pos = target_bottom + (uv_value + 5000) / 10000 * (target_top - target_bottom)
            scale_positions.append(y_pos)
            scale_values.append(uv_value)
    
    return np.array(scale_positions), np.array(scale_values)

def create_eeg_video(mat_file, out_file='eeg_video.mp4', window_size=6, fps=30, adaptive_scale=True):
    print(f"Loading data from {mat_file}...")
    
    # Try to determine if file is v7.3 format
    try:
        # First try with scipy (for older .mat files)
        data = io.loadmat(mat_file)
        if 'raw_data' in data:
            raw_signal = data['raw_data']
            if raw_signal.ndim > 1:
                raw_signal = raw_signal[0, :]  # Extract first row/channel
            raw_signal = np.squeeze(raw_signal)
        else:
            raise ValueError(f"Variable 'raw_data' not found in {mat_file}")
    except Exception as e:
        if 'Please use HDF reader' in str(e):
            print(f"File is in MATLAB v7.3 format, using h5py instead...")
            try:
                # Try with h5py for newer v7.3 .mat files
                with h5py.File(mat_file, 'r') as f:
                    # Look for raw_data in the file structure
                    if 'raw_data' in f:
                        # h5py reads MATLAB arrays transposed
                        raw_signal = f['raw_data'][()]
                        if raw_signal.ndim > 1:
                            # Take first channel if multiple channels exist
                            raw_signal = raw_signal[:, 0]  # Note: transposed indexing
                        raw_signal = np.squeeze(raw_signal)
                    else:
                        # If raw_data not found, look for other dataset names
                        datasets = list(f.keys())
                        print(f"Available datasets: {datasets}")
                        # Try to use the first dataset that looks like a signal
                        for key in datasets:
                            if isinstance(f[key], h5py.Dataset) and len(f[key].shape) > 0:
                                raw_signal = f[key][()]
                                if raw_signal.ndim > 1:
                                    raw_signal = raw_signal[:, 0]  # Note: transposed indexing
                                raw_signal = np.squeeze(raw_signal)
                                print(f"Using dataset '{key}' instead of 'raw_data'")
                                break
                        else:
                            raise ValueError(f"No suitable dataset found in {mat_file}")
            except Exception as h5_error:
                print(f"Error loading MAT file with h5py: {h5_error}")
                return
        else:
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
    
    # Get actual signal range
    signal_min = np.min(filtered_signal)
    signal_max = np.max(filtered_signal)
    signal_range = signal_max - signal_min
    signal_mid = (signal_min + signal_max) / 2
    
    print(f"Signal range: {signal_min:.1f} to {signal_max:.1f} µV")
    print(f"Signal midpoint: {signal_mid:.1f} µV")
    print(f"Peak-to-peak amplitude: {signal_range:.1f} µV")
    print(f"Scale mode: {'Adaptive' if adaptive_scale else 'Fixed (-5000 to +5000 µV)'}")
    
    # Define target positions
    target_top = 0.8      # Position for signal maximum (2nd gridline from top)
    target_bottom = 0.2   # Position for signal minimum (4th gridline from bottom)
    
    # Signal scaling - different approaches for adaptive vs fixed
    if adaptive_scale:
        # ADAPTIVE SCALING: Map signal to use vertical space between 2nd and 4th gridlines
        # Target: signal_max at y=0.8 (2nd gridline), signal_min at y=0.2 (4th gridline)
        if signal_range > 0:
            scale_factor = (target_top - target_bottom) / signal_range
            # Map signal so that signal_min → target_bottom, signal_max → target_top
            signal = target_bottom + (filtered_signal - signal_min) * scale_factor
        else:
            # If signal is flat, center it
            signal = np.full_like(filtered_signal, 0.5)
        
        # The baseline position (where signal_mid maps to)
        baseline_position = target_bottom + (signal_mid - signal_min) * scale_factor if signal_range > 0 else 0.5
    else:
        # FIXED SCALING: Map signal so that -5000 µV → target_bottom, +5000 µV → target_top
        # This provides consistent scaling across all files
        fixed_range = 10000  # -5000 to +5000 µV
        scale_factor = (target_top - target_bottom) / fixed_range
        # Map signal so that -5000 µV maps to target_bottom and +5000 µV maps to target_top
        signal = target_bottom + (filtered_signal + 5000) * scale_factor
        
        # The baseline position (where 0 µV maps to)
        baseline_position = target_bottom + 5000 * scale_factor
    
    # Scale bar calculation
    scale_positions, scale_values = calculate_scale_positions(
        signal_min, signal_max, target_bottom, target_top, adaptive_scale)
    
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
    
    # Grid parameters - matching createGrid.m exactly
    major = 0.2   # 200ms
    minor = 0.04  # 40ms
    second = 1.0  # 1 second
    
    # Grid colors - matching createGrid.m exactly
    major_color = np.array([173, 172, 176]) / 255
    minor_color = np.array([236, 236, 236]) / 255
    second_color = np.array([140, 140, 140]) / 255
    
    # Line colors
    line_color = np.array([232, 60, 91]) / 255  # Signal line pink
    line_width = 1.0
    
    # Scale bar settings
    scale_text_color = np.array([100, 100, 100]) / 255  # Dark gray
    scale_font_size = 8
    
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
        
        # Create figure with slightly wider layout to accommodate scale bar
        fig, ax = plt.figure(figsize=(width/100, height/100), dpi=100, 
                           facecolor='white'), plt.axes([0, 0, 0.92, 1])  # Leave space for scale bar
        ax.set_xlim(visible_start, visible_end)
        ax.set_ylim(0, 1)
        ax.set_aspect('auto')
        ax.axis('off')
        
        # Draw grid exactly like createGrid.m
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
        
        # Add scale bar with calculated values
        scale_x_pos = 0.95  # Position relative to figure width
        
        # Draw scale bar marks at calculated positions
        min_scale_value = np.min(scale_values) if len(scale_values) > 0 else 0
        
        for i, (y_pos, value) in enumerate(zip(scale_positions, scale_values)):
            # Draw tick mark
            fig.text(scale_x_pos - 0.02, y_pos, '—', fontsize=scale_font_size, 
                    color=scale_text_color, ha='center', va='center', 
                    transform=fig.transFigure)
            
            # Format label with explicit + for positive values
            if value > 0:
                label = f'+{value:.0f}'
            elif value == 0:
                label = '0'
            else:
                label = f'{value:.0f}'
            
            # Add µV unit to the most negative (lowermost) visible label
            if value == min_scale_value:
                label += ' µV'
            
            # Draw label
            fig.text(scale_x_pos, y_pos, label, fontsize=scale_font_size, 
                    color=scale_text_color, ha='left', va='center', 
                    transform=fig.transFigure)
        
        # Add a thin vertical line connecting the scale marks
        if len(scale_positions) > 1:
            scale_line_x = scale_x_pos - 0.025
            fig.add_artist(plt.Line2D([scale_line_x, scale_line_x], 
                                     [np.min(scale_positions), np.max(scale_positions)], 
                                     color=scale_text_color, linewidth=0.5, alpha=0.5,
                                     transform=fig.transFigure))
        
        # Convert to image
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

def main():
    """Command line interface"""
    import argparse
    parser = argparse.ArgumentParser(description='Create professional-style EEG video')
    parser.add_argument('mat_file', nargs='?', help='Path to .mat file containing EEG data')
    parser.add_argument('--output', '-o', default=None, help='Output video file')
    parser.add_argument('--window', '-w', type=float, default=6, help='Visible window size in seconds')
    parser.add_argument('--fps', '-f', type=int, default=30, help='Frames per second')
    parser.add_argument('--adaptive-scale', action='store_true', help='Use adaptive scale bars (skip interactive prompt)')
    parser.add_argument('--fixed-scale', action='store_true', help='Use fixed -5000 to +5000 µV scale bars (skip interactive prompt)')
    parser.add_argument('--all', '-a', action='store_true', help='Process all .mat files in current directory')
    parser.add_argument('--pattern', '-p', default="*.mat", help='File pattern to match (e.g., "data/*.mat")')
    args = parser.parse_args()
    
    # Determine scaling mode
    if args.adaptive_scale and args.fixed_scale:
        print("Error: Cannot specify both --adaptive-scale and --fixed-scale")
        return
    elif args.adaptive_scale:
        adaptive_scale = True
    elif args.fixed_scale:
        adaptive_scale = False
    else:
        # Interactive prompt for scaling mode
        print("\nScale bar options:")
        print("1. Adaptive scaling - Scale bars adapt to the actual signal range (good for detailed analysis)")
        print("2. Fixed scaling - Always show -5000 to +5000 µV range (good for comparing multiple files)")
        
        while True:
            try:
                choice = input("\nChoose scaling mode (1 for adaptive, 2 for fixed) [default: 1]: ").strip()
                if choice == '' or choice == '1':
                    adaptive_scale = True
                    print("Using adaptive scaling")
                    break
                elif choice == '2':
                    adaptive_scale = False
                    print("Using fixed scaling (-5000 to +5000 µV)")
                    break
                else:
                    print("Please enter 1 or 2 (or press Enter for default)")
            except KeyboardInterrupt:
                print("\nOperation cancelled.")
                return
            except EOFError:
                # Handle case where input is not available (e.g., piped input)
                adaptive_scale = True
                print("Using default adaptive scaling")
                break
    
    # Process all .mat files mode
    if args.all:
        mat_files = glob.glob(args.pattern)
        if not mat_files:
            print(f"No .mat files found matching pattern: {args.pattern}")
            return
        
        print(f"Found {len(mat_files)} .mat files to process")
        scale_mode = "adaptive" if adaptive_scale else "fixed"
        print(f"Using {scale_mode} scale bars for all files")
        
        for mat_file in mat_files:
            # Remove extension and append _video
            base_name = os.path.splitext(mat_file)[0]
            output_video = base_name + f"_video_{scale_mode}.mp4"
            
            print(f"\nProcessing {mat_file} -> {output_video}")
            create_eeg_video(mat_file, out_file=output_video, window_size=args.window, 
                                 fps=args.fps, adaptive_scale=adaptive_scale)
    
    # Single file mode
    elif args.mat_file:
        scale_mode = "adaptive" if adaptive_scale else "fixed"
        default_output = os.path.splitext(args.mat_file)[0] + f"_video_{scale_mode}.mp4"
        output = args.output if args.output else default_output
        create_eeg_video(args.mat_file, output, args.window, args.fps, adaptive_scale)
    
    # No input provided
    else:
        parser.print_help()

if __name__ == '__main__':
    main() 
