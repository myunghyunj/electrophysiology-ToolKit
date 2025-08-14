import cv2
import numpy as np
import matplotlib.pyplot as plt
import os # Import os for path normalization

def measure_distance_from_video(video_path):
    """
    Opens the first frame of a video, allows the user to click two points,
    and calculates the pixel distance between them.
    """
    try:
        # 1. Open the video file
        # The os.path.normpath helps clean up the path format
        cap = cv2.VideoCapture(os.path.normpath(video_path))
        if not cap.isOpened():
            print(f"❌ Error: Could not open video file at {video_path}")
            return
            
        # 2. Read the first frame
        ret, frame = cap.read()
        if not ret:
            print("❌ Error: Could not read the first frame from the video.")
            cap.release()
            return

        # Convert color space for correct display
        frame_rgb = cv2.cvtColor(frame, cv2.COLOR_BGR2RGB)
        
        # 3. Display the frame
        fig, ax = plt.subplots(1)
        ax.imshow(frame_rgb)
        ax.set_title("Click on two points to measure. Close this window when done.", fontsize=10)
        
        print("\n✅ A window with your video frame has opened.")
        print("Please click two points on the image.")
        
        # 4. Get two points from user clicks
        points = plt.ginput(n=2, timeout=0)
        plt.close(fig)
        
        if len(points) < 2:
            print("\n⚠️ You did not select two points. Please run the script again.")
            return

        p1, p2 = points
        
        # 5. Calculate and print the distance
        pixel_distance = np.hypot(p2[0] - p1[0], p2[1] - p1[1])
        
        print("-" * 40)
        print("Measurement Complete:")
        print(f"Point 1: (x={p1[0]:.2f}, y={p1[1]:.2f})")
        print(f"Point 2: (x={p2[0]:.2f}, y={p2[1]:.2f})")
        print(f"📏 The direct pixel distance is: {pixel_distance:.2f} pixels")
        print("-" * 40)

    except Exception as e:
        print(f"❌ An unexpected error occurred: {e}")
    finally:
        if 'cap' in locals() and cap.isOpened():
            cap.release()

# --- Main part of the script ---
if __name__ == "__main__":
    video_file = input("Drag and drop your video file here and press Enter:\n> ")
    
    # --- THIS IS THE FIX ---
    # 1. Strip whitespace and quotes from the beginning/end
    cleaned_path = video_file.strip().strip("'\"")
    # 2. Remove the escaping backslashes that the terminal adds
    cleaned_path = cleaned_path.replace('\\', '')
    # ----------------------
    
    measure_distance_from_video(cleaned_path)
