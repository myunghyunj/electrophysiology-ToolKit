import time
from pynput import keyboard

filename = input("Enter Stopwatch Title or sample name for .txt file: ") + '.txt'

laps = []
running = False
paused = False
start_time = None
lap_counter = 1
elapsed_time = 0

def format_time(elapsed):
    hours, remainder = divmod(elapsed, 3600)
    mins, secs = divmod(remainder, 60)
    centisec = int((secs - int(secs)) * 100)
    return f"{int(hours):02}:{int(mins):02}:{int(secs):02}:{centisec:02}"

def on_press(key):
    global running, paused, start_time, lap_counter, elapsed_time

    if key == keyboard.Key.space:
        if not running and not paused:
            running = True
            start_time = time.time()
            print("Stopwatch started.")
        elif running:
            current_elapsed = elapsed_time + (time.time() - start_time)
            timestamp = format_time(current_elapsed)
            laps.append((lap_counter, timestamp))
            print(f"Lap #{lap_counter}: {timestamp}")
            lap_counter += 1

    elif key == keyboard.Key.enter:
        if running:
            paused = True
            running = False
            elapsed_time += time.time() - start_time
            print(f"Paused at {format_time(elapsed_time)}. Stop stopwatch? (y/n)")

def on_release(key):
    global running, paused, start_time

    if paused:
        if hasattr(key, 'char') and key.char:
            if key.char.lower() == 'y':
                with open(filename, 'w') as f:
                    for lap, stamp in laps:
                        f.write(f"Lap #{lap}: {stamp}\n")
                print("Timestamps saved. Exiting...")
                return False  # Stops the listener
            elif key.char.lower() == 'n':
                paused = False
                running = True
                start_time = time.time()
                print("Resumed.")

print("press SPACE to start and record laps, ENTER to pause, and ctrl+C to exit")
with keyboard.Listener(on_press=on_press, on_release=on_release) as listener:
    try:
        listener.join()
    except KeyboardInterrupt:
        print("\n Exiting.")
