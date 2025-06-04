## 250604 MATLAB Temporal Dynamics Visualization—Modifications

### Original Function Overview:
`plotTemporalPSDmap.m` creates a 2-panel time-resolved PSD (Power Spectral Density) map with spectrum overlay to visualize the temporal distribution of detected peaks in EEG data. The function compares STIM vs SHAM conditions with:
- Stacked bar charts showing z-score distributions over time bins
- Optional KDE (Kernel Density Estimation) smoothing curves
- Standard Error of Mean (SEM) error bars across multiple files
- Support for custom Adobe .ase color palettes
- Export capabilities to Excel, PNG, and MATLAB .fig formats

### Initial Modifications:

#### 1. **Z-score colorbar major ticks at 2 sigma intervals** (Line 469)
```matlab
% MODIFICATION 1: Set colorbar ticks to every 2.0 sigma for the finite range
tick_values = z_score_min:2:z_score_max;  % Changed from 1 to 2
```
Changes colorbar ticks from every 1σ (-4, -3, -2, -1, 0, 1, 2, 3, 4) to every 2σ (-4, -2, 0, 2, 4).

#### 2. **Inverted/flipped ASE color palette** (Line 265)
```matlab
% MODIFICATION 2: Flip/inverse the ASE color palette
ase_colors = flipud(ase_colors);  % Flip the color order
```
Flips the color palette vertically for reverse order application. Console output indicates "(INVERTED)" when custom palette loads.

#### 3. **Dynamic Y-axis tick intervals** (Lines 293-298)
```matlab
% MODIFICATION 3: Determine y-axis tick interval based on maxCount
if maxCount > 20    % Changed threshold from 100 to 20
    y_tick_interval = 20;  // Use 20 when max exceeds 20
else
    y_tick_interval = 10;  // Use 10 for smaller values
end
```
Y-axis major ticks appear every 10 units for smaller datasets, or every 20 units when maxCount exceeds 20.

### Additional Modifications (250604):

#### 4. **Fixed missing right axis border** (Multiple locations)
Added explicit box drawing commands to ensure complete rectangular borders:
```matlab
ax_right.Box = 'on';  // Ensures right axis border is drawn
set(gca, 'Box', 'on');  // Force complete box
drawnow;  // Force redraw
```

#### 5. **Single black border around stacked bars** (Lines ~310, ~420)
```matlab
% MODIFICATION 4: Add single black border around entire stacked bar
bar(cent, sum(stimStackedData, 2), 'FaceColor', 'none', 'EdgeColor', 'k', 'LineWidth', 0.5);
```
Replaced individual segment borders with a single outline around the entire stacked bar for cleaner visualization.

### Usage Notes:
- The function prompts for z-score increment, SEM display, and KDE display options
- Optional Adobe .ase color palette can be loaded for custom z-score coloring
- Exports include both raw data (Excel) and visualizations (PNG/FIG)
- Handles multiple files per condition with proper averaging and SEM calculations
