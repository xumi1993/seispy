# Visually pick depth of the discontinuities after CCP stacking

Before open an user interface, the file (`*.npz`) after CCP Stacking must be calculated, such as `stack_data_3D.npz`  

`pickdepth` command will pick the depths of the Moho, 410km and 660km discontinuities from the CCP Stacking files.

    usage: pickdepth [-h] [-d dep_min/dep_max] [-i index] [-s smooth] stack_data_path

    User interface for picking PRFs

    positional arguments:
        stack_data_path     Path to CCP stacked data

    options:
      -h, --help          show this help message and exit
      -d dep_min/dep_max  Depth range contain target interface
      -i index            Specify starting index of bins
      -s smooth           Smoothing scale in km

- `-i` Starting the picking from the specified bin.
- `-d` Depth range contain target interface, Such as 25/35, 370/450, 620/700 for Moho(30km), 410km and 660km discontinuities, respectively.
- `-s` Smoothing scale in km, such as `-s 5` for Mantle Transition Zone  
  
## Open the user interface

```
pickdepth -d370/450 -s 5 stack_data_3D.npz
```

The window will open as following image. Users can press <kbd>z</kbd> and <kbd>c</kbd> hotkey or click `back` and `next` button for loading data in previous and next bins, respectively.

```{figure} ../_static/pickdepth/pickdepth1.png
:alt: Main UI without any picks
:figwidth: 70%
Main UI for picking depth of d410
```

It shows the RF in the stacking Bin of NO. 577 at latitude: 49.14 and longitude: 102.48. 

- Bottom Middle Panel: The blue line is the mean stacked RF. The intervals between two dashed lines are the corresponding 95% confidence intervals after bootstapping.  
***<u>Green Line</u>*** show the picked depth at this bin. Orange line show another possible depth.

- Upper Panel shows RFs in the stacking bins along the Longitude.  
- Bottom Left Panel show RFs in the stacking bins along the Latitude.  
- Bottom Right Panel show number of RFs in the stacking bin at each depth.  

## Operations

- Directly click the line on the interface to set the depth and it will become green. Logs Panel shows the depths you picked.

```{figure} ../_static/pickdepth/pickdepth2b.png
:alt: Main UI without any picks
:figwidth: 70%
Main UI with picked depth (green line)
```

- Input the Bin Number and click `Load` in the Bin Location Panel

- Click the right mouse botton in the middel panel to set the depth as NAN if the data is bad and you need delete it.

```{figure} ../_static/pickdepth/pickdepth4b.png
:alt: Main UI without any picks
:figwidth: 70%
Main UI without any picks
```

- `Save` or <kbd>ctrl</kbd>+<kbd>s</kbd> to save the picked depths in a text file, which include following columns:

  :::{table} Content of picked depths in the text file
  :widths: auto

  | key          | Description |
  | ---          |         --- |
  | Latitude     | Latitude of each bin |
  | Longitude    | Longitude of each bin |
  | Depth        | Depth of the target discontinuity at each bin |
  | Low CI       | Lower bound of confidence interval |
  | Upper CI     | Upper bound of confidence interval |
  | Number       | Number of RFs stacked at each bin |
  :::

