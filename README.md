# BlenderSpike ⚡

<p align="center">
  <img src="assets/logo-light.png" width="550"
 </p>
 

An add-on for Blender, to create 3D activity animations of [NEURON](https://neuron.yale.edu/neuron/) models. 

- Easlily bring morphologies from NEURON to Blender 💻
- Animate dynamics of membare voltage in space and time ⚡
- Color by voltage using *matplotlib* and *seaborn* stunning colormaps 🌈

<p align="center">
  <img src="assets/overview.gif" width="1000"
 </p>



## Installation


- Install Blender 3.3 or newer: https://www.blender.org/download/ 
- Download the latest `BlenderSpike.zip` from the Releases section 
- In Blender go to Edit > Preferences... > Add-ons > Install... open the downloaded ZIP file
- Enable the addon by going to Edit > Preferences... > Add-ons and enabling "Add Mesh: BlenderSpike"

If you wish to run your own simulations, you will need to install a copy of NEURON (see [instructions](https://nrn.readthedocs.io/en/8.2.2/install/install.html)). After this, install the Python companion module by running:

```python
pip install .
```

(install it to the same environment that contains NEURON)

## 📄 Basic usage

> **Note**
> For more detailed examples see `demos` directory

BlenderSpike constists of 2 parts.

**I)** **blenderspike_py** is a Python companion module which makes it easy to record voltage across all segments in your NEURON model and save the data as a `.pickle` file.

**II)** **BlenderSpike addon** adds a set of UI panels, which let you load the resulting `.pickle` file and bring your morphology into Blender along with voltage animation. 



<p align="center">
  <img src="assets/How to use.png" width="1000"
 </p>

### NEURON → pickle

BlenderSpike pulls geometry and 3D location data from NEURON objects and operates with `nmn.Section` instances (see [documentation](https://nrn.readthedocs.io/en/latest/python/modelspec/programmatic/topology/geometry.html)). This implies that you use Python to initialize and run NEURON models using the `neuron.h` module.

1) Set up your model morphology and biophysics in NEURON (**but don't run the simulation yet**)

2) Create an instance of `CellRecorder` class by passing a list of all Sections, which constitute the modelled cell

3) Run the NEURON model (using `h.continuerun` or any other way)

4) Export the model and its activity by calling `CellRecorder.save_pickle()`

### pickle → Blender

1) In the Blender 3D view window navigate to the BlenderSpike panel (press `n` to show the menu)

2) In the "Neuron Bulder" tab click on the folder icon (Path to file), select the `.pickle` you created and click "Build a neuron"

3) To add a voltage material, navigate to "Shading Manager" panel, choose a colormap and voltage color limits 

4) Select the neuron parent object and click "Create a voltage coloring"

<video src="assets/BlenderSpike walkthrough.mp4",controls autoplay loop> </video>



## 🎨 Customizing neurons

### Coordinates 🛰️

BlenderSpike uses internal coordinates of NEURON (eg. `nrn.Section.x3d`) to position neurons in the scene. These often be way off (especially for reconstructed morphologies), so two options are provided during the neuron creation in Blender.


1) To make sure that soma is located in the center of the scene `(0,0,0)` check the "Center at origin" option 
2) To adjust for different scales (for example NEURON units can be much larger / smaller than Blender's) change the Downscale factor slider.

<p align="center">
  <img src="assets/Customization – coordinates.png" width="500"
 </p>

```
> **Note**
> Coordinates from `.pickle` will be divided (not multiplied) by the downscale factor
```

### Geometry  📐

#### Segmentation

BlenderSpike allows you to control spatial resolution of both the morphology and voltage animation with a single **segmentation** parameter. More a more detailed explation of what it stands for, see [Segmentations and linear interpolation](###Segmentations and linear interpolation).

On a high level – bigger **segmentation** values make the morphology and the voltage dynamics look more detailed and visually pleasing, but require more time to render and may cause more crashes.

<p align="center">
  <img src="assets/Customization – segmentation.png" width="700"
 </p>


#### Branch thickness 

Thickness factor of branches is controlled by the **branch thickness** slider. Depending on the relationship between units used in NEURON and Blender optimal values may vary.

<p align="center">
  <img src="assets/Customization – branch thickness.png" width="700"
 </p>



#### Thickness homogeneity

For the sake of visualization, BlenderSpike allows you to control the displayed "homogeneity" of branch thickness. This does not affect the simulation result and is used only to see thin branches better.

Be careful, since increasing the homogeneity will probably require to tweak the thickness factor. Just play around with thickness homogeneity and baseline branch thickness to find the optimal balance ;)

<p align="center">
  <img src="assets/Customization – homogeneity.png" width="700"
 </p>



### Coloring 🌈

In order to color the neuron according to the values of membrane voltage you need to choose a [colormap](https://matplotlib.org/stable/tutorials/colors/colormaps.html) and set voltage limits. To set a colormap type its name in the text field and press "Create a voltage coloring" (make sure that the parent empty object is selected).

**Available colormaps:** Colormaps are parsed by name using `seaborn.color_palette()`, so please refer to the [documentation](https://seaborn.pydata.org/generated/seaborn.color_palette.html)

There is also an option to create a sub-map from a portion of an existing colormap by specifying normalized `start color` and `stop color` (from a range between 0 and 1)

<p align="center">
  <img src="assets/Customization – colormaps.jpg" width="700"
 </p>



## ❗Things to keep in mind❗

### ❗Reloading data

Note, that when you save a Blender file as `.blend`, meshes and metadata are saved, but the Python objects and `frame_change_post` handlers are discarded when the application is closed. This is why after each opening the `.blend` file you will need to manually reload the animation handlers by selecting the NEURON parent Null object and clicking on "Reload animation data" in the BlenderSpike panel.

### ❗Blender crashes

For some reason during rendering long animations (1000s of frames) and/or handling scenes with large number of segments **Blender often crashes**. At the moment I'm not sure what causes this (maybe not enough video memory or something). Just keep this in mind and remember to save your files frequently 🙂

<p align="center">
  <img src="assets/Blender crash.png" width="600"
 </p>


## ⚙️ How it works 

### Recording voltage
In Python, `CellRecorder` class initializes an array of `h.Vector` instances to record the voltage of each segment along the neuron tree over the entire simulation.

### Dumping into pickle

When you call `save_pickle()`, the `CellRecorder` iterates over all its sections and for each, constructs the following dictionary:

- **`ID`**: an integer number serving as a unique identifier used by Blender. Typically, it is equal to the index of the particular section in the list of all sections.
- **`type`**: a String specifying the section's type (Soma, Dendrite, Axon etc.). Only used in Blender for display purposes. The data for the type is parsed from the `name` attribute of the corresponding NEURON section object.
- **`X`**: an array with length `section.n3d()`, containing X-coordinates of points, specifying section's shape
- **`Y`**: an array with length `section.n3d()`, containing Y-coordinates of points, specifying section's shape
- **`Z`**: an array with length `section.n3d()`, containing Z-coordinates of points, specifying section's shape
- **`DIAM`**: an array with length `section.n3d()`, containing the diameters of anchor points, specifying section's shape
- **`Voltage`**: a dictionary with frame-wise animation data of the form `{FRAME: VOLTAGE_ARRAY}`, where `VOLTAGE_ARRAY ` is an array with `section.nseg` points, specifying the voltage profile along the segments for a given frame. Note, that the maximum number of frames is given by `FRAME_NUM` argument, when calling `.save_pickle()`

A list of such section dictionaries is what is dumped into the pickle. 

Frame-wise voltage is obtained by resampling source voltage array with `FRAME_NUM` points. 

### Constructing the mesh

- A neuron is built branch-by-branch using Blender's [Bezier splines](https://docs.blender.org/manual/en/latest/modeling/curves/structure.html#splines), where the `.co` and `.radius` attributes of `spline.bezier_points` is controlled by the morphology data, exported as `X`,`Y`,`Z` and `DIAM` arrays in `.pickle` as well as the  [**segmentation**](###Segmentations and linear interpolation).
- Soma is built as a sphere if "**Simplify soma**" is checked. If unchecked – soma has a cylidrical shape as is built similarly to the branches.
- Branches are converted from splines into `Mesh` objects
- Both soma and branches are nested inside a parent **Empty** object (an empty axes in Blender), which holds metadata (such as path to the pickle) is custom object properties.

### Updating voltage attributes

Inside Blender the spatial profile of voltage is represented as a custom **Vertex attribute** (you could see it in Blender's *Spreadsheet* tab)

<p align="center">
  <img src="assets/Blender interface - spreadsheet.png" width="400"
 </p>
Updating the voltage profile in time relies on using `frame_change_post` [handler](https://docs.blender.org/api/current/bpy.app.handlers.html).

Every time you build a neuron, internally an instance of `BlenderNeuron` class is created. It has a `voltage_handler` method that loops over all the sections and sets the voltage attribute according to a current frame (as a lookup from the .pickle file). To be automatically called every time the frame changes, the `voltage_handler` method of a given neuron should be appended to the list of Blender's `frame_change_post` handlers.

Because Blender doesn't save the handler objects in the `.blend` file, **every time you close and open Blender, you need to re-create the handlers** (see [Reloading data](###reloading-data))

### Segmentations and linear interpolation

Under the hood, NEURON typically ignores how the morphology is laid out in 3D space (as specified by `.x3d()`, `.y3d()` etc). It only uses the section lengths, radii and the topology (which sections are connected to which) to run the simulations. This results in **2 types of discrectization**: 

1) Spatial discrectization of the **actual morphology** (points to specify 3D positions and angles of sections)  $N=\text{n3d()}$ 

2) Spatial discrectization for **numerical simulations** (each `Section` is broken up into `.nseg` segments that are treated as cylinders to run biophysical calculations) – $N=\text{nseg}$

To bring the model into Blender interpolation is performed to contstuct a one-to-one mapping between the two discrectization for each Section.

- Array of coordinates (containing $N=\text{n3d()}$ elements) is linearly interpolated to $N_\text{segmentation}$ elements
- Voltage array for each section ($N=\text{nseg}$) is interpolated to $N_\text{segmentation}$ elements

Resulting interpolated coordinates are using during branches contruction, while interpolated voltage arrays are used to set point attributes of branch vertices on every frame in a consistent manner. 

Thus, the `segmentation` parameter **controls both the morphological detailization and the detailization of voltage profile** simultaneously.



<p align="center">
  <img src="assets/Segmentation diagram.png" width="800"
 </p>



> **Warning**
> 
> Be aware that increasing the `segmentation` parameter in Blender **will not** make the result more detailed than what is stored in `.pickle`. So make sure you simulate the model with enough `nseg` and with the detailed enough `.swc` morphology to begin with.



## Licence

MIT

