import bpy
import importlib
import subprocess
import sys

bl_info = {
    "name" : "BlenderSpike",
    "author" : "Artem Kirsanov",
    "description" : "Bring NEURON animations to Blender",
    "blender" : (3, 3, 0),
    "version" : (1, 0, 0),
    "location" : "View3D > BlenderSpike",
    "warning" : "",
    "category" : "3D View"
}

def check_and_install_modules():
    '''
        Automatically install required Python modules
    '''
    required_modules_import_names = ["matplotlib","seaborn", "cmasher", "numpy", "scipy"]  # Required Python modules
    required_modules_install_names = ["matplotlib","seaborn", "cmasher", "numpy", "scipy"]


    missing_modules = []
    for k,module_name in enumerate(required_modules_import_names):
        try:
            importlib.import_module(module_name)
            print(f"Module {module_name} is already installed. Importing...")
        except ImportError:
            missing_modules.append(required_modules_install_names[k])
            
    if missing_modules:
        print("Found missing modules: ", missing_modules)
        for module in missing_modules: 
            try:
                subprocess.check_call([sys.executable, "-m", "pip", "install", module])
                print(f"{module} installed successfully.")
            except subprocess.CalledProcessError:
                print(f"Failed to install {module}.")


check_and_install_modules() # This is called before any imports from the submodules

from .neuron_builder import NeuronBuilderProps, BLENDERSPIKE_OT_NeuronBuilder
from .animation_manager import BLENDERSPIKE_OT_HandlerRemover,BLENDERSPIKE_OT_AnimationLoader
from .materials import VoltageMaterialProps, BLENDERSPIKE_OT_MaterialCreator, BLENDERSPIKE_OT_RemoveMatertials,BLENDERSPIKE_OT_SetupWorld
from .UI_panels import BLENDERSPIKE_PT_NeuronBuilder,BLENDERSPIKE_PT_MaterialCreator, BLENDERSPIKE_PT_AnimationManager

ordered_classes = [
    # Property Groups
    NeuronBuilderProps,
    VoltageMaterialProps,

    # Operators
    BLENDERSPIKE_OT_NeuronBuilder,
    BLENDERSPIKE_OT_HandlerRemover,
    BLENDERSPIKE_OT_MaterialCreator,
    BLENDERSPIKE_OT_RemoveMatertials,
    BLENDERSPIKE_OT_AnimationLoader,
    BLENDERSPIKE_OT_SetupWorld,

    # UI Panels
    BLENDERSPIKE_PT_NeuronBuilder,
    BLENDERSPIKE_PT_AnimationManager,
    BLENDERSPIKE_PT_MaterialCreator
]

def register():
    for cl in ordered_classes:
        bpy.utils.register_class(cl)

    bpy.types.Scene.blenderspike_neuronbuild = bpy.props.PointerProperty(type = NeuronBuilderProps)
    bpy.types.Scene.blenderspike_materials = bpy.props.PointerProperty(type =VoltageMaterialProps )

def unregister():
    for cl in reversed(ordered_classes):
        bpy.utils.unregister_class(cl)
    del bpy.types.Scene.blenderspike_neuronbuild
    del bpy.types.Scene.blenderspike_materials

if __name__ == "__main__":
    register()
    
