import bpy
import numpy as np
import pickle
from .utils import linear_interpolation

from .neuron_builder import BlenderNeuron, BlenderSection

class BLENDERSPIKE_OT_HandlerRemover(bpy.types.Operator):   
    '''
    Operator to clear frame_change_post handers.

    Note: this should remove a handler linked to a selected object, but I don't know how to implement this. 
    So, as a temporary band aid, it just removes all handlers. 
    '''
    bl_idname = 'blenderspike.remove_handlers'
    bl_label =  'Remove all voltage handlers'

    def execute(self, context):
        bpy.app.handlers.frame_change_post.clear()
        return {"FINISHED"}
    
class BLENDERSPIKE_OT_AnimationLoader(bpy.types.Operator):   
    '''Reload animation data for selected neurons''' 

    bl_idname = 'blenderspike.reload_animations'
    bl_label =  'Reload animation data'

    def execute(self, context):
        '''
            The operator reinstantiates BlenderNeuron and nested BlenderSection objects from the metadata of selected object
        '''
        for ob in context.selected_objects:
            neuron = BlenderNeuron(
                    filepath=ob["filepath"],
                    with_caps=bool(ob["with_caps"]),
                    simplify_soma=bool(ob["simplify_soma"]),
                    segmentation=ob["segmentation"],
                    parent_ob=ob,
                    DOWNSCALE_FACTOR=ob["DOWNSCALE_FACTOR"],
                    branch_base_thickness=ob["branch_base_thickness"],
                    branch_thickness_homogeneity = ob["branch_thickness_homogeneity"]
                )
                
            neuron.reinstantiate_sections_from_childen()
            neuron.add_voltage_handler()
        return {"FINISHED"}