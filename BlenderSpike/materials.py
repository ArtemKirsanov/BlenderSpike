import bpy
import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib
import seaborn as sns
import cmasher

def srgb2lin(s):
    if s <= 0.0404482362771082:
        lin = s / 12.92
    else:
        lin = pow(((s + 0.055) / 1.055), 2.4)
    return lin

def lin2srgb(lin):
    if lin > 0.0031308:
        s = 1.055 * (pow(lin, (1.0 / 2.4))) - 0.055
    else:
        s = 12.92 * lin
    return s

def to_blender_color(rgba):
    return np.array([
        srgb2lin(rgba[0]),
        srgb2lin(rgba[1]),
        srgb2lin(rgba[2]),
        rgba[3]])

def get_cmap_by_name(cmap_name):
    return sns.color_palette(cmap_name, as_cmap=True)

def get_enum_items():
    cmap_ids = sorted(plt.colormaps())
    return [(i,i,"") for i in cmap_ids]

def create_material(name = "SectionMaterial",
                     min_voltage_value = -70,
                     max_voltage_value = 20,
                     cmap_name="plasma",
                     cmap_start=0,
                     cmap_end=1,
                     emission_strength = 2,
                     colormap_steps = 10
                    ):

    cmap = get_cmap_by_name(cmap_name)
    cmap = cmasher.get_sub_cmap(cmap,cmap_start,cmap_end)

    mat = bpy.data.materials.new(name) 
    mat.use_nodes = True
    nodes = mat.node_tree.nodes
    
    nodes.remove(nodes.get("Principled BSDF"))
    output_node = nodes.get("Material Output")
    output_node.location = (150,0)

    emission_node = nodes.new("ShaderNodeEmission")
    emission_node.location = (0, 0)
    emission_node.inputs[1].default_value = emission_strength
    
    # Add a ColorRamp and set its location:    
    ramp_node = nodes.new('ShaderNodeValToRGB')
    ramp_node.color_ramp.elements[0].color = to_blender_color(cmap(0))
    ramp_node.color_ramp.elements.remove(ramp_node.color_ramp.elements[-1])
    
    for k,p in enumerate(np.linspace(0,1,colormap_steps)[1::]):
        ramp_node.color_ramp.elements.new(p)
        ramp_node.color_ramp.elements[k+1].color = to_blender_color(cmap(p))
    ramp_node.location = (-350,0)
    
    # Voltage attribute node
    attribute_node = nodes.new("ShaderNodeAttribute")
    attribute_node.location=(-900, 150)
    attribute_node.attribute_name="Voltage"
    
    
    # Color limits nodes
    min_clim_node = nodes.new("ShaderNodeValue")
    min_clim_node.name = "Min value"
    min_clim_node.label = "Min value"
    min_clim_node.outputs[0].default_value = min_voltage_value
    min_clim_node.location = (-900,-50)
    
    max_clim_node = nodes.new("ShaderNodeValue")
    max_clim_node.name = "Max value"
    max_clim_node.label = "Max value"
    max_clim_node.outputs[0].default_value = max_voltage_value
    max_clim_node.location = (-900,-200)
    
    # Math nodes
    subtract_node_1 = nodes.new("ShaderNodeMath")
    subtract_node_1.location = (-700, 70)
    subtract_node_1.operation = "SUBTRACT"
    
    subtract_node_2 = nodes.new("ShaderNodeMath")
    subtract_node_2.location = (-700, -150)
    subtract_node_2.operation = "SUBTRACT"
    
    divide_node = nodes.new("ShaderNodeMath")
    divide_node.location = (-500,0)
    divide_node.operation="DIVIDE"
    
    # Connecting nodes
    mat.node_tree.links.new(emission_node.outputs[0], output_node.inputs[0])
    mat.node_tree.links.new(ramp_node.outputs[0], emission_node.inputs[0])
    mat.node_tree.links.new(divide_node.outputs[0], ramp_node.inputs[0])


    mat.node_tree.links.new(subtract_node_1.outputs[0], divide_node.inputs[0])
    mat.node_tree.links.new(subtract_node_2.outputs[0], divide_node.inputs[1])
    
    mat.node_tree.links.new(attribute_node.outputs[2], subtract_node_1.inputs[0])
    mat.node_tree.links.new(min_clim_node.outputs[0], subtract_node_1.inputs[1])
    
    mat.node_tree.links.new(min_clim_node.outputs[0], subtract_node_2.inputs[1])
    mat.node_tree.links.new(max_clim_node.outputs[0], subtract_node_2.inputs[0])


    return mat


class VoltageMaterialProps(bpy.types.PropertyGroup):
    min_value : bpy.props.FloatProperty(
        name="Min voltage",
        default = -70
    )
    
    max_value : bpy.props.FloatProperty(
        name="Max voltage",
        default = 20
    )

    colormap : bpy.props.StringProperty(
        name = "Colormap",
        default = "plasma"
    )

    cmap_start : bpy.props.FloatProperty(
        name = "Start color",
        default = 0,
        min = 0,
        max = 1
    )

    cmap_end : bpy.props.FloatProperty(
        name = "End color",
        default = 1,
        min = 0,
        max = 1
    )

    emission_strength : bpy.props.FloatProperty(
        name = "Emission Strength",
        default = 2,
        min = 0,
        soft_max = 50
    )

    colormap_steps : bpy.props.IntProperty(
        name = "Colormap steps",
        default = 10,
        min=2,
        max=30
    )
    
class BLENDERSPIKE_OT_MaterialCreator(bpy.types.Operator):
    '''
        Create a material to color-code for voltage in set limits
    '''

    bl_idname = 'blenderspike.create_material'
    bl_label = 'Create a voltage coloring'

    def execute(self, context):

        props = context.scene.blenderspike_materials

        

        mat = create_material(
            min_voltage_value=props.min_value,
            max_voltage_value=props.max_value,
            cmap_name=props.colormap,
            cmap_start = props.cmap_start,
            cmap_end = props.cmap_end,
            emission_strength= props.emission_strength,
            colormap_steps= props.colormap_steps
        )

        # Cleaning all children materials and assigning a new one

        for ob in context.selected_objects:
            for sec in ob.children:
                sec.data.materials.clear()
                sec.data.materials.append(mat)

        return {'FINISHED'}


class BLENDERSPIKE_OT_RemoveMatertials(bpy.types.Operator):
    bl_idname = 'blenderspike.remove_materials'
    bl_label = 'Remove materials'

    def execute(self, context):
        for sec in context.object.children:
            sec.data.materials.clear()
        return {'FINISHED'} 

class BLENDERSPIKE_OT_SetupWorld(bpy.types.Operator):
    bl_idname = 'blenderspike.setup_world'
    bl_label = 'Auto world setting'

    def execute(self, context):
        bpy.data.worlds["World"].node_tree.nodes["Background"].inputs[0].default_value = (0, 0, 0, 1)
        bpy.context.scene.eevee.use_bloom = True
        return {'FINISHED'}
