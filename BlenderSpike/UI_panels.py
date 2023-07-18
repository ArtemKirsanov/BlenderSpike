import bpy

# ----------------------- NEURON BUILDER UI -----------------------


class BLENDERSPIKE_PT_NeuronBuilder(bpy.types.Panel):
    
    bl_label =  'Neuron Builder'
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "BlenderSpike"
    
    def draw(self, context):
        layout = self.layout
        props = context.scene.blenderspike_neuronbuild
        col = layout.column()
        col.prop(props, "filepath")
        col.label(text="Coordinates", icon="GRID")
        col.prop(props, "center_at_origin")
        col.prop(props, "downscale_factor")

        col.separator()
        col.label(text="Morphology", icon="MESH_UVSPHERE")
        col.prop(props, "segmentation")
        col.prop(props, "branch_base_thickness")
        col.prop(props, "branch_thickness_homogeneity")
        col.prop(props, "simplify_soma")
        col.prop(props, "with_caps")
        
    

        row = layout.row()
        row.operator("blenderspike.build_neuron")


# ----------------------- ANIMATION MANAGER UI -----------------------
class BLENDERSPIKE_PT_AnimationManager(bpy.types.Panel):
    bl_label = 'Animation manager'
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "BlenderSpike"
    
    def draw(self, context):
        layout = self.layout
        
        row = layout.row()
        row.operator("blenderspike.reload_animations") 

        row = layout.row()
        row.operator("blenderspike.remove_handlers")



# ----------------------- SHADING UI -----------------------

class BLENDERSPIKE_PT_MaterialCreator(bpy.types.Panel):
    bl_label = 'Shading manager'
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"
    bl_category = "BlenderSpike"
    
    def draw(self, context):
        layout = self.layout
        
        props = context.scene.blenderspike_materials
        col = layout.column()
        col.prop(props, "min_value")
        col.prop(props, "max_value")
        col.prop(props, "colormap")
        col.prop(props, "cmap_start")
        col.prop(props, "cmap_end")
        col.prop(props, "colormap_steps")
        col.prop(props, "emission_strength")
        


        row = layout.row()
        row.operator("blenderspike.create_material")
        col = layout.column()
        col.operator("blenderspike.remove_materials")
        col.operator("blenderspike.setup_world")

