import bpy
import numpy as np
import pickle
from .utils import linear_interpolation


class BlenderSection():
    '''
        A Class for representing one NEURON Section as a single branch object in Blender.
    '''
    def cast_segment_data_to_verts(self, data):
        ''' Function to repeat segment data (Nseg elemements) in a point-by-point basis for mesh '''
        if self.mesh_Npoints is None:
             self.calculate_mesh_points()
        if self.with_caps:
            Nseg = self.Nseg+2
            padded_data = np.pad(data, ((1,1)), mode="edge")
            return np.array([[i]*(self.mesh_Npoints//Nseg) for i in padded_data]).reshape(-1)
        return np.array([[i]*(self.mesh_Npoints//self.Nseg) for i in data]).reshape(-1)

               
    def __init__(self, X, Y, Z, DIAM,
                branch_ID, type,
                parent_ob=None, 
                with_caps=False, 
                simplify_soma=True,
                ):
        '''
            X,Y,Z - segment-wise coordinates of the branch
            DIAM - 
            parent_ob (bpy.types.Object) - a Blender parent object (usually a an EMPTY axis, which is created automatically with BlenderNeuron)
            with_caps (Bool) - whether to fill the caps on sections. May look a bit nicer, but the processing is slower
            simplify_soma (Bool) - whether to represent a soma as a sphere with homogeneous voltage. 
                    If true, all points have voltage as the mean across all soma segments in a given frame.
                    If false, treats soma as any other section
        '''
        self.X = X
        self.Y = Y
        self.Z = Z
        self.DIAM = DIAM
        self.ID = branch_ID
        self.type = type
        self.parent_ob = parent_ob
        self.mesh_Npoints = None            # Number of vertices in the actual Blender mesh. Set automatically while building
        self.with_caps = with_caps
        self.simplify_soma = simplify_soma
        self.Nseg = len(X)
        
        self.assign_UV = False # Whether to assign U value as a Vertex attribute
        self.attr_name = "Voltage" # Name of the custom attribute 
        
    def build_soma(self):
        ''' Build a simplified soma as a sphere'''
        center_seg_id = self.Nseg//2
        soma_coords = np.array([self.X[center_seg_id], self.Y[center_seg_id], self.Z[center_seg_id]])
        soma_radius = (self.DIAM[center_seg_id]/2)*1.5 # Soma is rendered a bit bigger for vizualization purposes
        bpy.ops.mesh.primitive_uv_sphere_add(radius = soma_radius, location = soma_coords)
        self.ob = bpy.context.active_object
        self.ob.parent = self.parent_ob
        self.ob.name = "soma"
        return self.ob

    def build_curves(self, resolution_u, bevel_depth):
        '''
            Builds a Section in Blender with Bezier curves 

            resolution_u - resolution of the crosssection circle
            bevel_depth - thickness of the curve when converted to mesh
        '''
        object_name = "{}_{}".format(self.type, self.ID)
        # Creating BLENDER curves
        tracer = bpy.data.curves.new(self.type,'CURVE')
        tracer.dimensions = '3D'
        spline = tracer.splines.new('BEZIER')

        branch_ob = bpy.data.objects.new(object_name,tracer)
        bpy.context.collection.objects.link(branch_ob) # Adding to collection
        
        # Adding thickness to convert to mesh
        tracer.resolution_u = resolution_u
        tracer.bevel_resolution = 5
        tracer.fill_mode = 'FULL'
        tracer.bevel_depth = bevel_depth
        tracer.use_fill_caps = self.with_caps
      
        spline.bezier_points.add(self.Nseg-1)
        for num in range(self.Nseg):
            spline.bezier_points[num].co = [self.X[num], self.Y[num], self.Z[num]]
            spline.bezier_points[num].radius = self.DIAM[num]/2

            spline.bezier_points[num].handle_right_type='VECTOR'
            spline.bezier_points[num].handle_left_type='VECTOR'
        
        if self.parent_ob is not None:
            branch_ob.parent = self.parent_ob
        
        self.ob = branch_ob
        return self.ob

    def build(self, resolution_u=20, bevel_depth = 2):
        '''Build the section'''
        if self.type=="soma" and self.simplify_soma:
            self.build_soma()
        else:
            self.build_curves(resolution_u,bevel_depth)


    def convert_to_mesh(self):
        ''' Convert Bezier curves to mesh'''
        if self.type=="soma" and self.simplify_soma: # No need to a simplified soma, which is already a mesh
            self.mesh_Npoints = len(self.ob.data.vertices)
            return 

        # If not a soma
        bpy.context.view_layer.objects.active = self.ob
        self.ob.select_set(True)
        bpy.ops.object.convert(target="MESH")
    
        if self.assign_UV: # Assign UV values if necessary
            UVvalues =  self.cast_segment_data_to_verts(np.linspace(0,1,self.Nseg))
            self.ob.data.attributes.new(name="Uvalue",  type="FLOAT", domain="POINT")
            self.ob.data.attributes['Uvalue'].data.foreach_set("value",UVvalues)

    def create_voltage_attribute(self):
        '''Create a custom Vertex mesh attribute'''
        self.voltage_attr = self.ob.data.attributes.new(name=self.attr_name,  type="FLOAT", domain="POINT")

    def calculate_mesh_points(self):
        self.mesh_Npoints = len(self.ob.data.vertices)

    def set_voltage_data(self,data):
        voltage_attr = self.ob.data.attributes[self.attr_name] # Getting Vertex attribute

        if self.mesh_Npoints is None:
            self.calculate_mesh_points()
            
        # --- If soma - is siplified, its voltage is set to the mean across the section (homogeneous voltage)
        if self.type=="soma" and self.simplify_soma:
            soma_voltage = np.mean(data)
            for point in range(self.mesh_Npoints):
                voltage_attr.data[point].value = soma_voltage
            return 

        # --- If not a soma - cast segment data to mesh points
        mesh_animation_data = self.cast_segment_data_to_verts(data)
        for point in range(self.mesh_Npoints):
            voltage_attr.data[point].value = mesh_animation_data[point]

    def set_metadata_custom_properties(self):
        '''Sets section ID as a custom property of the object to be saved in .blend file'''
        self.ob["ID"] = self.ID


## ------------------------------ Blender Neuron container -----------------------------------

class BlenderNeuron():
    '''
        Container class for storing a generated a neuron object and metadata
    '''

    def load_sections_dicts(self):
        ''' Load the dictionary of Sections data (exported from neuron) into the sections_dicts attribute'''
        with open(self.filepath, "rb") as f:
            self.sections_dicts = pickle.load(f)

    def __init__(self, 
                filepath,
                name="NEURON",
                center_at_origin = False,
                with_caps=False, 
                simplify_soma=True,
                segmentation=5,
                parent_ob = None,
                DOWNSCALE_FACTOR=25,
                branch_base_thickness=2,
                branch_thickness_homogeneity=0
                ):
        self.filepath = filepath
        self.name = name
        self.with_caps = with_caps
        self.center_at_origin = center_at_origin
        self.simplify_soma = simplify_soma
        self.segmentation = segmentation
        self.DOWNSCALE_FACTOR = DOWNSCALE_FACTOR
        self.branch_base_thickness = branch_base_thickness
        self.branch_thickness_homogeneity = branch_thickness_homogeneity


        self.ALL_SECTIONS = []


        self.load_sections_dicts() # Loading sections dictionary

        if parent_ob is None:
            self.create_parent_empty()
            self.set_parent_metadata()
        else:
            self.parent_ob = parent_ob
            self.name = parent_ob.name
        
        self.calculate_mean_branch_thickness()
        self.calculate_center_of_mass()

    def calculate_mean_branch_thickness(self):
        self.mean_branch_thickness = np.mean(np.array([np.mean(self.sections_dicts[k]["DIAM"]) for k in range(len(self.sections_dicts))]))

    def calculate_center_of_mass(self):
        self.center_of_mass = [0,0,0]
        # if self.center_at_origin:
        #     self.center_of_mass = [np.mean([np.mean(self.sections_dicts[k][coord]) for k in range(len(self.sections_dicts))]) for coord in ["X","Y","Z"]]
        # print("Center of mass", self.center_of_mass)


    def set_parent_metadata(self):
        ''' Store metadata in a custom properties of the parent EMPTY object '''
        attrs_to_save = ["filepath", "center_at_origin", "with_caps", "simplify_soma", "segmentation", "DOWNSCALE_FACTOR", "branch_base_thickness", "branch_thickness_homogeneity"]
        for attr in attrs_to_save:
            self.parent_ob[attr] = getattr(self, attr)

    def create_parent_empty(self):
        ''' Create a parent EMPTY Blender object, which holds metadata'''

        print("Creating parent object")
        bpy.ops.object.empty_add(type='ARROWS',location=(self.sections_dicts[0]["X"][0],self.sections_dicts[0]["Y"][0],self.sections_dicts[0]["Z"][0]), rotation=(0, 0, 0))
        self.parent_ob = bpy.context.selected_objects[0]
        self.parent_ob.name = self.name

    def get_voltage_data(self,branch_ID, frame):
        return self.sections_dicts[branch_ID]["Voltage"][frame]

    def get_branch_coordinates(self,branch_ID):
        branch_dict = self.sections_dicts[branch_ID]
        # Blender coordinates for Bezier points are constructed by interpolating the sourse NEURON array of coordinates with specified resolution (segmentation)

        X = linear_interpolation((branch_dict["X"] - self.center_of_mass[0])/ self.DOWNSCALE_FACTOR, self.segmentation) # X coordinates of segments
        Y = linear_interpolation((branch_dict["Y"] - self.center_of_mass[1]) / self.DOWNSCALE_FACTOR, self.segmentation)  # Y coordinates of segments
        Z = linear_interpolation((branch_dict["Z"] - self.center_of_mass[2])/ self.DOWNSCALE_FACTOR , self.segmentation) # Z coordinates of segments
        return [X,Y,Z]

    def get_branch_diam(self, branch_ID):
        branch_dict = self.sections_dicts[branch_ID]
        raw_diam = branch_dict["DIAM"] / self.DOWNSCALE_FACTOR
        scaled_diam = self.branch_thickness_homogeneity*self.mean_branch_thickness + (1-self.branch_thickness_homogeneity)*raw_diam
        return linear_interpolation(scaled_diam, self.segmentation) # segment diameters

    def get_branch_type(self, branch_ID):
        return self.sections_dicts[branch_ID]["type"]

    def build_branches(self):
        for i in range(len(self.sections_dicts)):
            X,Y,Z = self.get_branch_coordinates(i)
            DIAM = self.get_branch_diam(i)
            section = BlenderSection(
                                    X,Y,Z,DIAM,
                                    branch_ID=i,
                                    type=self.get_branch_type(i),
                                    parent_ob=self.parent_ob, 
                                    with_caps=self.with_caps,
                                    simplify_soma=self.simplify_soma)

            section.build(bevel_depth=self.branch_base_thickness)
            section.convert_to_mesh()
            section.create_voltage_attribute()
            section.set_metadata_custom_properties()
            self.ALL_SECTIONS.append(section)

    def voltage_handler(self,scene,*args):
        frame = scene.frame_current
        for k,sec in enumerate(self.ALL_SECTIONS):
            voltage_data = linear_interpolation(self.get_voltage_data(k, frame), self.segmentation) # Interpolating from source voltage data depending on segmentation
            try:
                sec.set_voltage_data(voltage_data)
            except:
                continue # In case the section object was deleted

    def add_voltage_handler(self):
        bpy.app.handlers.frame_change_post.append(self.voltage_handler)
        bpy.context.scene.render.use_lock_interface = True # This is to ensure render doesn't crash

    def remove_voltage_handlder(self):
        raise NotImplementedError
    
    def reinstantiate_sections_from_childen(self):
        self.ALL_SECTIONS = [0]*len(self.sections_dicts)

        for child_ob in self.parent_ob.children:
            section_ID = child_ob["ID"]
            X,Y,Z = self.get_branch_coordinates(section_ID)
            DIAM = self.get_branch_diam(section_ID)
            section = BlenderSection(
                                    X,Y,Z,DIAM,
                                    branch_ID=section_ID,
                                    type=self.get_branch_type(section_ID),
                                    parent_ob=self.parent_ob, 
                                    with_caps=self.with_caps,
                                    simplify_soma=self.simplify_soma)

            section.ob = child_ob
            self.ALL_SECTIONS[section_ID] = section



## ------------------------------ OPERATORS -----------------------------------

class NeuronBuilderProps(bpy.types.PropertyGroup):
    '''
        Property group for holding neuron builder parameters
    '''

    center_at_origin : bpy.props.BoolProperty(
        name = "Center at origin",
        default = True
    )

    segmentation : bpy.props.IntProperty(
        name="Segmentation",
        min=3,
        soft_max=101,
        default = 5
    )
    
    simplify_soma : bpy.props.BoolProperty(
        name = "Simplify soma",
        default = True
    )
    
    with_caps: bpy.props.BoolProperty(
        name = "Fill caps",
        default = False
    )

    downscale_factor : bpy.props.FloatProperty(
        name = "Downscaling factor",
        default = 25
    )

    branch_base_thickness : bpy.props.FloatProperty(
        name = "Branch thickness",
        default = 2,
        min=0,
        soft_min=1,
        soft_max=10
    )

    branch_thickness_homogeneity : bpy.props.FloatProperty(
        name = "Thickness homogeneity",
        default = 0,
        min=0,
        max=1
    )

    filepath: bpy.props.StringProperty(
        name="Path to .pickle",
        subtype="FILE_PATH"
    )
    
class BLENDERSPIKE_OT_NeuronBuilder(bpy.types.Operator):
    '''
       Operator to load the NEURON dictionary and create the mesh
    '''
    
    bl_idname = 'blenderspike.build_neuron'
    bl_label =  'Build a neuron'
    
    def execute(self, context):

        props = context.scene.blenderspike_neuronbuild


        neuron = BlenderNeuron(filepath=props.filepath, 
                                with_caps=props.with_caps,
                                center_at_origin = props.center_at_origin,
                                simplify_soma=props.simplify_soma,
                                segmentation = props.segmentation,
                                DOWNSCALE_FACTOR=props.downscale_factor,
                                branch_base_thickness=props.branch_base_thickness,
                                branch_thickness_homogeneity=props.branch_thickness_homogeneity
                                )
        neuron.build_branches()
        neuron.add_voltage_handler()

        if props.center_at_origin:
            neuron.parent_ob.location[0]=0
            neuron.parent_ob.location[1]=0
            neuron.parent_ob.location[2]=0
        

        print("Built a neuron from {}".format(props.filepath))
        return {"FINISHED"}

