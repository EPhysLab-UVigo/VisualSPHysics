#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# VisualSPHysics
# Copyright (C) 2020 Orlando Garcia-Feal orlando@uvigo.es

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

bl_info = {
    "name": "VisualSPHysics Blender Addon",
    "category": "Object",
    "blender": (2, 80, 0),
}

import bpy, sys, math, re, os, time, array, ctypes, bmesh, mathutils
import xml.etree.ElementTree as ET
import traceback
 
import vtkimporter 
import diffuseparticles

#
#   Auxiliary functions
# 

# @brief Deletes a mesh.
# This function deletes the mesh pointed by parameter, freeing the used memory.
# @param meshName Name of the mesh that has to be deleted.
#
def deleteMesh (meshName):
    mesh = bpy.data.meshes.get(meshName)
    try:
        mesh.user_clear()
        bpy.data.meshes.remove(mesh)
    except:
        print("Cannot remove mesh from memory")

def magnitude(v):
    return math.sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2])

def normalize(v):
    m = magnitude(v)
    if m == 0 :
        return v    
    return [v[0]/m, v[1]/m, v[2]/m]

## @brief Generates a rope type object.
#  This function generates an object of rope type, useful for mooring simulations
#  where the rope is given by a list of points corresponding to the line traced by the rope.
#  A geometry is created for the rope to give volume to the rope.
#  @param vertices List of points that defines the rope.
#  @return Tuple of points and triangles.
#
def makeRope (vertices, lines):
    rsize=.250
    pts = []
    tris = []
    start = 0
    
    for line in lines:
        first = True
        prev = [0,0,0]
        for lp in line:
            v = vertices[lp]

            if prev == v :
                prev = [0,0,0]            

            vec = mathutils.Vector((v[0] - prev[0], v[1] - prev[1], v[2] - prev[2]))
            mrot = mathutils.Matrix.Rotation(math.radians(90.0), 4, 'X')
            
            nvec = vec * mrot
            nvec2 = nvec * mrot                                
            
            vec = normalize(vec)
            nvec = normalize(nvec)
            nvec2 = normalize(nvec2)                

            pt0 = [v[0], v[1], v[2]]                
            pt1 = [v[0]+nvec[0]*rsize, v[1]+nvec[1]*rsize, v[2]+nvec[2]*rsize]                
            pt2 = [v[0]+nvec2[0]*rsize, v[1]+nvec2[1]*rsize, v[2]+nvec2[2]*rsize]                

            pts.append(pt0)
            pts.append(pt1)
            pts.append(pt2)
                
            prev = v
            
        for i in range(start, len(pts)-3):
            m = i % 3
            r = i - m
            tris.append([(i), (r+((m+1)%3)) , ((r+((m+1)%3)) + 3) , (i+3)])
            
        start = len(pts)                   
       
    return (pts,tris)

## @brief Creates a 3D object.
# This function generates a new 3D object.
def createObject (objName, fileName, pathName, baseName, extension, objType, smooth, validate, uv, blur, startFrame, endFrame):
    filePath = bpy.path.abspath(os.path.join(pathName, fileName))
    print("Loading "+fileName+" ...")

    pts = []
    tris = []
    vels = []

    if objType == "ROPE":
        t_pts, t_lines = vtkimporter.loadrope(filePath)
        pts, tris = makeRope(t_pts, t_lines)
    elif objType == "FOAM":
        pts, tris, vels = vtkimporter.loaddiffuse(filePath)
    else:
        pts, tris = vtkimporter.load(filePath)

    newMesh = bpy.data.meshes.new(objName+'Mesh') # New mesh
    newMesh.from_pydata(pts,[],tris)    # edges or faces should be [], or you ask for problems
    newMesh.update(calc_edges=False)    # Update mesh with new data

    obj = bpy.data.objects.new(objName, newMesh)       # Create an object with that mesh
    coll = bpy.context.view_layer.active_layer_collection.collection
    coll.objects.link(obj) # Link object to scene
    
    # Smooth shading
    if smooth :
        obj = bpy.data.objects[objName]
        for poly in obj.data.polygons:
            poly.use_smooth = smooth
            
    # Set properties
    obj["DsphPathName"] = pathName
    obj["DsphBaseName"] = baseName
    obj["DsphExtension"] = extension
    obj["DsphObjType"] = objType
    obj["DsphSmooth"] = smooth
    obj["DsphValidate"] = validate
    obj["DsphUV"] = uv
    obj["DsphBlur"] = blur
    obj["DsphStartFrame"] = startFrame
    obj["DsphEndFrame"] = endFrame


## Parse DualSPHysics XML configuration
def parseXML(context):
    fpath = bpy.path.abspath(context.scene.DsphFoamXML)
    
    if os.path.isfile(fpath) :
        with open(fpath) as f:
          xmlString = f.readlines()
          # Remove the first comment from xml files generated with DesignSPHysics
          if (xmlString[0].count("<!--") > 0) :
            xmlString = xmlString[1:]
          xmlString = '\n'.join(xmlString)
          
        try:
            #tree = ET.parse(fpath)
            #root = tree.getroot()
            root = ET.fromstring(xmlString)            
            
            bpy.context.scene.DsphFoamH = float(root.find("./execution/constants/h").get("value"))
            context.scene.DsphFoamMass = float(root.find("./execution/constants/massfluid").get("value"))
            context.scene.DsphFoamTimeStep = float(root.find("./execution/parameters/parameter[@key='TimeOut']").get("value"))
                    
            if not context.scene.DsphFoamCustomDomain:
                #Get domain from xml file
                pmin = root.find("./casedef/geometry/definition/pointmin")
                context.scene.DsphFoamMinX = float(pmin.get("x"))
                context.scene.DsphFoamMinY = float(pmin.get("y"))
                context.scene.DsphFoamMinZ = float(pmin.get("z"))
                pmax = root.find("./casedef/geometry/definition/pointmax")
                context.scene.DsphFoamMaxX = float(pmax.get("x"))
                context.scene.DsphFoamMaxY = float(pmax.get("y"))
                context.scene.DsphFoamMaxZ = float(pmax.get("z"))
            return True
        
        except Exception:
            print(traceback.format_exc())
            return False                  
    else:
        return False

    
#
#   User interface
#
 
from bpy.props import *

def showPopup(text = "", title = "Info", icon = 'INFO'):
    def draw(self, context):
        self.layout.label(text=text)
    bpy.context.window_manager.popup_menu(draw, title = title, icon = icon)


def updateObjectProperty(self, context):    
    if self['DsphObjType'] == 0 :
        context.active_object['DsphObjType'] = 'FLUID'
    elif self['DsphObjType'] == 1 :
        context.active_object['DsphObjType'] = 'FOAM'
    elif self['DsphObjType'] == 2 :
        context.active_object['DsphObjType'] = 'ROPE'
    elif self['DsphObjType'] == 3 :
        context.active_object['DsphObjType'] = 'OTHER'
    
    #context.active_object['DsphSmooth'] = self['DsphSmooth']
    context.active_object['DsphValidate'] = self['DsphValidate']
    context.active_object['DsphUV'] = self['DsphUV']
    context.active_object['DsphBlur'] = self['DsphBlur']
    context.active_object['DsphStartFrame'] = self['DsphStartFrame']
    context.active_object['DsphEndFrame'] = self['DsphEndFrame']
    

class VIEW3D_PT_DsphObjectPanel(bpy.types.Panel):
    bl_label = "DualSPHysics Object Properties"
    bl_space_type = "VIEW_3D"
    bl_region_type = "UI"

    # Getters and setters for properties. With this we can change
    # the object id props and mantain checkboxes without using rna props
    # for every object.

    def getType(self):
        #return bpy.context.active_object['DsphObjType']
        val = bpy.context.active_object['DsphObjType']
        if val == 'FLUID' :
            return 0
        elif val == 'FOAM' :
            return 1
        elif val == 'ROPE' :
            return 2
        else :
            return 3
        
    def setType(self,value):
        if value == 0 :
            bpy.context.active_object['DsphObjType'] = 'FLUID'
        elif value == 1 :
            bpy.context.active_object['DsphObjType'] = 'FOAM'
        elif value == 2 :
            bpy.context.active_object['DsphObjType'] = 'ROPE'
        elif value == 3 :
            bpy.context.active_object['DsphObjType'] = 'OTHER'        

    
    def getSmooth(self):
        return bpy.context.active_object['DsphSmooth'] == True
        
    def setSmooth(self,value):
        bpy.context.active_object['DsphSmooth'] = value


    def getValidate(self):
        return bpy.context.active_object['DsphValidate'] == True
        
    def setValidate(self,value):
        bpy.context.active_object['DsphValidate'] = value


    def getUV(self):
        return bpy.context.active_object['DsphUV'] == True
        
    def setUV(self,value):
        bpy.context.active_object['DsphUV'] = value


    def getBlur(self):
        return bpy.context.active_object['DsphBlur'] == True
        
    def setBlur(self,value):
        bpy.context.active_object['DsphBlur'] = value


    # Properties to be shown in the panel
        
    bpy.types.Scene.DsphObjType = bpy.props.EnumProperty(
            items = [('FLUID', 'Fluid Object', "", "MOD_WAVE", 0),
                    ('FOAM', 'Foam Object', "", "MOD_PARTICLES", 1), 
                    ('ROPE', 'Rope Object', "", "MOD_CURVE", 2), 
                    ('OTHER', 'Other', "", "MESH_CUBE", 3)],
            name = "Object Type",
            description = "Select the type of object from a DualSPHysics simulation.",
            get = getType,
            set = setType)

    bpy.types.Scene.DsphSmooth = BoolProperty(
        name = "Smooth Shading", 
        description = "Can lead to artifacts on glossy objects.",
        get = getSmooth,
        set = setSmooth)

    bpy.types.Scene.DsphValidate = BoolProperty(
        name = "Validate mesh", 
        description = "Check it in order to avoid artifacts on fluid objects.",
        get = getValidate,
        set = setValidate)

    bpy.types.Scene.DsphUV = BoolProperty(
        name = "Transfer UV Maps", 
        description = "Needed for textured objects.",
        get = getUV,
        set = setUV)
        
    bpy.types.Scene.DsphBlur = BoolProperty(
        name = "Motion Blur", 
        description = "Enable motion blur for this object.",
        get = getBlur,
        set = setBlur)        
        
        
    @classmethod
    def poll(self, context):
        return context.object and 'DsphObjType' in context.object
 
    def draw(self, context):
        layout = self.layout
               
        layout.prop(context.scene, "DsphObjType")
        layout.prop(context.scene, "DsphSmooth")
        
        if context.object["DsphStartFrame"] != context.object["DsphEndFrame"] :
            layout.prop(context.scene, "DsphValidate")
            layout.prop(context.scene, "DsphUV")                
            layout.prop(context.scene, "DsphBlur")            
            layout.prop(context.object, '["DsphStartFrame"]')
            layout.prop(context.object, '["DsphEndFrame"]')

class FILE_BROWSER_OT_BrowseFileOperator(bpy.types.Operator):
    bl_idname = "dialog.browsefile"
    bl_label = "Select DualSPHysics VTK file"
    
    bl_space_type = "FILE_BROWSER"
    bl_region_type = "CHANNELS"
    
    filter_glob: bpy.props.StringProperty(default="*.vtk")

    filename: bpy.props.StringProperty()

    directory: bpy.props.StringProperty(subtype='DIR_PATH')
    
    DsphObjEnum: bpy.props.EnumProperty(
            items = [('FLUID', 'Fluid Object', "", "MOD_WAVE", 0),
                    ('FOAM', 'Foam Object', "", "MOD_PARTICLES", 1), 
                    ('ROPE', 'Rope Object', "", "MOD_CURVE", 2), 
                    ('OTHER', 'Other', "", "MESH_CUBE", 3)],
            name = "Object Type",
            description = "Select the type of object from a DualSPHysics simulation.")
    
    DsphSmooth: bpy.props.BoolProperty(
        name = "Smooth Shading", 
        description = "Can lead to artifacts on glossy objects.")

    DsphValidate: bpy.props.BoolProperty(
        name = "Validate mesh", 
        description = "Check it in order to avoid artifacts on fluid objects.")

    DsphUV: bpy.props.BoolProperty(
        name = "Transfer UV Maps", 
        description = "Needed for textured objects.")
        
    DsphBlur: bpy.props.BoolProperty(
        name = "Motion Blur", 
        description = "Enable motion blur for this object.") 
    
    def draw(self, context):
        layout = self.layout

        layout.prop(self, "DsphObjEnum")
        layout.prop(self, "DsphSmooth")
        layout.prop(self, "DsphValidate")
        layout.prop(self, "DsphUV")
        layout.prop(self, "DsphBlur")
        
    def execute(self, context):
        print(self.filename)        
        print(self.directory)
        #Let's detect sequence numbers
        p = re.compile('(.*)(\d{4,4})(\.vtk)',re.IGNORECASE)
        m = p.match(self.filename)
        
        if (m):
            #Seems a sequence
            fileBaseName = m.group(1)
            seqStr = m.group(2)
            fileExtension = m.group(3)

            startFrame = int(seqStr)
            endFrame = int(seqStr)

            #Look for starting frame
            while os.path.exists( os.path.join(self.directory,  fileBaseName + str(startFrame - 1).zfill(4) + fileExtension )) :
                startFrame = startFrame - 1            
            
            #Look for ending frame
            while os.path.exists( os.path.join(self.directory,  fileBaseName + str(endFrame + 1).zfill(4) + fileExtension )) :
                endFrame = endFrame + 1
                
            print("Start frame: "+str(startFrame)+" End frame: "+str(endFrame))
            print(fileBaseName + self.DsphObjEnum)
            
            createObject (fileBaseName + "_" + self.DsphObjEnum,
                         self.filename, self.directory, fileBaseName, fileExtension,
                         self.DsphObjEnum, self.DsphSmooth, self.DsphValidate,
                         self.DsphUV, self.DsphBlur, startFrame, endFrame)
            
        else:
            #This doesn't seems a sequence
            print("Not a sequence!!")
            
            #Maybe it is a lonely vtk
            p = re.compile('(.*)(\.vtk)',re.IGNORECASE)
            m = p.match(self.filename)
            
            if (m):
                createObject (m.group(1) + "_" + self.DsphObjEnum,
                             self.filename, self.directory, m.group(1), m.group(2),
                             self.DsphObjEnum, self.DsphSmooth, self.DsphValidate,
                             self.DsphUV, self.DsphBlur, 0, 0)
        
        return {'FINISHED'}
 
    def invoke(self, context, event):
        context.window_manager.fileselect_add(self)
        return {'RUNNING_MODAL'}

class MESH_OT_dsph_object_add(bpy.types.Operator):
    '''Add a DualSPHysics Object'''
    bl_idname = "mesh.dsph_object_add"
    bl_label = "Add a DualSPHysics Object"
    bl_options = {'REGISTER', 'UNDO'}
    
    def execute(self, context):
        bpy.ops.dialog.browsefile('INVOKE_DEFAULT')
        return {'FINISHED'}

class OBJECT_OT_RunFoamSimulation(bpy.types.Operator):
    bl_idname = "foam.runsimulation"
    bl_label = "Run Foam Simulation"
     
    def execute(self, context):
        if parseXML(context):
            try:
                diffuseparticles.run(
                    bpy.path.abspath(context.scene.DsphFoamInputPath),  # Input path
                    context.scene.DsphFoamInputPrefix,                  # Input files prefix
                    bpy.path.abspath(context.scene.DsphFoamPath),       # Output path
                    context.scene.DsphFoamPrefix,                       # Output files prefix
                    "",                                                 # Exclusion zone file
                    context.scene.DsphFoamStart,                        # Starting timestep
                    context.scene.DsphFoamEnd,                          # Ending timestep
                    4,                                                  # Zero padding
                    False,                                              # Output text files
                    True,                                               # Output vtk files
                    False,                                              # Output vtk extra info files
                    False,                                              # Output vtk fluid info files
                    context.scene.DsphFoamH,                            # H value
                    context.scene.DsphFoamMass,                         # Fluid particle mass value
                    context.scene.DsphFoamTimeStep,                     # Timestep duration
                    context.scene.DsphFoamMinX,                         # Domain limits
                    context.scene.DsphFoamMinY,
                    context.scene.DsphFoamMinZ,
                    context.scene.DsphFoamMaxX,
                    context.scene.DsphFoamMaxY,
                    context.scene.DsphFoamMaxZ,
                    context.scene.DsphFoamMinTrappedAir,
                    context.scene.DsphFoamMaxTrappedAir,
                    context.scene.DsphFoamMinWaveCrests,
                    context.scene.DsphFoamMaxWaveCrests,
                    context.scene.DsphFoamMinKinetic,
                    context.scene.DsphFoamMaxKinetic,
                    context.scene.DsphFoamTAMult,
                    context.scene.DsphFoamWCMult,
                    context.scene.DsphFoamSprayDensity,
                    context.scene.DsphFoamBubblesDensity,
                    context.scene.DsphFoamLifetime,
                    context.scene.DsphFoamBuoyancy,
                    context.scene.DsphFoamDrag)
                    
                fileName = context.scene.DsphFoamPrefix + str(context.scene.DsphFoamStart).zfill(4) + ".vtk"

                createObject (context.scene.DsphFoamPrefix + "_FOAM",
                    fileName,
                    context.scene.DsphFoamPath,
                    context.scene.DsphFoamPrefix,
                    ".vtk",
                    "FOAM",
                    True,
                    False,
                    False,
                    True,
                    context.scene.DsphFoamStart,
                    context.scene.DsphFoamEnd)                
                        
                return{'FINISHED'}  
            except:
                showPopup("Something went wrong with the simulation. Take a look to the system console to get more information.", "Error", "ERROR")
                return{'CANCELLED'}
        else:
            showPopup("Error: the XML cannot be loaded. Make sure to choose the XML file generated by the gencase tool.", "Error", "ERROR")
            return{'CANCELLED'}


# TODO: scene properties to register/unregister functions 
# http://blender.stackexchange.com/questions/2382/how-to-add-a-select-path-input-in-a-ui-addon-script
# Foam simulation panel
class PROPERTIES_PT_FoamSimulationPanel(bpy.types.Panel):
    bl_label = "DualSPHysics foam simulation"
    bl_space_type = "PROPERTIES"
    bl_region_type = "WINDOW"
    bl_context = "physics"

    bpy.types.Scene.DsphFoamInputPath = bpy.props.StringProperty(name = "Input Path",
                                                                 description = "Input path of particle data.",
                                                                 subtype = 'DIR_PATH')
    
    bpy.types.Scene.DsphFoamInputPrefix = bpy.props.StringProperty(name = "Input Prefix",
                                                                   description = "Input files preffix")
    
    bpy.types.Scene.DsphFoamPath = bpy.props.StringProperty(name = "Output Path",
                                                            description = "Output path for foam simulation",
                                                            subtype = 'DIR_PATH')
    
    bpy.types.Scene.DsphFoamPrefix = bpy.props.StringProperty(name = "Output Prefix",
                                                              description = "Output files preffix")
    
    bpy.types.Scene.DsphFoamXML = bpy.props.StringProperty(name = "XML",
                                                           description = "DualSPHysics XML file",
                                                           subtype = 'FILE_PATH')
    
    bpy.types.Scene.DsphFoamStart = bpy.props.IntProperty(name = "Starting step",
                                                          description = "First step of the foam simulation")
    
    bpy.types.Scene.DsphFoamEnd = bpy.props.IntProperty(name = "Ending step",
                                                        description = "Last step of the foam simulation")
    
    bpy.types.Scene.DsphFoamMinTrappedAir = bpy.props.FloatProperty(name = "Min Trapped Air Threshold",
                                                                    description = "Min Trapped Air Threshold",
                                                                    min = 0,
                                                                    default = 5.)

    bpy.types.Scene.DsphFoamMaxTrappedAir = bpy.props.FloatProperty(name = "Max Trapped Air Threshold",
                                                                    description = "Max Trapped Air Threshold",
                                                                    min = 0,
                                                                    default = 20.)
    
    bpy.types.Scene.DsphFoamMinWaveCrests = bpy.props.FloatProperty(name = "Min Wave Crests Threshold",
                                                                    description = "Min Wave Crests Threshold",
                                                                    min = 0,
                                                                    default = 2.)

    bpy.types.Scene.DsphFoamMaxWaveCrests = bpy.props.FloatProperty(name = "Max Wave Crests Threshold",
                                                                    description = "Max Wave Crests Threshold",
                                                                    min = 0,
                                                                    default = 8.)
    
    bpy.types.Scene.DsphFoamMinKinetic = bpy.props.FloatProperty(name = "Min Kinetic Energy Threshold",
                                                                 description = "Min Kinetic Energy Threshold",
                                                                 min = 0,
                                                                 default = 5.)

    bpy.types.Scene.DsphFoamMaxKinetic = bpy.props.FloatProperty(name = "Max Kinetic Energy Threshold",
                                                                 description = "Max Kinetic Energy Threshold",
                                                                 min = 0,
                                                                 default = 50.)
    
    bpy.types.Scene.DsphFoamTAMult = bpy.props.FloatProperty(name = "Trapper air multiplier",
                                                             description = "Amount of diffuse material created by trapper air",
                                                             min = 0,
                                                             default = 40.)
    
    bpy.types.Scene.DsphFoamWCMult = bpy.props.FloatProperty(name = "Wave crests multiplier",
                                                             description = "Amount of diffuse material created by wave crests",
                                                             min = 0,
                                                             default = 40.)
    
    bpy.types.Scene.DsphFoamSprayDensity = bpy.props.FloatProperty(name = "Max spray density factor",
                                                                   description = "Max spray density factor",
                                                                   min = 0,
                                                                   default = 6.0)
    
    bpy.types.Scene.DsphFoamBubblesDensity = bpy.props.FloatProperty(name = "Min bubbles density factor",
                                                                     description = "Min bubbles density factor",
                                                                     min = 0,
                                                                     default= 9.0)
    
    bpy.types.Scene.DsphFoamLifetime = bpy.props.FloatProperty(name = "Foam lifetime factor",
                                                               description = "Foam particles lifetime factor",
                                                               min = 0,
                                                               default= 10.0)
    
    bpy.types.Scene.DsphFoamBuoyancy = bpy.props.FloatProperty(name = "Buoyancy Control",
                                                               description = "Bubble particles buoyancy control",
                                                               min = 0,
                                                               max = 1,
                                                               default= 0.8)
    
    bpy.types.Scene.DsphFoamDrag = bpy.props.FloatProperty(name = "Drag control",
                                                           description = "Bubble particles drag control",
                                                           min = 0,
                                                           max = 1,
                                                           default= 0.5)
    
    bpy.types.Scene.DsphFoamCustomDomain = bpy.props.BoolProperty(name = "Enable custom domain limits", 
                                                                  description = "Enable custom domain limits",
                                                                  default = False)
    
    bpy.types.Scene.DsphFoamMinX = bpy.props.FloatProperty(name = "Min X",
                                                           description = "Min X")

    bpy.types.Scene.DsphFoamMinY = bpy.props.FloatProperty(name = "Min Y",
                                                           description = "Min Y")

    bpy.types.Scene.DsphFoamMinZ = bpy.props.FloatProperty(name = "Min Z",
                                                           description = "Min Z")

    bpy.types.Scene.DsphFoamMaxX = bpy.props.FloatProperty(name = "Max X",
                                                           description = "Max X")

    bpy.types.Scene.DsphFoamMaxY = bpy.props.FloatProperty(name = "Max Y",
                                                           description = "Max Y")

    bpy.types.Scene.DsphFoamMaxZ = bpy.props.FloatProperty(name = "Max Z",
                                                           description = "Max Z")
        
    # Hidden properties
    bpy.types.Scene.DsphFoamH = bpy.props.FloatProperty(name = "H")
    bpy.types.Scene.DsphFoamMass = bpy.props.FloatProperty(name = "Mass")
    bpy.types.Scene.DsphFoamTimeStep = bpy.props.FloatProperty(name = "TimeStep")

    
    @classmethod
    def poll(self, context):
        return context.object and 'DsphObjType' in context.object and context.object['DsphObjType'] == "FLUID"
    
    def draw(self, context):        
        layout = self.layout
        
        ds = False
        
        layout.label(text="File paths:")

        row = layout.row()
        row.alert = not os.path.isdir(bpy.path.abspath(context.scene.DsphFoamPath))
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamInputPath")
        
        layout.prop(context.scene, "DsphFoamInputPrefix")
        
        row = layout.row()
        row.alert = not os.path.isdir(bpy.path.abspath(context.scene.DsphFoamPath))
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamPath")
        
        layout.prop(context.scene, "DsphFoamPrefix")
        
        row = layout.row()
        row.alert = (not os.path.isfile(bpy.path.abspath(context.scene.DsphFoamXML)))
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamXML")   
        
        layout.label(text="Time range:")
        row = layout.row()
        row.alert = context.scene.DsphFoamStart > context.scene.DsphFoamEnd
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamStart")
        row.prop(context.scene, "DsphFoamEnd")
        
        layout.label(text="Configurable thresholds:")
        row = layout.row()
        row.alert = context.scene.DsphFoamMinTrappedAir > context.scene.DsphFoamMaxTrappedAir
        ds = ds or row.alert              
        row.prop(context.scene, "DsphFoamMinTrappedAir")
        row.prop(context.scene, "DsphFoamMaxTrappedAir")   
        
        row = layout.row()
        row.alert = context.scene.DsphFoamMinWaveCrests > context.scene.DsphFoamMaxWaveCrests
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamMinWaveCrests")
        row.prop(context.scene, "DsphFoamMaxWaveCrests")
        
        row = layout.row()
        row.alert = context.scene.DsphFoamMinKinetic > context.scene.DsphFoamMaxKinetic
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamMinKinetic")
        row.prop(context.scene, "DsphFoamMaxKinetic")
        
        layout.label(text="Diffuse material volume:")
        layout.prop(context.scene, "DsphFoamTAMult")
        layout.prop(context.scene, "DsphFoamWCMult")
        
        layout.label(text="Particle category:")
        layout.prop(context.scene, "DsphFoamSprayDensity")
        layout.prop(context.scene, "DsphFoamBubblesDensity")
        
        layout.label(text="Foam lifetime:")
        layout.prop(context.scene, "DsphFoamLifetime")
        
        layout.label(text="Bubbles advection:")
        layout.prop(context.scene, "DsphFoamBuoyancy")
        layout.prop(context.scene, "DsphFoamDrag")
        
        layout.label(text="Domain limits:")
        layout.prop(context.scene, "DsphFoamCustomDomain")
        
        domaincol = layout.column()
        row = domaincol.row()
        row.alert = context.scene.DsphFoamMinX > context.scene.DsphFoamMaxX        
        ds = ds or row.alert        
        row.prop(context.scene, "DsphFoamMinX")
        row.prop(context.scene, "DsphFoamMaxX")
        
        row = domaincol.row()
        row.alert = context.scene.DsphFoamMinY > context.scene.DsphFoamMaxY        
        ds = ds or row.alert
        row.prop(context.scene, "DsphFoamMinY")
        row.prop(context.scene, "DsphFoamMaxY")        

        row = domaincol.row()
        row.alert = context.scene.DsphFoamMinZ > context.scene.DsphFoamMaxZ
        ds = ds or row.alert                
        row.prop(context.scene, "DsphFoamMinZ")
        row.prop(context.scene, "DsphFoamMaxZ")
        
        domaincol.enabled = context.scene.DsphFoamCustomDomain            
        
        buttoncol = layout.column()
        buttoncol.enabled = not ds
        buttoncol.operator("foam.runsimulation")
        
#
#   Handlers
#

## @brief Handler for the frame change.
# This is a Blender handler. It updates the objects of the simulation each time a frame change is produced.
# @param scene Current scene object.
def DsphHandler (scene):
    nFrame = scene.frame_current 
        
    print("Frame Change", nFrame)
    
    for do in bpy.context.scene.objects:
        if 'DsphObjType' in do and nFrame >= do["DsphStartFrame"] and nFrame <= do["DsphEndFrame"] :
            
            fileName = bpy.path.abspath(os.path.join(do["DsphPathName"],  do["DsphBaseName"] + str(nFrame).zfill(4) + do["DsphExtension"]))
            
            if not os.path.exists(fileName):
                print("Error: The following file does not exists: " + fileName)
                return

            print("Reading file " + fileName + " ...")

            pts = []
            tris = []
            vels = []

            try:
                if do["DsphObjType"] == "ROPE":
                    t_pts, t_tris = vtkimporter.load(fileName)
                    pts, tris = makeRope(t_pts)
                elif do["DsphObjType"] == "FOAM":
                    pts, tris, vels = vtkimporter.loaddiffuse(fileName)
                elif do["DsphObjType"] == "FLUID" and do["DsphBlur"]:
                    pts, tris, vels = vtkimporter.loadvel(fileName)
                else:
                    pts, tris = vtkimporter.load(fileName)                                    
            except:
                print("Error: cant read "+fileName)
                return
            
            print("Reading finished!")    

            newMesh = bpy.data.meshes.new(do.name + 'Mesh' + str(nFrame)) # New mesh
            newMesh.from_pydata(pts,[],tris)   # edges or faces should be [], or you ask for problems
            newMesh.update(calc_edges=False)    # Update mesh with new data

            if do["DsphValidate"] :
                newMesh.validate() # This is necessary to avoid lines along the mesh
                newMesh.update(calc_edges=False)
             
            oldMesh = do.data              
             
            #Transfer UVMaps. Based on source code from Blender's object.py file.
            if do["DsphUV"] and oldMesh.uv_layers :
                nbr_loops = len(oldMesh.loops)

                # seems to be the fastest way to create an array
                uv_array = array.array('f', [0.0] * 2) * nbr_loops
                oldMesh.uv_layers.active.data.foreach_get("uv", uv_array)

                uv_other = newMesh.uv_layers.active
                if not uv_other:
                    newMesh.uv_layeras.new()
                    uv_other = newMesh.uv_layers.active
                    if not uv_other:
                        self.report({'ERROR'}, "Could not add a new UV map to object "
                                    f"'{obj.name}' (Mesh '{newMesh.name}')\n")

                    # finally do the copy
                    uv_other.data.foreach_set("uv", uv_array)

                newMesh.update(calc_edges=False)

            # Copy materials           
            for tempMaterial in oldMesh.materials:
                if tempMaterial != None:
                    newMesh.materials.append(tempMaterial)
            newMesh.update(calc_edges=False)
                    
            # This is to copy materials for floating objects
            # TODO: debug this
            if(do["DsphObjType"] == "OTHER" and 
               oldMesh.polygons.__len__() == newMesh.polygons.__len__()):
                for i in range(0, oldMesh.polygons.__len__()) :
                    newMesh.polygons[i].material_index = oldMesh.polygons[i].material_index
                newMesh.update(calc_edges=False)

            # Smooth shading
            if do["DsphSmooth"]:
                for poly in newMesh.polygons:
                    poly.use_smooth = (do["DsphSmooth"] == True)
                newMesh.update(calc_edges=False)

            if do["DsphValidate"]:    
                newMesh.validate()
                
            newMesh.update(calc_edges=False)
            do.data=newMesh

            # Motion blur code
            if do["DsphBlur"]:                    
                if do["DsphObjType"] == "FLUID" or do["DsphObjType"] == "FOAM":
                    do.shape_key_add(name="Base",from_mix=False) 
                    bm = bmesh.new()
                    bm.from_mesh(do.data)
                    key_shape = bm.verts.layers.shape.new("key_blur")

                    c=0;
                    for i in bm.verts:
                        i[key_shape] = vels[c]
                        c = c+1

                    bm.to_mesh(do.data)
       
                    do.data.shape_keys.animation_data_clear()
                    kb = do.data.shape_keys.key_blocks["key_blur"]
                    kb.value = 1.
                    kb.keyframe_insert("value",frame=nFrame+1)
                    kb.value = 0.
                    kb.keyframe_insert("value",frame=nFrame-1)                  

                    scene.view_layers.update()
                    bm.free()
                    # =========================
                elif do["DsphObjType"] == "OTHER":
                    do.shape_key_add(name="Base",from_mix=False) 
                    bm = bmesh.new()
                    bm.from_mesh(do.data)
                    key_shape = bm.verts.layers.shape.new("key_blur")

                    fileName = os.path.join(do["DsphPathName"],  do["DsphBaseName"] + str(nFrame + 1).zfill(4) + do["DsphExtension"] )
                    pts2, tris2 = vtkimporter.load(fileName)

                    c=0;
                    for i in bm.verts:
                        i[key_shape] = pts2[c]
                        c = c+1

                    bm.to_mesh(do.data)
       
                    do.data.shape_keys.animation_data_clear()
                    kb = do.data.shape_keys.key_blocks["key_blur"]
                    kb.value = 1.
                    kb.keyframe_insert("value",frame=nFrame+1)
                    kb.value = 0.
                    kb.keyframe_insert("value",frame=nFrame-1)                  

                    scene.view_layers.update()
                    bm.free()                    
            
            deleteMesh(oldMesh.name)  
            print("Object updated!")  

def preRenderHandler (scene):
    bpy.app.handlers.frame_change_pre.remove(DsphHandler)
    bpy.app.handlers.render_pre.append(DsphHandler)
    
def postRenderHandler (scene):
    bpy.app.handlers.render_pre.remove(DsphHandler)
    bpy.app.handlers.frame_change_pre.append(DsphHandler)    
 
#
#    Registration
#    Makes it possible to access the script from the Add > Mesh menu
#
 
def menu_func(self, context):
    self.layout.operator("mesh.dsph_object_add", 
        text="DualSPHysics Object", 
        icon='MOD_FLUIDSIM')

classes = (
    VIEW3D_PT_DsphObjectPanel, 
    FILE_BROWSER_OT_BrowseFileOperator, 
    MESH_OT_dsph_object_add,
    OBJECT_OT_RunFoamSimulation, 
    PROPERTIES_PT_FoamSimulationPanel
)

def register():
    print("Register")

    for cls in classes:
        bpy.utils.register_class(cls)

    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)

    bpy.app.handlers.persistent(DsphHandler)
    bpy.app.handlers.persistent(preRenderHandler)
    bpy.app.handlers.persistent(postRenderHandler)
    
    bpy.app.handlers.frame_change_pre.append(DsphHandler)
    bpy.app.handlers.render_init.append(preRenderHandler)
    bpy.app.handlers.render_cancel.append(postRenderHandler)
    bpy.app.handlers.render_complete.append(postRenderHandler)
 
def unregister():
    print("Unregister")

    for cls in reversed(classes):
        bpy.utils.unregister_class(cls)

    bpy.types.VIEW3D_MT_mesh_add.append(menu_func)
    
    bpy.app.handlers.frame_change_pre.remove(DsphHandler)
    bpy.app.handlers.render_init.remove(preRenderHandler)
    bpy.app.handlers.render_cancel.remove(postRenderHandler)
    bpy.app.handlers.render_complete.remove(postRenderHandler)    
 
if __name__ == "__main__":
    register()
