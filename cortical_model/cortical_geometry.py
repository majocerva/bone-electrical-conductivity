#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  3 11:17:08 2025

#modificado Julio 2025 - falta red

@author: mjcervantes
"""

import os
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random
from scipy.stats import qmc
from scipy.spatial import distance_matrix
import math

import sys
import gmsh

import subprocess
import time

from mpi4py import MPI
import ufl
from dolfinx import mesh
from dolfinx import fem
from dolfinx import default_scalar_type
from dolfinx.io import gmshio
from dolfinx.fem.petsc import LinearProblem
from dolfinx import io



#%%


# ==============================================================================
# PARAMETROS (default)
# ==============================================================================

# - Clases de parámetros del la osteona, los canalaes de Havers (verticales) y los canales de Volkmanns (horizontales).
#
class SAMPLE_parameters:
  
   # Diameters Osteon (mm)
   cube_size = 1
   density = 15
  
class OSTEON_parameters:
  
   # Diameters Osteon (mm)
   d_min = 0.1
   d_max = 0.25
   
class HAVERSIAN_CANALS_parameters:
    
   # Diameters Haversian canals (mm)
   d_min = 0.04#0.05
   d_max = 0.09#0.1
   # Inclination angle Haversian canals
   a_min = 0
   a_max = 15

class VOLKMANNS_CANALS_parameters:
   
      # Diameters Volkmanns canals (mm)
   d_min = 0.04#0.05
   d_max = 0.05#0.08
   # Inclination angle Volkmanns canals
   a_min = 0
   a_max = 15
   # Distance Volkmanns canals (mm)
   #Volk_dist_min = 0.5
   #Volk_dist_max = 0.15

class DIELECTRIC_parameters:
    
    sigma_matrix = 3.841  # en mS/m
    sigma_marrow = 267 # en mS/m
   
class Seed:   
      semilla = 2



#%%
# ==============================================================================
# 1-Genera Modelo
# ==============================================================================

def generate_model(name, folder ):
      
      rng = np.random.default_rng(Seed.semilla)
      
      # Funciones -----------
      
      def find_coord(center, direction, z_z):

          x = center[0]
          y = center[1]
          z = center[2]

          dx = direction[0]/SAMPLE_parameters.cube_size
          dy = direction[1]/SAMPLE_parameters.cube_size
          dz = direction[2]/SAMPLE_parameters.cube_size

          t = (z_z-z)/dz

          x_z = x + dx*t
          y_z = y + dy*t
          z_z = z_z

          center = x_z, y_z, z_z

          return center

      def generate_random_XY_Vertical(num_Vertical, area_size, min_distance):
          XY_Vertical = []

          l_bounds = [-area_size/2]
          u_bounds = [area_size/2]
          sampler = qmc.Sobol(2, scramble=True, seed=rng)
          X = sampler.random_base2(22)
          X = qmc.scale(X, l_bounds, u_bounds)

          i = -1
          while len(XY_Vertical) < num_Vertical:

              i = i+1

              # Generar coordenadas aleatorias para un nodo
              x = float(X[i][0])
              y = float(X[i][1])

              # Verificar si el nodo se solapa con algún nodo existente
              is_overlapping = False
              for (n_x, n_y) in XY_Vertical:
                  distance = math.sqrt((x - n_x) ** 2 + (y - n_y) ** 2)
                  if min_distance > distance:
                      is_overlapping = True
                      break

                  # Si no se solapa, añadir el nodo a la lista
              if not is_overlapping:

                  XY_Vertical.append([x, y])

          return XY_Vertical  # lista

      def connection(XY_Vertical):

          G = nx.Graph()
          for i in range(int(num_Vertical)):
              G.add_node(i, pos=XY_Vertical)
          
          conexion = []
          for i in range(len(XY_Vertical)):

              distancias = distance_matrix(XY_Vertical, XY_Vertical)
              # Ordenar distancias y obtener los índices de los cilindros más cercanos
              indices_cercanos = np.argsort(distancias[i])

              # [1] = más cercano, [2] = segundo más cercano
              vecinos_cercanos = indices_cercanos[1:10]

              conectados = 0  # Contador de conexiones exitosas

              for j in vecinos_cercanos:
                  if not G.has_edge(i, j):
                      G.add_edge(i, j)
                      conexion.append((i, int(j)))
                      conectados += 1

                      if conectados >= 4:  # Si ya conectó con 2, pasar al siguiente nodo
                       break
                   
          #pos=XY_Vertical  
          #nx.draw(G, pos, node_size=1, node_color='red')
          return conexion  # lista de tuplas

      def generate_vertical(XY_Vertical):

          sampler = qmc.Sobol(1, scramble=True, seed=rng)
          R = sampler.random_base2(14)

          R_O = qmc.scale(R, OSTEON_parameters.d_min /
                          2, OSTEON_parameters.d_max/2)
          R_H = qmc.scale(R, HAVERSIAN_CANALS_parameters.d_min /
                          2, HAVERSIAN_CANALS_parameters.d_max/2)

          n_z = -cube_size2/2

          Cyl_V = []
          k = 0
          for (n_x, n_y) in XY_Vertical:

              ang_z = rng.uniform(-HAVERSIAN_CANALS_parameters.a_max,
                                  HAVERSIAN_CANALS_parameters.a_max)
              ang_fi = rng.uniform(0, 360)

              x = n_x
              y = n_y
              z = n_z

              center = (x, y, z)

              dx = math.sin(math.pi*ang_z/180) * \
                  math.cos(math.pi*ang_fi/180) * cube_size2
              dy = math.sin(math.pi*ang_z/180) * \
                  math.sin(math.pi*ang_fi/180) * cube_size2
              dz = math.cos(math.pi*ang_z/180) * cube_size2

              direction = (dx, dy, dz)

              radio_O = R_O[k][0]
              radio_H = R_H[k][0]

              radio = (radio_O, radio_H)

              Cyl_V.append([center, direction, radio])  # para el gmsh

              k = k+1

          return Cyl_V

      def generate_horizontal(conexion):

          Cyl_H = []
          for (n1, n2) in conexion:

              radio = rng.uniform(VOLKMANNS_CANALS_parameters.d_min/2,
                                  VOLKMANNS_CANALS_parameters.d_max/2)
              z_1 = rng.uniform(-SAMPLE_parameters.cube_size/2,
                                SAMPLE_parameters.cube_size/2)
              ang_z = rng.uniform(0,  VOLKMANNS_CANALS_parameters.a_max)

          # suponemos que arranca en cyl 1
              x1 = Cyl_V[n1][0][0]
              y1 = Cyl_V[n1][0][1]
              z1 = Cyl_V[n1][0][2]

              center1 = (x1, y1, z1)

              dx1 = Cyl_V[n1][1][0]
              dy1 = Cyl_V[n1][1][1]
              dz1 = Cyl_V[n1][1][2]

              direction1 = (dx1, dy1, dz1)

              # busco el centro del cilinro en z_1
              center1_z = find_coord(center1,  direction1, z_1)

              x2 = Cyl_V[n2][0][0]
              y2 = Cyl_V[n2][0][1]
              z2 = Cyl_V[n2][0][2]

              center2 = (x2, y2, z2)  # en z negativo

              dx2 = Cyl_V[n2][1][0]
              dy2 = Cyl_V[n2][1][1]
              dz2 = Cyl_V[n2][1][2]

              direction2 = (dx2, dy2, dz2)

              center2_z = find_coord(center2,  direction2, z_1)  # En z1

              dist_x = center2_z[0] - center1_z[0]
              dist_y = center2_z[1] - center1_z[1]
              dist_horizontal = math.sqrt(dist_x**2 + dist_y**2)

              z_2 = math.sin(ang_z*math.pi/180) * dist_horizontal + z_1
              center3_z = find_coord(center2,  direction2, z_2)

              dx = (center3_z[0]-center1_z[0])
              dy = (center3_z[1]-center1_z[1])
              dz = z_2 - z_1

              direction = (dx, dy, dz)

              Cyl_H.append([center1_z, direction, radio])

          return Cyl_H
      
      
      # Creación de goemtria -----------
      
      area_size = SAMPLE_parameters.cube_size**2
      cube_size2 = SAMPLE_parameters.cube_size + 2 * HAVERSIAN_CANALS_parameters.d_max #cubo mas grande que la muestra para segruar ostonas y conexiones afuera
      area_size2 = cube_size2 ** 2
      
      num_Vertical = int(SAMPLE_parameters.density * area_size2)
      min_distance = (OSTEON_parameters.d_max + OSTEON_parameters.d_min)/ 2
      
      XY_Vertical = generate_random_XY_Vertical(num_Vertical, area_size2, min_distance)
      conexion = connection(XY_Vertical)
        
      Cyl_V = generate_vertical(XY_Vertical)
      Cyl_H = generate_horizontal(conexion)

      
      # Plot -------
      # Creo un grafo  2D para plotear
      G = nx.Graph()
      # Generar nodos
      for i in range (int( num_Vertical)):
          G.add_node(i, pos=XY_Vertical)
      pos=XY_Vertical
      # Generar aristas
      conexiones = conexion
      for i, j in conexiones:
              G.add_edge(i, j)
        
      fig, axs = plt.subplots(figsize=(8, 8))
      # Dibujar el grafo usando las posiciones especificadas
      nx.draw(G, pos, node_size=1, node_color='red')
     
      # Visualizar nodos
      rect = plt.Rectangle((-area_size/2, -area_size/2), area_size, area_size, linewidth=1, edgecolor='black', facecolor='azure', alpha = 0.6)
      axs.add_patch(rect)

      for i in range(len(XY_Vertical)):
          x = XY_Vertical[i][0]
          y = XY_Vertical[i][1]
          
          radius_Ost = Cyl_V[i][2][0]
          radius_Hav = Cyl_V[i][2][1]
      
          circle=plt.Circle((x,y),radius_Ost,  facecolor='wheat')
          axs.add_patch(circle)
          
          circle2=plt.Circle((x,y),radius_Hav,color='red')
          axs.add_patch(circle2)
          
        
          # axs.plot(x, y, 'ko', markersize=10)  # 
          # axs.text(x, y, str(i), fontsize=12, ha='right')  # Etiqueta con el índice del nodo
      rect = plt.Rectangle((-area_size/2, -area_size/2), area_size, area_size, linewidth=1, edgecolor='black', facecolor='azure', alpha = 0.3)
      axs.add_patch(rect)
      axs.set_xlim(-area_size/2, area_size/2)
      axs.set_ylim(-area_size/2, area_size/2)
      axs.grid(True)
      
      
      #save 
      # plt.savefig('Diag_Matriz_cortical.jpg', bbox_inches='tight')  #jpg
      # plt.savefig('Diag_Matriz_cortical.pdf', bbox_inches='tight')  #PDF
      # plt.show()
      #  -------
    
      # Genereter Gmsh Script  (.geo) -----------
    
      gmsh_script_introduction = """//                                Cube with cylindrical inclusions\n
      //----------------------------------------------------------------------------------------------------------------------------
      \n\n//Geometry\n\n
        
      SetFactory("OpenCASCADE");\n"""
        
      gmsh_script_environment = """\n\n// Environment\n\n"""
      gmsh_script_environment +="Box(1)"+" = {{{}, {}, {}, {}, {}, {}}};\n".format(-SAMPLE_parameters.cube_size/2,-SAMPLE_parameters.cube_size/2,-SAMPLE_parameters.cube_size/2,SAMPLE_parameters.cube_size,SAMPLE_parameters.cube_size,SAMPLE_parameters.cube_size) 
      gmsh_script_inclusion ="""\n\n// Inclusions\n\n"""
        
        
      # Osteons-Havers Canals
      for i, (center,direction,radio) in enumerate(Cyl_V):
        gmsh_script_inclusion += "Cylinder({}) = {{{}, {}, {},{}, {}, {},{}}} ;\n".format(i+2, center[0], center[1], center[2],direction[0],direction[1],direction[2], radio[1])
 
      # Volkmann Canals  
      for  j, (center,direction,radio) in enumerate(Cyl_H): 
          gmsh_script_inclusion += "Cylinder({}) = {{{}, {}, {},{}, {}, {},{}}} ;\n".format(j+3+i,center[0], center[1], center[2],direction[0],direction[1],direction[2], radio)
        
      k=j+3+i
      # Boolean Operation
      gmsh_script_Bool = """\n\n// Boolean Operation\n\n"""
      gmsh_script_Bool += "volUnion() = BooleanUnion{Volume{2:"+"{}".format(i+2)+"}; Delete;}{Volume{"+"{}".format(i+2+1)+":"+"{}".format(k)+"}; Delete;};\n"
      gmsh_script_Bool += "g() = BooleanIntersection{ Volume{volUnion()};Delete; }{ Volume{1};  };\n"
      gmsh_script_Bool += "f() = BooleanDifference{ Volume{1}; Delete;  }{ Volume{g()};};\n"
        
      # Physical Volume
      n=500
      gmsh_script_physical = """\n\n//Physical Volume\n\n"""
     
      gmsh_script_physical += "Physical Volume({})".format(n)+"={g()}"+";\n"
      gmsh_script_physical +="Physical Volume({})".format(n+100)+"={f()};\n"
        
        
      # # # Physical Surface
      # # n=10
      # gmsh_script_physical_s = """\n\n//Dummy Physical Surface\n\n"""
      # gmsh_script_physical_s += "Physical Surface(10)={160};\n"
      # gmsh_script_physical_s += "Physical Surface(20)={163};\n"
        
      # Characteristic Length
      gmsh_script_characteristic = "q() = PointsOf{Volume{g()}; }; Characteristic Length{q()} = "+"{}".format(VOLKMANNS_CANALS_parameters.d_min/2)+";\n"
      gmsh_script_characteristic += "Coherence;\n"
      gmsh_script_characteristic += "\nMesh.CharacteristicLengthMax = {}".format(HAVERSIAN_CANALS_parameters.d_min/2)+";\n"
      #Escala
      gmsh_script_characteristic += "\nMesh.ScalingFactor = {}".format(1e-3)+";"
        
        
      # generate .geo
      path_save = os.path.join(folder, name + ".geo") 
      with open(path_save, "w") as f:
        
          f.write(gmsh_script_introduction)
          f.write(gmsh_script_environment)
          f.write(gmsh_script_inclusion)
          # f.write(gmsh_script_plane_s)
          f.write(gmsh_script_Bool)
          f.write(gmsh_script_physical)
          #f.write(gmsh_script_physical_s)
          f.write(gmsh_script_characteristic)
            
          print(f"--{name} generado con éxito--")
          
      return len(Cyl_H), len(Cyl_V)
  
#%%
# ==============================================================================
# 2- Calcula porosidad
# ==============================================================================
def calculate_porosity (name, folder):
      
    geo_file =  os.path.join(folder, f"{name}.geo")
    
    # Inicializar y cargar el archivo .geo
    gmsh.initialize()
    
   
    gmsh.open(geo_file)

    # Sincronizar la geometría
    gmsh.model.occ.synchronize()
    
       
    hueso_tag = gmsh.model.getEntitiesForPhysicalGroup(3, 600)
    poros_tags = gmsh.model.getEntitiesForPhysicalGroup(3, 500)
    
    # volumen total del hueso
    vol_hueso = sum(gmsh.model.occ.getMass(3, tag) for tag in hueso_tag)
    
    # volumen de los poros
    vol_poros = sum(gmsh.model.occ.getMass(3, tag) for tag in poros_tags)
    
    # volumen total (hueso + poros)
    vol_total = vol_hueso + vol_poros
    #print(vol_total)
    # porosidad
    phi = (vol_poros / vol_total)
    
    #print(f"--Porosidad de {modelo}  = {phi:.4f}--\n")
     # Imprime solo con dos decimales
    
    # # Finalizar Gmsh
    gmsh.finalize()
    print(f"--Porosidad calculada con éxito--")  # Confirmación breve
    
    return phi 

#%%
# ==============================================================================
# 3-Genera Mallado
# ==============================================================================

import os
import subprocess
import time

def generate_mesh(name, folder):
    msh_file = os.path.join(folder, f"{name}.msh")
    geo_file = os.path.join(folder, f"{name}.geo")

    def run_gmsh(geo_file, timeout=30):
        """Ejecuta Gmsh para mallar un archivo .geo con límite de tiempo."""
        try:
            subprocess.run(
                ["gmsh", "-3", geo_file],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                timeout=timeout
            )
        except subprocess.TimeoutExpired:
            raise TimeoutError(f"Gmsh tardó más de {timeout} segundos")
        except subprocess.CalledProcessError as e:
            error_msg = f"Error en Gmsh: {e.stderr.decode().strip()}"
            raise RuntimeError(error_msg) from e

    try:
        if not os.path.exists(geo_file):
            raise FileNotFoundError(f"No existe el archivo .geo: {geo_file}")
        
        # Primer intento
        try:
            run_gmsh(geo_file, timeout=30)
        except TimeoutError as e:
            print(f"[Aviso] Primer intento de mallado falló por tiempo: {e}")
            print("[Info] Reintentando generar la malla...")
            time.sleep(2)  # Esperar un momento antes de reintentar
            run_gmsh(geo_file, timeout=30)  # Segundo intento

        # Verificar si se generó el .msh
        if not os.path.exists(msh_file):
            raise RuntimeError(f"No se generó el archivo .msh: {msh_file}")

        print(f"--Malla generada con éxito--")
        return msh_file

    except Exception as e:
        print(f"Error en generate_mesh: {str(e)}")
        raise  # Relanza la excepción para que la capture el llamador

