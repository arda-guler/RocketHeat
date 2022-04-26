import OpenGL
from OpenGL.GL import *
from OpenGL.GLU import *
import glfw

import keyboard

def import_model():
    import_file = open("3d_model.txt", "r")
    import_lines = import_file.readlines()

    axials = [ [ [],[] ] ]
    x_index = 0
    outer = 0

    for line in import_lines:
        
        if line == "OUTERSHELL\n":
            outer = 1

        elif line == "NEWX\n":
            axials.append([[],[]])
            x_index += 1
            outer = 0
            
        else:
            line = line[1:-2]
            line = line.split(",")
            x = float(line[0])
            y = float(line[1])
            z = float(line[2])

            axials[x_index][outer].append([x, y, z])

    return axials

def main():
    print("Reading 3D model...")
    model_data = import_model()

    glfw.init()

    window = glfw.create_window(800, 600, "Engine Viewer", None, None)
    glfw.set_window_pos(window,200,200)
    glfw.make_context_current(window)
    
        
    gluPerspective(70, 800/600, 0.0005, 250.0)
    glEnable(GL_CULL_FACE)
    glPolygonMode(GL_FRONT_AND_BACK, GL_POINT)
    glTranslate(0,0,-0.2)
    
    while not glfw.window_should_close(window):
        glfw.poll_events()

        glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT)

        # controls
        if keyboard.is_pressed("w"):
            glTranslate(0, 0, 0.001)
        if keyboard.is_pressed("s"):
            glTranslate(0, 0, -0.001)
        if keyboard.is_pressed("d"):
            glTranslate(-0.001, 0, 0)
        if keyboard.is_pressed("a"):
            glTranslate(0.001, 0, 0)
        if keyboard.is_pressed("f"):
            glTranslate(0, 0.001, 0)
        if keyboard.is_pressed("r"):
            glTranslate(0, -0.001, 0)

        if keyboard.is_pressed("i"):
            glRotate(1, 1,0,0)
        if keyboard.is_pressed("k"):
            glRotate(1, -1,0,0)
        if keyboard.is_pressed("l"):
            glRotate(1, 0,1,0)
        if keyboard.is_pressed("j"):
            glRotate(1, 0,-1,0)
        if keyboard.is_pressed("o"):
            glRotate(1, 0,0,1)
        if keyboard.is_pressed("u"):
            glRotate(1, 0,0,-1)

        # rendering
        for axial in model_data[::20]:

            glColor(0.8, 0.2, 0.2)
            # inner wall vertices
            glPushMatrix()
            glBegin(GL_LINES)
            for ivertex in axial[0][::15]:
                glVertex3f(ivertex[0], ivertex[1], ivertex[2])
            glEnd()
            glPopMatrix()

            glColor(0.2, 0.5, 0.8)
            # outer wall vertex
            glPushMatrix()
            glBegin(GL_LINES)
            for overtex in axial[1]:
                glVertex3d(overtex[0], overtex[1], overtex[2])
            glEnd()
            glPopMatrix()
                
        glfw.swap_buffers(window)

    glfw.terminate()

main()
