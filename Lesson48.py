# -*- coding: utf-8 -*-
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
import math
import numpy
import PIL.Image

from ArcBall import * 				# ArcBallT and this tutorials set of points/vectors/matrix types
from geometry import *
from graph2 import *
# from matrix import *

PI2 = 2.0*3.1415926535			# 2 * PI (not squared!) 		// PI Squared

# *********************** Globals *********************** 
# Python 2.2 defines these directly
try:
	True
except NameError:
	True = 1==1
	False = 1==0

g_Transform = Matrix4fT ()
g_LastRot = Matrix3fT ()
g_ThisRot = Matrix3fT ()

g_ArcBall = ArcBallT (640, 480)
g_isDragging = False
g_quadratic = None

boundingBox = None;

solidFaces = [];
objectsToDraw = [];

graph = Graph();

angularSpeed = 0.6;

#          blue    green     cyan     red      pink     yellow   white
colors = ([0,0,1], [0,1,0], [0,1,1], [1,0,0], [1,0,1], [1,1,0], [1,1,1])

color_dict = {
	(0,0,1): "blue",
	(0,1,0): "green",
	(0,1,1): "cyan",
	(1,0,0): "red",
	(1,0,1): "pink",
	(1,1,0): "yellow",
	(1,1,1): "white",
	(0.502, 0.502, 0.502): "gray"
}

imageID = 0;
hasTexture = False;

# A general OpenGL initialization function.  Sets all of the initial parameters. 
def Initialize (Width, Height, argv):				# We call this right after our OpenGL window is created.
	global g_quadratic, solidFaces, objectsToDraw, graph, imageID;

	if len(argv) < 2:
		print "Entre com o nome do arquivo na linha de comando"
		return False

	glClearColor(0.0, 0.0, 0.0, 1.0)					# This Will Clear The Background Color To Black
	glClearDepth(1.0)									# Enables Clearing Of The Depth Buffer
	glDepthFunc(GL_LEQUAL)								# The Type Of Depth Test To Do
	glEnable(GL_DEPTH_TEST)								# Enables Depth Testing
	glShadeModel (GL_FLAT);								# Select Flat Shading (Nice Definition Of Objects)
	glHint (GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST) 	# Really Nice Perspective Calculations

	g_quadratic = gluNewQuadric();
	gluQuadricNormals(g_quadratic, GLU_SMOOTH);
	gluQuadricDrawStyle(g_quadratic, GLU_FILL); 
	# Why? this tutorial never maps any textures?! ? 
	# gluQuadricTexture(g_quadratic, GL_TRUE);			# // Create Texture Coords

	# glEnable (GL_LIGHT0)
	# glEnable (GL_LIGHTING)

	# glEnable (GL_COLOR_MATERIAL)

	vertices, faces = GetSolidFromFile(argv);
	solidFaces = BuildSolidStructure(vertices, faces);
	objectsToDraw.append(solidFaces);

	graph = BuildGraph(solidFaces);

	imageID = LoadImage();

	return True

def BuildGraph(polygons):
	graph = Graph();

	for polygon in polygons:
		graph.add_vertex(polygon);

	for polygon1 in polygons:
		for polygon2 in polygons:
			if polygon1.points != polygon2.points and DoPolygonsHaveAnEdgeInCommon(polygon1, polygon2):
				graph.add_edge({polygon1, polygon2});
	return graph;

def DoPolygonsHaveAnEdgeInCommon(poly1, poly2):
	edges1 = GeneratePolygonEdges(poly1);
	edges2 = GeneratePolygonEdges(poly2);

	for line1 in edges1:
		for line2 in edges2:
			if (line1.p1 == line2.p1 and line1.p2 == line2.p2) or (line1.p1 == line2.p2 and line1.p2 == line2.p1):
				return True;
	return False;

def GetEdgesInCommon(poly1, poly2):
	edges1 = GeneratePolygonEdges(poly1);
	edges2 = GeneratePolygonEdges(poly2);
	edgesInCommon = [];

	for line1 in edges1:
		for line2 in edges2:
			if (line1.p1 == line2.p1 and line1.p2 == line2.p2) or (line1.p1 == line2.p2 and line1.p2 == line2.p1):
				edgesInCommon.append(line1);
	return edgesInCommon;

def GetSolidFromFile(argv):
	global solidFaces, objectsToDraw;
	if len(argv) < 2:
		return None;

	f = open(argv[1]);
	lines = f.readlines();

	line = lines[3].split(" ");
	nVertices = int(line[2]);

	line = lines[7].split(" ");
	nFaces = int(line[2]);

	vertexIndex = 0;
	vertices = [];
	while vertexIndex < nVertices:
		v = lines[10 + vertexIndex].split(" ");
		vertices.append(Point(float(v[0]), float(v[1]), float(v[2])));
		vertexIndex += 1;

	faceIndex = 0;
	faces = [];
	while faceIndex < nFaces:
		face = lines[10 + vertexIndex + faceIndex].split(" ");
		faces.append([int(fac) for fac in face[1:] if str.isdigit(fac)]);
		faceIndex += 1;

	f.close();
	return vertices, faces;

def LoadImage(imageName = "Imagens/Velazquez.jpg"):
	"""Load an image file as a 2D texture using PIL"""

	# PIL defines an "open" method which is Image specific!
	im = PIL.Image.open(imageName)
	try:
		ix, iy, image = im.size[0], im.size[1], im.tobytes("raw", "RGBA", 0, -1)
	except (SystemError, ValueError):
		ix, iy, image = im.size[0], im.size[1], im.tobytes("raw", "RGBX", 0, -1)
	except AttributeError:
		ix, iy, image = im.size[0], im.size[1], im.tostring("raw", "RGBX", 0, -1)

	# Generate a texture ID
	ID = glGenTextures(1)

	# Make our new texture ID the current 2D texture
	glBindTexture(GL_TEXTURE_2D, ID)
	glPixelStorei(GL_UNPACK_ALIGNMENT,1)

	# Copy the texture data into the current texture ID
	glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, ix, iy, 0, GL_RGBA, GL_UNSIGNED_BYTE, image)

	# Note that only the ID is returned, no reference to the image object or the 
	# string data is stored in user space. 
	# The data is only present within the GL after this call exits.
	return ID

def SetUpTexture():
	"""Render-time texture environment setup"""

	# Configure the texture rendering parameters
	glEnable(GL_TEXTURE_2D)

	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST)
	glTexParameterf(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST)

	glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_DECAL)
	# Re-select our texture, could use other generated textures if we had generated them earlier...
	glBindTexture(GL_TEXTURE_2D, imageID)

# retorna o array de polígonos que forma o cubo, com base nas variáveis globais acima: vertices e surfaces
def BuildSolidStructure(vertices, surfaces):
	i = 0;
	polygons = [];
	for listOfVertexIndexes in surfaces:
		arrayOfPoints = [];
		for vertexIndex in listOfVertexIndexes:
			arrayOfPoints.append(vertices[vertexIndex]);
		polygons.append(Polygon(arrayOfPoints, colors[i % len(colors)]));
		i += 1;

	return polygons;

def Upon_Drag (cursor_x, cursor_y):
	""" Mouse cursor is moving
		Glut calls this function (when mouse button is down)
		and pases the mouse cursor postion in window coords as the mouse moves.
	"""
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot

	if (g_isDragging):
		mouse_pt = Point2fT (cursor_x, cursor_y)
		ThisQuat = g_ArcBall.drag (mouse_pt)						# // Update End Vector And Get Rotation As Quaternion
		g_ThisRot = Matrix3fSetRotationFromQuat4f (ThisQuat)		# // Convert Quaternion Into Matrix3fT
		# Use correct Linear Algebra matrix multiplication C = A * B
		g_ThisRot = Matrix3fMulMatrix3f (g_LastRot, g_ThisRot)		# // Accumulate Last Rotation Into This One
		g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot)	# // Set Our Final Transform's Rotation From This One
	return

def Upon_Click (button, button_state, cursor_x, cursor_y):
	""" Mouse button clicked.
		Glut calls this function when a mouse button is
		clicked or released.
	"""
	global g_isDragging, g_LastRot, g_Transform, g_ThisRot

	g_isDragging = False
	if (button == GLUT_RIGHT_BUTTON and button_state == GLUT_UP):
		# Right button click
		g_LastRot = Matrix3fSetIdentity ();							# // Reset Rotation
		g_ThisRot = Matrix3fSetIdentity ();							# // Reset Rotation
		g_Transform = Matrix4fSetRotationFromMatrix3f (g_Transform, g_ThisRot);	# // Reset Rotation
	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_UP):
		# Left button released
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One
	elif (button == GLUT_LEFT_BUTTON and button_state == GLUT_DOWN):
		# Left button clicked down
		g_LastRot = copy.copy (g_ThisRot);							# // Set Last Static Rotation To Last Dynamic One
		g_isDragging = True											# // Prepare For Dragging
		mouse_pt = Point2fT (cursor_x, cursor_y)
		g_ArcBall.click (mouse_pt);								# // Update Start Vector And Prepare For Dragging
		x, y = ScreenToOGLCoords(cursor_x, cursor_y);
		pickedFace = PickSurface(x, y, solidFaces);

		if pickedFace != None:
			visited = bfs_keeping_track_of_parents(graph, pickedFace);

	return

def bfs_keeping_track_of_parents(graph, start):
	parent = {};
	queue = [start];
	visited = set();
	parent[start] = start;
	while queue:
		node = queue.pop(0);
		if node not in visited:
			visited.add(node);

			for adjacent in graph.vertex_neighbours(node): # <<<<< record its parent
				parent[adjacent] = node

			Visit(node, parent[node]);
			queue.extend(set(graph.vertex_neighbours(node)) - visited);

def Visit(thisPoly, prevPoly):
	if thisPoly != prevPoly:
		n1 = [thisPoly.normal[0], thisPoly.normal[1], thisPoly.normal[2]];
		n2 = [prevPoly.normal[0], prevPoly.normal[1], prevPoly.normal[2]];
		angle = Angle(n1, n2);
		angle = numpy.rad2deg(angle);
		edges1 = GeneratePolygonEdges(thisPoly);
		edges2 = GeneratePolygonEdges(prevPoly);
		commonEdges = GetEdgesInCommon(thisPoly, prevPoly);
		if commonEdges == []:
			return
		fixedPoint = commonEdges[0].midpoint();
		axis = numpy.cross(n1, n2);
		thisPoly.transAngle = 1;
		thisPoly.transPoint = fixedPoint;
		thisPoly.transAxis = axis;
		thisPoly.transMaxAngle = angle;
		thisPoly.transParent = prevPoly;
		trans = translateAndRotate(thisPoly.transAngle, thisPoly.transPoint, thisPoly.transAxis, thisPoly.transParent);
		thisPoly.transform = trans;
	else:
		thisPoly.transParent = None;
	return

def ResetColors(polygons):
	i = 0;
	for poly in polygons:
		poly.color = colors[i % len(colors)];
		i += 1;
	return polygons;

# Recebe o ponto onde o mouse clicou, e um array de polígonos, e retorna o polígono que foi clicado
def PickSurface(mouseX, mouseY, polygons):
	global objectsToDraw;

	# passar os pontos mouseX e mouseY para a rotação certa
	p1 = [mouseX, mouseY, 2];
	p2 = [mouseX, mouseY, -2];
	p1Certo = Matrix3fMulMatrix3f(g_ThisRot, p1);
	p2Certo = Matrix3fMulMatrix3f(g_ThisRot, p2);
	pontoOrigem = Point(p1Certo[0], p1Certo[1], p1Certo[2]);
	pontoDestino = Point(p2Certo[0], p2Certo[1], p2Certo[2]);
	line = Line(pontoOrigem, pontoDestino);

	interceptedPolygons = WhichPolygonsDoesLineCross(line, polygons);
	pickedPolygon = WhichPolygonIsCloserToScreen(interceptedPolygons, line);

	if pickedPolygon != None:
		objectsToDraw[0] = ResetColors(objectsToDraw[0]);
		objectsToDraw[0][objectsToDraw[0].index(pickedPolygon)].color = [0.502, 0.502, 0.502] # gray
	return pickedPolygon;

# Recebe uma linha e um array de polígonos, e retorna quais polígonos são interceptados pela linha
def WhichPolygonsDoesLineCross(line, polygons):
	polyArray = [];
	for poly in polygons:
		contains, p, t = poly.doesLineCrossPolygon(line);
		if (contains):
			polyArray.append(poly);
	return polyArray;

# Recebe um array de polígonos, e retorna aquele que estiver mais próximo da tela
def WhichPolygonIsCloserToScreen(polygons, clickLine):
	closerPoly = None;
	distanceToOrigin = 2;

	for poly in polygons:
		intersection = clickLine.intersectToPlane(poly);

		if intersection[1] < distanceToOrigin:
			distanceToOrigin = intersection[1];
			closerPoly = poly;

	return closerPoly;

def ScreenToOGLCoords(cursor_x, cursor_y):
	viewport = glGetDoublev(GL_VIEWPORT);

	cursor_x = float (cursor_x);
	cursor_y = float (viewport[3]) - float (cursor_y);

	cursor_z = glReadPixels(cursor_x, int(cursor_y), 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT);

	posX, posY, posZ = gluUnProject(cursor_x, cursor_y, cursor_z, None, None, None);

	return posX, posY;

def translateAndRotate(ang, p, axis, parent):
	glPushMatrix();
	glLoadIdentity();
	glMultMatrixf(parent.transform);
	glTranslate(p[0],p[1],p[2]);
	glRotate(ang, axis[0], axis[1], axis[2]);
	glTranslate(-p[0],-p[1],-p[2]);
	T = glGetDoublev ( GL_MODELVIEW_MATRIX );
	glPopMatrix();
	return T;

def GeneratePolygonEdges(polygon):
	edges = [];
	for i in range(len(polygon.points) - 1):
		edges.append(Line(polygon.points[i], polygon.points[i+1]));
	edges.append(Line(polygon.points[-1], polygon.points[0]));
	return edges;

# modelview x g_transform x translaterotate
def BuildBoundingBox():
	global solidFaces;
	# Adiciona todos os pontos transformados numa bounding box
	box = Box();
	for polygon in solidFaces:
		for point in polygon.points:
			auxPoint = [point.x, point.y, point.z, 1];
			transformedPoint = 0;
			if polygon.transParent != None:
				# translaterotate = translateAndRotate(polygon.transMaxAngle, polygon.transPoint, polygon.transAxis, polygon.transParent);
				# transform = Matrix3fMulMatrix3f(auxPoint, modelview);
				auxPointR = Matrix3fMulMatrix3f(auxPoint, polygon.transform); # Aplica a transformação de cada polígono a cada um de seus pontos, para que eu possa trabalhar com o sólido aberto
				transformedPoint = Point(auxPointR[0], auxPointR[1], auxPointR[2]);
				# print GetColor(solidFaces[polyIndex].color), auxPoint, " ~~> ", auxPointR
			else:
				transformedPoint = point;
				# print GetColor(solidFaces[polyIndex].color), auxPoint, " ~~> ", auxPoint
			box.add(transformedPoint);
		# print ""

	return box

def MakeTextureCoords(point, polygon, box, translaterotate):
	transformedPoint = 0;
	if polygon.transParent == None:
		transformedPoint = point;
	else:
		auxPoint = [point[0], point[1], point[2], 1];
		transformedPoint = Matrix3fMulMatrix3f(auxPoint, translaterotate);
		transformedPoint = Point(transformedPoint[0], transformedPoint[1], transformedPoint[2]);
	return box.normalize(transformedPoint);

debugCount = 0;
pointArray = []
texArray = []
arrayIndex = 0;

def DrawSolid():
	global boundingBox, solidFaces, angularSpeed, hasTexture, debugCount, pointArray, texArray, arrayIndex
	i = 0;
	
	for polygon in solidFaces:
		
		translaterotate = 0;
		if hasTexture and polygon.transParent != None:
			translaterotate = translateAndRotate(polygon.transMaxAngle, polygon.transPoint, polygon.transAxis, polygon.transParent);

		# ANIMAÇÃO:
		if angularSpeed != 0:
			if polygon.transMaxAngle != None:
				polygon.transform = translateAndRotate(polygon.transAngle, polygon.transPoint, polygon.transAxis, polygon.transParent);
				polygon.transAngle += angularSpeed;
				if polygon.transAngle >= polygon.transMaxAngle and angularSpeed > 0:
					angularSpeed *= -1;
					if not hasTexture:
						hasTexture = True;
						boundingBox = BuildBoundingBox();
						# Pausa a animação
						# angularSpeed = 0;

						for polygon2 in solidFaces:
							translaterotate = 0;
							if polygon2.transParent != None:
								translaterotate = translateAndRotate(polygon2.transMaxAngle, polygon2.transPoint, polygon2.transAxis, polygon2.transParent);
							for point in polygon2.points:
								pointArray.append([GetColor(polygon2.color), point]);
								texArray.append(MakeTextureCoords(point, polygon2, boundingBox, translaterotate));
								# print GetColor(polygon2.color), point, arrayIndex, " --> ", texArray[arrayIndex]
								arrayIndex += 1;
						return;

				elif polygon.transAngle <= 0 and angularSpeed < 0:
					angularSpeed *= -1;


		glMultMatrixf(polygon.transform);


		glBegin(GL_POLYGON);
		if polygon.color == None:
			glColor3fv(colors[i % len(colors)]);
			i += 1;
		else:
			glColor3fv(polygon.color);
		for point in polygon.points:
			if hasTexture:
				index = pointArray.index([GetColor(polygon.color), point])
				glTexCoord2f(texArray[index][0], texArray[index][1]);

				# if debugCount < 24:
					# print GetColor(polygon.color), point, index, "-->", texArray[index]
					# debugCount += 1

			glVertex3f(point[0], point[1], point[2]);
		glEnd();

		glLoadIdentity();
		glTranslatef(0.0,0.0,-6.0);
		glMultMatrixf(g_Transform);

	# for polygon in solidFaces:
	# 	glBegin(GL_LINES);
	# 	glColor3f(1,1,1);

	# 	for i in range(len(polygon.points) - 1):
	# 		glVertex3f(polygon.points[i][0], polygon.points[i][1], polygon.points[i][2]);
	# 		glVertex3f(polygon.points[i+1][0], polygon.points[i+1][1], polygon.points[i+1][2]);

	# 	glVertex3f(polygon.points[-1][0], polygon.points[-1][1], polygon.points[-1][2]);
	# 	glVertex3f(polygon.points[0][0], polygon.points[0][1], polygon.points[0][2]);
	# 	glEnd();

	return

def Draw ():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
	glLoadIdentity();												# // Reset The Current Modelview Matrix
	glTranslatef(0.0,0.0,-6.0);										# // Move  Into The Screen 6.0

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform

	if hasTexture:
		SetUpTexture();

	# n = len(objectsToDraw);
	# for i in range(n):
	# 	DrawSolid(objectsToDraw[i]);
	DrawSolid();

	# drawCube(); # DEBUG

	glPopMatrix();													# // NEW: Unapply Dynamic Transform

	glFlush ();														# // Flush The GL Rendering Pipeline
	glutSwapBuffers()
	return

# DEBUG
def drawCube():
		"""Draw a cube with texture coordinates"""

		SetUpTexture()

		glBegin(GL_QUADS)
		glColor3f(1.0,0.0,0.0)
		glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f( 1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f( 1.0,  1.0,  1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f(-1.0,  1.0,  1.0);
		glColor3f(0.0,1.0,0.0)
		glTexCoord2f(1.0, 0.0); glVertex3f(-1.0, -1.0, -1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f(-1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f( 1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f( 1.0, -1.0, -1.0);
		glColor3f(0.0,0.0,1.0)
		glTexCoord2f(0.0, 1.0); glVertex3f(-1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f(-1.0,  1.0,  1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f( 1.0,  1.0,  1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f( 1.0,  1.0, -1.0);
		glColor3f(1.0,1.0,0.0)
		glTexCoord2f(1.0, 1.0); glVertex3f(-1.0, -1.0, -1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f( 1.0, -1.0, -1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f( 1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f(-1.0, -1.0,  1.0);
		glColor3f(0.0,1.0,1.0)
		glTexCoord2f(1.0, 0.0); glVertex3f( 1.0, -1.0, -1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f( 1.0,  1.0, -1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f( 1.0,  1.0,  1.0);
		glTexCoord2f(0.0, 0.0); glVertex3f( 1.0, -1.0,  1.0);
		glColor3f(1.0,0.0,1.0)
		glTexCoord2f(0.0, 0.0); glVertex3f(-1.0, -1.0, -1.0);
		glTexCoord2f(1.0, 0.0); glVertex3f(-1.0, -1.0,  1.0);
		glTexCoord2f(1.0, 1.0); glVertex3f(-1.0,  1.0,  1.0);
		glTexCoord2f(0.0, 1.0); glVertex3f(-1.0,  1.0, -1.0);
		glEnd()

		glDisable(GL_TEXTURE_2D)

def GetColor(color):
	return color_dict[tuple(color)];

def Length(v):
	return math.sqrt(Vector3fDot(v, v))

def Angle(v1, v2):
	return math.acos(Vector3fDot(v1, v2) / (Length(v1) * Length(v2)))