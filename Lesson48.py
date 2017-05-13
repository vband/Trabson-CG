# -*- coding: utf-8 -*-
from OpenGL.GL import *
from OpenGL.GLUT import *
from OpenGL.GLU import *
import sys
import copy
import math
import numpy

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

solidFaces = [];
objectsToDraw = [];

graph = Graph();

angularSpeed = 0.1;

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

# A general OpenGL initialization function.  Sets all of the initial parameters. 
def Initialize (Width, Height, argv):				# We call this right after our OpenGL window is created.
	global g_quadratic, solidFaces, objectsToDraw, graph;

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

	# for fa in solidFaces:
	# 	print color_dict[tuple(fa.color)], ":"
	# 	print "tem aresta em comum com:"
	# 	for face in solidFaces:
	# 		if fa != face and DoPolygonsHaveAnEdgeInCommon(fa, face):
	# 			print color_dict[tuple(face.color)]
	# 	print ""

	# print graph.is_connected()

	return True

def BuildGraph(polygons):
	graph = Graph();

	for polygon in polygons:
		graph.add_vertex(polygon);

	for polygon1 in polygons:
		# print color_dict[tuple(polygon1.color)], ":"
		# print "tem aresta em comum com:"
		for polygon2 in polygons:
			if polygon1.points != polygon2.points and DoPolygonsHaveAnEdgeInCommon(polygon1, polygon2):
				# print color_dict[tuple(polygon2.color)]
				graph.add_edge({polygon1, polygon2});
				# graph.add_edge({polygon2, polygon1});
		# print ""





	# for polygon in polygons:
	# 	graph.add_vertex(polygon);
	# 	# print color_dict[tuple(polygon.color)]
	
	# for i in range(len(polygons)):
	# 	for j in range(len(polygons)):
	# 		if j != i:
	# 			if DoPolygonsHaveAnEdgeInCommon(polygons[i], polygons[j]):
	# 				# if {polygons[i], polygons[j]} not in graph.edges():
	# 				graph.add_edge({polygons[i], polygons[j]});
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
			visited = dfs_v3(graph, pickedFace);
			# for v in visited:
			# 	print GetColor(v.color)

			# for nei in graph.vertex_neighbours(pickedFace):
			# 	Visit(nei, pickedFace);

	return

def dfs_v3(graph, start):
	visited, stack = set(), [start]
	previous = start;
	while stack:
		vertex = stack.pop();
		if vertex not in visited:
			visited.add(vertex);
			Visit(vertex, previous);
			stack.extend(set(graph.vertex_neighbours(vertex)) - visited);
			# stack.extend(set(graph.vertex_children(vertex)) - visited);
			# stack.extend(set(graph.vertex_parents(vertex)) - visited);
			previous = vertex;
	return visited

def Visit(thisPoly, prevPoly):
	# print "thisPoly: ", GetColor(thisPoly.color);
	# print "prevPoly: ", GetColor(prevPoly.color);
	if thisPoly != prevPoly:
		n1 = [thisPoly.normal[0], thisPoly.normal[1], thisPoly.normal[2]];
		n2 = [prevPoly.normal[0], prevPoly.normal[1], prevPoly.normal[2]];
		angle = Angle(n1, n2);
		angle = numpy.rad2deg(angle);
		# print "angle: ", angle
		edges1 = GeneratePolygonEdges(thisPoly);
		# print "edges do this: ", edges1
		edges2 = GeneratePolygonEdges(prevPoly);
		# print "edges do prev: ", edges2
		commonEdges = GetEdgesInCommon(thisPoly, prevPoly);
		if commonEdges == []:
			return
		# print "commonEdges: ", commonEdges
		fixedPoint = commonEdges[0].midpoint();
		# fixedPoint = [fixedPoint[0], fixedPoint[1], fixedPoint[2]];
		# print "fixedPoint: ", fixedPoint
		# axis = Vector3fCross(n1, n2);
		axis = numpy.cross(n1, n2);
		# print "axis: ", axis
		# fixedPoint = numpy.array(fixedPoint);
		# axis = numpy.array(axis);
		thisPoly.transAngle = 1;
		thisPoly.transPoint = fixedPoint;
		thisPoly.transAxis = axis;
		thisPoly.transMaxAngle = angle;
		thisPoly.ancestralTransform = prevPoly.transform;
		# thisPoly.transform = dot(thisPoly.transform, translateAndRotate(1, fixedPoint, axis));
		trans = translateAndRotate(thisPoly.transMaxAngle, thisPoly.transPoint, thisPoly.transAxis, thisPoly.ancestralTransform);
		thisPoly.transform = trans;

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

# Desenha na tela uma linha entrando na tela, representando graficamente o clique do mouse
# def DrawClickLine(mouseX, mouseY):
# 	# mouseX, mouseY = ScreenToOGLCoords(mouseX, mouseY);
# 	glBegin(GL_LINES);
# 	glColor3f(0,0,1);
# 	glVertex3f(mouseX, mouseY, 2);
# 	glVertex3f(mouseX, mouseY, -2);
# 	glEnd();
# 	return

def translateAndRotate(ang, p, axis, transform):
	# T = translate(p[0],p[1],p[2]) * \
        # rotate(ang, axis[0], axis[1], axis[2]) * \
        # translate(-p[0],-p[1],-p[2])

	glPushMatrix();
	glLoadIdentity();
	glMultMatrixf(transform);
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

def DrawSolid(solidFaces):
	global angularSpeed
	i = 0;
	for polygon in solidFaces:

		# Aplica a transformação
		# if polygon.transform != None:
			# print "vai aplicar a transform do", GetColor(polygon.color)
			# t = polygon.transform;


			# ANIMAÇÃO:
		# if polygon.transMaxAngle != None:
		# 	polygon.transform = translateAndRotate(polygon.transAngle, polygon.transPoint, polygon.transAxis, polygon.ancestralTransform);
		# 	# polygon.transform = rotate(polygon.transAngle, polygon.transAxis[0], polygon.transAxis[1], polygon.transAxis[2]);
		# 	polygon.transAngle += angularSpeed;
		# 	if polygon.transAngle >= polygon.transMaxAngle and angularSpeed > 0:
		# 		angularSpeed *= -1;
		# 	elif polygon.transAngle <= 0 and angularSpeed < 0:
		# 		angularSpeed *= -1;


		# print solidFaces[0].transAngle, solidFaces[0].transMaxAngle
		glMultMatrixf(polygon.transform);

		glBegin(GL_POLYGON);
		if polygon.color == None:
			glColor3fv(colors[i % len(colors)]);
			i += 1;
		else:
			glColor3fv(polygon.color);
		for point in polygon.points:
			glVertex3f(point[0], point[1], point[2]);
		glEnd();

		# Desfaz a transformação
		# if polygon.transform != None:
		glLoadIdentity();
		glTranslatef(0.0,0.0,-6.0);
		glMultMatrixf(g_Transform);

	for polygon in solidFaces:
		glBegin(GL_LINES);
		glColor3f(1,1,1);

		for i in range(len(polygon.points) - 1):
			glVertex3f(polygon.points[i][0], polygon.points[i][1], polygon.points[i][2]);
			glVertex3f(polygon.points[i+1][0], polygon.points[i+1][1], polygon.points[i+1][2]);

		glVertex3f(polygon.points[-1][0], polygon.points[-1][1], polygon.points[-1][2]);
		glVertex3f(polygon.points[0][0], polygon.points[0][1], polygon.points[0][2]);
		glEnd();

	return

def Draw ():
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);				# // Clear Screen And Depth Buffer
	glLoadIdentity();												# // Reset The Current Modelview Matrix
	glTranslatef(0.0,0.0,-6.0);										# // Move  Into The Screen 6.0

	glPushMatrix();													# // NEW: Prepare Dynamic Transform
	glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
	# glColor3f(0.75,0.75,1.0);
	# Torus(0.30,1.00);
	n = len(objectsToDraw);
	for i in range(n):
		DrawSolid(objectsToDraw[i]);

	

	glPopMatrix();													# // NEW: Unapply Dynamic Transform

	# if hasClicked:
	# 	DrawClickLine(lastMouseClickX, lastMouseClickY);
	

	# glLoadIdentity();												# // Reset The Current Modelview Matrix
	# glTranslatef(1.5,0.0,-6.0);										# // Move Right 1.5 Units And Into The Screen 7.0

	# glPushMatrix();													# // NEW: Prepare Dynamic Transform
	# glMultMatrixf(g_Transform);										# // NEW: Apply Dynamic Transform
	# glColor3f(1.0,0.75,0.75);
	# gluSphere(g_quadratic,1.3,20,20);
	# glPopMatrix();													# // NEW: Unapply Dynamic Transform

	glFlush ();														# // Flush The GL Rendering Pipeline
	glutSwapBuffers()
	return

def GetColor(color):
	return color_dict[tuple(color)];

# def Dotproduct(v1, v2):
# 	return sum((a*b) for a, b in zip(v1, v2))

def Length(v):
	return math.sqrt(Vector3fDot(v, v))

def Angle(v1, v2):
	return math.acos(Vector3fDot(v1, v2) / (Length(v1) * Length(v2)))