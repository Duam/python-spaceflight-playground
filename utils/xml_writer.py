#!/usr/bin/python3

##
# @file xml_writer.py
# @author Paul Daum
##

import xml.etree.cElementTree as etree


##
# @brief Writes data to XML. Contains:
#        - Time horizon
#        - Number of samples
#        - Orbital elements of target orbit
#        - Position and velocity trajectory in cartesian coordinates
#        - Mass of the spacecraft
#        - Forces on the spacecraft
#        It's currently very limited and should be extended later
def write_to_xml(
    filename,
    T, N,
    e, h,
    xPoses, yPoses,
    xVelos, yVelos,
    masses,
    xForces, yForces
):

    # Check if the sizes are correct
    if (
        len(xPoses) != N+1 or
        len(yPoses) != N+1 or
        len(xVelos) != N+1 or
        len(yVelos) != N+1 or
        len(masses) != N+1 or
        len(xForces) != N or
        len(yForces) != N
    ):
        print("One or more sizes are incorrect!")

    # Add a root
    root = etree.Element("Trajectory_Data")

    # Parameters
    params = etree.SubElement(root, "Parameters")
    etree.SubElement(params, "Time_Horizon", name="T").text = str(T)
    etree.SubElement(params, "Samples_Number", name="N").text = str(N)

    # Target orbit
    target_orbit = etree.SubElement(root, "Target_Orbit")
    etree.SubElement(target_orbit, "Eccentricity", name="e_x").text = str(e[0])
    etree.SubElement(target_orbit, "Eccentricity", name="e_y").text = str(e[1])
    etree.SubElement(target_orbit, "Angular_Momentum", name="h").text = str(h)
    
    # Initial state
    x0 = etree.SubElement(root, "x0")
    etree.SubElement(x0, "xPos", name="xPos").text = str(xPoses[0])
    etree.SubElement(x0, "yPos", name="yPos").text = str(yPoses[0])
    etree.SubElement(x0, "xVel", name="xVel").text = str(xVelos[0])
    etree.SubElement(x0, "yVel", name="yVel").text = str(yVelos[0])

    # States
    xs = etree.SubElement(root, "xs")
    for k in range(1,N+1):
        etree.SubElement(xs, "x"+str(k), name="xPos").text = str(xPoses[k])
        etree.SubElement(xs, "x"+str(k), name="yPos").text = str(yPoses[k])
        etree.SubElement(xs, "x"+str(k), name="xVel").text = str(xVelos[k])
        etree.SubElement(xs, "x"+str(k), name="yVel").text = str(yVelos[k])
        etree.SubElement(xs, "x"+str(k), name="mass").text = str(masses[k])

    # Controls
    us = etree.SubElement(root, "us")
    for k in range(N):
        etree.SubElement(us, "u"+str(k), name="xFor").text = str(xForces[k])
        etree.SubElement(us, "u"+str(k), name="yFor").text = str(yForces[k])

    # Create tree and write to file
    tree = etree.ElementTree(root)
    tree.write(filename)


##
# Execute this script to test the xml writer
##
if __name__ == '__main__':
    
    write_to_xml(
        filename = 'test.xml', 
        T = 2,
        N = 1,
        e = [0.5, 0.4],
        h = 9,
        xPoses = [1,1],
        yPoses = [1,2],
        xVelos = [1,3],
        yVelos = [1,4],
        masses = [0,5],
        xForces = [0.2],
        yForces = [0.1]
    )