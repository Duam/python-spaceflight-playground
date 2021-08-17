#!/usr/bin/python3

##
# @file xml_writer.py
# @author Paul Daum
##

import xml.etree.cElementTree as etree
import numpy as np

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
        len(masses) != N+1):
        print("One or more state trajectory lengths are incorrect!")

    if (len(xForces) != N or
        len(yForces) != N):
        print("One or more control trajectory lenghts are incorrect!")

    # Add a root
    root = etree.Element("Trajectory_Data")

    # Parameters
    params = etree.SubElement(root, "Parameters")
    etree.SubElement(params, "T").text = str(T)
    etree.SubElement(params, "N").text = str(N)

    # Target orbit
    target_orbit = etree.SubElement(root, "Target_Orbit")
    etree.SubElement(target_orbit, "e_x").text = str(e[0])
    etree.SubElement(target_orbit, "e_y").text = str(e[1])
    etree.SubElement(target_orbit, "h").text = str(h)
    
    # States
    xs = etree.SubElement(root, "xs")
    for k in range(0,N+1):
        etree.SubElement(xs, "state",
            k=str(k),
            xPos=str(xPoses[k]),
            yPos=str(yPoses[k]),
            xVel=str(xVelos[k]),
            yVel=str(yVelos[k]),
            mass=str(masses[k])
        )

    # Controls
    us = etree.SubElement(root, "us")
    for k in range(N):
        etree.SubElement(us, "control",
            k=str(k),
            xFor=str(xForces[k]),
            yFor=str(yForces[k])
        )

    # Create tree and write to file
    tree = etree.ElementTree(root)
    tree.write(filename)


def read_from_xml(filename):

    # Grab the root of the element tree
    root = etree.parse(filename).getroot()
    
    # Get the parameters, x0, xs, us
    parameters = root.find("Parameters")
    target_orbit = root.find("Target_Orbit")
    xs = root.find("xs")
    us = root.find("us")

    # Grab parameters
    params = {
        'T': float(parameters.find('T').text),
        'N': int(parameters.find('N').text)
    }

    # Grab target orbit
    target_orbit = {
        'e_x' : float(target_orbit.find('e_x').text),
        'e_y' : float(target_orbit.find('e_y').text),
        'h' : float(target_orbit.find('h').text)
    }

    # Grab states
    xs_out = np.zeros((5,params['N']+1)) # MAGICNUMBER
    for x in xs.findall('state'):
        k = int(x.get('k'))
        xPos = float(x.get('xPos'))
        yPos = float(x.get('yPos'))
        xVel = float(x.get('xVel'))
        yVel = float(x.get('yVel'))
        mass = float(x.get('mass'))
        xs_out[:,k] = np.array([xPos, yPos, xVel, yVel, mass])

    # Grab controls
    us_out = np.zeros((2,params['N'])) # MAGICNUMBER
    for u in us.findall('control'):
        k = int(u.get('k'))
        xFor = float(u.get('xFor'))
        yFor = float(u.get('yFor'))
        us_out[:,k] = np.array([xFor, yFor])

    return params, target_orbit, xs_out, us_out


##
# Execute this script to test the xml writer
##
if __name__ == '__main__':
    
    # Write something to xml
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

    # Read something from xml
    params, x0, xs, us = read_from_xml('test.xml')

    print(params)
    print(xs)
    print(us)