<Please submit this file with your solution.>

CSCI 520, Assignment 1

Isabella Rocha
================

<Description of what you have accomplished>

1. Modeled a jello cube as a mass-spring system
- Created a cube of 512 discrete mass points and connected using structural, shear, and bend springs
- These 3 sets of springs the jello cube to stretch, contract, oscillate, change velocity, and bounce off of walls based on the laws of physics

2. Correctly respond to forces
- Elastic and damping forces are applied to all springs based on its neighbors and type of spring
- Deforms based on where it hits the bounding box walls, abides collision penalty laws
  - A collision spring is generated at the point the jello cube hits the wall
- Deforms based on force fields imposed on the jello cube
