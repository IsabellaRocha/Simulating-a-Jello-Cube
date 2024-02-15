/*

  USC/Viterbi/Computer Science
  "Jello Cube" Assignment 1 starter code

*/

#include "jello.h"
#include "physics.h"
#include <vector>

point hookForce(point* L, double restLength, double kh) {
	double lenOfVector = 0;
	lenOfVector = pLENGTH(*L);
	if (lenOfVector == 0) 
		return { 0, 0, 0 };
	//Apply hooks law: fh = -kh(|L|-R)(L/|L|), L is vector pointing from B to A, R is spring at rest length
	point hookForce = { 0, 0, 0 };
	pMULTIPLY(*L, -kh * (lenOfVector - restLength) / lenOfVector, hookForce);
	return hookForce;
}

point dampForce(point* L, point* velocity, double kd) {
	double lenOfVector = 0;
	lenOfVector = pLENGTH(*L);
	if (lenOfVector == 0) 
		return { 0, 0, 0 };
	//Apply damping in 3D: -kd * (((vA-vB) dot L)/|L|) * (L/|L|)
	point dampingForce = { 0, 0, 0 };
	pMULTIPLY(*L, -kd * (pDOTPRODUCT(*velocity, *L)) / lenOfVector / lenOfVector, dampingForce);
	return dampingForce;
}

/*
	According to Wikipedia
	Interpolate along x:
	c00 = c000(1-xd) + c100*xd
	c01 = c001(1-xd) + c101*xd
	c10 = c010(1-xd) + c110*xd
	c11 = c011(1-xd) + c111*xd

	c0 = c00(1-yd) + c10*yd
	c1 = c01(1-yd) + c11*yd

	c = c0(1-zd) + c1zd
*/
void interpolate3D(point cubeCoors[8], double x, double y, double z, point* forceFieldForce) {
	// Interpolate along x-axis between cubeCoors[0] and cubeCoors[1] based on x
	point c00;
	c00.x = cubeCoors[0].x * (1 - x) + cubeCoors[1].x * x;
	c00.y = cubeCoors[0].y * (1 - x) + cubeCoors[1].y * x;
	c00.z = cubeCoors[0].z * (1 - x) + cubeCoors[1].z * x;

	// Interpolate along x-axis between cubeCoors[2] and cubeCoors[3] based on x
	point c01;
	c01.x = cubeCoors[2].x * (1 - x) + cubeCoors[3].x * x;
	c01.y = cubeCoors[2].y * (1 - x) + cubeCoors[3].y * x;
	c01.z = cubeCoors[2].z * (1 - x) + cubeCoors[3].z * x;

	// Interpolate along x-axis between cubeCoors[4] and cubeCoors[5] based on x
	point c10;
	c10.x = cubeCoors[4].x * (1 - x) + cubeCoors[5].x * x;
	c10.y = cubeCoors[4].y * (1 - x) + cubeCoors[5].y * x;
	c10.z = cubeCoors[4].z * (1 - x) + cubeCoors[5].z * x;

	// Interpolate along x-axis between cubeCoors[6] and cubeCoors[7] based on x
	point c11;
	c11.x = cubeCoors[6].x * (1 - x) + cubeCoors[7].x * x;
	c11.y = cubeCoors[6].y * (1 - x) + cubeCoors[7].y * x;
	c11.z = cubeCoors[6].z * (1 - x) + cubeCoors[7].z * x;

	// Interpolate along y-axis between c00 and c01 based on y
	point c0;
	c0.x = c00.x * (1 - y) + c01.x * y;
	c0.y = c00.y * (1 - y) + c01.y * y;
	c0.z = c00.z * (1 - y) + c01.z * y;

	// Interpolate along y-axis between c10 and c11 based on y
	point c1;
	c1.x = c10.x * (1 - y) + c11.x * y;
	c1.y = c10.y * (1 - y) + c11.y * y;
	c1.z = c10.z * (1 - y) + c11.z * y;

	// Interpolate along z-axis between c0 and c1 based on z
	forceFieldForce->x = c0.x * (1 - z) + c1.x * z;
	forceFieldForce->y = c0.y * (1 - z) + c1.y * z;
	forceFieldForce->z = c0.z * (1 - z) + c1.z * z;
}

void forceFieldForce(struct world* jello, struct point currPoint, point* forceFieldForce) {
	struct point cell;
	struct point cubeCoors[8];

	//Multiply everything by 1/h, and then take the floor, this gets you what cell you are currently in
	const double hinv = 1.0 / 4.0;
    pMULTIPLY(currPoint, hinv, currPoint);
	//Correct if out of bounds (caused error when didn't do this)
	if (currPoint.x > 1) currPoint.x = 1;
	else if (currPoint.x < 0) currPoint.x = 0;
	if (currPoint.y > 1) currPoint.y = 1;
	else if (currPoint.y < 0) currPoint.y = 0;
	if (currPoint.z > 1) currPoint.z = 1;
	else if (currPoint.z < 0) currPoint.z = 0;

	//Identify which cell will be the bottom left of the cube you're currently in
    cell.x = floor(currPoint.x);
    cell.y = floor(currPoint.y);
    cell.z = floor(currPoint.z);

    // Get the 8 mass points that form the cube with the first point being the bottom left
	cubeCoors[0] = jello->forceField[int(cell.x * jello->resolution * jello->resolution + cell.y * jello->resolution + cell.z * jello->resolution)];
	cubeCoors[1] = jello->forceField[int((cell.x+1) * jello->resolution * jello->resolution + cell.y * jello->resolution + cell.z * jello->resolution)];
	cubeCoors[2] = jello->forceField[int((cell.x+1) * jello->resolution * jello->resolution + (cell.y+1) * jello->resolution + cell.z * jello->resolution)];
	cubeCoors[3] = jello->forceField[int((cell.x+1) * jello->resolution * jello->resolution + (cell.y+1) * jello->resolution + (cell.z+1) * jello->resolution)];
	cubeCoors[4] = jello->forceField[int(cell.x * jello->resolution * jello->resolution + (cell.y+1) * jello->resolution + cell.z * jello->resolution)];
	cubeCoors[5] = jello->forceField[int(cell.x * jello->resolution * jello->resolution + (cell.y+1) * jello->resolution + (cell.z+1) * jello->resolution)];
	cubeCoors[6] = jello->forceField[int(cell.x * jello->resolution * jello->resolution + cell.y * jello->resolution + (cell.z+1) * jello->resolution)];
	cubeCoors[7] = jello->forceField[int((cell.x+1) * jello->resolution * jello->resolution + cell.y * jello->resolution + (cell.z+1) * jello->resolution)];

	/*
	* Do pDifference because
	* xd = (x-x0)/(x1-x0)
	* yd = (y-y0)/(y1-y0)
	* zd = (z-z0)/(z1-z0)
	*/
	pDIFFERENCE(currPoint, cell, currPoint);
	interpolate3D(cubeCoors, currPoint.x, currPoint.y, currPoint.z, forceFieldForce);
}


// Function to check for collision and calculate forces
void checkForCollision(int direction, double boxCoor, double jelloCoor, double* lengthPoint, double kh, double kd, point* lengthVector, point* velocity, point* totalForce) {
	memset(lengthVector, 0.0, sizeof(point));
	*lengthPoint = jelloCoor - boxCoor;

	// Check if the difference between the wall and the object is less than 0, indicating a collision
	if (*lengthPoint * direction < 0) {
		point fh = hookForce(lengthVector, 0, kh);
		point fd = dampForce(lengthVector, velocity, kd);
		point hookAndDamp = { 0, 0, 0 };
		pSUM(fh, fd, hookAndDamp);
		pSUM(*totalForce, hookAndDamp, *totalForce); // Add up forces for all coordinates
	}
}


void collisionForce(struct point* p, point velocity, double kh, double kd, point* totalForce) {
	//Bounding box is set to -2 to 2
	const double boxCoor[] = { -2, 2 };
	//Check both directions
	const int directions[] = { 1, -1 };
	point L = { 0, 0, 0 };
	// Make sure to get rid of previous forces
	memset(totalForce, 0, sizeof(point)); 

	// Check for collisions in all directions for x, y, and z coordinates
	for (int i = 0; i < 2; ++i) {
		checkForCollision(directions[i], boxCoor[i], p->x, &L.x, kh, kd, &L, &velocity, totalForce);
		checkForCollision(directions[i], boxCoor[i], p->y, &L.y, kh, kd, &L, &velocity, totalForce);
		checkForCollision(directions[i], boxCoor[i], p->z, &L.z, kh, kd, &L, &velocity, totalForce);
	}
}

//Calculate the forces for structural
void calculateStructuralForce(struct world* jello, int x, int y, int z, point* velocityVector, point* L, point* structuralSpringForce, point* finalForce) {
	double restLength = 0;
	double restLengthCoefficient = 1.0 / 7.0;

	int neighbors[6][3] = { {1, 0, 0}, {-1, 0, 0}, {0, 1, 0}, {0, -1, 0}, {0, 0, 1}, {0, 0, -1} };
	for (int idx = 0; idx < 6; idx++) {
		int i = neighbors[idx][0];
		int j = neighbors[idx][1];
		int k = neighbors[idx][2];
		// Calculate rest length
		restLength = sqrt(i * i + j * j + k * k);
		// Check bounds
		if (x + i >= 0 && x + i <= 7 && y + j >= 0 && y + j <= 7 && z + k >= 0 && z + k <= 7) {
			// Calculate length and velocity of the spring
			pDIFFERENCE(jello->p[x][y][z], jello->p[x + i][y + j][z + k], *L);
			pDIFFERENCE(jello->v[x][y][z], jello->v[x + i][y + j][z + k], *velocityVector);
			// Calculate forces
			struct point fh = hookForce(L, restLengthCoefficient * restLength, jello->kElastic);
			struct point fd = dampForce(L, velocityVector, jello->dElastic);
			pSUM(fh, fd, *finalForce);
			// Add to structural spring force
			pSUM(*structuralSpringForce, *finalForce, *structuralSpringForce);
		}
	}
}

//Calculate the forces for structural and shear springs
void calculateShearForce(struct world* jello, int x, int y, int z, point* velocityVector, point* L, point* shearSpringForce, point* finalForce) {
	double restLength = 0;
	double restLengthCoefficient = 1.0 / 7.0;
	int neighbors[20][3] = { {1, 1, 1}, {-1, 1, 1}, {-1, -1, 1}, {1, -1, 1}, {1, 1, -1}, {-1, 1, -1}, {-1, -1, -1}, {1, -1, -1}, 
		{1, 1, 0}, {-1, 1, 0}, {-1, -1, 0}, {1, -1, 0}, {0, 1, 1}, {0, -1, 1}, {0, -1, -1}, {0, 1, -1}, {1, 0, 1}, {-1, 0, 1}, {-1, 0, -1}, {1, 0, -1} };
	// For loop to look at neighboring springs
	for (int idx = 0; idx < 20; idx++) {
		int i = neighbors[idx][0];
		int j = neighbors[idx][1];
		int k = neighbors[idx][2];
		// Calculate rest length
		restLength = sqrt(i * i + j * j + k * k);
		// Check bounds
		if (x + i >= 0 && x + i <= 7 && y + j >= 0 && y + j <= 7 && z + k >= 0 && z + k <= 7) {
			// Calculate length and velocity of the spring from its neighboring springs
			pDIFFERENCE(jello->p[x][y][z], jello->p[x + i][y + j][z + k], *L);
			pDIFFERENCE(jello->v[x][y][z], jello->v[x + i][y + j][z + k], *velocityVector);
			// Calculate forces
			struct point fh = hookForce(L, restLengthCoefficient * restLength, jello->kElastic);
			struct point fd = dampForce(L, velocityVector, jello->dElastic);
			pSUM(fh, fd, *finalForce);
			// Add to shear spring force
			pSUM(*shearSpringForce, *finalForce, *shearSpringForce);
		}
	}
}

//Calculate the forces for bend springs
void calculateBendForces(struct world* jello, int x, int y, int z, point* velocityVector, point* L, point* bendSpringForce, point* finalForce) {
	//Rest length is 2 since bend springs are connected to their second neighbor
	double restLength = 2;
	double restLengthCoefficient = 1.0 / 7.0;
	int neighbors[6][3] = { {2, 0, 0}, {0, 2, 0}, {0, 0, 2}, {-2, 0, 0}, {0, -2, 0}, {0, 0, -2} };
	for (int idx = 0; idx < 6; idx++) {
		int i = neighbors[idx][0];
		int j = neighbors[idx][1];
		int k = neighbors[idx][2];
		// Check bounds
		if (x + i >= 0 && x + i <= 7  && y + j >= 0 && y + j <= 7 && z + k >= 0 && z + k <= 7) {
			// Calculate length and velocity of the spring from its neighboring springs
			pDIFFERENCE(jello->v[x][y][z], jello->v[x + i][y + j][z + k], *velocityVector);
			pDIFFERENCE(jello->p[x][y][z], jello->p[x + i][y + j][z + k], *L);
			// Calculate forces
			struct point fh = hookForce(L, restLengthCoefficient * restLength, jello->kElastic);
			struct point fd = dampForce(L, velocityVector, jello->dElastic);
			pSUM(fh, fd, *finalForce);
			// Add to bend spring force
			pSUM(*bendSpringForce, *finalForce, *bendSpringForce);
		}
	}
}

void computeAcceleration(struct world* jello, struct point a[8][8][8]) {
    struct point L = { 0, 0, 0 };
    struct point velocityVector = { 0, 0, 0 };
    struct point finalForce = { 0, 0, 0 };

	double massInv = 1.0 / jello->mass;

    for (int x = 0; x < 8; x++) {
        for (int y = 0; y < 8; y++) {
            for (int z = 0; z < 8; z++) {
                struct point* acceleration = &(a[x][y][z]);
                memset(acceleration, 0, sizeof(struct point));

                struct point position = jello->p[x][y][z];
                struct point velocity = jello->v[x][y][z];

                // Set up forces for all three types of springs
                struct point structuralSpringForce = { 0, 0, 0 };
                struct point shearSpringForce = { 0, 0, 0 };
                struct point bendSpringForce = { 0, 0, 0 };

				calculateStructuralForce(jello, x, y, z, &velocityVector, &L, &structuralSpringForce, &finalForce);
				calculateShearForce(jello, x, y, z, &velocityVector, &L, &shearSpringForce, &finalForce);
				calculateBendForces(jello, x, y, z, &velocityVector, &L, &bendSpringForce, &finalForce);
               
                // Add structural, shear, and bend force to accumulated force
                pSUM(*acceleration, structuralSpringForce, *acceleration);
                pSUM(*acceleration, shearSpringForce, *acceleration);
                pSUM(*acceleration, bendSpringForce, *acceleration);

                // Calculate the penalty force from collision
                collisionForce(&position, velocity, jello->kCollision, jello->dCollision, &finalForce);
                pSUM(*acceleration, finalForce, *acceleration);

                // Calculate the forcefield force
                forceFieldForce(jello, position, &finalForce);
                pSUM(*acceleration, finalForce, *acceleration);

                pMULTIPLY(*acceleration, massInv, *acceleration);
            }
        }
    }
}



/* performs one step of Euler Integration */
/* as a result, updates the jello structure */
void Euler(struct world* jello)
{
	int i, j, k;
	point a[8][8][8];

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				jello->p[i][j][k].x += jello->dt * jello->v[i][j][k].x;
				jello->p[i][j][k].y += jello->dt * jello->v[i][j][k].y;
				jello->p[i][j][k].z += jello->dt * jello->v[i][j][k].z;
				jello->v[i][j][k].x += jello->dt * a[i][j][k].x;
				jello->v[i][j][k].y += jello->dt * a[i][j][k].y;
				jello->v[i][j][k].z += jello->dt * a[i][j][k].z;

			}
}

/* performs one step of RK4 Integration */
/* as a result, updates the jello structure */
void RK4(struct world* jello)
{
	point F1p[8][8][8], F1v[8][8][8],
		F2p[8][8][8], F2v[8][8][8],
		F3p[8][8][8], F3v[8][8][8],
		F4p[8][8][8], F4v[8][8][8];

	point a[8][8][8];


	struct world buffer;

	int i, j, k;

	buffer = *jello; // make a copy of jello

	computeAcceleration(jello, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				pMULTIPLY(jello->v[i][j][k], jello->dt, F1p[i][j][k]);
				pMULTIPLY(a[i][j][k], jello->dt, F1v[i][j][k]);
				pMULTIPLY(F1p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F1v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F2p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F2p[i][j][k]);
				// F2v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F2v[i][j][k]);
				pMULTIPLY(F2p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F2v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);

	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F3p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F3v[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 0.5, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 0.5, buffer.v[i][j][k]);
				pSUM(jello->p[i][j][k], buffer.p[i][j][k], buffer.p[i][j][k]);
				pSUM(jello->v[i][j][k], buffer.v[i][j][k], buffer.v[i][j][k]);
			}

	computeAcceleration(&buffer, a);


	for (i = 0; i <= 7; i++)
		for (j = 0; j <= 7; j++)
			for (k = 0; k <= 7; k++)
			{
				// F3p = dt * buffer.v;
				pMULTIPLY(buffer.v[i][j][k], jello->dt, F4p[i][j][k]);
				// F3v = dt * a(buffer.p,buffer.v);     
				pMULTIPLY(a[i][j][k], jello->dt, F4v[i][j][k]);

				pMULTIPLY(F2p[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3p[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1p[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4p[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->p[i][j][k], jello->p[i][j][k]);

				pMULTIPLY(F2v[i][j][k], 2, buffer.p[i][j][k]);
				pMULTIPLY(F3v[i][j][k], 2, buffer.v[i][j][k]);
				pSUM(buffer.p[i][j][k], buffer.v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F1v[i][j][k], buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], F4v[i][j][k], buffer.p[i][j][k]);
				pMULTIPLY(buffer.p[i][j][k], 1.0 / 6, buffer.p[i][j][k]);
				pSUM(buffer.p[i][j][k], jello->v[i][j][k], jello->v[i][j][k]);
			}

	return;
}