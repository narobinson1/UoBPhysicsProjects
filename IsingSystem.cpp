//
//  IsingSystem.cpp
//

#include "IsingSystem.h"

// colors
namespace colours {
	// note the f's here avoid warnings by telling C++ to read the numbers as floats and not doubles 
	GLfloat blue[] = { 0.1f, 0.3f, 0.9f, 1.0f };   // blue
	GLfloat red[] = { 1.0f, 0.2f, 0.1f, 0.2f };   // red
	GLfloat green[] = { 0.3f, 0.6f, 0.3f, 1.0f };   // green
}


// constructor
IsingSystem::IsingSystem(Window *set_win) {
	cout << "creating system, gridSize " << gridSize << endl;
	win = set_win;

	inverseTemperatureBeta = 1/3.0;
	slowNotFast = 1;
	isActive = 0;

	// Allocate memory for the grid, remember to free the memory in destructor
	//   the point here is that each row of the grid is an array
	//   the grid itself is a an array of pointers, one for each row
	// Here we allocate the array of pointers
	grid = new int*[gridSize];
	// Now allocate the indidual rows
	for (int i = 0; i<gridSize; i++) {
		grid[i] = new int[gridSize];
	}

	logfile.open("5Beta0.74.txt");
	// this sets the temperatre and initialises the spins grid
	Reset();
}

void IsingSystem::Reset() {

	double initialTemp = 3.0;

	setTemperature(initialTemp);

	// set the grid to -1
	for (int i = 0; i<gridSize; i++) {
		for (int j = 0; j<gridSize; j++) {
			// position is (i,j)
			int pos[2] = { i,j };
			// set this spin to state -1
			setGrid(pos, -1);
		}
	}
}


// destructor
IsingSystem::~IsingSystem() {
	// Close the file (if open)
	if (logfile.is_open())
		logfile.close();

	// Delete the window
	delete win;

	// Delete the grid
	// First we delete the individual rows
	for (int i = 0; i<gridSize; i++)
		delete[] grid[i];
	// Finally delete the array of pointers
	delete[] grid;
}

// this draws the system
void IsingSystem::DrawSquares() {

	double drawScale = 2.0 / (gridSize * 1.1);

	// draw the particles
	double halfSize = 0.5;
	int halfGrid = gridSize / 2;
	for (int x = 0; x<gridSize; x++) {
		for (int y = 0; y<gridSize; y++) {

			double vec[2];
			vec[0] = x - halfGrid;
			vec[1] = y - halfGrid;

			// openGL magic
			glPushMatrix();
			// choose a color
			if (grid[x][y] == -1)
				glColor4fv(colours::green);
			else
				glColor4fv(colours::blue);
			// draw a rectangle for the particle
			glRectd(drawScale*(vec[0] - halfSize),
				drawScale*(vec[1] - halfSize),
				drawScale*(vec[0] + halfSize),
				drawScale*(vec[1] + halfSize));
			// openGL magic
			glPopMatrix();
		}
	}

	// print some information (at top left)
	// this ostringstream is a way to make a string with numbers and words (similar to cout << ... )
	ostringstream str;
	str << "beta " << inverseTemperatureBeta << " size " << gridSize;
	win->displayString(str, -0.9, 0.94, colours::red);

}


// attempt N spin flips, where N is the number of spins
void IsingSystem::MCsweep() {
	for (int i = 0; i<gridSize*gridSize; i++)
		attemptSpinFlip();
}

// here we attempt to flip a spin and accept/reject with Metropolis rule
void IsingSystem::attemptSpinFlip() {
	int pos[2];

	// random site
	pos[0] = rgen.randomInt(gridSize);
	pos[1] = rgen.randomInt(gridSize);

	double hloc = computeLocalField(pos);
	
	double dE = 2.0 * hloc * readGrid(pos);
	if (dE < 0)
		flipSpin(pos);
	else if (rgen.random01() < exp(-dE))
		flipSpin(pos);
}

int IsingSystem::Magnetisation() {
	int spin_sum = 0;
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			spin_sum += grid[i][j];
		}
	}
	return spin_sum;
}

double IsingSystem::computeLocalField(int pos[]) {
	double result = 0.0;
	for (int i = 0; i < 4; i++) {
		int nborPos[2];
		setPosNeighbour(nborPos, pos, i);
		result += readGrid(nborPos);
	}
	result *= inverseTemperatureBeta;
	return result;
}


double IsingSystem::computeTotalEnergy() {
	float E_total = 0;
	for (int i = 0; i < gridSize; i++) {
		for (int j = 0; j < gridSize; j++) {
			// position is (i,j)
			int pos[2] = { i,j };
			E_total = E_total - computeLocalField(pos)*readGrid(pos);			
		}
	}
	E_total = 0.5 * E_total;

	return E_total;
}

// NOTE: this returns the local field *divided by the temperature* (dimensionless quantity)

// set the value of a grid cell for a particular position
void IsingSystem::setGrid(int pos[], int val) {
	grid[pos[0]][pos[1]] = val;
}

// read the grid cell for a given position
int IsingSystem::readGrid(int pos[]) {
	return grid[pos[0]][pos[1]];
}

// read the grid cell for a given position
void IsingSystem::flipSpin(int pos[]) {
	grid[pos[0]][pos[1]] = -grid[pos[0]][pos[1]];
}


// send back the position of a neighbour of a given grid cell
// NOTE: we take care of periodic boundary conditions, also positions are integers now not doubles
void IsingSystem::setPosNeighbour(int setpos[], int pos[], int val) {
	switch (val) {
	case 0:
		setpos[0] = (pos[0] + 1) % gridSize;
		setpos[1] = pos[1];
		break;
	case 1:
		setpos[0] = (pos[0] - 1 + gridSize) % gridSize;
		setpos[1] = pos[1];
		break;
	case 2:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] + 1) % gridSize;
		break;
	case 3:
		setpos[0] = pos[0];
		setpos[1] = (pos[1] - 1 + gridSize) % gridSize;
		break;
	}
}

// this is the update function which at the moment just does one mc sweep
void IsingSystem::Update() {
	number_of_sweeps += 1;
	if (logfile.is_open()) {
		n += 1;
		//sum_magnetisation = sum_magnetisation + Magnetisation();
		//sum_totalenergy = sum_totalenergy + computeTotalEnergy();
		//logfile << "Magnetisation:" << sum_magnetisation/n << " | Total energy:" << computeTotalEnergy()/n << endl;

		logfile << Magnetisation() << ", " << computeTotalEnergy() << endl;
	}


	// transient stage taken as 1000 sweeps, number of sweeps between measurements m



	if (number_of_sweeps > 10 & number_of_sweeps % 50 == 0) {

		// number of data points used in average

		m += 1;

		sum_magnetisation = sum_magnetisation + Magnetisation();

		// sum_magnetisation_squared used in Var(M)

		sum_magnetisation_squared = sum_magnetisation_squared + (Magnetisation() * Magnetisation());

		sum_totalenergy = sum_totalenergy + computeTotalEnergy();

		// sum_totalenergy_squared used in Var(M)

		sum_totalenergy_squared = sum_totalenergy_squared + (computeTotalEnergy() * computeTotalEnergy());

	}



	if (m == 50) {

		cout << "Average magnetisation:" << sum_magnetisation / m << " |  Average total energy:" << sum_totalenergy / (m * gridSize * gridSize) << "Magnetisation Modulus:" << abs(sum_magnetisation / m) << endl;

		cout << "Specific heat:" << (k_b / (gridSize * gridSize * pow((1/inverseTemperatureBeta),2))) * (sum_totalenergy_squared/m - pow(sum_totalenergy/m, 2)) << endl;

		cout << "Magnetic susceptibility:" << (gridSize * gridSize / (1/inverseTemperatureBeta)) * (sum_magnetisation_squared/m - pow(sum_magnetisation/m, 2)) << endl;

	}

	MCsweep();

}