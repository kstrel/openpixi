package org.openpixi.pixi.physics.fields;

import edu.emory.mathcs.jtransforms.fft.*;

import org.openpixi.pixi.parallel.cellaccess.CellAction;
import org.openpixi.pixi.physics.grid.Cell;
import org.openpixi.pixi.physics.grid.Grid;

public class PoissonSolverFFTPeriodic extends FieldSolver {
	
	private SolveForFields solveForFields;
	
	/**Solves the electrostatic Poisson equation with FFT assuming periodic boundaries.
	 * 
	 * <p>This method should be called every time when new particles
	 * are loaded into the simulation area (i.e. a new charge
	 * distribution is introduced) It calculates the electrostatic
	 * potential caused by this distribution, calculates the electric
	 * fields by applying the negative nabla operator and saves them
	 * in the field variables of the Grid class. Note that periodic
	 * boundaries are assumed by the transformation itself AND by the
	 * derivative of the potential!</p>
	 * @param grid Grid on which the calculation should be performed
	 */
	public void step(Grid grid, double timeStep) {
		
		//size of the array to be transformed
		int columns = grid.getNumCellsX();
		int rows = grid.getNumCellsY();
		double cellArea = grid.getCellWidth() * grid.getCellHeight();
		//JTransforms saves the imaginary part as a second row entry
		//therefore there must be twice as many rows
		double[][] trho = new double[columns][2*rows];
		double[][] phi = new double[columns][2*rows];

		DoubleFFT_2D fft = new DoubleFFT_2D(columns, rows);
		
		//prepare input for fft
		for(int i = 0; i < columns; i++) {
			for(int j = 0; j < rows; j++) {
				trho[i][2*j] = grid.getRho(i,j);
				trho[i][2*j+1] = 0;
			}
		}
		
		//perform Fourier transformation
		fft.complexForward(trho);
		
		//Solve Poisson equation in Fourier space
		//We omit the term with i,j=0 where d would become 0. This term only contributes a constant term
		//to the potential and can therefore be chosen arbitrarily.
		for(int i = 1; i < columns; i++) {
			for(int j = 0; j < rows; j++) {
				double d = (4 - 2 * Math.cos((2 * Math.PI * i) / grid.getNumCellsX()) - 2 * Math.cos((2 * Math.PI * j) / grid.getNumCellsY()));
				phi[i][2*j] = (cellArea * trho[i][2*j]) / d;
				phi[i][2*j+1] = (cellArea * trho[i][2*j+1]) / d;
			}
		}
		//i=0 but j!=0
		for(int j = 1; j < rows; j++) {
			double d = (2 - 2 * Math.cos((2 * Math.PI * j) / grid.getNumCellsY()));
			phi[0][2*j] = (cellArea * trho[0][2*j]) / d;
			phi[0][2*j+1] = (cellArea * trho[0][2*j+1]) / d;		
		}
		
		//perform inverse Fourier transform
		fft.complexInverse(phi, true);
		
		solveForFields = new SolveForFields(phi);
		cellIterator.execute(grid, solveForFields);
		
		solveForFields = null;
		phi = null;
		trho = null;
		
		return;
	}
	
	private class SolveForFields implements CellAction {

		private int CORRECTION;
		private double[][] phi;
		
		SolveForFields(double[][] phi) {
			this.phi = phi;
			
		}
		
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid grid, int x, int y) {
	
			//the electric field in x direction is equal to the negative derivative of the 
			//potential in x direction, analogous for y direction
			//using central difference, omitting imaginary part since it should be 0 anyway
			grid.setEx(x, y, -(phi[x+1][2*y] - phi[x-1][2*y]) / (2 * grid.getCellWidth()));
			grid.setEy(x, y, -(phi[x][2*(y+1)] - phi[x][2*(y-1)]) / (2 * grid.getCellHeight()));
			
		}
		
		private int index(int clientIdx) {
			return CORRECTION + clientIdx;
		}
	}
	
}
