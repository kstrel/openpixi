package org.openpixi.pixi.ui;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;

import org.openpixi.pixi.parallel.cellaccess.CellAction;
import org.openpixi.pixi.parallel.cellaccess.CellIterator;
import org.openpixi.pixi.physics.GeneralBoundaryType;
import org.openpixi.pixi.physics.Settings;
import org.openpixi.pixi.physics.grid.Grid;
import org.openpixi.pixi.physics.grid.Cell;

public class FieldPropagation {

	private static CellIterator cellIterator;
	private static BufferedWriter output;
	SetInitialFields setInitialFields = new SetInitialFields();
	SetIncidentWave setIncidentWave = new SetIncidentWave();
	OutputAction outputAction = new OutputAction();
	int interval;

	
	public static void setParameters(Settings stt) {
		
		stt.setIterations(400);
		stt.setSimulationWidth(100);
		stt.setSimulationHeight(100);
		stt.setGridCellsX(100);
		stt.setGridCellsY(100);
		
		stt.setTimeStep(0.6);
		stt.setSpeedOfLight(1);
		
		stt.setNumOfParticles(0);
		
		stt.setBoundary(GeneralBoundaryType.Periodic);
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {
		
		FieldPropagation instance = new FieldPropagation();
		
		try{
			output = writerFactory("/home/kirill/fieldpropagation.dat");
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		// Creates a settings class with the default parameters
		Settings settings = new Settings();
		setParameters(settings);
		
		Grid grid = new Grid(settings);
		
		//Accessible replica of the FieldSolver cellIterator
		cellIterator = settings.getCellIterator();
		cellIterator.setNormalMode(grid.getNumCellsX(), grid.getNumCellsY());
		
		//write header of file
		instance.interval = settings.getIterations() / 200;
		output.write("# " + settings.getGridCellsX() + "\t" + settings.getGridCellsY() + "\t" + settings.getIterations() + "\t" + instance.interval + "\n");
		
		//Set initial conditions and write them to file
		instance.setTimeStep(settings);
		cellIterator.execute(grid, instance.setIncidentWave);
		cellIterator.execute(grid, instance.outputAction);
		output.write("\n" + "\n");
		
		int next = instance.interval;
		for(int i = 0; i < settings.getIterations(); i++) {
			grid.updateGrid(settings.getTimeStep());
			if ( i == next ) {
				output.write("# " + i + "\n");
				cellIterator.execute(grid, instance.outputAction);
				output.write("\n" + "\n");
				next += instance.interval;
			}
		}
		
		closeStreams();
	}
	
	private class SetInitialFields implements CellAction {
		
		private double timeStep;
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid grid, int x, int y) {
			double kx = 2*Math.PI* 1 / (grid.getNumCellsX()*grid.getCellWidth());
			double ky = 0*2*Math.PI* 1 / (grid.getNumCellsX()*grid.getCellWidth());
			grid.setEx(x, y, 0);
			grid.setEy(x, y, Math.sin(-kx * (x + grid.getCellWidth()/2) + ky*y) + Math.sin(-10*kx * (x + grid.getCellWidth()/2) + ky*y));
			grid.setBzo(x, y, Math.sin(-kx * x + ky*y + kx * timeStep / 2) + Math.sin(-10*kx * x + ky*y + kx * timeStep / 2));
			grid.setBz(x, y, Math.sin(-kx * x + ky*y - kx * timeStep / 2) + Math.sin(- 10*kx * x + ky*y - kx * timeStep / 2));
		}
	}
	
	private class SetIncidentWave implements CellAction {
		
		private double timeStep;
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid grid, int x, int y) {
			if( ((int) (grid.getNumCellsX()/2 - 10) < x && x < (int) (grid.getNumCellsX()/2 + 10)) && 
					((int) (grid.getNumCellsY()/2 - 10) < y && y < (int) (grid.getNumCellsY()/2 + 10))) {
				double kx = 2*Math.PI* 1 / (grid.getNumCellsX()*grid.getCellWidth());
				grid.setEx(x, y, 0);
				grid.setEy(x, y, 10);
				grid.setBzo(x, y, 10);
				grid.setBz(x, y, 10);
			};
		}
	}
	
	private class OutputAction implements CellAction {
		
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid g, int x, int y) {
			try {
				double Ex = (g.getEx(x, y) + g.getEx(x, y-1))/2;
				double Ey = (g.getEy(x,y) + g.getEy(x-1, y))/2;
				double Bz = (g.getBzo(x, y) + g.getBz(x, y))/2;
				output.write(x + "\t" +  y + "\t" + Ex + "\t" + Ey + "\t" + Bz + "\t" + (Ex*Ex + Ey*Ey + Bz*Bz) + "\n");
//				output.write(x + "\t" +  y + "\t" + g.getEx(x, y) + "\t" + g.getEy(x,y) + "\t" + g.getBzo(x, y) + "\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private void setTimeStep(Settings stt) {
		setInitialFields.timeStep = stt.getTimeStep();
		setIncidentWave.timeStep = stt.getTimeStep();
	}
	
	public static void closeStreams() {
		try {
			output.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	
	public static BufferedWriter writerFactory(String path) throws IOException {
		File file;
		FileWriter fstream;
		BufferedWriter out;
	
		file = new File(path);
		fstream = new FileWriter(file);
		out = new BufferedWriter(fstream);
		return out;
	}
}


