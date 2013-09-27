package org.openpixi.pixi.ui;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;

import org.openpixi.pixi.parallel.cellaccess.CellAction;
import org.openpixi.pixi.parallel.cellaccess.CellIterator;
import org.openpixi.pixi.physics.GeneralBoundaryType;
import org.openpixi.pixi.physics.Settings;
import org.openpixi.pixi.physics.grid.Grid;
import org.openpixi.pixi.physics.grid.Cell;

public class FieldPropagation {

	private static CellIterator cellIterator;
	private static BufferedWriter output;
	private static BufferedWriter slicedOutput;
	private static BufferedWriter speedOutput;


	SetInitialFields setInitialFields = new SetInitialFields();
	SetIncidentWave setIncidentWave = new SetIncidentWave();
	SetSuperpositionOfFields setSuperpositionOfFields = new SetSuperpositionOfFields();
	OutputAction outputAction = new OutputAction();
	SlicedOutputAction slicedOutputAction = new SlicedOutputAction();
	int interval;
	ArrayList<Integer> speed = new ArrayList<Integer>();
	int n;

	
	public static void setParameters(Settings stt) {
		
		stt.setIterations(500);
		stt.setSimulationWidth(200);
		stt.setSimulationHeight(10);
		stt.setGridCellsX(200);
		stt.setGridCellsY(10);
			
		stt.setTimeStep(0.1);
		stt.setSpeedOfLight(1);
		
		stt.setNumOfParticles(0);
		
		stt.setBoundary(GeneralBoundaryType.Periodic);
	}
	
	public static void main(String[] args) throws FileNotFoundException, IOException, InterruptedException {

		// Creates a settings class with the default parameters
		Settings settings = new Settings();
		setParameters(settings);
		
		try{
			speedOutput = writerFactory("/home/kirill/speed.dat");
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		for(int i = 197; i < 198 ; i++) {
			speedOutput.write(i + "\t" + calculation(i, settings) + "\n");
		}
		
		speedOutput.close();
		return;
	}
	
	public static double calculation(int number, Settings settings) throws FileNotFoundException, IOException, InterruptedException {
		
		FieldPropagation instance = new FieldPropagation();
		
		instance.n = number;
		
		try{
			output = writerFactory("/home/kirill/fieldpropagation.dat");
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		try{
			slicedOutput = writerFactory("/home/kirill/timeslices.dat");
		} catch(IOException e) {
			e.printStackTrace();
		}
		
		Grid grid = new Grid(settings);
		
		//Accessible replica of the FieldSolver cellIterator
		cellIterator = settings.getCellIterator();
		cellIterator.setNormalMode(grid.getNumCellsX(), grid.getNumCellsY());
		
		//write header of file
		instance.interval = settings.getIterations() / 500;
		output.write("# " + settings.getGridCellsX() + "\t" + settings.getGridCellsY() + "\t" + settings.getIterations() + "\t" + instance.interval + "\n");
		
		
		//Set initial conditions and write them to two files
		instance.setTimeStep(settings);
		cellIterator.execute(grid, instance.setInitialFields);
		int index = instance.getZeroCrossing(grid, 0);
		instance.speed.add(index);
		instance.speed.add(0);
		
		
		cellIterator.execute(grid, instance.outputAction);
		output.write("\n" + "\n");
		cellIterator.execute(grid, instance.slicedOutputAction);
		instance.slicedOutputAction.i++;;
		
		int next = instance.interval;
		for(int i = 0; i < settings.getIterations(); i++) {
			grid.updateGrid(settings.getTimeStep());
			
			index = instance.getZeroCrossing(grid, index);
			instance.speed.add(index);
			instance.speed.add(i);
			
			if ( i == next ) {
				output.write("# " + i + "\n");
				cellIterator.execute(grid, instance.outputAction);
				output.write("\n" + "\n");
				
				cellIterator.execute(grid, instance.slicedOutputAction);
				instance.slicedOutputAction.i++;
				
				next += instance.interval;
			}
		}
		
		double v = instance.calculateSpeed(settings.getTimeStep(), grid.getCellWidth(), grid.getNumCellsX(), instance.speed);
		
		closeStreams();
		
		return v;
	}
	
	private class SetInitialFields implements CellAction {
		
		private double timeStep;
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();

		}
		public void execute(Grid grid, int x, int y) {
			
			double kx = n*2*Math.PI / (grid.getCellWidth()*grid.getNumCellsX());
			double ky = 0*2*Math.PI / (grid.getNumCellsY()*grid.getCellHeight());
			grid.setEx(x, y, 0);
			grid.setEy(x, y, Math.sin((kx * (x + grid.getCellWidth()/2) + ky*y)));
			grid.setBzo(x, y, Math.sin((kx * x + ky*y + kx * timeStep / 2)));
			grid.setBz(x, y, Math.sin((kx * x + ky*y - kx * timeStep / 2)));
			
//			grid.setEy(x, y, Math.sin(kx * (x + grid.getCellWidth()/2) + ky*y) + Math.sin(20 * (kx * (x + grid.getCellWidth()/2) + ky*y)));
//			grid.setBzo(x, y, Math.sin(kx * x + ky*y + kx * timeStep / 2) + Math.sin(20 * (kx * x + ky*y + kx * timeStep / 2)));
//			grid.setBz(x, y, Math.sin(kx * x + ky*y - kx * timeStep / 2) + Math.sin(20 * (kx * x + ky*y - kx * timeStep / 2)));
		}
	}
	
	private class SetSuperpositionOfFields implements CellAction {
		
		private double timeStep;
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid grid, int x, int y) {
			
			double Bz = 0;
			for (int i = 1; i < 100; i++) {
				double kx = 2*Math.PI* 1 / (grid.getCellWidth()*grid.getNumCellsX());
				double ky = 0*2*Math.PI* 1 / (grid.getNumCellsY()*grid.getCellHeight());
				grid.setEx(x, y, 0);
				grid.addEy(x, y, (Math.sin(12*kx*i)/i)*Math.sin(i*(kx * (x + grid.getCellWidth()/2) + ky*y)));
				Bz += (Math.sin(12*kx*i)/i)*Math.sin(i*(kx * x + ky*y + kx * timeStep / 2));
				grid.addBz(x, y, (Math.sin(12*kx*i)/i)*Math.sin(i*(kx * x + ky*y - kx * timeStep / 2)));
			}
			grid.setBzo(x, y, Bz);

//			grid.setEy(x, y, Math.sin(kx * (x + grid.getCellWidth()/2) + ky*y) + Math.sin(20 * (kx * (x + grid.getCellWidth()/2) + ky*y)));
//			grid.setBzo(x, y, Math.sin(kx * x + ky*y + kx * timeStep / 2) + Math.sin(20 * (kx * x + ky*y + kx * timeStep / 2)));
//			grid.setBz(x, y, Math.sin(kx * x + ky*y - kx * timeStep / 2) + Math.sin(20 * (kx * x + ky*y - kx * timeStep / 2)));
		}
	}
	
	private class SetIncidentWave implements CellAction {
		
		private double timeStep;
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid grid, int x, int y) {
			if( ((int) (grid.getNumCellsX()/2 - 10) < x && x < (int) (grid.getNumCellsX()/2 + 10))) { //&&					((int) (grid.getNumCellsY()/2 -10) < y && y < (int) (grid.getNumCellsY()/2 + 10))) {
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
//				output.write(x + "\t" +  y + "\t" + g.getEx(x, y) + "\t" + g.getEy(x,y) + "\t" + g.getBz(x, y) + "\t" + (g.getEx(x, y)*g.getEx(x, y)+g.getEy(x,y)*g.getEy(x,y)+ g.getBz(x, y)* g.getBz(x, y)) + "\n");
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
	
	private class SlicedOutputAction implements CellAction {
		
		int i = 0;
		public void execute(Cell cell) {
			throw new UnsupportedOperationException();
		}

		public void execute(Grid g, int x, int y) {
			if (y==0) {
				try {
					double Ex = (g.getEx(x, y) + g.getEx(x, y-1))/2;
					double Ey = (g.getEy(x,y) + g.getEy(x-1, y))/2;
					double Bz = (g.getBzo(x, y) + g.getBz(x, y))/2;
					slicedOutput.write(x + "\t" +  i + "\t" + (Ex*Ex + Ey*Ey + Bz*Bz) + "\n");
//					slicedOutput.write(x + "\t" +  i + "\t" + (g.getEx(x, y)*g.getEx(x, y)+g.getEy(x,y)*g.getEy(x,y)+ g.getBz(x, y)* g.getBz(x, y)) + "\n");

				} catch (IOException e) {
					e.printStackTrace();
				}
			}
		}
	}
	
	private double calculateSpeed(double timeStep, double cellWidth, int n, ArrayList<Integer> speed) {
		
		int location = 0;
		int time = 0;
		int iterations = 0;
		
		for(int i = 0; i < speed.size(); i += 2) {
			if (speed.get(i) > (n-10)) {
				break;
			}
			
			location = speed.get(i);
			time = speed.get(i+1);
			
			iterations++;
		}
		
		double v = cellWidth*(location - speed.get(0))/( (double) time*timeStep);
				
		return v;
	}

	private int getZeroCrossing(Grid g, int j) {
		
		for(int i = j; i < g.getNumCellsX(); i++) {
			
			if ((int) Math.signum(g.getBz(i,0)) != (int) Math.signum(g.getBz(i+1, 0))) {
					return i;
				}
		}
		
		return 0;
	}
	
	private void setTimeStep(Settings stt) {
		setInitialFields.timeStep = stt.getTimeStep();
		setIncidentWave.timeStep = stt.getTimeStep();
	}
	
	public static void closeStreams() {
		try {
			output.close();
			slicedOutput.close();
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


