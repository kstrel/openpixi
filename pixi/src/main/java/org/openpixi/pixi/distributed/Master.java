package org.openpixi.pixi.distributed;

import org.openpixi.pixi.distributed.ibis.MasterCommunicator;
import org.openpixi.pixi.distributed.ibis.IbisRegistry;
import org.openpixi.pixi.distributed.partitioning.Partitioner;
import org.openpixi.pixi.distributed.partitioning.SimplePartitioner;
import org.openpixi.pixi.physics.Particle;
import org.openpixi.pixi.physics.ParticleGridInitializer;
import org.openpixi.pixi.physics.Settings;
import org.openpixi.pixi.physics.grid.Cell;
import org.openpixi.pixi.physics.grid.Grid;
import org.openpixi.pixi.physics.util.IntBox;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

/**
 * Distributes the problem and collects the results.
 * To each worker we assign exactly one problem.
 */
public class Master {

	private MasterCommunicator communicator;
	private Settings settings;

	private Grid initialGrid;
	private List<Particle> initialParticles;

	private Grid finalGrid;
	private List<Particle> finalParticles;

	/**
	 * Problem decomposition table.
	 * N-th element in the table belongs to the n-th worker.
	 */
	IntBox[] partitions;

	Grid getInitialGrid() {
		return initialGrid;
	}

	List<Particle> getInitialParticles() {
		return initialParticles;
	}

	public Grid getFinalGrid() {
		return finalGrid;
	}

	public List<Particle> getFinalParticles() {
		return finalParticles;
	}


	public Master(IbisRegistry registry, Settings settings) throws Exception {
		this.settings = settings;
		communicator = new MasterCommunicator(registry);

		initialGrid = new Grid(settings);
		initialParticles = settings.getParticles();

		ParticleGridInitializer pgi = new ParticleGridInitializer();
		pgi.initialize(
				settings.getInterpolator(), settings.getPoissonSolver(),
				initialParticles, initialGrid);
	}


	public void distributeProblem() throws IOException {
		// Partition the problem
		Partitioner partitioner = new SimplePartitioner();
		partitions = partitioner.partition(
				settings.getGridCellsX(), settings.getGridCellsY(), settings.getNumOfNodes());

		// Partition the data
		List<List<Particle>> particlePartitions = partitionParticles(
				partitions, initialParticles);
		Cell[][][] gridPartitions = partitionGrid(partitions, initialGrid);

		// Send to each worker
		for (int workerID = 0; workerID < partitions.length; workerID++) {
			communicator.sendProblem(
					workerID, partitions,
					particlePartitions.get(workerID), gridPartitions[workerID]);
		}
	}


	/**
	 * Divides cells according to partitions.
	 */
	private Cell[][][] partitionGrid(IntBox[] partitions, Grid grid) {
		Cell[][][] gridPartitions = new Cell[partitions.length][][];
		for (int i = 0; i < partitions.length; ++i) {
			gridPartitions[i] = getSubgrid(partitions[i], grid);
		}
		return gridPartitions;
	}


	/**
	 * At the borders returns also the extra cells!
	 * Since there is some initialization done on the grid
	 * (density charge interpolation + poisson solver)
	 * we need to distribute also the extra cells.
	 */
	private Cell[][] getSubgrid(IntBox partition, Grid grid) {
		int startX = partition.xmin();
		int endX = partition.xmax();
		int startY = partition.ymin();
		int endY = partition.ymax();
		if (startX == 0) {
			startX -= Grid.EXTRA_CELLS_BEFORE_GRID;
		}
		if (endX == settings.getGridCellsX() - 1) {
			endX += Grid.EXTRA_CELLS_AFTER_GRID;
		}
		if (startY == 0) {
			startY -= Grid.EXTRA_CELLS_BEFORE_GRID;
		}
		if (endY == settings.getGridCellsY() - 1) {
			endY += Grid.EXTRA_CELLS_AFTER_GRID;
		}

		Cell[][] subGrid = new Cell[endX - startX + 1][endY - startY + 1];
		for (int x = startX; x <= endX ; x++) {
			for (int y = startY; y <= endY; y++) {
				subGrid[x - startX][y - startY] = grid.getCell(x,y);
			}
		}
		return  subGrid;
	}


	/**
	 * Divides particles according to partitions they belong to.
	 */
	private List<List<Particle>> partitionParticles(IntBox[] partitions, List<Particle> particles) {
		List<List<Particle>> particlePartitions = new ArrayList<List<Particle>>();
		for (int i = 0; i < partitions.length; ++i) {
			particlePartitions.add(new ArrayList<Particle>());
		}

		for (Particle p: particles) {
			int partitionIndex = -1;
			int cellX = (int)Math.floor(p.getX() / initialGrid.getCellWidth());
			int cellY = (int)Math.floor(p.getY() / initialGrid.getCellHeight());

			for  (int i = 0; i < partitions.length; ++i) {
				if (partitions[i].contains(cellX, cellY)) {
					partitionIndex = i;
				}
			}

			assert partitionIndex != -1;
			particlePartitions.get(partitionIndex).add(p);
		}
		return particlePartitions;
	}


	public void collectResults() throws Exception {
		ResultsHolder results = communicator.collectResults();
		finalParticles = assembleParticles(results.particlePartitions);
		finalGrid = assembleGrid(results.gridPartitions);
	}


	/**
	 * Puts together the particle lists coming from workers.
	 */
	private List<Particle> assembleParticles(List<List<Particle>> particlePartitions) {
		List<Particle> assembledParticles = new ArrayList<Particle>();
		for (List<Particle> particles: particlePartitions) {
			assembledParticles.addAll(particles);
		}
		return assembledParticles;
	}


	/**
	 * Puts together the subgrids coming from workers.
	 */
	private Grid assembleGrid(Cell[][][] gridPartitions) {
		int totalXCells =
				Grid.EXTRA_CELLS_BEFORE_GRID +
				settings.getGridCellsX() +
				Grid.EXTRA_CELLS_AFTER_GRID;
		int totalYCells =
				Grid.EXTRA_CELLS_BEFORE_GRID +
				settings.getGridCellsY() +
				Grid.EXTRA_CELLS_AFTER_GRID;

		Cell[][] cells = new Cell[totalXCells][totalYCells];
		for (int workerID = 0; workerID < partitions.length; ++workerID) {
			fillSubgrid(partitions[workerID], gridPartitions[workerID], cells);
		}

		return new Grid(cells, settings.getSimulationWidth(), settings.getSimulationHeight());
	}


	/**
	 * We fill also the extra cells.
	 * In case of hardwall boundary we have additional information in them.
	 */
	private void fillSubgrid(IntBox partition, Cell[][] subgrid, Cell[][] cells) {
		int startX = partition.xmin();
		int endX = partition.xmax();
		int startY = partition.ymin();
		int endY = partition.ymax();
		if (startX == 0) {
			startX -= Grid.EXTRA_CELLS_BEFORE_GRID;
		}
		if (endX == settings.getGridCellsX() - 1) {
			endX += Grid.EXTRA_CELLS_AFTER_GRID;
		}
		if (startY == 0) {
			startY -= Grid.EXTRA_CELLS_BEFORE_GRID;
		}
		if (endY == settings.getGridCellsY() - 1) {
			endY += Grid.EXTRA_CELLS_AFTER_GRID;
		}

		for (int x = startX; x <= endX ; x++) {
			for (int y = startY; y <= endY; y++) {
				cells[x + Grid.EXTRA_CELLS_BEFORE_GRID][y + Grid.EXTRA_CELLS_BEFORE_GRID] =
						subgrid[x - startX][y - startY];
			}
		}
	}
}
