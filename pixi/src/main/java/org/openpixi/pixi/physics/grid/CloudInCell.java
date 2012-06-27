package org.openpixi.pixi.physics.grid;

import java.util.ArrayList;
import org.openpixi.pixi.physics.Debug;
import org.openpixi.pixi.physics.Particle;

public class CloudInCell extends Interpolator {

	@Override
	public void interpolateToGrid(ArrayList<Particle> particles, Grid g) {
		g.resetCurrentAndCharge();

		for(Particle p : particles)
		{
			int xCellPosition = (int) (Math.floor((p.getX() / g.getCellWidth())));
			int yCellPosition = (int) (Math.floor((p.getY() / g.getCellHeight())));

			int xCellPosition2 = xCellPosition;
			int yCellPosition2 = yCellPosition;

			if (Debug.asserts) {
				// Assert conditions for interpolation
				assert xCellPosition2 * g.getCellWidth() > p.getX() : p.getX();
				assert p.getX() > (xCellPosition2 - 1) * g.getCellWidth() : p.getX();
				assert yCellPosition2 * g.getCellHeight() > p.getY() : p.getY();
				assert p.getY() > (yCellPosition2 - 1) * g.getCellHeight() : p.getY();
			}

			g.addJx(xCellPosition, yCellPosition,
					p.getCharge() * p.getVx() *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) /
					(g.getCellWidth() * g.getCellHeight()));
			g.addJx(xCellPosition + 1, yCellPosition,
					p.getCharge() * p.getVx() *
					(p.getX() - (xCellPosition2-1) * g.getCellWidth()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) /
					(g.getCellWidth() * g.getCellHeight()));
			g.addJx(xCellPosition,yCellPosition + 1,
					p.getCharge() * p.getVx() *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(p.getY() - (yCellPosition2-1) * g.getCellHeight()) /
					(g.getCellWidth() * g.getCellHeight()));
			g.addJx(xCellPosition + 1,yCellPosition + 1,
					p.getCharge() * p.getVx() *
					(p.getX() - (xCellPosition2-1) * g.getCellWidth()) *
					(p.getY() - (yCellPosition2-1) * g.getCellHeight()) /
					(g.getCellWidth() * g.getCellHeight()));

			g.addJy(xCellPosition, yCellPosition,
					p.getCharge() * p.getVy() *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) /
					(g.getCellWidth() * g.getCellHeight()));
			g.addJy(xCellPosition + 1, yCellPosition,
					p.getCharge() * p.getVy() *
					(p.getX() - (xCellPosition2-1) * g.getCellWidth()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) /
					(g.getCellWidth() * g.getCellHeight()));
			g.addJy(xCellPosition, yCellPosition + 1,
					p.getCharge() * p.getVy() *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(p.getY() - (yCellPosition2-1) * g.getCellHeight()) /
					(g.getCellWidth() * g.getCellHeight()));
			g.addJy(xCellPosition + 1, yCellPosition + 1,
					p.getCharge() * p.getVy() *
					(p.getX() - (xCellPosition2-1) * g.getCellWidth()) *
					(p.getY() - (yCellPosition2-1) * g.getCellHeight()) /
					(g.getCellWidth() * g.getCellHeight()));
		}

	}

	public void interpolateToParticle(ArrayList<Particle> particles, Grid g) {

		for (int i = 0; i < particles.size(); i++) {

		Particle p = g.simulation.particles.get(i);
		int xCellPosition = (int) Math.floor(p.getX() / g.getCellWidth());
		int yCellPosition = (int) Math.floor(p.getY() / g.getCellHeight());

		int xCellPosition2 = xCellPosition;
		int yCellPosition2 = yCellPosition;

		if (Debug.asserts) {
			// Assert conditions for interpolation
			assert xCellPosition2 * g.getCellWidth() > p.getX() : p.getX();
			assert p.getX() > (xCellPosition2 - 1) * g.getCellWidth() : p.getX();
			assert yCellPosition2 * g.getCellHeight() > p.getY() : p.getY();
			assert p.getY() > (yCellPosition2 - 1) * g.getCellHeight() : p.getY();
		}

			particles.get(i).setEx((g.getEx(xCellPosition, yCellPosition) *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) +
					g.getEx(xCellPosition + 1, yCellPosition) *
					(p.getX() - (xCellPosition2 - 1) * g.getCellWidth()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) +
					g.getEx(xCellPosition, yCellPosition + 1) *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(p.getY() - (yCellPosition2 - 1) * g.getCellHeight()) +
					g.getEx(xCellPosition + 1, yCellPosition + 1) *
					(p.getX() - (xCellPosition2 - 1) * g.getCellWidth()) *
					(p.getY() - (yCellPosition2 - 1) * g.getCellHeight())) /
					(g.getCellWidth() * g.getCellHeight()));

			particles.get(i).setEy((g.getEy(xCellPosition, yCellPosition) *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) +
					g.getEy(xCellPosition + 1, yCellPosition) *
					(p.getX() - (xCellPosition2 - 1) * g.getCellWidth()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) +
					g.getEy(xCellPosition, yCellPosition + 1) *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(p.getY() - (yCellPosition2 - 1) * g.getCellHeight()) +
					g.getEy(xCellPosition + 1, yCellPosition + 1) *
					(p.getX() - (xCellPosition2 - 1) * g.getCellWidth()) *
					(p.getY() - (yCellPosition2 - 1) * g.getCellHeight())) /
					(g.getCellWidth() * g.getCellHeight()));

			particles.get(i).setBz((g.getBz(xCellPosition, yCellPosition) *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) +
					g.getBz(xCellPosition + 1, yCellPosition) *
					(p.getX() - (xCellPosition2 - 1) * g.getCellWidth()) *
					(yCellPosition2 * g.getCellHeight() - p.getY()) +
					g.getBz(xCellPosition, yCellPosition + 1) *
					(xCellPosition2 * g.getCellWidth() - p.getX()) *
					(p.getY() - (yCellPosition2 - 1) * g.getCellHeight()) +
					g.getBz(xCellPosition + 1, yCellPosition + 1) *
					(p.getX() - (xCellPosition2 -1) * g.getCellWidth()) *
					(p.getY() - (yCellPosition2 -1) * g.getCellHeight())) /
					(g.getCellWidth() * g.getCellHeight()));

		}

	}


}
